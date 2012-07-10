/*
 *
 */
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  /*--------------------------------------------------------------*/
  /*- process options  -------------------------------------------*/
  /*--------------------------------------------------------------*/
  char *Geometry;
  double LBV[2][2]={ {1.0, 0.0}, {0.0, 1.0} };
  double P[2]={0.0, 0.0};
  cdouble Omega;
  /* name, type, #args, max_instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"Geometry",  PA_STRING,  1, 1, (void *)&GeoFileName, 0,        ".scuffgeo file"},
     {"LBV1",      PA_DOUBLE,  2, 1, (void *)LBV[0],       0,        "lattice basis vector 1"},
     {"LBV2",      PA_DOUBLE,  2, 1, (void *)LBV[1],       0,        "lattice basis vector 2"},
     {"P",         PA_DOUBLE,  2, 1, (void *)P,            0,        "bloch wavevector"},
     {"Omega",     PA_CDOUBLE, 1, 1, (void *)&Omega,       0,        "angular frequency"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (Geometry==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");
  SetLogFileName("tPBC.log");
  RWGGeometry::SetLogLevel(SCUFF_VERBOSELOGGING);

  RWGGeometry *G=new RWGGeometry(GeoFileName);
  PBCGeometry *PG=new PBCGeometry(G, LBV);
  
  HMatrix *M = G->AllocateBEMMatrix();
  PG->AssembleBEMMatrix_PBC(Omega, P, M);
  M->LUFactorize();

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  cdouble E0[3]  = {1.0, 0.0,  0.0};
  double nHat[3] = {0.0, 0.0, -1.0};
  PlaneWave PW(E0, nHat);
  HMatrix *KN = G->AllocateRHSVector();
  AssembleRHSVector(Omega, &PW, KN);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  M->LUSolve(KN);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double XMin = -LBV[0][0]; 
  double XMax =  2.0*LBV[0][0]; 
  double DeltaX = (XMax-XMin) / 50.0;

  double YMin = -LBV[1][1]; 
  double YMax =  2.0*LBV[1][1]; 
  double DeltaX = (YMax-YMin) / 5.0;
 
  double X[3];
  cdouble EH[6];
  FILE *f=vfopen("%s.out","w",GetFileBase(GeoFileName));
  X[2]=1.0;
  for(X[1] = YMin; X[1]<=YMax; X[1]+=DeltaY );
   for(X[0] = XMin; X[0]<=XMax; X[0]+=DeltaX );
     { 
       fprintf(f,"%e %e %e ",X[0],X[1],X[2]);

       PG->GetScatteredFields(X, EH);
       EH[0] += 1.0;
       EH[1] -= 1.0/ZVAC;
       fprintf(f,"%e %e %e %e %e %e ", real(X[0]),imag(X[0]),
                                       real(X[1]),imag(X[1]),
                                       real(X[2]),imag(X[2]));

       fprintf(f,"\n");
       if( float(X[0]) == float(XMax) )
        fprintf(f,"\n\n");
       
     };
  fclose(f);

}
