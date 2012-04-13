/*
 *
 */
#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#include <libscuff.h>

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  char *GeoFile=0;
  cdouble Omega=1.0;
  int PlotCurrents=0;
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { {"geometry",       PA_STRING,  1, 1, (void *)&GeoFile,      0, "geometry file"},
     {"Omega",          PA_CDOUBLE, 1, 1, (void *)&Omega,        0, "angular frequency"},
     {"PlotCurrents",   PA_BOOL,    0, 1, (void *)&PlotCurrents, 0, "plot surface currents"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (GeoFile==0)
   OSUsage(argv[0], OSArray, "--geometry option is mandatory");

  SetLogFileName("Fresnel.log");

  char *CacheFileName=vstrdup("%s.scuffcache",GetFileBase(GeoFile));

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  PreloadGlobalFIPPICache(CacheFileName);
  RWGGeometry *G=new RWGGeometry(GeoFile);
  G->SetLogLevel(SCUFF_VERBOSELOGGING);

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  HMatrix *M = G->AllocateBEMMatrix();

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  cdouble E0[3]   = { 1.0, 0.0,  0.0 };
  double nHat[3] = { 0.0, 0.0, -1.0 };
  PlaneWave *PW = new PlaneWave( E0, nHat );
  HVector *RHS=G->AssembleRHSVector(Omega, PW);

  HVector *KN=G->AllocateRHSVector();

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  double XAbove[3]  = {0.0, 0.0,  1.0 }; // eval point above the slab
  double XInside[3] = {0.0, 0.0, -1.0 }; // eval point inside the slab 
  cdouble EHAbove[6], EHInside[6];
  
  FILE *f=vfopen("%s.out","w",GetFileBase(GeoFile));

  cdouble EpsSlab(0.0, 5.0);
  //double EpsMult = exp( log(20.0/2.0) / 20.0 );
  double EpsStep = 1.0;
  int StoredCache=0;
  for( real(EpsSlab)=2.0; real(EpsSlab)<=20.0; real(EpsSlab)+=EpsStep )
   { 
     G->Objects[0]->MP->SetEps(EpsSlab); 

     G->AssembleBEMMatrix(Omega, M);
     M->LUFactorize();
     KN->Copy(RHS);
     M->LUSolve(KN);

     G->GetFields( XAbove, Omega, KN, EHAbove);
     G->GetFields( XInside, Omega, KN, EHInside);

     fprintf(f,"%e %e %e %e %e \n", real(EpsSlab), 
                abs(EHAbove[0]),  abs(EHAbove[4]),
                abs(EHInside[0]),  abs(EHInside[4]));

     if (!StoredCache)
      { StoreGlobalFIPPICache(CacheFileName);
        StoredCache=1;
      };

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     if (PlotCurrents)
      { double XYMin=-2.0;
        double XYMax=+2.0;
        double XYStep=(XYMax - XYMin) / (50.0);
        double X[3]={0.0, 0.0, 0.0};
        cdouble KNX[3];
        FILE *f2=vfopen("%s.currents","w",GetFileBase(GeoFile));
        for(X[0]=XYMin; X[0]<=XYMax; X[0]+=XYStep)
         for(X[1]=XYMin; X[1]<=XYMax; X[1]+=XYStep)
          {
            G->EvalCurrentDistribution(X, KN, KNX);
            fprintf(f2,"%e %e %e %e %e %e %e %e %e \n",
                        X[0],X[1],X[2],
                        abs(KNX[0]),abs(KNX[1]),abs(KNX[2]),
                        abs(KNX[3]),abs(KNX[4]),abs(KNX[5]));
            
          };
        fclose(f2);
        exit(1);
      };

   };

  fclose(f);

}
