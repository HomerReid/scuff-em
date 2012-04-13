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
  cdouble Omega=0.01;
  cdouble EpsSlab=2.0;
  int PlotCurrents=0;
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { {"geometry",      PA_STRING,  1, 1, (void *)&GeoFile,      0, "geometry file"},
     {"Omega",         PA_CDOUBLE, 1, 1, (void *)&Omega,        0, "angular frequency"},
     {"EpsSlab",       PA_CDOUBLE, 1, 1, (void *)&EpsSlab,      0, "slab permittivity"},
     {"PlotCurrents",  PA_BOOL,    0, 1, (void *)&PlotCurrents, 0, "plot surface currents"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (GeoFile==0)
   OSUsage(argv[0], OSArray, "--geometry option is mandatory");

  SetLogFileName("SlabTransmission.log");

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
  G->Objects[0]->MP->SetEps(EpsSlab);
  HMatrix *M=G->AssembleBEMMatrix(Omega, M);

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  cdouble E0[3]  = { 1.0, 0.0, 0.0 };
  double nHat[3] = { 0.0, 0.0, 1.0 };
  PlaneWave *PW = new PlaneWave( E0, nHat );
  HVector *KN=G->AssembleRHSVector(Omega, PW);

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  M->LUFactorize();
  M->LUSolve(KN);

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  double XAbove[3]  = {0.0, 0.0,  1.5 }; // eval point above the slab
  double XInside[3] = {0.0, 0.0,  0.5 }; // eval point inside the slab 
  double XBelow[3]  = {0.0, 0.0, -0.5 }; // eval point below the slab 
  cdouble EH[6], EH2[6];
  
  FILE *f=vfopen("%s.out","a",GetFileBase(GeoFile));
  fprintf(f,"%e %e %e ",real(Omega),real(EpsVal),imag(EpsVal));

  G->GetFields( XAbove, Omega, KN, EH);
  PW->GetFields( XBelow,  Omega, KN, EH2);
  EH[0]+=EH2[0];
  EH[4]+=EH2[4];
  fprintf(f,"%e %e %e %e ",abs(EH[0]),arg(EH[0])
                           abs(EH[4]),arg(EH[4]));

  G->GetFields( XInside, Omega, KN, EH);
  fprintf(f,"%e %e %e %e ",abs(EH[0]),arg(EH[0])
                           abs(EH[4]),arg(EH[4]));


  G->GetFields( XBelow,  Omega, KN, EH);
  PW->GetFields( XBelow,  Omega, KN, EH2);
  EH[0]+=EH2[0];
  EH[4]+=EH2[4];
  fprintf(f,"%e %e %e %e ",abs(EH[0]),arg(EH[0])
                           abs(EH[4]),arg(EH[4]));
  fprintf(f,"\n");
  fclose(f);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (PlotCurrents)
   { double XYMin=-2.0;
     double XYMax=+2.0;
     double XYStep=(XYMax - XYMin) / (50.0);
     double X[3];
     cdouble KNX[3];

     FILE *f2=vfopen("%s.currents","w",GetFileBase(GeoFile));

     X[2]=0.0;
     for(X[0]=XYMin; X[0]<=XYMax; X[0]+=XYStep)
      for(X[1]=XYMin; X[1]<=XYMax; X[1]+=XYStep)
       {
         G->EvalCurrentDistribution(X, KN, KNX);
         fprintf(f2,"%e %e %e %e %e %e %e %e %e \n",
                     X[0],X[1],X[2],
                     abs(KNX[0]),abs(KNX[1]),abs(KNX[2]),
                     abs(KNX[3]),abs(KNX[4]),abs(KNX[5]));
            
       };

     X[2]=1.0;
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

   };

}
