/*
 * tInterpND
 *
 * homer reid      -- 3/14/2011
 */
#include <stdio.h>
#include <stdlib.h>

#include <libhrutil.h>
#include <libhmat.h>

#include "libMDInterp.h"

#include <readline/readline.h>
#include <readline/history.h>

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define NFUN 2
#define X10 0.2
#define X20 0.3
#define X30 0.4
void Phi(double *X, void *UserData, double *PhiVD)
{
  int D = *((int *)UserData);
  static double Alpha[4]={0.2,0.3,0.4,0.5};
  static double X0[4]={-0.1,0.2,0.35,0.5};
  double XBar[4];

  double ExpArg[2] = {0.0, 0.0};
  for(int d=0; d<D; d++)
   { XBar[d]=X[d]-X0[d];
     ExpArg[0] += Alpha[d]*XBar[d];
     ExpArg[1] += Alpha[d]*XBar[d]*XBar[d];
   };

  int NVD=(1<<D);

  double ExpFac[2];
  ExpFac[0] = exp(-ExpArg[0]);
  ExpFac[1] = exp(-ExpArg[1]);

  iVec Twos(D,2);
  LOOP_OVER_IVECS(nVD, sigmaVec, Twos)
   { double Factor[2]={1.0, 1.0};
     for(int d=0; d<D; d++)
      if (sigmaVec[d]==1)
       { Factor[0]*=-1.0*Alpha[d];
         Factor[1]*=-2.0*Alpha[d]*XBar[d];
       };
     PhiVD[0*NVD + nVD] = Factor[0]*ExpFac[0];
     PhiVD[1*NVD + nVD] = Factor[1]*ExpFac[1];
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /*--------------------------------------------------------------*/
  /*- process command-line arguments -----------------------------*/
  /*--------------------------------------------------------------*/
  char *XXNFile=0;
  bool Console=false;
  /* name, type, #args, max_instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"XXNFile",  PA_STRING,  1, 1, (void *)&XXNFile,  0,  ""},
     {"console",  PA_BOOL,    0, 1, (void *)&Console,  0,  ""},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  //if (GeoFile==0)
  // OSUsage(argv[0],OSArray,"--geometry option is mandatory");

  // default grid
  int D=2;
  dVec XMin(D), XMax(D);
  iVec NVec(D);
  XMin[0] =  0.0;  XMax[0] = 2.0;  NVec[0] = 5;
  XMin[1] = -1.0;  XMax[1] = 3.0;  NVec[1] = 3;

  if (XXNFile)
   { HMatrix *M=new HMatrix(XXNFile);
     D=M->NR;
     XMin.resize(D);
     XMax.resize(D);
     NVec.resize(D);
     for(int d=0; d<D; d++)
      { XMin[d] = M->GetEntryD(d,0);
        XMax[d] = M->GetEntryD(d,1);
        NVec[d] = (int)(M->GetEntryD(d,2));
      };
   };

  dVec DX(D);
  FILE *f=fopen("/tmp/InterpError.out","w");
  double *MeanRelError = new double[D*NFUN];
  double *MeanAbsError = new double[D*NFUN];
  for(double dx=1.0; dx>=1.0e-3; dx*=0.5)
   for(double dy=1.0; dy>=1.0e-3; dy*=0.5)
    { DX[0]=dx;
      DX[1]=dy;
      HMatrix PhiVEMatrix(NFUN,2);
      double Err = GetInterpolationError(Phi, (void *)&D, NFUN, XMin, DX,
                                         MeanRelError, MeanAbsError);
      fprintf(f,"%e %e %e ",dx,dy,Err);
      for(int nfd=0; nfd<D*NFUN; nfd++)
       fprintf(f,"%e %e ",MeanRelError[nfd], MeanAbsError[nfd]);
      fprintf(f,"\n");
    }
  fclose(f);

  InterpND Interp(XMin, XMax, NVec, NFUN, Phi, (void *)&D);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (Console)
   { srand48(time(0));
     char *p;
     do
      { 
        p=readline("xVec: "); 
      } while(!p);
     add_history(p);
     write_history(0);
     int D=Interp.D;
     dVec xVec(D);
     for(int d=0; d<D; d++)
      xVec[d] = XMin[d] + drand48()*(XMax[d]-XMin[d]);

     switch(D)
      { case 1: sscanf(p,"%le",&(xVec[0])); 
                break;
        case 2: sscanf(p,"%le %le",&(xVec[0]), &(xVec[1])); 
                break;
        case 3: sscanf(p,"%le %le %le",&(xVec[0]), &(xVec[1]), &(xVec[2]));
                break;
      };

     printf("at x={");
     for(int d=0; d<D; d++)
      printf("%g ",xVec[d]);
     printf("}: \n");
     //printf("Grid cell %i ",Interp->GetCellIndex(

     int NVD = Interp.NVD;
     double PhiVDExact[NFUN*NVD], PhiExact[NFUN], PhiHR[NFUN];
     Phi(&(xVec[0]), (void *)&D, PhiVDExact);
     for(int nf=0; nf<NFUN; nf++)
      PhiExact[nf]=PhiVDExact[nf*NVD];
     
     Interp.Evaluate(&(xVec[0]),PhiHR);
     Compare(PhiExact, PhiHR, NFUN, "Exact", "Interp");
   }

}
