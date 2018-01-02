/*
 * tInterp1D.cc   -- second test program for Interp1D class
 *
 * homer reid     -- 3/14/2011
 */
#include <stdio.h>
#include <stdlib.h>

#include <libhrutil.h>
#include <libhmat.h>
#include "libMDInterp.h"

#define NFUN 2

/***************************************************************/
/* the function to be interpolated   ***************************/
/***************************************************************/
#define A 0.5
#define B 0.1
void Phi(double X, void *UserData, double *PhiVD)
{

  PhiVD[0] = exp( -A*X*X  );      /* Phi_1       */
  PhiVD[1] = -2.0*A*X*PhiVD[0];   /* dPhi_1 / dX */
  PhiVD[2] = exp( -B*X*X  );      /* Phi_2       */
  PhiVD[3] = -2.0*B*X*PhiVD[0];   /* dPhi_2 / dX */

}

/***************************************************************/ 
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  Interp1D *I1D;

  /*--------------------------------------------------------------*/
  /*- process command-line arguments -----------------------------*/
  /*--------------------------------------------------------------*/
  int nThread;
  double XMin, XMax;
  int N;
  char *DataFile;
  int xCol, yCol;
  ArgStruct ASArray[]=
   { {"XMin",     PA_DOUBLE, (void *)&XMin,     "-5.0", "XMin"},
     {"XMax",     PA_DOUBLE, (void *)&XMax,     "+5.0", "XMax"},
     {"N",        PA_INT,    (void *)&N,         "100", "N"},
     {"DataFile", PA_STRING, (void *)&DataFile,      0, "Data file "},
     {"xCol",     PA_INT,    (void *)&xCol,        "1", "x column index (one-based)"},
     {"yCol",     PA_INT,    (void *)&yCol,        "2", "y column index (one-based)"},
     {"nThread",  PA_INT,    (void *)&nThread,    "-1", "number of threads"},
     {0,0,0,0,0}
   };
  ProcessArguments(argc, argv, ASArray);
  if (nThread==-1)
   nThread=GetNumThreads();

  /*--------------------------------------------------------------*/ 
  /*--------------------------------------------------------------*/ 
  /*--------------------------------------------------------------*/ 
  if (DataFile)
   { 
     HMatrix *M=new HMatrix(DataFile,LHM_TEXT);
     if (M->ErrMsg)
      ErrExit(M->ErrMsg);
     if (xCol<1 || xCol>M->NC)
      ErrExit("invalid x column specification %i (file has %i columns)",xCol,M->NC);
     if (yCol<1 || yCol>M->NC)
      ErrExit("invalid y column specification %i (file has %i columns)",yCol,M->NC);

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     int NumPoints=M->NR;
     double *XData = new double[NumPoints], *YData = new double[NumPoints];
     int nr;
     for(nr=0; nr<NumPoints; nr++)
      { XData[nr]=M->GetEntryD(nr,xCol-1);
        YData[nr]=M->GetEntryD(nr,yCol-1);
      };
     I1D=new Interp1D(XData, YData, NumPoints, 1);

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     FILE *f=fopen("tInterp1D.out","w");
     for(nr=0; nr<NumPoints; nr++)
      fprintf(f,"%e %e \n",XData[nr],YData[nr]);
     fprintf(f,"\n\n");
     double X;
     // double DeltaX=XData[NumPoints-1] - XData[0];
    // for(X=XData[0] - 0.1*DeltaX; X<XData[0]+1.1*DeltaX; X+=0.001*DeltaX)
 //     fprintf(f,"%e %e \n",X,I1D->Evaluate(X));
     double XMult=pow( 10.0, 1.0/100.0);
printf("NumPoints=%i.\n",NumPoints);
     for(X=1e-6; X<10*XData[NumPoints-1]; X*=XMult)
      fprintf(f,"%e %e \n",X,I1D->Evaluate(X));
     fclose(f);

     f=popen("gnuplot -persist","w");
     fprintf(f,"set logscale xy\n");
     fprintf(f,"plot 'tInterp1D.out' i 0 u 1:2 w p pt 7 ps 1.5,");
     fprintf(f," '' i 1 u 1:2 w l\n");
     fclose(f);
   
   }
  else
   {
     /*--------------------------------------------------------------*/ 
     /*--------------------------------------------------------------*/ 
     /*--------------------------------------------------------------*/ 
     I1D=new Interp1D(XMin, XMax, N, NFUN, nThread, Phi, 0);
     FILE *f;
     double X;
     double PhiExact[4], PhiInterp[2];
     double DX=I1D->DX;
  
     f=CreateUniqueFile("tInterp1D.out",1);
     fprintf(f,"# columns: \n");
     fprintf(f,"# 1 X\n");
     fprintf(f,"# 2 F1(x) (exact) \n");
     fprintf(f,"# 3 F2(x) (exact) \n");
     fprintf(f,"# 4 F1(x) (interpolated) \n");
     fprintf(f,"# 5 F2(x) (interpolated) \n");
     for(X=XMin - 2.0*DX; X<XMax+2.0*DX; X+=0.1*DX)
      { 
        Phi(X, 0, PhiExact);
        I1D->Evaluate(X, PhiInterp);
  
        fprintf(f,"%e %e %e %e %e \n",X,PhiExact[0],PhiExact[2],PhiInterp[0],PhiInterp[1]);
  
       };
     fclose(f);

   };
  
}
