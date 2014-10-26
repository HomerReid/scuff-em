/*
 * tlibInterp3D.cc -- test program for libInterp3D 
 *
 * homer reid      -- 3/14/2011
 */
#include <stdio.h>
#include <stdlib.h>

#include <libhrutil.h>
#include <libhmat.h>
#include "libMDInterp.h"

/***************************************************************/
/* the functions to be interpolated: just two gaussians with   */
/* different widths.                                           */
/***************************************************************/
#define NFUN 2
#define A1 0.2
#define A2 0.3
#define A3 0.4
#define B1 0.5
#define B2 0.6
#define B3 0.7
#define X10 0.2
#define X20 0.3
#define X30 0.4
void Phi(double X1, double X2, double X3, void *UserData, double *PhiVD)
{
  double X1Bar, X2Bar, X3Bar;

  X1Bar=X1-X10;
  X2Bar=X2-X20;
  X3Bar=X3-X30;

  /* first function and its first derivatives */
  PhiVD[0] = exp( -A1*X1Bar*X1Bar - A2*X2Bar*X2Bar - A3*X3Bar*X3Bar );
  PhiVD[1] = -2.0*A1*X1Bar*PhiVD[0];
  PhiVD[2] = -2.0*A2*X2Bar*PhiVD[0];
  PhiVD[3] = -2.0*A3*X3Bar*PhiVD[0];
  PhiVD[4] = 4.0*A1*A2*X1Bar*X2Bar*PhiVD[0];
  PhiVD[5] = 4.0*A1*A3*X1Bar*X3Bar*PhiVD[0];
  PhiVD[6] = 4.0*A2*A3*X2Bar*X3Bar*PhiVD[0];
  PhiVD[7] = -8.0*A1*A2*A3*X1Bar*X2Bar*X3Bar*PhiVD[0];

  /* second function and its first derivatives */
  PhiVD[8]  = exp( -B1*X1Bar*X1Bar - B2*X2Bar*X2Bar - B3*X3Bar*X3Bar );
  PhiVD[9]  = -2.0*B1*X1Bar*PhiVD[8];
  PhiVD[10] = -2.0*B2*X2Bar*PhiVD[8];
  PhiVD[11] = -2.0*B3*X3Bar*PhiVD[8];
  PhiVD[12] = 4.0*B1*B2*X1Bar*X2Bar*PhiVD[8];
  PhiVD[13] = 4.0*B1*B3*X1Bar*X3Bar*PhiVD[8];
  PhiVD[14] = 4.0*B2*B3*X2Bar*X3Bar*PhiVD[8];
  PhiVD[15] = -8.0*B1*B2*B3*X1Bar*X2Bar*X3Bar*PhiVD[8];

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  int nThread;
  char *X1File, *X2File, *X3File, *DataFile;

  double *X1Points, *X2Points, *X3Points;
  int N1, N2, N3;
  Interp3D *I3D;

  /*--------------------------------------------------------------*/
  /*- process command-line arguments -----------------------------*/
  /*--------------------------------------------------------------*/
  ArgStruct ASArray[]=
   { {"X1Points",   PA_STRING, (void *)&X1File,       0,  "list of X1Points"},
     {"X2Points",   PA_STRING, (void *)&X2File,       0,  "list of X2Points"},
     {"X3Points",   PA_STRING, (void *)&X3File,       0,  "list of X3Points"},
     {"DataFile",   PA_STRING, (void *)&DataFile,     0,  "binary data file"},
     {"nThread",    PA_INT,    (void *)&nThread,    "-1", "number of threads"},
     {0,0,0,0,0}
   };
  ProcessArguments(argc, argv, ASArray);
  if (nThread==-1)
   nThread=GetNumThreads();

  /*--------------------------------------------------------------*/
  /*- if the user specified a binary data file, attempt to read --*/
  /*- in the Interp3D object from the data file                 --*/
  /*--------------------------------------------------------------*/
  if (DataFile)
   { 
     I3D=new Interp3D(DataFile);

     X1Points=I3D->X1Points; N1=I3D->N1;
     X2Points=I3D->X2Points; N2=I3D->N2;
     X3Points=I3D->X3Points; N3=I3D->N3;
   }
  else
    {
      /*--------------------------------------------------------------*/
      /*- otherwise initialize the Interp3D object from scratch and   */
      /*- export it to the binary data file                           */
      /*--------------------------------------------------------------*/
      if( !X1File )
       ErrExit("--X1Points option is mandatory");
      HVector *X1Vec=new HVector(X1File,LHM_TEXT);
      if (X1Vec->ErrMsg)
       ErrExit("%s",X1Vec->ErrMsg);
      X1Points=X1Vec->DV; N1=X1Vec->N;
 
      if( !X2File )
       ErrExit("--X2Points option is mandatory");
      HVector *X2Vec=new HVector(X2File,LHM_TEXT);
      if (X2Vec->ErrMsg)
       ErrExit("%s",X2Vec->ErrMsg);
      X2Points=X2Vec->DV; N2=X2Vec->N;

      if( !X3File )
       ErrExit("--X3Points option is mandatory");
      HVector *X3Vec=new HVector(X3File,LHM_TEXT);
      if (X3Vec->ErrMsg)
       ErrExit("%s",X3Vec->ErrMsg);
      X3Points=X3Vec->DV; N3=X3Vec->N;

      I3D=new Interp3D(X1Points, N1, X2Points, N2, X3Points, N3, 
                       NFUN, nThread, Phi, 0);

      I3D->WriteToFile("tInterp3D.dat");
    };

  /*--------------------------------------------------------------*/ 
  /*- test the interpolating function at the center of each grid -*/ 
  /*- cell                                                       -*/ 
  /*--------------------------------------------------------------*/ 
  double RD1, MaxRD1=-1.0e9, MaxRDPE1, MaxRDPI1, MaxRDC1[3];
  double RD2, MaxRD2=-1.0e9, MaxRDPE2, MaxRDPI2, MaxRDC2[3];
  double PhiExact[16], PhiInterp[2];
  double X1, X2, X3;
  int n1, n2, n3;

  FILE *f=CreateUniqueFile("tInterp3D.out",1);
  fprintf(f,"# columns: \n");
  fprintf(f,"# 1 Delta\n");
  fprintf(f,"# 2 Lambda\n");
  fprintf(f,"# 3 Gamma\n");
  fprintf(f,"# 4 W(exact)\n");
  fprintf(f,"# 5 W(interpolated)\n");
  fprintf(f,"# 6 |W(e) - W(i)| / |W(e)|\n");

  for(n1=0; n1<N1-1; n1++)
   for(n2=0; n2<N2-1; n2++)
    for(n3=0; n3<N3-1; n3++)
     { 
       X1=X1Points[n1] + 0.5*(X1Points[n1+1]-X1Points[n1]);
       X2=X2Points[n2] + 0.5*(X2Points[n2+1]-X2Points[n2]);
       X3=X3Points[n3] + 0.5*(X3Points[n3+1]-X3Points[n3]);

       Phi(X1, X2, X3, 0, PhiExact);
       I3D->Evaluate(X1, X2, X3, PhiInterp);

       RD1=fabs(PhiExact[0]-PhiInterp[0])/fabs(PhiExact[0]);
       if (RD1>MaxRD1)
        { MaxRD1=RD1;
          MaxRDPE1=PhiExact[0];
          MaxRDPI1=PhiInterp[0];
          MaxRDC1[0]=X1;
          MaxRDC1[1]=X2;
          MaxRDC1[2]=X3;
        };

       RD2=fabs(PhiExact[8]-PhiInterp[1])/fabs(PhiExact[8]);
       if (RD2>MaxRD2)
        { MaxRD2=RD2;
          MaxRDPE2=PhiExact[8];
          MaxRDPI2=PhiInterp[1];
          MaxRDC2[0]=X1;
          MaxRDC2[1]=X2;
          MaxRDC2[2]=X3;
        };

       fprintf(f,"%e %e %e %.12e %.12e %e %.12e %.12e %e\n",
                  X1,X2,X3,
                  PhiExact[0],PhiInterp[0],RD1,
                  PhiExact[8],PhiInterp[1],RD2);

       if ( n3==N3-2 )
        fprintf(f,"\n\n");
     };
  fclose(f);

  /*--------------------------------------------------------------*/ 
  /*--------------------------------------------------------------*/ 
  /*--------------------------------------------------------------*/ 
  printf("Function 1: \n");
  printf("Maximum rel. diff. is %5.2e \n",MaxRD1);
  printf("at (X1,X2,X3)=(%+10.3e,%+10.3e,%+10.3e)\n",MaxRDC1[0],MaxRDC1[1],MaxRDC1[2]);
  printf(" Exact  = %+17.10e\n",MaxRDPE1);
  printf(" Interp = %+17.10e\n",MaxRDPI1);
  printf("\n");

  printf("Function 2: \n");
  printf("Maximum rel. diff. is %5.2e \n",MaxRD2);
  printf("at (X1,X2,X3)=(%+10.3e,%+10.3e,%+10.3e)\n",MaxRDC2[0],MaxRDC2[1],MaxRDC2[2]);
  printf(" Exact  = %+17.10e\n",MaxRDPE2);
  printf(" Interp = %+17.10e\n",MaxRDPI2);
  printf("\n");

  printf("\n");
  printf("Thank you for your support.\n");

}
