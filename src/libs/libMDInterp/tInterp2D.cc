/*
 * tInterp2D.cc -- test program for Interp2D class 
 *
 * homer reid      -- 3/14/2011
 */
#include <stdio.h>
#include <stdlib.h>

#include <libhrutil.h>
#include "libMDInterp.h"

/***************************************************************/
/* the function to be interpolated   ***************************/
/***************************************************************/
#define A1 0.5
#define A2 2.0
void Phi(double X1, double X2, void *data, double *PhiVD)
{

  PhiVD[0] = exp( -A1*X1*X1 -A2*X2*X2 );
  PhiVD[1] = -2.0*A1*X1*PhiVD[0];                 /* dPhi / dX1 */
  PhiVD[2] = -2.0*A2*X2*PhiVD[0];                 /* dPhi / dX2 */
  PhiVD[3] = +4.0*A1*A2*X1*X2*PhiVD[0];           /* d^2Phi / dX1dX2 */

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  Interp2D *I2D;
  int nThread;
  char *DeltaList, *LambdaList;

  /*--------------------------------------------------------------*/
  /*- process command-line arguments -----------------------------*/
  /*--------------------------------------------------------------*/
  ArgStruct ASArray[]=
   { {"DeltaList",  PA_STRING, (void *)&DeltaList,     0, "list of Delta points"},
     {"LambdaList", PA_STRING, (void *)&LambdaList,    0, "list of Lambda points"},
     {"nThread",    PA_INT,    (void *)&nThread,    "-1", "number of threads"},
     {0,0,0,0,0}
   };
  ProcessArguments(argc, argv, ASArray);
  if (DeltaList==0)
   ASUsage(argv[0],ASArray,"--DeltaList option is mandatory");
  if (LambdaList==0)
   ASUsage(argv[0],ASArray,"--LambdaList option is mandatory");
  if (nThread==-1)
   nThread=GetNumThreads();

  /*--------------------------------------------------------------*/
  /*- read lists of points in from user-specified files ----------*/
  /*--------------------------------------------------------------*/
  FILE *f;
  int n1, n2;
  int N1, N2;
  double *X1Points, *X2Points;
  char buffer[100];

  if ( !(f=fopen(DeltaList,"r")) )
   ErrExit("could not open file %s",DeltaList);
  N1=0;
  while(fgets(buffer,100,f))
   N1++;
  rewind(f);
  X1Points=(double *)mallocEC(N1*sizeof(double));
  n1=0;
  while(fgets(buffer,100,f))
   sscanf(buffer,"%le",X1Points + (n1++) );
  fclose(f);
 
  if ( !(f=fopen(LambdaList,"r")) )
   ErrExit("could not open file %s",LambdaList);
  N2=0;
  while(fgets(buffer,100,f))
   N2++;
  rewind(f); 
  X2Points=(double *)mallocEC(N2*sizeof(double));
  n2=0;
  while(fgets(buffer,100,f))
   sscanf(buffer,"%le",X2Points + (n2++) );
  fclose(f);

  /*--------------------------------------------------------------*/ 
  /*- create the interpolation table -----------------------------*/ 
  /*--------------------------------------------------------------*/ 
  I2D=new Interp2D(X1Points, N1, X2Points, N2, 1, nThread, Phi, NULL);

  /*--------------------------------------------------------------*/ 
  /*- test the interpolating function at the center of each grid -*/ 
  /*- cell                                                       -*/ 
  /*--------------------------------------------------------------*/ 
  double RD, MaxRD=-1.0e9, MaxRDPE, MaxRDPI, MaxRDC[2];
  double PhiExact[4], PhiInterp;
  double X1, X2;

  f=CreateUniqueFile("tInterp2D.out",1);
  fprintf(f,"# columns: \n");
  fprintf(f,"# 1 Delta\n");
  fprintf(f,"# 2 Lambda\n");
  fprintf(f,"# 3 W(exact)\n");
  fprintf(f,"# 4 W(interpolated)\n");
  fprintf(f,"# 5 |W(e) - W(i)| / |W(e)|\n");

  for(n1=0; n1<N1-1; n1++)
   for(n2=0; n2<N2-1; n2++)
    { 
      X1=X1Points[n1] + 0.5*(X1Points[n1+1]-X1Points[n1]);
      X2=X2Points[n2] + 0.5*(X2Points[n2+1]-X2Points[n2]);

      Phi(X1, X2, NULL, PhiExact);
      I2D->Evaluate(X1, X2, &PhiInterp);

      RD=fabs(PhiExact[0]-PhiInterp)/fabs(PhiExact[0]);
      if (RD>MaxRD)
       { MaxRD=RD;
         MaxRDPE=PhiExact[0];
         MaxRDPI=PhiInterp;
         MaxRDC[0]=X1;
         MaxRDC[1]=X2;
       };

      fprintf(f,"%e %e %.12e %.12e %e \n",X1,X2,PhiExact[0],PhiInterp,RD);

      if ( n2==N2-2 )
       fprintf(f,"\n\n");
    };
 fclose(f);

  /*--------------------------------------------------------------*/ 
  /*--------------------------------------------------------------*/ 
  /*--------------------------------------------------------------*/ 
  printf("Maximum rel. diff. is %5.2e \n",MaxRD);
  printf("at (D,L)=(%+10.3e,%+10.3e)\n",MaxRDC[0],MaxRDC[1]);
  printf(" Exact  = %+17.10e\n",MaxRDPE);
  printf(" Interp = %+17.10e\n",MaxRDPI);

  printf("\n");
  printf("Thank you for your support.\n");

}
