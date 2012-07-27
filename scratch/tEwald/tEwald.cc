/*
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#include <readline/readline.h>
#include <readline/history.h>

#include <libhrutil.h>
#include <libscuff.h>

#define ABSTOL 1.0e-4
#define RELTOL 1.0e-2

#define II cdouble(0,1)

#define NSUM 8

int ExcludeFirst9;

/***************************************************************/
/* brute-force evaluation of sum *******************************/
/***************************************************************/
void GBarVDBF(double *R, cdouble k, double *P, double **LBV, int ExcludeFirst9,
              double AbsTol, double RelTol, int *nCells, cdouble *Sum);

/***************************************************************/
/* ewald evaluation of sum *************************************/
/***************************************************************/
namespace scuff{

void ComputeG1(double *R, cdouble k, double *P, double **LBV,
               double E, int *pnCells, cdouble *Sum);

void ComputeG2(double *R, cdouble k, double *P, double **LBV,
               double E, int *pnCells, cdouble *Sum);

void AddGBFContribution(double R[3], cdouble k, double P[2],
                        double Lx, double Ly, cdouble *Sum);

void GBarVDEwald(double *R, cdouble k, double *P, double *LBV[2],
                 double E, int ExcludeFirst9, cdouble *GBarVD);

void ComputeGBFFirst9(double *R, cdouble k, double *P, double *LBV[2], cdouble *Sum);
}
using namespace scuff;


/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  double LBV1[2]={1.0, 0.0};
  double LBV2[2]={0.0, 1.0};
  double *LBV[2]={LBV1, LBV2};

  double E;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  char *p;
  char *Tokens[100];
  int nt, NumTokens;
  double kR, kI, P[2], R[3];
  double AbsTol, RelTol;
  int nCells1, nCells2;
  cdouble k;
  cdouble GBF[NSUM], GEwald[NSUM]; 
  cdouble G1[NSUM], G2[NSUM], GBFFirst9[NSUM];
  using_history();
  read_history(0);
  srand48(time(0));
  while(1)
   {
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     kR=2.0*drand48();
     kI=2.0*drand48();
     P[0]=M_PI*drand48();
     P[1]=M_PI*drand48();
     R[0]=0.5 - drand48();
     R[1]=0.5 - drand48();
     R[2]=0.5 - drand48();
     E=sqrt( M_PI / (LBV[0][0]*LBV[1][1]) );
     AbsTol=1.0e-8;
     RelTol=1.0e-2;

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     printf("options: \n");
     printf(" --kR xx --kI xx \n");
     printf(" --Px xx --Py xx \n");
     printf(" --x xx --y xx --z xx ");
     printf(" --E xx \n");
     printf(" --abstol xx --reltol xx \n");
     printf(" --RetainFirst9 \n");
     p=readline("Enter options: ");
     add_history(p); 
     write_history(0); 
     NumTokens=Tokenize(p, Tokens, 100);
     ExcludeFirst9=1;
     for(nt=0; nt<NumTokens; nt++)
      {  
        if( !strcasecmp(Tokens[nt],"--kR") )
         sscanf(Tokens[++nt],"%le",&kR);
        if( !strcasecmp(Tokens[nt],"--kI") )
         sscanf(Tokens[++nt],"%le",&kI);
        else if( !strcasecmp(Tokens[nt],"--Px") )
         sscanf(Tokens[++nt],"%le",P+0);
        else if( !strcasecmp(Tokens[nt],"--Py") )
         sscanf(Tokens[++nt],"%le",P+1);
        else if( !strcasecmp(Tokens[nt],"--x") )
         sscanf(Tokens[++nt],"%le",R);
        else if( !strcasecmp(Tokens[nt],"--y") )
         sscanf(Tokens[++nt],"%le",R+1);
        else if( !strcasecmp(Tokens[nt],"--z") )
         sscanf(Tokens[++nt],"%le",R+2);
        else if( !strcasecmp(Tokens[nt],"--E") )
         sscanf(Tokens[++nt],"%le",&E);
        else if( !strcasecmp(Tokens[nt],"--AbsTol") )
         sscanf(Tokens[++nt],"%le",&AbsTol);
        else if( !strcasecmp(Tokens[nt],"--RelTol") )
         sscanf(Tokens[++nt],"%le",&RelTol);
        else if( !strcasecmp(Tokens[nt],"--RetainFirst9") )
         ExcludeFirst9=0;
      };
     k=cdouble(kR, kI);

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     printf("\n ** Parameters: \n");
     printf(" --kR %g --kI %g  \n",kR, kI);
     printf(" --Px %g --Py %g \n",P[0],P[1]);
     printf(" --x %g --y %g --z %g \n",R[0],R[1],R[2]);
     printf(" --E %g\n",E);
     printf(" --AbsTol %g --RelTol %g\n",AbsTol,RelTol);
     printf(" --ExcludeFirst9: %i\n",ExcludeFirst9);
     printf("\n");

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     GBarVDBF(R, k, P, LBV, ExcludeFirst9, AbsTol, RelTol, &nCells1, GBF);
     printf("**Direct method: %i cells summed.\n",nCells1);
#if 0
     ComputeG1(R, k, P, LBV, E, &nCells1, G1);
     ComputeG2(R, k, P, LBV, E, &nCells2, G2);

if (ExcludeFirst9) 
 ComputeGBFFirst9(R, k, P, LBV, GBFFirst9);
else
 memset(GBFFirst9,0,NSUM*sizeof(cdouble));
 
     for(int ns=0; ns<NSUM; ns++)
      GEwald[ns]=G1[ns] + G2[ns] - GBFFirst9[ns];

     SetDefaultCD2SFormat("(%+20.12e, %+20.12e)");
     printf("**Ewald method: %i momentum-space / %i real-space cells summed.\n",
             nCells1, nCells2);
     printf("\n");
  
     printf("G (sum):   %s \n",CD2S(GBF[0]));
     printf("G (Ewald): %s (%e) \n",CD2S(GEwald[0]),RD(GBF[0],GEwald[0]));
     printf("\nContributions to Ewald sum: \n");
     printf(" Momentum space: %s\n",CD2S(G1[0]));
     printf(" Real space    : %s\n",CD2S(G2[0]));
     printf(" \n\n");
#endif

     GBarVDEwald(R, k, P, LBV, -1.0, ExcludeFirst9, GEwald);

     SetDefaultCD2SFormat("(%+20.12e, %+20.12e)");
     printf("**Ewald method: %i momentum-space / %i real-space cells summed.\n",
             nCells1, nCells2);
     printf("\n");
 
     SetDefaultCD2SFormat("(%+12.4e, %+12.4e)");
     for(int n=0; n<8; n++)
      printf("%i: | %s | %s | %.1e \n",n,CD2S(GEwald[n]),CD2S(GBF[n]),RD(GEwald[n],GBF[n]));
    
   };

}
