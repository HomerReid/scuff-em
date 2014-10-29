/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * tOverlap.cc -- test of libscuff's routines for computing
 *                overlap integrals between RWG basis functions 
 * 
 * homer reid  -- 3/2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <config.h>

#ifdef HAVE_LIBREADLINE
 #include <readline/readline.h>
 #include <readline/history.h>
#else
 #include "readlineReplacement.h"
#endif

#include "libhrutil.h"
#include "libSGJC.h"
#include "libscuff.h"
#include "libscuffInternals.h"

#if defined(_WIN32)
#  define srand48 srand
#  define drand48 my_drand48
static double my_drand48(void) {
  return rand() * 1.0 / RAND_MAX;
}
static long int lrand48() { return rand(); }
#endif

#define MAXSTR 1000
#define MAXEVAL 10000
#define ABSTOL 1.0e-8
#define RELTOL 1.0e-4

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define NUMOVERLAPS 20

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct GOBFData 
 {
   double *V0, *V1, *V2, *QBeta, *ZHat;
   double PreFac;
 } GOBFData;

void GetOverlapBFIntegrand(unsigned ndim, const double *x, void *params,
			   unsigned fdim, double *fval)
{
  (void) ndim;
  (void) fdim;

  GOBFData *GOBFD = (GOBFData *)params;
  double u=x[0], v=u*x[1], Jacobian=u;

  double *V0    = GOBFD->V0;
  double *V1    = GOBFD->V1;
  double *V2    = GOBFD->V2;
  double *QBeta = GOBFD->QBeta;
  double *ZHat  = GOBFD->ZHat;
  double PreFac = GOBFD->PreFac;

  double X[3], bAlpha[3], bBeta[3]; 

  int Mu;
  for(Mu=0; Mu<3; Mu++)
   { X[Mu]      = V0[Mu] + u*(V1[Mu]-V0[Mu]) + v*(V2[Mu]-V1[Mu]);
     bAlpha[Mu] = X[Mu] - V0[Mu];
     bBeta[Mu]  = X[Mu] - QBeta[Mu];
   };

  double nxbAlpha[3];
  nxbAlpha[0] = ZHat[1]*bAlpha[2] - ZHat[2]*bAlpha[1];
  nxbAlpha[1] = ZHat[2]*bAlpha[0] - ZHat[0]*bAlpha[2];
  nxbAlpha[2] = ZHat[0]*bAlpha[1] - ZHat[1]*bAlpha[0];

  double nxbBeta[3];
  nxbBeta[0] = ZHat[1]*bBeta[2] - ZHat[2]*bBeta[1];
  nxbBeta[1] = ZHat[2]*bBeta[0] - ZHat[0]*bBeta[2];
  nxbBeta[2] = ZHat[0]*bBeta[1] - ZHat[1]*bBeta[0];

  double DotProd = (bAlpha[0]*bBeta[0] + bAlpha[1]*bBeta[1] + bAlpha[2]*bBeta[2]);

  double RxN[3];
  RxN[0] = X[1]*ZHat[2] - X[2]*ZHat[1];
  RxN[1] = X[2]*ZHat[0] - X[0]*ZHat[2];
  RxN[2] = X[0]*ZHat[1] - X[1]*ZHat[0];

  double RxNxbAlpha[3];
  RxNxbAlpha[0] = X[1]*nxbAlpha[2] - X[2]*nxbAlpha[1];
  RxNxbAlpha[1] = X[2]*nxbAlpha[0] - X[0]*nxbAlpha[2];
  RxNxbAlpha[2] = X[0]*nxbAlpha[1] - X[1]*nxbAlpha[0];

  fval[0] = Jacobian*PreFac*DotProd;
  fval[1] = Jacobian*PreFac*(bAlpha[0]*nxbBeta[0] + bAlpha[1]*nxbBeta[1] + bAlpha[2]*nxbBeta[2]);

  fval[2] =     Jacobian*PreFac*ZHat[0]*DotProd;
  fval[3] =     Jacobian*PreFac*ZHat[0]*4.0;
  fval[4] = Jacobian*PreFac*nxbAlpha[0]*2.0;

  fval[5] =     Jacobian*PreFac*ZHat[1]*DotProd;
  fval[6] =     Jacobian*PreFac*ZHat[1]*4.0;
  fval[7] = Jacobian*PreFac*nxbAlpha[1]*2.0;

  fval[8] =     Jacobian*PreFac*ZHat[2]*DotProd;
  fval[9] =     Jacobian*PreFac*ZHat[2]*4.0;
  fval[10]= Jacobian*PreFac*nxbAlpha[2]*2.0;

  fval[11] =        Jacobian*PreFac*RxN[0]*DotProd;
  fval[12] =        Jacobian*PreFac*RxN[0]*4.0;
  fval[13] = Jacobian*PreFac*RxNxbAlpha[0]*2.0;

  fval[14] =        Jacobian*PreFac*RxN[1]*DotProd;
  fval[15] =        Jacobian*PreFac*RxN[1]*4.0;
  fval[16] = Jacobian*PreFac*RxNxbAlpha[1]*2.0;

  fval[17] =        Jacobian*PreFac*RxN[2]*DotProd;
  fval[18] =        Jacobian*PreFac*RxN[2]*4.0;
  fval[19] = Jacobian*PreFac*RxNxbAlpha[2]*2.0;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetOverlapBF(RWGSurface *S, int nea, int neb, double OValues[NUMOVERLAPS])
{
  RWGEdge *Ea=S->Edges[nea], *Eb=S->Edges[neb];

  double *QP = S->Vertices + 3*Ea->iQP;
  double *V1 = S->Vertices + 3*Ea->iV1;
  double *V2 = S->Vertices + 3*Ea->iV2;
  double *QM = S->Vertices + 3*Ea->iQM;

  double *QPP = S->Vertices + 3*Eb->iQP;
  double *QMP = S->Vertices + 3*Eb->iQM;

  GOBFData MyGOBFD, *GOBFD = &MyGOBFD;
  GOBFD->V1   = V1;
  GOBFD->V2   = V2;

  RWGPanel *PP= S->Panels[Ea->iPPanel];
  double PArea = PP->Area;

  RWGPanel *PM= S->Panels[Ea->iMPanel];
  double MArea = PM->Area;

  double Lower[2]={0.0, 0.0};
  double Upper[2]={1.0, 1.0};
  double I[NUMOVERLAPS], E[NUMOVERLAPS];

  int n;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  memset(OValues,0,NUMOVERLAPS*sizeof(double));
  if ( Ea->iPPanel == Eb->iPPanel )
   { 
     GOBFD->V0     = QP; 
     GOBFD->QBeta  = QPP;
     GOBFD->ZHat   = PP->ZHat;
     GOBFD->PreFac = (Ea->Length * Eb->Length)/(2.0*PArea);

     adapt_integrate(NUMOVERLAPS, GetOverlapBFIntegrand, (void *)GOBFD, 2, 
		     Lower, Upper, MAXEVAL, ABSTOL, RELTOL, I, E);
     for(n=0; n<NUMOVERLAPS; n++)
      OValues[n] += I[n];
   };

  if ( Ea->iPPanel == Eb->iMPanel )
   { 
     GOBFD->V0     = QP; 
     GOBFD->QBeta  = QMP;
     GOBFD->ZHat   = PP->ZHat;
     GOBFD->PreFac = (Ea->Length * Eb->Length)/(2.0*PArea);

     adapt_integrate(NUMOVERLAPS, GetOverlapBFIntegrand, (void *)GOBFD, 2, 
		     Lower, Upper, MAXEVAL, ABSTOL, RELTOL, I, E);
     for(n=0; n<NUMOVERLAPS; n++)
      OValues[n] -= I[n];

   };

  if ( Ea->iMPanel == Eb->iPPanel )
   { 
     GOBFD->V0     = QM; 
     GOBFD->QBeta  = QPP;
     GOBFD->ZHat   = PM->ZHat;
     GOBFD->PreFac = (Ea->Length * Eb->Length)/(2.0*MArea);

     adapt_integrate(NUMOVERLAPS, GetOverlapBFIntegrand, (void *)GOBFD, 2, 
		     Lower, Upper, MAXEVAL, ABSTOL, RELTOL, I, E);
     for(n=0; n<NUMOVERLAPS; n++)
      OValues[n] -= I[n];

   };

  if ( Ea->iMPanel == Eb->iMPanel )
   { 
     GOBFD->V0     = QM; 
     GOBFD->QBeta  = QMP;
     GOBFD->ZHat   = PM->ZHat;
     GOBFD->PreFac = (Ea->Length * Eb->Length)/(2.0*MArea);

     adapt_integrate(NUMOVERLAPS, GetOverlapBFIntegrand, (void *)GOBFD, 2, 
		     Lower, Upper, MAXEVAL, ABSTOL, RELTOL, I, E);
     for(n=0; n<NUMOVERLAPS; n++)
      OValues[n] += I[n];
   };

}

/***************************************************************/
/* print a console prompt, then get and parse a line of        */
/* input data                                                  */
/***************************************************************/
typedef struct Request
 {
   int ne1, ne2;
 } Request;

void GetRequest(RWGGeometry *G, Request *R)
{
  static int init=0;

  if (init==0)
   { init=1;
     using_history();
     read_history(0);
     srand48(time(0));
   };

  /*--------------------------------------------------------------*/
  /*- print prompt and get input string --------------------------*/
  /*--------------------------------------------------------------*/
  printf("\n");
  printf(" options: --ne1 xx \n");
  printf("          --ne2 xx \n");
  char *p;
  do
   { 
     p=readline("enter options: "); 
   } while(!p);
  add_history(p);
  write_history(0);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int nt, NumTokens;
  int ncv;
  char *Tokens[50];

  R->ne1   =  -1;
  R->ne2   =  -1;
  ncv      = -1;

  NumTokens=Tokenize(p,Tokens,50);
  for(nt=0; nt<NumTokens; nt++)
   if ( !strcasecmp(Tokens[nt],"--ne1") )
    sscanf(Tokens[nt+1],"%i",&(R->ne1));
  for(nt=0; nt<NumTokens; nt++)
   if ( !strcasecmp(Tokens[nt],"--ne2") )
    sscanf(Tokens[nt+1],"%i",&(R->ne2));
  for(nt=0; nt<NumTokens; nt++)
   if ( !strcasecmp(Tokens[nt],"--ncv") )
    sscanf(Tokens[nt+1],"%i",&ncv);
  for(nt=0; nt<NumTokens; nt++)
   if ( !strcasecmp(Tokens[nt],"--quit") )
    exit(1);

  free(p);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGSurface *S = G->Surfaces[0];

  if ( R->ne1==-1 && R->ne2==-1 && ncv==-1) 
   { R->ne1=lrand48() % S->NumEdges;
     ncv = (drand48() > 0.5) ? 3 : 4;
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (ncv!=-1)
   { for(R->ne2=0; R->ne2 < (S->NumEdges-1); R->ne2++)
      if ( ncv == NumCommonBFVertices(S, R->ne1, S, R->ne2) )
       break;     
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /***************************************************************/
  /* process command-line arguments                              */
  /***************************************************************/
  char *GeoFileName;
  cdouble Omega=1.0;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { 
     {"geometry", PA_STRING,  1, 1, (void *)&GeoFileName,  0,  ".scuffgeo file"},
     {"Omega",    PA_CDOUBLE, 3, 3, (void *)&Omega,        0,  "angular frequency"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFileName==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  RWGGeometry *G=new RWGGeometry(GeoFileName);
  RWGSurface *S=G->Surfaces[0];

  /***************************************************************/
  /* enter command loop ******************************************/
  /***************************************************************/
  for(;;)
   { 
     Request MyRequest, *R=&MyRequest;
     GetRequest(G,R);
     int ne1 = R->ne1;
     int ne2 = R->ne2;

     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     double OverlapBF[NUMOVERLAPS];
     GetOverlapBF(S, ne1, ne2, OverlapBF);

     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     double OverlapHR[NUMOVERLAPS];
     S->GetOverlaps(ne1, ne2, OverlapHR);

     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     printf("\n*\n* --ne1 %i --ne2 %i (--ncv %i) \n",ne1,ne2,NumCommonBFVertices(S, ne1, S, ne2)); 
     for(int no=0; no<NUMOVERLAPS; no++)
      printf("%2i | %+12.5e | %+12.5e | %.1e %.5e\n",no, OverlapHR[no], OverlapBF[no],
              RD(OverlapHR[no],OverlapBF[no]), OverlapHR[no]/OverlapBF[no] );

   };

}
