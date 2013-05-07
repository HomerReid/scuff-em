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
 * scuff-test-EEIs.cc -- a test program for libscuff's routines for
 *                    -- computing edge-edge interactions
 * 
 * homer reid         -- 11/2005 -- 11/2011
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <readline/readline.h>
#include <readline/history.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <libhrutil.h>
#include <PanelCubature.h>

#include "scuff-test-EEIs.h"

/***************************************************************/
/* integrand for brute-force edge-edge integration            **/
/***************************************************************/
void Integrand(double *X, double *XP, PPCData *PPCD,
               void *UserData, double *Result)
{ 
  cdouble *F1    = PPCD->K1;    // RWG basis function 1 
  cdouble DivF1  = PPCD->DivF1;
  cdouble *F2    = PPCD->K2;    // RWG basis function 2
  cdouble DivF2  = PPCD->DivF2;

  double R[3];
  VecSub(X, XP, R);
  double r2=(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
  double r=sqrt(r2);

  double DotProduct    = real(F1[0]*F2[0] + F1[1]*F2[1] + F1[2]*F2[2]);
  double ScalarProduct = real(DivF1 * DivF2);
  double TripleProduct =  F1[0] * (R[1]*F2[2]-R[2]*F2[1])
                         +F1[1] * (R[2]*F2[0]-R[0]*F2[2])
                         +F1[2] * (R[0]*F2[1]-R[1]*F2[0]);

  cdouble ik  = II * PPCD->Omega;
  cdouble ikr = ik*r;
  cdouble Phi = exp(ikr)/(4.0*M_PI*r);
  cdouble Psi = (ikr-1.0)*Phi/r2;

  cdouble *zf = (cdouble *)Result;
  zf[0] = ( DotProduct + ScalarProduct/(ik*ik) ) * Phi;
  zf[1] = TripleProduct * Psi;

}

/***************************************************************/
/* compute edge-edge integrals by brute force                  */
/***************************************************************/
void GetEEIs_BruteForce(RWGGeometry *G, GetEEIArgStruct *Args)
{
  GetBFBFCubature(G, Args->Sa->Index, Args->nea,
                  Args->Sb->Index, Args->neb,
                  Args->Displacement,
                  Integrand, 0, 4, 0, 1.0e-10, 0.0,
                  Args->k, 0, (cdouble *)Args->GC);
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetEEIs_CartesianMultipole(GetEEIArgStruct *Args)
{
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
#if 0
int AssessEdgePair(RWGObject *Oa, int nea, RWGObject *Ob, int neb,
                   double *rRel)
{
  /***************************************************************/
  /* get relative distance ***************************************/
  /***************************************************************/
  RWGEdge *Ea=Oa->Edges[nea];
  RWGEdge *Eb=Ob->Edges[neb];
  *rRel=VecDistance(Ea->Centroid, Eb->Centroid) / fmax(Ea->Radius, Eb->Radius);

  /***************************************************************/
  /* count common vertices ***************************************/
  /***************************************************************/
  double *Va[4], *Vb[4];

  Va[0] = Oa->Vertices + 3*(Ea->iQP);
  Va[1] = Oa->Vertices + 3*(Ea->iV1);
  Va[2] = Oa->Vertices + 3*(Ea->iV2);
  Va[3] = Oa->Vertices + 3*(Ea->iQM);

  Vb[0] = Ob->Vertices + 3*(Eb->iQP);
  Vb[1] = Ob->Vertices + 3*(Eb->iV1);
  Vb[2] = Ob->Vertices + 3*(Eb->iV2);
  Vb[3] = Ob->Vertices + 3*(Eb->iQM);

  int i, j, ncv;
  ncv=0;
  for(i=0; i<4; i++)
   for(j=0; j<4; j++)
    if (Va[i]==Vb[j]) ncv++;

  return ncv;
 
} 
#endif

/***************************************************************/
/* print a console prompt, then get an parse a command         */
/***************************************************************/
typedef struct Request
 {
   int nsa, nea, nsb, neb;
   cdouble k;
   bool Gradient;
   double DX[3];
 } Request;

void GetRequest(RWGGeometry *G, Request *R)
{
  static int init=0;

  if (init==0)
   { init=1;
     using_history();
     read_history(0);
   };

  /*--------------------------------------------------------------*/
  /*- print prompt and get input string --------------------------*/
  /*--------------------------------------------------------------*/
  printf("\n");
  printf(" options: --nsa  xx \n");
  printf("          --nea  xx \n");
  printf("          --nsb  xx \n");
  printf("          --neb  xx \n");
  printf("          --rRel xx \n");
  printf("          --DX   xx xx xx \n");
  printf("          --k       \n");
  printf("          --gradient\n");
  printf("          --quit\n");
  char *p;
  do
   { 
     p=readline("enter options: "); 
   } while(!p);
  add_history(p);
  write_history(0);

  /*--------------------------------------------------------------*/
  /*- parse input string to construct a Request structure        -*/
  /*--------------------------------------------------------------*/
  int nt, NumTokens;
  char *Tokens[50];

  R->nsa=R->nsb=0;
  R->nea=R->neb=-1;
  double rRel=-1.0;
  R->DX[0]=R->DX[1]=R->DX[2]=0.0;
  R->Gradient=true;
  R->k=0.1;

  NumTokens=Tokenize(p,Tokens,50);
  for(nt=0; nt<NumTokens; nt++)
   if ( !strcasecmp(Tokens[nt],"--nsa") )
    sscanf(Tokens[nt+1],"%i",&(R->nsa));
  for(nt=0; nt<NumTokens; nt++)
   if ( !strcasecmp(Tokens[nt],"--nea") )
    sscanf(Tokens[nt+1],"%i",&(R->nea));
  for(nt=0; nt<NumTokens; nt++)
   if ( !strcasecmp(Tokens[nt],"--nsb") )
    sscanf(Tokens[nt+1],"%i",&(R->nsb));
  for(nt=0; nt<NumTokens; nt++)
   if ( !strcasecmp(Tokens[nt],"--neb") )
    sscanf(Tokens[nt+1],"%i",&(R->neb));
  for(nt=0; nt<NumTokens; nt++)
   if ( !strcasecmp(Tokens[nt],"--rRel") )
    sscanf(Tokens[nt+1],"%le",&rRel);
  for(nt=0; nt<NumTokens; nt++)
   if ( !strcasecmp(Tokens[nt],"--k") )
    S2CD(Tokens[nt+1],&(R->k));
  for(nt=0; nt<NumTokens; nt++)
   if ( !strcasecmp(Tokens[nt],"--dx") )
    { sscanf(Tokens[nt+1],"%le",(R->DX)+0); 
      sscanf(Tokens[nt+1],"%le",(R->DX)+1); 
      sscanf(Tokens[nt+1],"%le",(R->DX)+2);
    };
  for(nt=0; nt<NumTokens; nt++)
   if ( !strcasecmp(Tokens[nt],"--quit") )
    exit(1);

  free(p);

  /*--------------------------------------------------------------*/
  /*- if the user didn't specify anything about the panels then --*/
  /*- choose a random pair of panels                            --*/
  /*--------------------------------------------------------------*/
  RWGSurface *Sa = G->Surfaces[R->nsa], *Sb = G->Surfaces[R->nsb];
  if (R->nea==-1) R->nea=lrand48() % Sa->NumEdges;
  if (R->neb==-1) R->neb=lrand48() % Sb->NumEdges;

  /*--------------------------------------------------------------*/
  /* if the user specified a relative distance then try to replace*/
  /* neb with an edge whose relative distance to nea is within    */
  /* 20% of that request                                          */
  /*--------------------------------------------------------------*/
  RWGEdge *Ea=G->Surfaces[R->nsa]->Edges[R->nea];
  if( rRel != -1.0 )
   { 
     while(int neb=0; int neb<Sb->NumEdges; int neb++)
      { RWGEdge *Eb=G->Surfaces[R->nsb]->Edges[R->neb];
        double MaxRadius=fmax(Ea->Radius, Eb->Radius);
        double ThisrRel = VecDistance(Ea->Centroid, Eb->Centroid) / MaxRadius;
        if ( (0.8*rRel <= ThisrRel) && (ThisrRel <= 1.2*rRel) ) 
         { R->neb=neb;
           break;
         }l
      };
   };

#if 0
  if ( 0<=ncv && ncv<=3 )
   { R->nsb=R->nsa;
     R->npa=lrand48() % Sa->NumPanels;
     do
      { npb=lrand48() % Sa->NumPanels;
      } while( NumCommonVertices(Sa,R->npa,Sa,R->npb)!=ncv );
     R->iQa=lrand48() % 3;
     R->iQb=lrand48() % 3;
   };
#endif

}


/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[]) 
{ 
  /***************************************************************/
  /* process command-line arguments ******************************/
  /***************************************************************/
  char *GeoFileName=0;
  int Visualize=0;
  ArgStruct ASArray[]=
   { {"geometry",  PA_STRING, (void *)&GeoFileName, 0, ".rwggeo file"},
     {0,0,0,0,0}
   };
  ProcessArguments(argc, argv, ASArray);
  if (GeoFileName==0)
   ASUsage(argv[0],ASArray,"--geometry option is mandatory");

  srand48(time(0));
  SetDefaultCD2SFormat("(%+10.3e,%+10.3e)");

  /***************************************************************/
  /* create the geometry *****************************************/
  /***************************************************************/
  RWGGeometry *G = new RWGGeometry(GeoFileName);

  /***************************************************************/
  /* preinitialize an argument structure for the edge-edge       */
  /* interaction routines                                        */
  /***************************************************************/
  GetEEIArgStruct MyArgs, *Args=&MyArgs;

  /***************************************************************/
  /* enter command loop ******************************************/
  /***************************************************************/
  Request MyRequest, *R=&MyRequest;
  cdouble GCLS[2], GradGCLS[6]; // G,C integrals by libscuff method
  cdouble GCCM[2], GradGCCM[6]; // G,C integrals by cartesian multipole
  cdouble GCBF[2], GradGCBF[6]; // G,C integrals by brute force
  double LSTime, CMTime;
  for(;;)
   {
     GetRequest(R);

     /*--------------------------------------------------------------------*/
     /* print a little summary of the edge pair we will be considering     */
     /*--------------------------------------------------------------------*/
     //ncv=AssessEdgePair(Sa, nea, Ob, neb, &rRel);
     printf("*\n");
     printf("* --nsa %i --nea %i ",R->nsa, R->nea);
     printf("* --nsb %i --neb %i ",R->nsb, R->neb);
     printf("* --k %s",CD2S(R->k));
     printf("* relative distance: %+7.3e (DBFThreshold=10.0)\n",rRel);
     if (DZ!=0.0)
      printf("*  DZ:                %e\n",DZ);
     printf("*\n\n");

     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     Args->Sa                    = G->Surfaces[R->nsa];
     Args->Sb                    = G->Surfaces[R->nsb];
     Args->nea                   = R->nea;
     Args->neb                   = R->neb;
     Args->k                     = R->k;
     Args->Displacement          = R->DX;
     Args->NumGradientComponents = R->Gradient ? 3 : 0;
     Args->NumTorqueAxes         = 0;

     /*--------------------------------------------------------------------*/
     /* get edge-edge interactions by libscuff methods                     */
     /*--------------------------------------------------------------------*/
     Tic();
     for(int Times=0; Times<NUMTIMES; Times++)
      GetEdgeEdgeInteractions(Args);
     LSTime=Toc() / NUMTIMES;
     memcpy(GCLS, Args->GC, 2*sizeof(cdouble));
     memcpy(GradGCLS, Args->GradGC, 6*sizeof(cdouble));

     /*--------------------------------------------------------------------*/
     /* get edge-edge interactions by cartesian multipole method           */
     /*--------------------------------------------------------------------*/
     Tic();
     for(int Times=0; Times<NUMTIMES; Times++)
      GetEEIs_CartesianMultipole(Args);
     Toc();
     CMTime=Toc() / NUMTIMES;
     memcpy(GCCM, Args->GC, 2*sizeof(cdouble));
     memcpy(GradGCCM, Args->GradGC, 6*sizeof(cdouble));

     /*--------------------------------------------------------------------*/
     /* get panel-panel integrals by brute-force methods                   */
     /*--------------------------------------------------------------------*/
     GetEEIs_BruteForce(Args);
     memcpy(GCBF, Args->GC, 2*sizeof(cdouble));
     memcpy(GradGCBF, Args->GradGC, 6*sizeof(cdouble));


     /*--------------------------------------------------------------------*/
     /*- print results ----------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     printf("libscuff:  %e us\n",LSTime*1.0e6);
     printf("Cartesian: %e us\n",CMTime*1.0e6);
     printf("Quantity          | %23s | %23s | %23s |%s|%s\n",
            "     libscuff          ",
            "     C multipole       ",
            "     brute force       ", 
            "rdLS",
            "rdCM");
     printf("------------------|-%23s-|-%23s-|-%23s-|%s|%s\n", 
            "-----------------------", "-----------------------",
            "-----------------------", "----","----");
     printf("<fa|G|fb>         | %s | %s | %s |%4.1e|%4.1e\n", 
             CD2S(GCLS[0]), CD2S(GCCM[0]), CD2S(GCBF[0]), 
             RD(GCLS[0],GCBF[0]), RD(GCCM[0],GCBF[0]));
     printf("<fa|C|fb>         | %s | %s | %s |%4.1e|%4.1e\n", 
             CD2S(GCLS[1]), CD2S(GCCM[1]), CD2S(GCBF[1]), 
             RD(GCLS[1],GCBF[1]), RD(GCCM[1],GCBF[1]));

     if (Gradient)
      { 
         printf("\n");
         printf("(d/dx) <fa|G|fb>  | %s | %s | %s |%4.1e|%4.1e\n", 
                 CD2S(GradGCLS[0]), CD2S(GradGCCM[0]), CD2S(GradGCBF[0]), 
                 RD(GradGCLS[0],GradGCBF[0]), RD(GradGCCM[0],GradGCBF[0]));
         printf("(d/dx) <fa|C|fb>  | %s | %s | %s |%4.1e|%4.1e\n", 
                 CD2S(GradGCLS[1]), CD2S(GradGCCM[1]), CD2S(GradGCBF[1]), 
                 RD(GradGCLS[1],GradGCBF[1]), RD(GradGCCM[1],GradGCBF[1]));
      };
     printf("\n");

   }; // end of main command loop [ for(;;) ... ]

}
