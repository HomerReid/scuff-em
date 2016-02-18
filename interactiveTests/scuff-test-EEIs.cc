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

#include <math.h>
#include <readline/readline.h>
#include <readline/history.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <libhrutil.h>
#include <PanelCubature.h>

#ifdef HAVE_LIBPOLYFIT
#include <libPolyFit.h>
#endif

#include "scuff-test-EEIs.h"

#define II cdouble(0.0,1.0)
#define NUMTIMES 10

/***************************************************************/
/* integrand for brute-force edge-edge integration            **/
/***************************************************************/
void Integrand(double *X, double *XP, PPCData *PPCD,
               void *UserData, double *Result)
{ 
  double Eta;
  if (UserData)
   Eta = *((double *)UserData);
  else
   Eta = 0.0;

  cdouble *F1    = PPCD->K1;    // RWG basis function 1 
  cdouble DivF1  = PPCD->DivK1;
  cdouble *F2    = PPCD->K2;    // RWG basis function 2
  cdouble DivF2  = PPCD->DivK2;

  double R[3];
  VecSub(X, XP, R);
  double r2=(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
  double rE2=r2+Eta*Eta;
  double r=sqrt(r2);
  double rE=sqrt(rE2);

  cdouble DotProduct    = F1[0]*F2[0] + F1[1]*F2[1] + F1[2]*F2[2];
  cdouble ScalarProduct = DivF1 * DivF2;
  cdouble TripleProduct =  F1[0] * (R[1]*F2[2]-R[2]*F2[1])
                          +F1[1] * (R[2]*F2[0]-R[0]*F2[2])
                          +F1[2] * (R[0]*F2[1]-R[1]*F2[0]);

  cdouble ik  = II * PPCD->Omega;
  cdouble ikr = ik*r;
  cdouble Phi = exp(ikr)/(4.0*M_PI*rE);
  cdouble Psi = (ikr-1.0)*Phi/rE2;

  cdouble *zf = (cdouble *)Result;
  zf[0] = ( DotProduct + ScalarProduct/(ik*ik) ) * Phi;
  zf[1] = TripleProduct * Psi;

}

/***************************************************************/
/* compute edge-edge integrals by brute force                  */
/***************************************************************/
void GetEEIs_BruteForce(RWGGeometry *G, GetEEIArgStruct *Args)
{
  int nea = Args->nea;
  int neb = Args->neb;

  int ncv = AssessBFPair(Args->Sa,nea,Args->Sb,neb);

  if (ncv==0)
   { 
#if 0
     GetBFBFCubature(G, 
                     Args->Sa->Index, Args->nea,
                     Args->Sb->Index, Args->neb,
                     Args->Displacement,
                     Integrand, 0, 4, 100000, 1.0e-10, 0.0,
                     Args->k, 0, (double *)Args->GC);
#endif
   }
  else
   {
#ifdef HAVE_LIBPOLYFIT
     #define NUMETAS 5
     double EtaValues[NUMETAS]={0.002,0.004,0.006,0.008,0.01};
     double R1Values[NUMETAS], I1Values[NUMETAS], R2Values[NUMETAS], I2Values[NUMETAS];

     for(int neta=0; neta<NUMETAS; neta++)
      { 
        double Eta = EtaValues[neta];
        double Result[4];
        GetBFBFCubature(G, Args->Sa->Index, Args->nea, 
                           Args->Sb->Index, Args->neb, 
                           0, Integrand, (void *)&Eta, 4, 
                           1000000, 1.0e-10, 0.0, Args->k, 0, Result);
        R1Values[neta] = Result[0];
        I1Values[neta] = Result[1];
        R2Values[neta] = Result[2];
        I2Values[neta] = Result[3];
      };

     PolyFit *PF;
     PF = new PolyFit(EtaValues, R1Values, NUMETAS, 2);
     real(Args->GC[0]) = PF->f(0.0);
     PF->PlotFit(EtaValues, R1Values, NUMETAS);
     delete PF;

     PF = new PolyFit(EtaValues, I1Values, NUMETAS, 2);
     imag(Args->GC[0]) = PF->f(0.0);
     PF->PlotFit(EtaValues, I1Values, NUMETAS);
     delete PF;

     PF = new PolyFit(EtaValues, R2Values, NUMETAS, 2);
     real(Args->GC[1]) = PF->f(0.0);
     PF->PlotFit(EtaValues, R2Values, NUMETAS);
     delete PF;

     PF = new PolyFit(EtaValues, I2Values, NUMETAS, 2);
     imag(Args->GC[1]) = PF->f(0.0);
     PF->PlotFit(EtaValues, I2Values, NUMETAS);
     delete PF;
#endif

   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetEEIs_CartesianMultipole(GetEEIArgStruct *Args)
{ (void *)Args;
}

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
  int ncv=-1;
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
   if ( !strcasecmp(Tokens[nt],"--ncv") )
    sscanf(Tokens[nt+1],"%i",&ncv);
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
     for(int neb=0; neb<Sb->NumEdges; neb++)
      { RWGEdge *Eb=G->Surfaces[R->nsb]->Edges[R->neb];
        double MaxRadius=fmax(Ea->Radius, Eb->Radius);
        double ThisrRel = VecDistance(Ea->Centroid, Eb->Centroid) / MaxRadius;
        if ( (0.8*rRel <= ThisrRel) && (ThisrRel <= 1.2*rRel) ) 
         { R->neb=neb;
           break;
         };
      };
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if ( 0<=ncv && ncv<=4 )
   { R->nsb=R->nsa;
     R->nea=lrand48() % Sa->NumEdges;
     do
      { R->neb=lrand48() % Sa->NumEdges;
      } while( AssessBFPair(Sa,R->nea,Sa,R->neb)!=ncv );
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void DoVisualization(RWGGeometry *G, Request *R)
{ 
  int nea        = R->nea;
  int neb        = R->neb;

  RWGSurface *Sa = G->Surfaces[R->nsa];
  RWGEdge    *Ea = Sa->Edges[nea];
  double *QPa    = Sa->Vertices + 3*(Ea->iQP);
  double *V1a    = Sa->Vertices + 3*(Ea->iV1);
  double *V2a    = Sa->Vertices + 3*(Ea->iV2);
  double *QMa    = (Ea->iQM==-1) ? 0 : Sa->Vertices + 3*(Ea->iQM);

  RWGSurface *Sb = G->Surfaces[R->nsb];
  RWGEdge    *Eb = Sb->Edges[neb];
  double *QPb    = Sb->Vertices + 3*(Eb->iQP);
  double *V1b    = Sb->Vertices + 3*(Eb->iV1);
  double *V2b    = Sb->Vertices + 3*(Eb->iV2);
  double *QMb    = (Eb->iQM==-1) ? 0 : Sb->Vertices + 3*(Eb->iQM);

  char FileName[100];
  snprintf(FileName,100,"/tmp/E%i_%i.log",nea,neb);
  FILE *f=fopen(FileName,"w");

  fprintf(f,"%e %e %e \n",QPa[0],QPa[1],QPa[2]);
  fprintf(f,"%e %e %e \n",V1a[0],V1a[1],V1a[2]);
  fprintf(f,"%e %e %e \n",V2a[0],V2a[1],V2a[2]);
  fprintf(f,"%e %e %e \n",QPa[0],QPa[1],QPa[2]);
  fprintf(f,"\n\n");

  if (QMa)
   {
     fprintf(f,"%e %e %e \n",QMa[0],QMa[1],QMa[2]);
     fprintf(f,"%e %e %e \n",V1a[0],V1a[1],V1a[2]);
     fprintf(f,"%e %e %e \n",V2a[0],V2a[1],V2a[2]);
     fprintf(f,"%e %e %e \n",QMa[0],QMa[1],QMa[2]);
     fprintf(f,"\n\n");
   };

  fprintf(f,"%e %e %e \n",QPb[0],QPb[1],QPb[2]);
  fprintf(f,"%e %e %e \n",V1b[0],V1b[1],V1b[2]);
  fprintf(f,"%e %e %e \n",V2b[0],V2b[1],V2b[2]);
  fprintf(f,"%e %e %e \n",QPb[0],QPb[1],QPb[2]);
  fprintf(f,"\n\n");

  if (QMb)
   {
     fprintf(f,"%e %e %e \n",QMb[0],QMb[1],QMb[2]);
     fprintf(f,"%e %e %e \n",V1b[0],V1b[1],V1b[2]);
     fprintf(f,"%e %e %e \n",V2b[0],V2b[1],V2b[2]);
     fprintf(f,"%e %e %e \n",QMb[0],QMb[1],QMb[2]);
     fprintf(f,"\n\n");
   };
  fclose(f);

  f=vfopen("/tmp/E%i_%i.gp","w",nea,neb);
  fprintf(f,"splot '%s' i 0 t 'E%i_P' w lp ",FileName,nea);
  int Index=1;
  if (QMa)
   fprintf(f,", '%s' i %i t 'E%i_M' w lp ",FileName,Index++,neb);
  fprintf(f,",  '%s' i %i t 'E%i_P' w lp ",FileName,Index++,neb);
  if (QMb)
   fprintf(f,", '%s' i %i t 'E%i_M' w lp ",FileName,Index++,neb);
  fprintf(f,"\n\n");
  fclose(f);

  //system("gnuplot -persist /tmp/Edges.gp & ");

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
     {"visualize", PA_BOOL,   (void *)&Visualize,   0, "visualize"},
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
  InitGetEEIArgs(Args);

  /***************************************************************/
  /* enter command loop ******************************************/
  /***************************************************************/
  for(;;)
   {
     Request MyRequest, *R=&MyRequest;
     GetRequest(G,R);

     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     if (Visualize)
      DoVisualization(G,R);

     /*--------------------------------------------------------------------*/
     /* print a little summary of the edge pair we will be considering     */
     /*--------------------------------------------------------------------*/
     //ncv=AssessEdgePair(Sa, nea, Ob, neb, &rRel);
     printf("*\n");
     printf("* --nsa %i --nea %i ",R->nsa, R->nea);
     printf("* --nsb %i --neb %i ",R->nsb, R->neb);
     printf("* --k %s",CD2S(R->k));
     printf("* --ncv %i",AssessBFPair(G->Surfaces[R->nsa], R->nea, G->Surfaces[R->nsb], R->neb));
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
     cdouble GCLS[2], GradGCLS[6]; // G,C integrals by libscuff method
     Tic();
     for(int Times=0; Times<NUMTIMES; Times++)
      GetEdgeEdgeInteractions(Args);
     double LSTime=Toc() / NUMTIMES;
     memcpy(GCLS, Args->GC, 2*sizeof(cdouble));
     memcpy(GradGCLS, Args->GradGC, 6*sizeof(cdouble));

     /*--------------------------------------------------------------------*/
     /* get edge-edge interactions by cartesian multipole method           */
     /*--------------------------------------------------------------------*/
     cdouble GCCM[2], GradGCCM[6]; // G,C integrals by cartesian multipole
     double CMTime=0.0;
#if 0
     Tic();
     for(int Times=0; Times<NUMTIMES; Times++)
      GetEEIs_CartesianMultipole(Args);
     Toc();
     CMTime=Toc() / NUMTIMES;
     memcpy(GCCM, Args->GC, 2*sizeof(cdouble));
     memcpy(GradGCCM, Args->GradGC, 6*sizeof(cdouble));
#endif

     /*--------------------------------------------------------------------*/
     /* get panel-panel integrals by brute-force methods                   */
     /*--------------------------------------------------------------------*/
     cdouble GCBF[2], GradGCBF[6]; // G,C integrals by brute force
#if 0
     GetEEIs_BruteForce(G, Args);
     memcpy(GCBF, Args->GC, 2*sizeof(cdouble));
     memcpy(GradGCBF, Args->GradGC, 6*sizeof(cdouble));
#endif

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

     if (R->Gradient)
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
