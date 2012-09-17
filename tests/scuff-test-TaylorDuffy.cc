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
 * scuff-test-TaylorDuffy.cc -- a test program for libscuff's implementation
 *                           -- of the Taylor-Duffy method for evaluating 
 *                           -- panel-panel integrals between pairs of panels
 *                           -- with common vertices
 * 
 * homer reid                -- 11/2005 -- 1/2012
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "config.h"

#ifdef HAVE_LIBREADLINE
 #include <readline/readline.h>
 #include <readline/history.h>
#else
 // dummy functions
 void using_history() {} 
 void read_history(const char *) {}
 void add_history(const char *) {}
 void write_history(const char *) {}
 char *readline(const char *prompt) 
  { char str[1000];
    fputs(prompt, stdout);
    fgets(str,1000,stdin);
    return strdup(str);
  }
#endif

#include <libhrutil.h>
#include <libscuff.h>
#include <libscuffInternals.h>
#include <TaylorDuffy.h>

using namespace scuff;

#define ABSURD -128.0

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PrintParameters(TaylorDuffyArgStruct *Args)
{
  int WhichCase   = Args->WhichCase;
  int WhichG      = Args->WhichG;
  int WhichH      = Args->WhichH;
  double *V1      = Args->V1;
  double *V2      = Args->V2;
  double *V3      = Args->V3;
  double *V2P     = Args->V2P;
  double *V3P     = Args->V3P;
  double *Q       = Args->Q;
  double *QP      = Args->QP;

  double A[3], AP[3], B[3], BP[3], L[3], D[3], DP[3], QmQP[3], QxQP[3], TV[3];
  TMWorkspace MyTMW, *TMW=&MyTMW;

  VecSub(V2,V1,A);
  VecSub(V3,V2,B);
  VecSub(V2P,V1,AP);
  VecSub(V3P,V2P,BP);
  VecSub(BP,B,L);

  TMW->A2    = VecDot(A,A);
  TMW->B2    = VecDot(B,B);
  TMW->AP2   = VecDot(AP,AP);
  TMW->BP2   = VecDot(BP,BP);
  TMW->L2    = VecDot(L,L);
  TMW->AdB   = VecDot(A,B);
  TMW->AdAP  = VecDot(A,AP);
  TMW->AdBP  = VecDot(A,BP);
  TMW->AdL   = VecDot(A,L);
  TMW->BdAP  = VecDot(B,AP);
  TMW->BdBP  = VecDot(B,BP);
  TMW->APdBP = VecDot(AP,BP);
  TMW->BPdL  = VecDot(BP,L);

  if (WhichH==TM_DOT || WhichH==TM_DOTPLUS )
   {
     VecSub(V1,Q,D);
     VecSub(V1,QP,DP);
     TMW->AdD   = VecDot(A,D);
     TMW->AdDP  = VecDot(A,DP);
     TMW->BdDP  = VecDot(B,DP);
     TMW->APdD  = VecDot(AP,D);
     TMW->BPdD  = VecDot(BP,D);
     TMW->DdDP  = VecDot(D,DP);
   }
  else if (WhichH==TM_CROSS)
   { 
     VecSub(Q,QP,QmQP);
     VecCross(Q,QP,QxQP);

     TMW->AdQxQP     = VecDot(A,  QxQP );
     TMW->APdQxQP    = VecDot(AP, QxQP );
     TMW->BdQxQP     = VecDot(B,  QxQP );
     TMW->BPdQxQP    = VecDot(BP, QxQP );
     TMW->LdQxQP     = VecDot(L,  QxQP );
     TMW->V1xAdQmQP  = VecDot( VecCross(V1,A,TV),  QmQP );
     TMW->V1xAPdQmQP = VecDot( VecCross(V1,AP,TV), QmQP );
     TMW->V1xBdQmQP  = VecDot( VecCross(V1,B,TV),  QmQP );
     TMW->V1xBPdQmQP = VecDot( VecCross(V1,BP,TV), QmQP );
     TMW->AxAPdQmQP  = VecDot( VecCross(A,AP,TV),  QmQP );
     TMW->AxBdQmQP   = VecDot( VecCross(A,B,TV),   QmQP );
     TMW->AxBPdQmQP  = VecDot( VecCross(A,BP,TV),  QmQP );
     TMW->BxAPdQmQP  = VecDot( VecCross(B,AP,TV),  QmQP );
     TMW->BxBPdQmQP  = VecDot( VecCross(B,BP,TV),  QmQP );

   }
  else if (WhichH==TM_ENORMAL)
   {
     /* note that the ENormal choice of 'h' function depends   */
     /* on a vector 'ZHat', which we pass via the 'Q' argument */
     /* which is otherwise not used for this h function        */
     TMW->APdZHat = VecDot(AP, Q);
     TMW->BPdZHat = VecDot(BP, Q);
   };

  TMW->A2    = VecDot(A,A);
  TMW->B2    = VecDot(B,B);
  TMW->AP2   = VecDot(AP,AP);
  TMW->BP2   = VecDot(BP,BP);
  TMW->L2    = VecDot(L,L);
  TMW->AdB   = VecDot(A,B);
  TMW->AdAP  = VecDot(A,AP);
  TMW->AdBP  = VecDot(A,BP);
  TMW->AdL   = VecDot(A,L);
  TMW->BdAP  = VecDot(B,AP);
  TMW->BdBP  = VecDot(B,BP);
  TMW->APdBP = VecDot(AP,BP);
  TMW->BPdL  = VecDot(BP,L);


  printf("%s = %.10f; \n","A2",    TMW->A2);
  printf("%s = %.10f; \n","B2",    TMW->B2);
  printf("%s = %.10f; \n","AP2",   TMW->AP2);
  printf("%s = %.10f; \n","BP2",   TMW->BP2);
  printf("%s = %.10f; \n","L2",    TMW->L2);
  printf("%s = %.10f; \n","AdB",   TMW->AdB);
  printf("%s = %.10f; \n","AdAP",  TMW->AdAP);
  printf("%s = %.10f; \n","AdBP",  TMW->AdBP);
  printf("%s = %.10f; \n","AdL",   TMW->AdL);
  printf("%s = %.10f; \n","BdAP",  TMW->BdAP);
  printf("%s = %.10f; \n","BdBP",  TMW->BdBP);
  printf("%s = %.10f; \n","APdBP", TMW->APdBP);
  printf("%s = %.10f; \n","BPdL",  TMW->BPdL);

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
  ArgStruct ASArray[]=
   { {"geometry",   PA_STRING, (void *)&GeoFileName, 0, ".scuffgeo file"},
     {0,0,0,0,0}
   };
  ProcessArguments(argc, argv, ASArray);
  if (GeoFileName==0)
   ASUsage(argv[0],ASArray,"--geometry option is mandatory");

  srand48(time(0));
  SetDefaultCD2SFormat("(%15.8e,%15.8e)");

  /***************************************************************/
  /* create the geometry *****************************************/
  /***************************************************************/
  RWGGeometry *G = new RWGGeometry(GeoFileName);
  RWGSurface *S=G->Surfaces[0];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  TaylorDuffyArgStruct Args;
  InitTaylorDuffyArgs(&Args);

  /***************************************************************/
  /* enter command loop ******************************************/
  /***************************************************************/
  int nt, NumTokens;
  char *Tokens[50];
  char *p;
  int npa, iQa, npb, iQb, ncv;
  double P;
  cdouble k, TDOne, TDDot;
  int Times, NumTimes;
  double ElapsedOne, ElapsedDot;
  double *Va[3], *OVa[3], *Vb[3], *OVb[3], *Qa, *Qb;
  for(;;)
   { 
     /*--------------------------------------------------------------*/
     /*- print prompt and get input string --------------------------*/
     /*--------------------------------------------------------------*/
     printf(" options: --npa xx \n");
     printf("          --iQa xx \n");
     printf("          --npb xx \n");
     printf("          --iQb xx \n");
     printf("          --ncv xx \n");
     printf("          --P (power of r)\n"); 
     printf("          --k (k parameter in e^{ikr}/r})\n");
     printf("          --same | --ns \n");
     printf("          --DZ \n");
     p=readline("enter options: ");
     if (!p) break;
     add_history(p);
     write_history(0);

     /*--------------------------------------------------------------*/
     /* parse input string                                          -*/
     /*--------------------------------------------------------------*/
     npa=iQa=npb=iQb=ncv=-1;
     NumTimes=10;
     P=ABSURD;
     k=1.0;
     NumTokens=Tokenize(p,Tokens,50);
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--npa") )
       sscanf(Tokens[nt+1],"%i",&npa);
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--iQa") )
       sscanf(Tokens[nt+1],"%i",&iQa);
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--npb") )
       sscanf(Tokens[nt+1],"%i",&npb);
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--iQb") )
       sscanf(Tokens[nt+1],"%i",&iQb);
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--ncv") )
       sscanf(Tokens[nt+1],"%i",&ncv);
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--NumTimes") )
       sscanf(Tokens[nt+1],"%i",&NumTimes);
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--P") )
       sscanf(Tokens[nt+1],"%le",&P);
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--k") )
       S2CD(Tokens[nt+1],&k);
     free(p);

     /*--------------------------------------------------------------*/
     /* if the user didn't specify a pair of panels, then find a pair*/
     /* with ncv common vertices, where ncv was either specified by  */
     /* the user or is chosen randomly by us between 1 and 3 inclusive*/
     /*--------------------------------------------------------------*/
     if ( npa==-1 || npb==-1 )
      { npa=lrand48() % S->NumPanels;
        if (ncv==-1) ncv=1+(lrand48()%3);
        printf("Looking for a panel pair with %i common vertices...\n",ncv);
        do
         { npb=lrand48() % S->NumPanels;
         } while( NumCommonVertices(S,npa,S,npb)!=ncv );
      };
     if ( iQa==-1 )
      iQa = lrand48() % 3;
     if ( iQb==-1 )
      iQb = lrand48() % 3;

     /*--------------------------------------------------------------*/
     /*- if the user specified specific panel indices then make sure */
     /*- they make sense                                             */
     /*--------------------------------------------------------------*/
     if ( npa >= S->NumPanels )
      { printf("whoops! Sject %s has only %i panels (you requested panel %i).\n",
                S->Label,S->NumPanels,npa);
        continue;
      };
     if ( npb >= S->NumPanels )
      { printf("whoops! Sject %s has only %i panels (you requested panel %i).\n",
                S->Label,S->NumPanels,npb);
        continue;
      };

     /*--------------------------------------------------------------------*/
     /* print a little summary of the panel pair we will be considering    */
     /*--------------------------------------------------------------------*/
     double rRel;
     ncv=AssessPanelPair(S, npa, S, npb, &rRel, Va, Vb);
     CanonicallyOrderVertices(Va, Vb, ncv, OVa, OVb);
     printf("*\n");
     printf("* --npa %i --iQa %i --npb %i --iQb %i \n",npa,npb,iQa,iQb);
     printf("*  common vertices: %i\n",ncv);
     if (P==ABSURD)
      printf("*  kernel: e^{ikr/4\\pi r} with k=%s\n",z2s(k));
     else
      printf("*  kernel: r^%g\n",P);
     printf("*\n\n");

     /*--------------------------------------------------------------------*/
     /*- print results ----------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     Args.V1  = OVa[0];
     Args.V2  = OVa[1];
     Args.V3  = OVa[2];
     Args.Q   = Va[iQa];
     Args.V2P = OVb[1];
     Args.V3P = OVb[2];
     Args.QP  = Vb[iQb];
     Args.WhichCase = ncv;
  
     if (P==ABSURD)
      { Args.WhichG = TM_EIKR_OVER_R;
        Args.GParam = k;
      }
     else
      { Args.WhichG = TM_RP;
        Args.GParam = P;
      };

     Args.WhichH    = TM_ONE;
     Tic();
     for(int Times=0; Times<NumTimes; Times++)
      TDOne = TaylorDuffy(&Args);
     ElapsedOne=Toc() / NumTimes;

     Args.WhichH    = TM_DOT;
     for(int Times=0; Times<NumTimes; Times++)
      TDDot = TaylorDuffy(&Args);
     ElapsedDot=Toc() / NumTimes;

     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     PrintParameters(&Args);
  
     /*--------------------------------------------------------------------*/
     /*- print results ----------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     printf("\n");
     printf("One (%10f ms) %s\n",1e3*ElapsedOne,CD2S(TDOne));
     printf("\n");
     printf("Dot (%10f ms) %s\n",1e3*ElapsedDot,CD2S(TDDot));

   }; // end of main command loop [ for(;;) ... ]

}
