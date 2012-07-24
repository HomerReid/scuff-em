/*
 * tSWPPI.cc    -- test short-wavelength panel-panel integrals 
 * 
 * homer reid   -- 11/2005 -- 11/2011
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#define HAVE_READLINE
#ifdef HAVE_READLINE
 #include <readline/readline.h>
 #include <readline/history.h>
#else
 #include "readlineReplacement.h"
#endif

#include <complex>

#include <libhrutil.h>

#include <libscuff.h>
#include <libMatProp.h>
#include "libscuffInternals.h"
#include "TaylorDuffy.h"

#define TDTIMES 10 

namespace scuff {

cdouble TaylorDuffy2(TaylorDuffyArgStruct *Args);
}

using namespace scuff;

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
   { {"geometry",  PA_STRING, (void *)&GeoFileName, 0, ".rwggeo file"},
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
  RWGObject *O=G->Objects[0];

  /***************************************************************/
  /* enter command loop ******************************************/
  /***************************************************************/
  using_history();
  read_history(0);
  int nt, NumTokens;
  char *Tokens[50];
  char *p;
  int npa, npb, iQa, iQb, ncv;
  cdouble K, Omega, Eps, Mu;
  cdouble HTD1[2], HTD2[2];
  RWGPanel *Pa, *Pb;
  double rRel;
  double *TVa[3], *Qa, *TVb[3], *Qb;
  double TD1GTime, TD1CTime, TD2GTime, TD2CTime;
  int nTimes;
  for(;;)
   { 
     /*--------------------------------------------------------------*/
     /*- print prompt and get input string --------------------------*/
     /*--------------------------------------------------------------*/
     printf("\n");
     printf(" options: --npa xx \n");
     printf("          --iQa xx \n");
     printf("          --npb xx \n");
     printf("          --iQb xx \n");
     printf("          --ncv xx \n");
     printf("          --kr     \n");
     printf("          --ki     \n");
     printf("          --quit\n");
     p=readline("enter options: ");
     if (!p) break;
     add_history(p);
     write_history(0);

     /*--------------------------------------------------------------*/
     /* parse input string                                          -*/
     /*--------------------------------------------------------------*/  
     NumTokens=Tokenize(p,Tokens,50);
     npa=npb=ncv=-1;
     iQa=iQb=0;
     Omega=1.0;
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
      if ( !strcasecmp(Tokens[nt],"--omega") )
       S2CD(Tokens[nt+1],&Omega);
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--quit") )
       exit(1);
       
     free(p);
  
     /*--------------------------------------------------------------*/
     /* if the user specified a number of common vertices then       */
     /* find a randomly-chosen pair of panels with that number of    */
     /* common vertices                                              */
     /*--------------------------------------------------------------*/
     if ( npa==-1 || npb==-1 ) 
      ncv=1 + lrand48()%2;
     if ( 0<ncv && ncv<=3 )
      { npa=lrand48() % O->NumPanels;
        do
         { npb=lrand48() % O->NumPanels;
         } while( NumCommonVertices(O,npa,O,npb)!=ncv );
      };

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     O->MP->GetEpsMu(Omega, &Eps, &Mu);
     K=csqrt2(Eps*Mu)*Omega;

     /*--------------------------------------------------------------------*/
     /* print a little summary of the panel pair we will be considering    */
     /*--------------------------------------------------------------------*/
     Pa=O->Panels[npa];
     Pb=O->Panels[npb];
     ncv=AssessPanelPair(O,npa,O,npb,&rRel,TVa,TVb);
     printf("*\n");
     printf("* --npa %i --iQa %i (V #%i) --npb %i --iQb %i (V #%i)\n",
            npa,iQa,Pa->VI[iQa],npb,iQb,Pb->VI[iQb]);
     printf("*  common vertices:   %i\n",ncv);
     printf("*  eps:               %s\n",CD2S(Eps));
     printf("*  K:                 %s\n",CD2S(K));
     printf("*  k*MaxRadius:       %.1e\n",abs(K*fmax(Pa->Radius, Pb->Radius)));
     printf("*\n\n");

     if (ncv==0) continue;

     /*--------------------------------------------------------------------*/
     /* get panel-panel integrals by taylor-duffy method                   */
     /*--------------------------------------------------------------------*/
     Qa = O->Vertices + 3*Pa->VI[iQa];
     Qb = O->Vertices + 3*Pb->VI[iQb];

     double *OVa[3], *OVb[3];
     int Flipped=CanonicallyOrderVertices(TVa, TVb, ncv, OVa, OVb);

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     TaylorDuffyArgStruct MyTDArgStruct, *TDArgs=&MyTDArgStruct;
     InitTaylorDuffyArgs(TDArgs);
     TDArgs->WhichCase = ncv;
     TDArgs->GParam    = K;
     TDArgs->V1        = OVa[0];
     TDArgs->V2        = OVa[1];
     TDArgs->V3        = OVa[2];
     TDArgs->V2P       = OVb[1];
     TDArgs->V3P       = OVb[2];
     TDArgs->Q         = Flipped ? Qb : Qa ;
     TDArgs->QP        = Flipped ? Qa : Qb ;

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     TDArgs->WhichG=TM_EIKR_OVER_R;
     TDArgs->WhichH=TM_DOTPLUS;
     Tic();
     for(nTimes=0; nTimes<TDTIMES; nTimes++)
      HTD1[0]=TaylorDuffy(TDArgs);
     TD1GTime=Toc() / TDTIMES;

     if (ncv==3)
      { HTD1[1]=0.0;
        TD1CTime=0.0;
      }
     else
      { TDArgs->WhichG=TM_GRADEIKR_OVER_R;
        TDArgs->WhichH=TM_CROSS;
        Tic();
        for(nTimes=0; nTimes<TDTIMES; nTimes++)
         HTD1[1]=TaylorDuffy(TDArgs);
        TD1CTime=Toc() / TDTIMES;
      }

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     TDArgs->WhichG=TM_EIKR_OVER_R;
     TDArgs->WhichH=TM_DOTPLUS;
     Tic();
     for(nTimes=0; nTimes<TDTIMES; nTimes++)
      HTD2[0]=TaylorDuffy2(TDArgs);
     TD2GTime=Toc() / TDTIMES;

     if (ncv==3)
      { HTD2[1]=0.0;
        TD2CTime=0.0;
      }
     else
      { TDArgs->WhichG=TM_GRADEIKR_OVER_R;
        TDArgs->WhichH=TM_CROSS;
        Tic();
        for(nTimes=0; nTimes<TDTIMES; nTimes++)
         HTD2[1]=TaylorDuffy2(TDArgs);
        TD2CTime=Toc() / TDTIMES;
      };

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     printf("\n\n");
     printf("libscuff:    G: %5.2e us    C: %5.2e us \n",1e6*TD1GTime,1e6*TD1CTime);
     printf("modified:    G: %5.2e us    C: %5.2e us \n",1e6*TD2GTime,1e6*TD2CTime);
     printf("\n");

     printf("G(scuff): %s \n",CD2S(HTD1[0]));
     printf("G(mod  ): %s \n",CD2S(HTD2[0]));
     printf("RD:       %5.2e \n",RD(HTD1[0],HTD2[0]));
     printf("\n");

     printf("C(scuff): %s \n",CD2S(HTD1[1]));
     printf("C(mod  ): %s \n",CD2S(HTD2[1]));
     printf("RD:       %5.2e \n",RD(HTD1[1],HTD2[1]));
     printf("\n");


   }; // end of main command loop [ for(;;) ... ]

}
