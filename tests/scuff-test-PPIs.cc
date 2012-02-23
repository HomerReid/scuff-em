/*
 * scuff-test-PPIs.cc -- a test program for libscuff's routines for
 *                    -- computing panel-panel integrals
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

#include <complex>

#include <libhrutil.h>

#include <libscuff.h>
#include <libscuffInternals.h>
#include "TaylorDuffy.h"

using namespace scuff;

//#define HRTIMES 100
#define HRTIMES 1
#define TDTIMES 1

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetPPIs_BruteForce(GetPPIArgStruct *Args, int PlotFits);

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PrintResults(cdouble *HLS, cdouble *GradHLS, 
                  cdouble *HBF, cdouble *GradHBF, 
                  int Gradient)
{ 
  printf("\n");
  printf("Quantity          |  %33s  |  %33s  | %s\n", 
         "            libscuff             ",
         "          brute force            ",
         "rel delta");
  printf("------------------|--%33s--|--%33s--|-%s\n", 
         "---------------------------------",
         "---------------------------------",
         "---------");
  printf("<fa|G|fb>         |  %s  |  %s  | %5.2e\n", CD2S(HLS[0]), CD2S(HBF[0]), RD(HLS[0],HBF[0]));
  printf("<fa|C|fb>         |  %s  |  %s  | %5.2e\n", CD2S(HLS[1]), CD2S(HBF[1]), RD(HLS[1],HBF[1]));
  if (Gradient)
   { 
      printf("\n");
      printf("(d/dx) <fa|G|fb>  |  %s  |  %s  | %5.2e\n", CD2S(GradHLS[0]), CD2S(GradHBF[0]), RD(GradHLS[0],GradHBF[0]));
      printf("(d/dy) <fa|G|fb>  |  %s  |  %s  | %5.2e\n", CD2S(GradHLS[2]), CD2S(GradHBF[2]), RD(GradHLS[2],GradHBF[2]));
      printf("(d/dz) <fa|G|fb>  |  %s  |  %s  | %5.2e\n", CD2S(GradHLS[4]), CD2S(GradHBF[4]), RD(GradHLS[4],GradHBF[4]));
      printf("(d/dx) <fa|C|fb>  |  %s  |  %s  | %5.2e\n", CD2S(GradHLS[1]), CD2S(GradHBF[1]), RD(GradHLS[1],GradHBF[1]));
      printf("(d/dy) <fa|C|fb>  |  %s  |  %s  | %5.2e\n", CD2S(GradHLS[3]), CD2S(GradHBF[3]), RD(GradHLS[3],GradHBF[3]));
      printf("(d/dz) <fa|C|fb>  |  %s  |  %s  | %5.2e\n", CD2S(GradHLS[5]), CD2S(GradHBF[5]), RD(GradHLS[5],GradHBF[5]));
   };
  printf("\n");

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void PrintResults(cdouble *HLS, cdouble *HTD, cdouble *HBF)
{ 
  printf("\n");
  printf("Quantity  | %33s | %33s | %33s |%s\n", 
         "            libscuff             ",
         "          taylor-duffy           ",
         "          brute force            ",
         "rd(LS-TD)");
  printf("----------|-%33s-|-%33s-|-%33s-|-%s\n",
         "---------------------------------",
         "---------------------------------",
         "---------------------------------",
         "---------");
  printf("<fa|G|fb> | %s | %s | %s | %5.2e\n", CD2S(HLS[0]), CD2S(HTD[0]), CD2S(HBF[0]), RD(HLS[0],HTD[0]));
  printf("<fa|C|fb> | %s | %s | %s | %5.2e\n", CD2S(HLS[1]), CD2S(HTD[1]), CD2S(HBF[1]), RD(HLS[1],HTD[1]));

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
     {"visualize", PA_BOOL,   (void *)&Visualize,   0, "write visualization files"},
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
  RWGObject *Oa=G->Objects[0];
  RWGObject *Ob=G->NumObjects>1 ? G->Objects[1] : Oa;

  /*--------------------------------------------------------------*/
  /*- write visualization files if requested ---------------------*/
  /*--------------------------------------------------------------*/
  if (Visualize)
   { char buffer[1000];
     char *Base=GetFileBase(G->GeoFileName);
     snprintf(buffer,1000,"%s.Visualization",Base);
     mkdir(buffer,0755);
     snprintf(buffer,1000,"%s.Visualization/%s.pp",Base,Base);
     G->WritePPMesh(buffer,"Geometry");
     G->WriteGPMeshPlus(buffer);
   };

  /***************************************************************/
  /* preinitialize an argument structure for the panel-panel     */
  /* integration routines                                        */
  /***************************************************************/
  GetPPIArgStruct MyArgs, *Args=&MyArgs;
  InitGetPPIArgs(Args);

  /***************************************************************/
  /* enter command loop ******************************************/
  /***************************************************************/
  using_history();
  read_history(0);
  int nt, NumTokens;
  char *Tokens[50];
  char *p;
  int npa, npb, iQa, iQb, ncv;
  int SameObject, Gradient, PlotFits;
  double rRel, rRelRequest, DZ;
  cdouble K;
  cdouble HLS[2], GradHLS[6]; // G,C integrals by libscuff 
  cdouble HTD[2];             // G,C integrals by taylor-duffy method
  cdouble HBF[2], GradHBF[6]; // G,C integrals by brute force
  RWGPanel *Pa, *Pb;
  double *TVa[3], *Qa, *TVb[3], *Qb;
  int nTimes;
  double HRTime, TDTime;
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
     printf("          --rRel xx \n");
     printf("          --fTD  (force taylor-duffy) \n");
     printf("          --same | --ns \n");
     printf("          --DZ \n");
     printf("          --kr     \n");
     printf("          --ki     \n");
     printf("          --gradient\n");
     printf("          --PlotFits\n");
     printf("          --SWPPITol xx\n");
     printf("          --quit\n");
     p=readline("enter options: ");
     if (!p) break;
     add_history(p);
     write_history(0);

     /*--------------------------------------------------------------*/
     /* parse input string                                          -*/
     /*--------------------------------------------------------------*/
     NumTokens=Tokenize(p,Tokens,50);
     npa=npb=iQa=iQb=ncv=SameObject=-1;
     Args->ForceTaylorDuffy=0;
     Gradient=PlotFits=0;
     DZ=rRelRequest=0.0;
     real(K) = imag(K) = INFINITY;
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
      if ( !strcasecmp(Tokens[nt],"--kr") )
       sscanf(Tokens[nt+1],"%le",&(real(K)));
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--ki") )
       sscanf(Tokens[nt+1],"%le",&(imag(K)));
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--same") )
       SameObject=1;
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--ns") )
       SameObject=0;
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--ftd") )
       Args->ForceTaylorDuffy=1;
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--DZ") )
       sscanf(Tokens[nt+1],"%le",&DZ);
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--rRel") )
       sscanf(Tokens[nt+1],"%le",&rRelRequest);
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--Gradient") )
       Gradient=1;
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--PlotFits") )
       PlotFits=1;
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--quit") )
       exit(1);
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--SWPPITol") )
       sscanf(Tokens[nt+1],"%le",&(RWGGeometry::SWPPITol));
       
     free(p);
  
     /*--------------------------------------------------------------*/
     /* if the user specified a number of common vertices then       */
     /* find a randomly-chosen pair of panels with that number of    */
     /* common vertices                                              */
     /*--------------------------------------------------------------*/
     if ( 0<=ncv && ncv<=3 )
      { SameObject=1;
        npa=lrand48() % Oa->NumPanels;
        do
         { npb=lrand48() % Oa->NumPanels;
         } while( NumCommonVertices(Oa,npa,Oa,npb)!=ncv );
        iQa=lrand48() % 3;
        iQb=lrand48() % 3;
      }
     /*--------------------------------------------------------------*/
     /*- if the user specified a value for rRel then try to find a --*/
     /*- panel pair whose relative distance is within 10% of rRel  --*/
     /*--------------------------------------------------------------*/
     else if ( rRelRequest!=0.0 )
      { 
        iQa=lrand48() % 3;
        iQb=lrand48() % 3;
        npa=lrand48() % Oa->NumPanels;

        /*--------------------------------------------------------------*/
        /*- first look on same object ----------------------------------*/
        /*--------------------------------------------------------------*/
        SameObject=1;
        for(npb=0; npb<Oa->NumPanels; npb++)
         { AssessPanelPair(Oa,npa,Oa,npb,&rRel);
           if ( 0.9*rRelRequest<rRel && rRel<1.1*rRelRequest )
            break;
         };

        /*--------------------------------------------------------------*/
        /*- look on second object if that didn't work ------------------*/
        /*--------------------------------------------------------------*/
        if (npb==Oa->NumPanels)
         { SameObject=0;
           for(npb=0; npb<Ob->NumPanels; npb++)
           { AssessPanelPair(Oa,npa,Ob,npb,&rRel);
             if ( 0.9*rRelRequest<rRel && rRel<1.1*rRelRequest )
              break;
           };
          if (npb==Ob->NumPanels)
           { printf("\n**\n** warning: could not find a panel pair with rRel=%e\n",rRelRequest);
             continue;
           };
         };
      }
     /*--------------------------------------------------------------*/
     /* otherwise choose a random pair of panels                     */
     /*--------------------------------------------------------------*/
     else 
      { if (SameObject==-1) SameObject=lrand48()%2;
        if (npa==-1) npa=lrand48() % Oa->NumPanels;
        if (npb==-1) npb=lrand48() % (SameObject ? Oa->NumPanels : Ob->NumPanels);
        if (iQa==-1) iQa=lrand48() % 3;
        if (iQb==-1) iQb=lrand48() % 3;
      };

     /*--------------------------------------------------------------*/
     /*- if the user specified only the real or imag part of k then -*/
     /*- assume he meant that the other part should be zero; OTOH if-*/
     /*- he didn't specify either part then choose random values for-*/
     /*- both parts                                                  */
     /*--------------------------------------------------------------*/
     if ( isinf(real(K)) && isinf(imag(K)) )
      K=cdouble( drand48(), drand48());
     else if ( isinf(real(K)) )
      real(K)=0.0;
     else if ( isinf(imag(K)) )
      imag(K)=0.0;

     /*--------------------------------------------------------------------*/
     /* print a little summary of the panel pair we will be considering    */
     /*--------------------------------------------------------------------*/
     Pa=Oa->Panels[npa];
     Pb=(SameObject ? Oa:Ob)->Panels[npb];
     ncv=AssessPanelPair(Oa,npa,(SameObject ? Oa : Ob),npb,&rRel,TVa,TVb);
     if (!SameObject) 
      ncv=0;
     printf("*\n");
     printf("* --npa %i --iQa %i (V #%i) --npb %i --iQb %i (V #%i) %s\n",
            npa,iQa,Pa->VI[iQa],npb,iQb,Pb->VI[iQb],SameObject ? "--same" : "--ns");
     printf("*  common vertices:   %i\n",ncv);
     printf("*  relative distance: %+7.3e\n",rRel);
     printf("*  wavevector:        %s\n",CD2S(K));
     printf("*  k*MaxRadius:       %.1e\n",abs(K*fmax(Pa->Radius, Pb->Radius)));
     if (DZ!=0.0)
      printf("*  DZ:                %e\n",DZ);
     printf("*\n\n");

     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     Args->Oa=Oa;
     Args->Ob=(SameObject) ? Oa : Ob;
     Args->npa=npa;
     Args->npb=npb;
     Args->iQa=iQa;
     Args->iQb=iQb;
     Args->k=K;
     Args->NumGradientComponents = Gradient ? 3 : 0;
     Args->NumTorqueAxes         = 0;

     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     if ( !SameObject && DZ!=0.0 )
      Ob->Transform("DISP 0 0 %e",DZ);

     /*--------------------------------------------------------------------*/
     /* get panel-panel integrals by libscuff method                       */
     /*--------------------------------------------------------------------*/
     Tic();
     for(nTimes=0; nTimes<HRTIMES; nTimes++)
      GetPanelPanelInteractions(Args);
     HRTime=Toc() / HRTIMES;
     memcpy(HLS, Args->H, 2*sizeof(cdouble));
     memcpy(GradHLS, Args->GradH, 6*sizeof(cdouble));

     /*--------------------------------------------------------------------*/
     /* get panel-panel integrals by taylor-duffy method if there are any  */
     /* common vertices                                                    */
     /*--------------------------------------------------------------------*/
     if (ncv>0)
      { 
        Qa = Oa->Vertices + 3*Pa->VI[iQa];
        Qb = Args->Ob->Vertices + 3*Pb->VI[iQb];

        TaylorDuffyArgStruct MyTDArgStruct, *TDArgs=&MyTDArgStruct;
        InitTaylorDuffyArgs(TDArgs);

        TDArgs->WhichCase = ncv;
        TDArgs->GParam    = K;
        TDArgs->V1        = TVa[0];
        TDArgs->V2        = TVa[1];
        TDArgs->V3        = TVa[2];
        TDArgs->V2P       = TVb[1];
        TDArgs->V3P       = TVb[2];
        TDArgs->Q         = Qa;
        TDArgs->QP        = Qb;

        Tic();
        for(nTimes=0; nTimes<TDTIMES; nTimes++)
         { 
           TDArgs->WhichG=TM_EIKR_OVER_R;
           TDArgs->WhichH=TM_DOTPLUS;
           HTD[0]=TaylorDuffy(TDArgs);

           TDArgs->WhichG=TM_GRADEIKR_OVER_R;
           TDArgs->WhichH=TM_CROSS;
           TDArgs->AbsTol = RWGGeometry::SWPPITol*abs(HTD[0]);
           HTD[1]=TaylorDuffy(TDArgs);
         };
        TDTime=Toc() / TDTIMES;

      };

     /*--------------------------------------------------------------------*/
     /* get panel-panel integrals by brute-force methods                   */
     /*--------------------------------------------------------------------*/
     GetPPIs_BruteForce(Args, PlotFits);
     memcpy(HBF, Args->H, 2*sizeof(cdouble));
     memcpy(GradHBF, Args->GradH, 6*sizeof(cdouble));

     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     if ( !SameObject && DZ!=0.0 )
      Ob->UnTransform();

     /*--------------------------------------------------------------------*/
     /*- print results ----------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     printf("libscuff computation time:     %.2f us.\n",HRTime*1.0e6); 
     if (ncv>0)
      { printf("Taylor-Duffy computation time: %.2f us.\n",TDTime*1.0e6); 
        PrintResults(HLS, HTD, HBF);
      }
     else
      PrintResults(HLS, GradHLS, HBF, GradHBF,  Gradient);

   }; // end of main command loop [ for(;;) ... ]

}
