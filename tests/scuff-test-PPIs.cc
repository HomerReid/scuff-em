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

#include <libhrutil.h>

#include "scuff-test-PPIs.h"

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

  /***************************************************************/
  /* enter command loop ******************************************/
  /***************************************************************/
  using_history();
  read_history(0);
  int nt, NumTokens;
  char *Tokens[50];
  char *p;
  int npa, npb, iQa, iQb, ncv, SameObject, Gradient;
  double rRel, DZ;
  cdouble K;
  cdouble HLS[2], GradHLS[6]; // G,C integrals by libscuff 
  cdouble HBF[2], GradHBF[6]; // G,C integrals by brute force
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
     printf("          --same | --ns \n");
     printf("          --DZ \n");
     printf("          --kr     \n");
     printf("          --ki     \n");
     printf("          --gradient\n");
     p=readline("enter options: ");
     if (!p) break;
     add_history(p);
     write_history(0);

     /*--------------------------------------------------------------*/
     /* parse input string                                          -*/
     /*--------------------------------------------------------------*/
     NumTokens=Tokenize(p,Tokens,50);
     npa=npb=iQa=iQb=ncv=SameObject=-1;
     Gradient=0;
     DZ=0.0;
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
      if ( !strcasecmp(Tokens[nt],"--DZ") )
       sscanf(Tokens[nt+1],"%le",&DZ);
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--ns") )
       SameObject=0;
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--Gradient") )
       Gradient=1;
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
     ncv=AssessPanelPair(Oa, npa, Ob, npb, &rRel);
     printf("*\n");
     printf("* --npa %i --iQa %i (V %i) --npb %i --iQb %i (V%i) %s\n",
            npa,iQa,Oa->Panels[npa]->VI[iQa],
            npb,iQb,(SameObject ? Oa:Ob)->Panels[npb]->VI[iQb],
            SameObject ? "--same" : "--ns");
     printf("*  common vertices:   %i\n",ncv);
     printf("*  relative distance: %+7.3e\n",rRel);
     printf("*  wavevector:        %s\n",CD2S(K));
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
     GetPanelPanelInteractions(Args);
     memcpy(HLS, Args->H, 2*sizeof(cdouble));
     memcpy(GradHLS, Args->GradH, 6*sizeof(cdouble));

     /*--------------------------------------------------------------------*/
     /* get panel-panel integrals by brute-force methods                   */
     /*--------------------------------------------------------------------*/
     GetPPIs_BruteForce(Args);
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

   }; // end of main command loop [ for(;;) ... ]

}
