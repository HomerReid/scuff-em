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
 * scuff-test-PPIs.cc -- a test program for libscuff's routines for
 *                    -- computing panel-panel integrals
 * 
 * homer reid         -- 11/2005 -- 11/2011
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

#define HAVE_LIBREADLINE
#ifdef HAVE_LIBREADLINE
 #include <readline/readline.h>
 #include <readline/history.h>
#else
 #include "readlineReplacement.h"
#endif

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
  RWGSurface *Sa=G->Surfaces[0];
  RWGSurface *Sb=G->NumSurfaces>1 ? G->Surfaces[1] : Sa;

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
  int SameSurface, Gradient, PlotFits;
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
     npa=npb=iQa=iQb=ncv=SameSurface=-1;
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
       SameSurface=1;
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--ns") )
       SameSurface=0;
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
       
     free(p);
  
     /*--------------------------------------------------------------*/
     /* if the user specified a number of common vertices then       */
     /* find a randomly-chosen pair of panels with that number of    */
     /* common vertices                                              */
     /*--------------------------------------------------------------*/
     if ( 0<=ncv && ncv<=3 )
      { SameSurface=1;
        npa=lrand48() % Sa->NumPanels;
        do
         { npb=lrand48() % Sa->NumPanels;
         } while( NumCommonVertices(Sa,npa,Sa,npb)!=ncv );
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
        npa=lrand48() % Sa->NumPanels;

        /*--------------------------------------------------------------*/
        /*- first look on same object ----------------------------------*/
        /*--------------------------------------------------------------*/
        SameSurface=1;
        for(npb=0; npb<Sa->NumPanels; npb++)
         { AssessPanelPair(Sa,npa,Sa,npb,&rRel);
           if ( 0.9*rRelRequest<rRel && rRel<1.1*rRelRequest )
            break;
         };

        /*--------------------------------------------------------------*/
        /*- look on second object if that didn't work ------------------*/
        /*--------------------------------------------------------------*/
        if (npb==Sa->NumPanels)
         { SameSurface=0;
           for(npb=0; npb<Sb->NumPanels; npb++)
           { AssessPanelPair(Sa,npa,Sb,npb,&rRel);
             if ( 0.9*rRelRequest<rRel && rRel<1.1*rRelRequest )
              break;
           };
          if (npb==Sb->NumPanels)
           { printf("\n**\n** warning: could not find a panel pair with rRel=%e\n",rRelRequest);
             continue;
           };
         };
      }
     /*--------------------------------------------------------------*/
     /* otherwise choose a random pair of panels                     */
     /*--------------------------------------------------------------*/
     else 
      { if (SameSurface==-1) SameSurface=lrand48()%2;
        if (npa==-1) npa=lrand48() % Sa->NumPanels;
        if (npb==-1) npb=lrand48() % (SameSurface ? Sa->NumPanels : Sb->NumPanels);
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
     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     if ( !SameSurface && DZ!=0.0 )
      Sb->Transform("DISP 0 0 %e",DZ);

     /*--------------------------------------------------------------------*/
     /* print a little summary of the panel pair we will be considering    */
     /*--------------------------------------------------------------------*/
     Pa=Sa->Panels[npa];
     Pb=(SameSurface ? Sa:Sb)->Panels[npb];
     ncv=AssessPanelPair(Sa,npa,(SameSurface ? Sa : Sb),npb,&rRel,TVa,TVb);
     printf("*\n");
     printf("* --npa %i --iQa %i (V #%i) --npb %i --iQb %i (V #%i) %s\n",
            npa,iQa,Pa->VI[iQa],npb,iQb,Pb->VI[iQb],SameSurface ? "--same" : "--ns");
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
     Args->Sa=Sa;
     Args->Sb=(SameSurface) ? Sa : Sb;
     Args->npa=npa;
     Args->npb=npb;
     Args->iQa=iQa;
     Args->iQb=iQb;
     Args->k=K;
     Args->NumGradientComponents = Gradient ? 3 : 0;
     Args->NumTorqueAxes         = 0;

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
        Qa = Sa->Vertices + 3*Pa->VI[iQa];
        Qb = Args->Sb->Vertices + 3*Pb->VI[iQb];

        double *OVa[3], *OVb[3];
        int Flipped=CanonicallyOrderVertices(TVa, TVb, ncv, OVa, OVb);
        int PIndex, KIndex;
        cdouble Result, Error;

        TaylorDuffyArgStruct MyTDArgStruct, *TDArgs=&MyTDArgStruct;
        InitTaylorDuffyArgs(TDArgs);

        TDArgs->WhichCase = ncv;
        TDArgs->NumPKs    = 1;
        TDArgs->KParam    = &K;
        TDArgs->PIndex    = &PIndex;
        TDArgs->KIndex    = &KIndex;
        TDArgs->Result    = &Result;
        TDArgs->Error     = &Error; 
        TDArgs->V1        = OVa[0];
        TDArgs->V2        = OVa[1];
        TDArgs->V3        = OVa[2];
        TDArgs->V2P       = OVb[1];
        TDArgs->V3P       = OVb[2];
        TDArgs->Q         = Flipped ? Qb : Qa ;
        TDArgs->QP        = Flipped ? Qa : Qb ;

        Tic();
        for(nTimes=0; nTimes<TDTIMES; nTimes++)
         { 
           PIndex=TM_DOTPLUS;
           KIndex=TM_HELMHOLTZ;
           TaylorDuffy(TDArgs);
           HTD[0]=TDArgs->Result[0];
          
           if (ncv==3)
            HTD[1]=0.0;
           else
            { PIndex=TM_CROSS;
              KIndex=TM_GRADHELMHOLTZ;
              TDArgs->AbsTol = 1.0e-8;
              TaylorDuffy(TDArgs);
              HTD[1]=TDArgs->Result[0];
            };
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
     if ( !SameSurface && DZ!=0.0 )
      Sb->UnTransform();

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
