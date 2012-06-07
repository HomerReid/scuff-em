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

#include "scuff-test-EEIs.h"

/***************************************************************/
/* get edge-edge interactions by brute force                  **/
/***************************************************************/
void GetEEIs_BruteForce(GetEEIArgStruct *Args)
{
  GetPPIArgStruct *GetPPIArgs;

  GetPPIs_BruteForce(GetPPIArgs, 0);
  
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
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
  SetDefaultCD2SFormat("(%+10.3e,%+10.3e)");

  /***************************************************************/
  /* create the geometry *****************************************/
  /***************************************************************/
  RWGGeometry *G = new RWGGeometry(GeoFileName);
  RWGObject *Oa=G->Objects[0];
  RWGObject *Ob=G->NumObjects>1 ? G->Objects[1] : Oa;

  /***************************************************************/
  /*- write visualization files if requested ---------------------*/
  /***************************************************************/
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
  /* preinitialize an argument structure for the edge-edge       */
  /* interaction routines                                        */
  /***************************************************************/
  GetEEIArgStruct MyArgs, *Args=&MyArgs;

  /***************************************************************/
  /* enter command loop ******************************************/
  /***************************************************************/
  using_history();
  read_history(0);
  int nt, NumTokens;
  char *Tokens[50];
  char *p;
  int nea, neb, SameObject, Gradient, ncv;
  double rRel, rRelRequest, DZ;
  cdouble K;
  cdouble GCPP[2], GradGCPP[6]; // G,C integrals by panel-panel integration
  cdouble GCSM[2], GradGCSM[6]; // G,C integrals by spherical multipole method 
  cdouble GCBF[2], GradGCBF[6]; // G,C integrals by brute force
  for(;;)
   { 
     /*--------------------------------------------------------------*/
     /*- print prompt and get input string --------------------------*/
     /*--------------------------------------------------------------*/
     printf(" options: --nea xx \n");
     printf("          --neb xx \n");
     printf("          --rRel xx \n");
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
     nea=neb=SameObject=-1;
     Gradient=0;
     DZ=rRelRequest=0.0;
     real(K) = imag(K) = INFINITY;
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--nea") )
       sscanf(Tokens[nt+1],"%i",&nea);
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--neb") )
       sscanf(Tokens[nt+1],"%i",&neb);
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
      if ( !strcasecmp(Tokens[nt],"--rRel") )
       sscanf(Tokens[nt+1],"%le",&rRelRequest);
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--ns") )
       SameObject=0;
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--Gradient") )
       Gradient=1;
     free(p);
  
     /*--------------------------------------------------------------*/
     /*- if the user specified a value for rRel then try to find an -*/
     /*- edge pair whose relative distance is within 10% of rRel    -*/
     /*--------------------------------------------------------------*/
     if ( rRelRequest!=0.0 )
      { 
        nea=lrand48() % Oa->NumEdges;

        /*--------------------------------------------------------------*/
        /*- first look on same object ----------------------------------*/
        /*--------------------------------------------------------------*/
        SameObject=1;
        for(neb=0; neb<Oa->NumEdges; neb++)
         { AssessEdgePair(Oa,nea,Oa,neb,&rRel);
           if ( 0.9*rRelRequest<rRel && rRel<1.1*rRelRequest )
            break;
         };

        /*--------------------------------------------------------------*/
        /*- look on second object if that didn't work ------------------*/
        /*--------------------------------------------------------------*/
        if (neb==Oa->NumEdges)
         { SameObject=0;
           for(neb=0; neb<Ob->NumEdges; neb++)
           { AssessEdgePair(Oa,nea,Ob,neb,&rRel);
             if ( 0.9*rRelRequest<rRel && rRel<1.1*rRelRequest )
              break;
           };
          if (neb==Ob->NumEdges)
           { printf("\n**\n** warning: could not find an edge pair with rRel=%e\n",rRelRequest);
             continue;
           };
         };
      };

     /*--------------------------------------------------------------*/
     /* choose random values for any quantities left unspecified     */
     /*--------------------------------------------------------------*/
     if (SameObject==-1) SameObject=lrand48()%2;
     if (nea==-1) nea=lrand48() % Oa->NumEdges;
     if (neb==-1) neb=lrand48() % (SameObject ? Oa->NumEdges: Ob->NumEdges);

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
     /* print a little summary of the edge pair we will be considering     */
     /*--------------------------------------------------------------------*/
     ncv=AssessEdgePair(Oa, nea, Ob, neb, &rRel);
     printf("*\n");
     printf("* --nea %i --neb %i %s\n",
            nea, neb, SameObject ? "--same" : "--ns");
     printf("*  common vertices:   %i\n",ncv);
     printf("*  relative distance: %+7.3e (DBFThreshold=10.0)\n",rRel);
     printf("*  wavevector:        %s\n",CD2S(K));
     if (DZ!=0.0)
      printf("*  DZ:                %e\n",DZ);
     printf("*\n\n");

     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     Args->Oa=Oa;
     Args->Ob=(SameObject) ? Oa : Ob;
     Args->nea=nea;
     Args->neb=neb;
     Args->k=K;
     Args->NumGradientComponents = Gradient ? 3 : 0;
     Args->NumTorqueAxes         = 0;

     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     if ( !SameObject && DZ!=0.0 )
      Ob->Transform("DISP 0 0 %e",DZ);

     /*--------------------------------------------------------------------*/
     /* get edge-edge interactions by panel-panel integration method       */
     /*--------------------------------------------------------------------*/
     Args->Force=EEI_FORCE_PP;
     GetEdgeEdgeInteractions(Args);
     memcpy(GCPP, Args->GC, 2*sizeof(cdouble));
     memcpy(GradGCPP, Args->GradGC, 6*sizeof(cdouble));

     /*--------------------------------------------------------------------*/
     /* get edge-edge interactions by spherical multipole method           */
     /*--------------------------------------------------------------------*/
     Args->Force=EEI_FORCE_SM;
     GetEdgeEdgeInteractions(Args);
     memcpy(GCSM, Args->GC, 2*sizeof(cdouble));
     memcpy(GradGCSM, Args->GradGC, 6*sizeof(cdouble));

     /*--------------------------------------------------------------------*/
     /* get panel-panel integrals by brute-force methods                   */
     /*--------------------------------------------------------------------*/
     GetEEIs_BruteForce(Args);
     memcpy(GCBF, Args->GC, 2*sizeof(cdouble));
     memcpy(GradGCBF, Args->GradGC, 6*sizeof(cdouble));

     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     if ( !SameObject && DZ!=0.0 )
      Ob->UnTransform();

     /*--------------------------------------------------------------------*/
     /*- print results ----------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     printf("Quantity          | %23s | %23s | %23s |%s|%s\n",
            "     panel-panel       ", "     S multipole       ",
            "     brute force       ", "rdPP","rdSM");
     printf("------------------|-%23s-|-%23s-|-%23s-|%s|%s\n", 
            "-----------------------", "-----------------------",
            "-----------------------", "----","----");
     printf("<fa|G|fb>         | %s | %s | %s |%4.1e|%4.1e\n", CD2S(GCPP[0]), CD2S(GCSM[0]), CD2S(GCBF[0]), RD(GCPP[0],GCBF[0]), RD(GCSM[0],GCBF[0]));
     printf("<fa|C|fb>         | %s | %s | %s |%4.1e|%4.1e\n", CD2S(GCPP[1]), CD2S(GCSM[1]), CD2S(GCBF[1]), RD(GCPP[1],GCBF[1]), RD(GCSM[1],GCBF[1]));
     if (Gradient)
      { 
         printf("\n");
         printf("(d/dx) <fa|G|fb>  | %s | %s | %s |%4.1e|%4.1e\n", CD2S(GradGCPP[0]), CD2S(GradGCSM[0]), CD2S(GradGCBF[0]), RD(GradGCPP[0],GradGCBF[0]), RD(GradGCSM[0],GradGCBF[0]));
         printf("(d/dy) <fa|G|fb>  | %s | %s | %s |%4.1e|%4.1e\n", CD2S(GradGCPP[2]), CD2S(GradGCSM[2]), CD2S(GradGCBF[2]), RD(GradGCPP[2],GradGCBF[2]), RD(GradGCSM[2],GradGCBF[2]));
         printf("(d/dz) <fa|G|fb>  | %s | %s | %s |%4.1e|%4.1e\n", CD2S(GradGCPP[4]), CD2S(GradGCSM[4]), CD2S(GradGCBF[4]), RD(GradGCPP[4],GradGCBF[4]), RD(GradGCSM[4],GradGCBF[4]));
         printf("(d/dx) <fa|C|fb>  | %s | %s | %s |%4.1e|%4.1e\n", CD2S(GradGCPP[1]), CD2S(GradGCSM[1]), CD2S(GradGCBF[1]), RD(GradGCPP[1],GradGCBF[1]), RD(GradGCSM[1],GradGCBF[1]));
         printf("(d/dy) <fa|C|fb>  | %s | %s | %s |%4.1e|%4.1e\n", CD2S(GradGCPP[3]), CD2S(GradGCSM[3]), CD2S(GradGCBF[3]), RD(GradGCPP[3],GradGCBF[3]), RD(GradGCSM[3],GradGCBF[3]));
         printf("(d/dz) <fa|C|fb>  | %s | %s | %s |%4.1e|%4.1e\n", CD2S(GradGCPP[5]), CD2S(GradGCSM[5]), CD2S(GradGCBF[5]), RD(GradGCPP[5],GradGCBF[5]), RD(GradGCSM[5],GradGCBF[5]));
      };
     printf("\n");

   }; // end of main command loop [ for(;;) ... ]

}
