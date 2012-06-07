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
 * scuff-test-FIPPIs.cc -- a test program for libscuff's routines for
 *                      -- computing frequency-independent panel-panel 
 *                      -- integrals
 * 
 * homer reid           -- 11/2005 -- 1/2012
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>

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

using namespace scuff;

#define HRTIMES 10
#define BFTIMES 1

/***************************************************************/
/***************************************************************/
/***************************************************************/
extern int MaxCalls;
extern double DeltaZFraction;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void ComputeQIFIPPIData(double **Va, double **Vb, QIFIPPIData *QIFD);
void ComputeQIFIPPIData_BruteForce(double **Va, double **Vb, QIFIPPIData *QIFD);
int CanonicallyOrderVertices(double **Va, double *Qa,
                             double **Vb, double *Qb,
                             double **OVa, double **OQa,
                             double **OVb, double **OQb);

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
   { {"geometry",   PA_STRING, (void *)&GeoFileName, 0, ".rwggeo file"},
     {"visualize",  PA_BOOL,   (void *)&Visualize,   0, "write visualization files"},
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
  RWGObject *Oa=G->Objects[0], *Ob;

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
  /* enter command loop ******************************************/
  /***************************************************************/
  int nt, NumTokens;
  char *Tokens[50];
  char *p;
  int npa, npb, ncv, SameObject;
  double rRel, rRelRequest, DZ;
  QIFIPPIData QIFDHRBuffer, *QIFDHR=&QIFDHRBuffer;
  QIFIPPIData QIFDBFBuffer, *QIFDBF=&QIFDBFBuffer;
  double *Va[3], *OVa[3], *Vb[3], *OVb[3], *Qa, *Qb;
  double HRTime, BFTime;
  for(;;)
   { 
     /*--------------------------------------------------------------*/
     /*- print prompt and get input string --------------------------*/
     /*--------------------------------------------------------------*/
     printf(" options: --npa xx \n");
     printf("          --npb xx \n");
     printf("          --ncv xx \n");
     printf("          --rRel xx \n");
     printf("          --same | --ns \n");
     printf("          --DZ \n");
     p=readline("enter options: ");
     if (!p) break;
     add_history(p);
     write_history(0);

     /*--------------------------------------------------------------*/
     /* parse input string                                          -*/
     /*--------------------------------------------------------------*/
     NumTokens=Tokenize(p,Tokens,50);
     npa=npb=ncv=SameObject=-1;
     DZ=rRelRequest=0.0;
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--npa") )
       sscanf(Tokens[nt+1],"%i",&npa);
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--npb") )
       sscanf(Tokens[nt+1],"%i",&npb);
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--ncv") )
       sscanf(Tokens[nt+1],"%i",&ncv);
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
      if ( !strcasecmp(Tokens[nt],"--DeltaZFraction") )
       sscanf(Tokens[nt+1],"%le",&DeltaZFraction);
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--MaxCalls") )
       sscanf(Tokens[nt+1],"%i",&MaxCalls);
     free(p);

     /*--------------------------------------------------------------*/
     /* if the user specified a number of common vertices then       */
     /* find a randomly-chosen pair of panels with that number of    */
     /* common vertices                                              */
     /*--------------------------------------------------------------*/
     if ( 0<=ncv && ncv<=3 )
      { SameObject=1;
        Ob=Oa;
        npa=lrand48() % Oa->NumPanels;
        printf("Looking for a panel pair with %i common vertices...\n",ncv);
        do
         { npb=lrand48() % Oa->NumPanels;
         } while( NumCommonVertices(Oa,npa,Oa,npb)!=ncv );
      }
     /*--------------------------------------------------------------*/
     /*- if the user specified a value for rRel then try to find a --*/
     /*- panel pair whose relative distance is within 10% of rRel  --*/
     /*--------------------------------------------------------------*/
     else if ( rRelRequest!=0.0 )
      { 
        npa=lrand48() % Oa->NumPanels;
        printf("Looking for a panel pair with rRel=%g...\n",rRel);

        /*--------------------------------------------------------------*/
        /*- first look on same object ----------------------------------*/
        /*--------------------------------------------------------------*/
        SameObject=1;
        Ob=Oa;
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
           Ob = (G->NumObjects>1 ? G->Objects[1] : Oa );
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
        if (SameObject)
         Ob=Oa;
        else
         Ob = (G->NumObjects>1 ? G->Objects[1] : Oa );
        if (npa==-1) npa=lrand48() % Oa->NumPanels;
        if (npb==-1) npb=lrand48() % (SameObject ? Oa->NumPanels : Ob->NumPanels);
      };

     /*--------------------------------------------------------------*/
     /*- if the user specified specific panel indices then make sure */
     /*- they make sense                                             */
     /*--------------------------------------------------------------*/
     if ( npa >= Oa->NumPanels )
      { printf("whoops! Object %s has only %i panels (you requested panel %i).\n",
                Oa->Label,Oa->NumPanels,npa);
        continue;
      };
     if ( npb >= Ob->NumPanels )
      { printf("whoops! Object %s has only %i panels (you requested panel %i).\n",
                Ob->Label,Ob->NumPanels,npb);
        continue;
      };

     /*--------------------------------------------------------------------*/
     /* print a little summary of the panel pair we will be considering    */
     /*--------------------------------------------------------------------*/
     ncv=AssessPanelPair(Oa, npa, Ob, npb, &rRel, Va, Vb);
     CanonicallyOrderVertices(Va, Vb, ncv, OVa, OVb);
     printf("*\n");
     printf("* --npa %i --npb %i %s\n",
            npa,npb,SameObject ? "--same" : "--ns");
     printf("*  common vertices:   %i\n",ncv);
     printf("*  relative distance: %+7.3e\n",rRel);
     if (DZ!=0.0)
      printf("*  DZ:                %e\n",DZ);
     printf("*\n\n");

     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     if ( !SameObject && DZ!=0.0 )
      Ob->Transform("DISP 0 0 %e",DZ);

     /*--------------------------------------------------------------------*/
     /* get FIPPIs by libscuff method                                      */
     /*--------------------------------------------------------------------*/
     Tic();
     for(nt=0; nt<HRTIMES; nt++)
      ComputeQIFIPPIData(OVa, OVb, ncv, QIFDHR);
     HRTime=Toc() / HRTIMES;
     
     /*--------------------------------------------------------------------*/
     /* get FIPPIs by brute-force method                                   */
     /*--------------------------------------------------------------------*/
     Tic();
     for(nt=0; nt<BFTIMES; nt++)
      ComputeQIFIPPIData_BruteForce(OVa, OVb, QIFDBF);
     BFTime=Toc() / BFTIMES;

     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     if ( !SameObject && DZ!=0.0 )
      Ob->UnTransform();

     /*--------------------------------------------------------------------*/
     /*- print results ----------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     printf(" ** HR time: %e ms \n",HRTime*1.0e3);
     printf(" ** BF time: %e ms \n",BFTime*1.0e3);
     printf("\n");
     printf("Quantity  |  %15s  |  %15s  | %s\n", 
            "   libscuff   ",
            " brute force  ",
            "rel delta");
     printf("----------|--%15s--|--%15s--|-%5s\n", 
            "---------------","---------------","-----");
     printf("RmRPx/r^3 |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->xMxpRM3[0], QIFDBF->xMxpRM3[0], RD(QIFDBF->xMxpRM3[0],QIFDHR->xMxpRM3[0]));
     printf("RmRPy/r^3 |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->xMxpRM3[1], QIFDBF->xMxpRM3[1], RD(QIFDBF->xMxpRM3[1],QIFDHR->xMxpRM3[1]));
     printf("RmRPz/r^3 |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->xMxpRM3[2], QIFDBF->xMxpRM3[2], RD(QIFDBF->xMxpRM3[2],QIFDHR->xMxpRM3[2]));
     printf("RxRPx/r^3 |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->xXxpRM3[0], QIFDBF->xXxpRM3[0], RD(QIFDBF->xXxpRM3[0],QIFDHR->xXxpRM3[0]));
     printf("RxRPy/r^3 |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->xXxpRM3[1], QIFDBF->xXxpRM3[1], RD(QIFDBF->xXxpRM3[1],QIFDHR->xXxpRM3[1]));
     printf("RxRPz/r^3 |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->xXxpRM3[2], QIFDBF->xXxpRM3[2], RD(QIFDBF->xXxpRM3[2],QIFDHR->xXxpRM3[2]));
     printf("\n");

     printf(" 1  / r^1 |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpRM1[0], QIFDBF->uvupvpRM1[0], RD(QIFDHR->uvupvpRM1[0],QIFDBF->uvupvpRM1[0]));
     printf(" up / r^1 |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpRM1[1], QIFDBF->uvupvpRM1[1], RD(QIFDHR->uvupvpRM1[1],QIFDBF->uvupvpRM1[1]));
     printf(" vp / r^1 |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpRM1[2], QIFDBF->uvupvpRM1[2], RD(QIFDHR->uvupvpRM1[2],QIFDBF->uvupvpRM1[2]));
     printf("  u / r^1 |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpRM1[3], QIFDBF->uvupvpRM1[3], RD(QIFDHR->uvupvpRM1[3],QIFDBF->uvupvpRM1[3]));
     printf("uup / r^1 |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpRM1[4], QIFDBF->uvupvpRM1[4], RD(QIFDHR->uvupvpRM1[4],QIFDBF->uvupvpRM1[4]));
     printf("uvp / r^1 |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpRM1[5], QIFDBF->uvupvpRM1[5], RD(QIFDHR->uvupvpRM1[5],QIFDBF->uvupvpRM1[5]));
     printf("  v / r^1 |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpRM1[6], QIFDBF->uvupvpRM1[6], RD(QIFDHR->uvupvpRM1[6],QIFDBF->uvupvpRM1[6]));
     printf("vup / r^1 |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpRM1[7], QIFDBF->uvupvpRM1[7], RD(QIFDHR->uvupvpRM1[7],QIFDBF->uvupvpRM1[7]));
     printf("vvp / r^1 |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpRM1[8], QIFDBF->uvupvpRM1[8], RD(QIFDHR->uvupvpRM1[8],QIFDBF->uvupvpRM1[8]));
     printf("\n");

     printf(" 1  * r   |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpR1[0], QIFDBF->uvupvpR1[0], RD(QIFDHR->uvupvpR1[0],QIFDBF->uvupvpR1[0]));
     printf(" up * r   |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpR1[1], QIFDBF->uvupvpR1[1], RD(QIFDHR->uvupvpR1[1],QIFDBF->uvupvpR1[1]));
     printf(" vp * r   |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpR1[2], QIFDBF->uvupvpR1[2], RD(QIFDHR->uvupvpR1[2],QIFDBF->uvupvpR1[2]));
     printf("  u * r   |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpR1[3], QIFDBF->uvupvpR1[3], RD(QIFDHR->uvupvpR1[3],QIFDBF->uvupvpR1[3]));
     printf("uup * r   |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpR1[4], QIFDBF->uvupvpR1[4], RD(QIFDHR->uvupvpR1[4],QIFDBF->uvupvpR1[4]));
     printf("uvp * r   |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpR1[5], QIFDBF->uvupvpR1[5], RD(QIFDHR->uvupvpR1[5],QIFDBF->uvupvpR1[5]));
     printf("  v * r   |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpR1[6], QIFDBF->uvupvpR1[6], RD(QIFDHR->uvupvpR1[6],QIFDBF->uvupvpR1[6]));
     printf("vup * r   |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpR1[7], QIFDBF->uvupvpR1[7], RD(QIFDHR->uvupvpR1[7],QIFDBF->uvupvpR1[7]));
     printf("vvp * r   |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpR1[8], QIFDBF->uvupvpR1[8], RD(QIFDHR->uvupvpR1[8],QIFDBF->uvupvpR1[8]));
     printf("\n");

     printf(" 1  * r2  |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpR2[0], QIFDBF->uvupvpR2[0], RD(QIFDHR->uvupvpR2[0],QIFDBF->uvupvpR2[0]));
     printf(" up * r2  |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpR2[1], QIFDBF->uvupvpR2[1], RD(QIFDHR->uvupvpR2[1],QIFDBF->uvupvpR2[1]));
     printf(" vp * r2  |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpR2[2], QIFDBF->uvupvpR2[2], RD(QIFDHR->uvupvpR2[2],QIFDBF->uvupvpR2[2]));
     printf("  u * r2  |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpR2[3], QIFDBF->uvupvpR2[3], RD(QIFDHR->uvupvpR2[3],QIFDBF->uvupvpR2[3]));
     printf("uup * r2  |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpR2[4], QIFDBF->uvupvpR2[4], RD(QIFDHR->uvupvpR2[4],QIFDBF->uvupvpR2[4]));
     printf("uvp * r2  |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpR2[5], QIFDBF->uvupvpR2[5], RD(QIFDHR->uvupvpR2[5],QIFDBF->uvupvpR2[5]));
     printf("  v * r2  |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpR2[6], QIFDBF->uvupvpR2[6], RD(QIFDHR->uvupvpR2[6],QIFDBF->uvupvpR2[6]));
     printf("vup * r2  |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpR2[7], QIFDBF->uvupvpR2[7], RD(QIFDHR->uvupvpR2[7],QIFDBF->uvupvpR2[7]));
     printf("vvp * r2  |  %15.8e  |  %15.8e  | %5.2e\n", QIFDHR->uvupvpR2[8], QIFDBF->uvupvpR2[8], RD(QIFDHR->uvupvpR2[8],QIFDBF->uvupvpR2[8]));
     printf("\n");

   }; // end of main command loop [ for(;;) ... ]

}
