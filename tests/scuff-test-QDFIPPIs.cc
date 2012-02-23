/*
 * scuff-test-QDFIPPIs.cc -- a test program for libscuff's routines for
 *                        -- computing Q-dependent frequency-independent 
 *                        -- panel-panel integrals
 * 
 * homer reid             -- 11/2005 -- 1/2012
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

#include <libscuff.h>
#include <libscuffInternals.h>

using namespace scuff;

#define HRTIMES 10
#define BFTIMES 1

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetQDFIPPIData(double **Va, double *Qa, double **Vb, double *Qb, 
                    void *opFIPPIDT, QDFIPPIData *QDFD);

void ComputeQDFIPPIData_BruteForce(double **Va, double *Qa, double **Vb, double *Qb, 
                                   QDFIPPIData *QDFD);

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
  /* enter command loop ******************************************/
  /***************************************************************/
  using_history();
  read_history(0);
  int nt, NumTokens;
  char *Tokens[50];
  char *p;
  int npa, iQa, npb, iQb, ncv, SameObject;
  double rRel, rRelRequest, DZ;
  QDFIPPIData QDFDHRBuffer, *QDFDHR=&QDFDHRBuffer;
  QDFIPPIData QDFDBFBuffer, *QDFDBF=&QDFDBFBuffer;
  double *Va[3], *Qa, *Vb[3], *Qb;
  double HRTime, BFTime;
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
     npa=iQa=npb=iQb=ncv=SameObject=-1;
     DZ=rRelRequest=0.0;
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
     free(p);
  
     /*--------------------------------------------------------------*/
     /* if the user specified a number of common vertices then       */
     /* find a randomly-chosen pair of panels with that number of    */
     /* common vertices                                              */
     /*--------------------------------------------------------------*/
     if ( 0<=ncv && ncv<=3 )
      { SameObject=1;
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
      };
  
     if (iQa<0 || iQa>2) iQa=lrand48() % 3;
     if (iQb<0 || iQb>2) iQb=lrand48() % 3;

     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     /*--------------------------------------------------------------------*/
     if ( !SameObject && DZ!=0.0 )
      Ob->Transform("DISP 0 0 %e",DZ);

     Va[0] = Oa->Vertices + 3*(Oa->Panels[npa]->VI[0]);
     Va[1] = Oa->Vertices + 3*(Oa->Panels[npa]->VI[1]);
     Va[2] = Oa->Vertices + 3*(Oa->Panels[npa]->VI[2]);
     Qa    = Va[iQa];

     Vb[0] = Ob->Vertices + 3*(Ob->Panels[npb]->VI[0]);
     Vb[1] = Ob->Vertices + 3*(Ob->Panels[npb]->VI[1]);
     Vb[2] = Ob->Vertices + 3*(Ob->Panels[npb]->VI[2]);
     Qb    = Vb[iQb];

     double *OVa[3], *OQa, *OVb[3], *OQb;
     int Flipped=CanonicallyOrderVertices(Va, Qa, Vb, Qb, OVa, &OQa, OVb, &OQb);

     /*--------------------------------------------------------------------*/
     /* print a little summary of the panel pair we will be considering    */
     /*--------------------------------------------------------------------*/
     ncv=AssessPanelPair(Oa, npa, Ob, npb, &rRel);
     printf("*\n");
     printf("* --npa %i --iQa %i --npb %i --iQb %i %s\n",
            npa,iQa,npb,iQb,SameObject ? "--same" : "--ns");
     printf("*  common vertices:   %i\n",ncv);
     printf("*  relative distance: %+7.3e\n",rRel);
     printf("*  (Flipped: %s)\n",Flipped ? "yes" : "no");
     if (DZ!=0.0)
      printf("*  DZ:                %e\n",DZ);
     printf("*\n\n");

     /*--------------------------------------------------------------------*/
     /* get FIPPIs by libscuff method                                      */
     /*--------------------------------------------------------------------*/
     Tic();
     for(nt=0; nt<HRTIMES; nt++)
      GetQDFIPPIData(Va, Qa, Vb, Qb, 0, QDFDHR);
     HRTime=Toc() / HRTIMES;
     
     /*--------------------------------------------------------------------*/
     /* get FIPPIs by brute-force method                                   */
     /*--------------------------------------------------------------------*/
     Tic();
     for(nt=0; nt<BFTIMES; nt++)
      ComputeQDFIPPIData_BruteForce(Va, Qa, Vb, Qb, QDFDBF);
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
            "            libscuff             ",
            "          brute force            ",
            "rel delta");
     printf("----------|--%15s--|--%15s--|-%5s\n", 
            "---------------","---------------","-----");
     printf("hTimes/r3 |  %15.8e  |  %15.8e  | %5.2e\n", QDFDHR->hTimesRM3, QDFDBF->hTimesRM3, RD(QDFDHR->hTimesRM3, QDFDBF->hTimesRM3));
     printf("\n");

     printf("hDot/r    |  %15.8e  |  %15.8e  | %5.2e\n", QDFDHR->hDotRM1, QDFDBF->hDotRM1, RD(QDFDHR->hDotRM1, QDFDBF->hDotRM1));
     printf("hNabla/r  |  %15.8e  |  %15.8e  | %5.2e\n", QDFDHR->hNablaRM1, QDFDBF->hNablaRM1, RD(QDFDHR->hNablaRM1, QDFDBF->hNablaRM1));
     printf("hTimes/r  |  %15.8e  |  %15.8e  | %5.2e\n", QDFDHR->hTimesRM1, QDFDBF->hTimesRM1, RD(QDFDHR->hTimesRM1, QDFDBF->hTimesRM1));
     printf("\n");

     printf("hDot      |  %15.8e  |  %15.8e  | %5.2e\n", QDFDHR->hDotR0, QDFDBF->hDotR0, RD(QDFDHR->hDotR0, QDFDBF->hDotR0));
     printf("hNabla    |  %15.8e  |  %15.8e  | %5.2e\n", QDFDHR->hNablaR0, QDFDBF->hNablaR0, RD(QDFDHR->hNablaR0, QDFDBF->hNablaR0));
     printf("hTimes    |  %15.8e  |  %15.8e  | %5.2e\n", QDFDHR->hTimesR0, QDFDBF->hTimesR0, RD(QDFDHR->hTimesR0, QDFDBF->hTimesR0));
     printf("\n");

     printf("hDot*r    |  %15.8e  |  %15.8e  | %5.2e\n", QDFDHR->hDotR1, QDFDBF->hDotR1, RD(QDFDHR->hDotR1, QDFDBF->hDotR1));
     printf("hNabla*r  |  %15.8e  |  %15.8e  | %5.2e\n", QDFDHR->hNablaR1, QDFDBF->hNablaR1, RD(QDFDHR->hNablaR1, QDFDBF->hNablaR1));
     printf("hTimes*r  |  %15.8e  |  %15.8e  | %5.2e\n", QDFDHR->hTimesR1, QDFDBF->hTimesR1, RD(QDFDHR->hTimesR1, QDFDBF->hTimesR1));
     printf("\n");

     printf("hDot*r2   |  %15.8e  |  %15.8e  | %5.2e\n", QDFDHR->hDotR2, QDFDBF->hDotR2, RD(QDFDHR->hDotR2, QDFDBF->hDotR2));
     printf("hNabla*r2 |  %15.8e  |  %15.8e  | %5.2e\n", QDFDHR->hNablaR2, QDFDBF->hNablaR2, RD(QDFDHR->hNablaR2, QDFDBF->hNablaR2));
     printf("\n");

   }; // end of main command loop [ for(;;) ... ]

}
