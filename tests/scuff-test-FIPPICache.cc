/*
 * scuff-test-FIPPITable.cc -- a test program for libscuff's routines for
 *                          -- computing panel-panel integrals
 * 
 * homer reid               -- 11/2005 -- 11/2011
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

#include "libscuff.h"
#include "libscuffInternals.h"

using namespace scuff;

extern int Found, NotFound;

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

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  RWGGeometry *G = new RWGGeometry(GeoFileName);

  FIPPITable *FT = new FIPPITable();
  FT->DoNotCompute=1;

  /***************************************************************/
  /* first pass to add the entries *******************************/
  /***************************************************************/
  int noa, nob;
  int npa, npb;
  RWGObject *Oa, *Ob;
  RWGPanel *Pa, *Pb;
  double *Va[3], *Vb[3];
  double *OVa[3], *OVb[3];
  double rRel;
  int ncv;
  int nNearby, nTotal;
  double AddTime, LookupTime;
  QIFIPPIData *QIFD;

  nNearby=0; 
  nTotal=0; 
  Tic();
printf("First pass:\n");
Found=NotFound=0;
  for(noa=0; noa<G->NumObjects; noa++)
   for(nob=noa; nob<G->NumObjects; nob++)
    {
      Oa=G->Objects[noa];
      Ob=G->Objects[nob];
      for(npa=0; npa<Oa->NumPanels; npa++)
       for(npb=(noa==nob ? npa : 0); npb<Ob->NumPanels; npb++)
        { 
          nTotal++;

          ncv=AssessPanelPair(Oa,npa,Ob,npb,&rRel,Va,Vb);
          if (rRel<4.0)
           { 
             nNearby++;
             CanonicallyOrderVertices(Va, Vb, ncv, OVa, OVb);
             QIFD=FT->GetQIFIPPIData(OVa, OVb, ncv);
           };
        };
    };
  AddTime=Toc();
printf(" %i found, %i not found\n",Found,NotFound);
  
  /***************************************************************/
  /* second pass to do the lookups *******************************/
  /***************************************************************/
printf("Second pass:\n");
nNearby=0;
Found=NotFound=0;
  Tic();
  for(noa=0; noa<G->NumObjects; noa++)
   for(nob=noa; nob<G->NumObjects; nob++)
    {
      Oa=G->Objects[noa];
      Ob=G->Objects[nob];
      for(npa=0; npa<Oa->NumPanels; npa++)
       for(npb=(noa==nob ? npa : 0); npb<Ob->NumPanels; npb++)
        { 
          ncv=AssessPanelPair(Oa,npa,Ob,npb,&rRel,Va,Vb);
          if (rRel<4.0)
           { 
//if (nNearby==3) break;
             nNearby++;
             CanonicallyOrderVertices(Va, Vb, ncv, OVa, OVb);
             QIFD=FT->GetQIFIPPIData(OVa, OVb, ncv);
           };
        };
    };
  LookupTime=Toc();
printf(" %i found, %i not found\n",Found,NotFound);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  printf("%i/%i nearby panel pairs \n",nNearby,nTotal);
  printf("Add time:    %e us / pair \n",AddTime*1e6/nNearby);
  printf("Lookup time: %e us / pair \n",LookupTime*1e6/nNearby);
}  
