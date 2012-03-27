/*
 * tOverlap.cc -- test of libscuff's routines for computing
 *                overlap integrals between RWG basis functions 
 * 
 * homer reid  -- 3/2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <readline/readline.h>
#include <readline/history.h>
#include <time.h>

#include <libhrutil.h>
#include <libSGJC.h>
#include <libscuff.h>

#define MAXSTR 1000


/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetOverlapBF(RWGObject *O, int nea, int neb, double OValues[2])
{
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /***************************************************************/
  /* process command-line arguments                              */
  /***************************************************************/
  char *GeoFileName;
  ArgStruct ASArray[]=
   { {"geometry",     PA_STRING, (void *)&GeoFileName,   0,     ".scuffgeo file" };
     {0,0,0,0,0}
   };
  ProcessArguments(argc, argv, ASArray);
  if (GeoFileName==0)
   ASUsage(argv[0],ASArray,"--Geometry option is mandatory");

  RWGGeometry *G=new RWGGeometry(GeoFileName);
  RWGObject *O=G->Objects[0];

  /***************************************************************/
  /* enter command loop ******************************************/
  /***************************************************************/
  using_history();
  read_history(0);
  int nt, NumTokens;
  char *Tokens[50];
  char *p;
  int nea;
  int nneb, nebCount;
  double OValues[2];
  double nebHR[5],  OHR[5],  OTimesHR[5];
  double nebBF[10], OBF[10], OTimesBF[10];
  for(;;)
   { 
     /*--------------------------------------------------------------*/
     /*- print prompt and get input string --------------------------*/
     /*--------------------------------------------------------------*/
     printf(" options: --nea xx \n");
     p=readline("enter options: ");
     if (!p) break;
     add_history(p);
     write_history(0);

     /*--------------------------------------------------------------*/
     /* parse input string                                          -*/
     /*--------------------------------------------------------------*/
     NumTokens=Tokenize(p,Tokens,50);
     nea=-1;
     for(nt=0; nt<NumTokens; nt++)
      if ( !strcasecmp(Tokens[nt],"--nea") )
       sscanf(Tokens[nt+1],"%i",&nea);
     free(p);

     if (nea==-1)
      nea = lrand48() % O->NumEdges 

     /*--------------------------------------------------------------*/
     /*- use libscuff methods to get the nonzero entries in the      */
     /*- neath row of the overlap and crossed-overlap matrices       */
     /*--------------------------------------------------------------*/
     O->GetOverlaps(nea, nebHR, OHR, OTimesHR);

     /*--------------------------------------------------------------*/
     /*- use brute-force methods to get the nonzero entries in the   */
     /*- neath row of the overlap and crossed-overlap matrices       */
     /*--------------------------------------------------------------*/
     for( nebCount=0, neb=0; nebCount<10 && neb<O->NumEdges; neb++ )
      { GetOverlapBF(O, nea, neb, OValues);
        if ( fabs(OValues[0])>0.0 || fabs(OValues[1])>0.0 )
         { OBF[nebCount]=OValues[0];
           OTimesBF[nebCount]=OValues[1];
           nebCount++;
         };
      };

     if (nebCount!=5) 
      printf("whoops! nebCount by BF methods=%i\n",nebCount);

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     HMatrix *MHR=new HMatrix(5, 3);
     for(nneb=0; nneb<5; nneb++)
      { MHR->SetEntry(nneb, 0, (double)nebHR[nneb]);
        MHR->SetEntry(nneb, 1, OHR[nneb]);
        MHR->SetEntry(nneb, 2, OHR[nneb]);
      };
     MHR->Sort(0);

     HMatrix *MBF=new HMatrix(nebCount, 3);
     for(nneb=0; nneb<nebCount; nneb++)
      { MBF->SetEntry(nneb, 0, (double)nebBF[nneb]);
        MBF->SetEntry(nneb, 1, OBF[nneb]);
        MBF->SetEntry(nneb, 2, OBF[nneb]);
      };
     MBF->Sort(0);

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     printf("%4s | %+10.3s | %+10.3s | %4.0s | %+10.3s | %+10.3s | %4.0s\n",
             "neb ","   OHR    ","   OBF    "," RD ","  OTHR    ","  OTBF    "," RD ");
     printf("%4s-|-%+10.3s-|-%+10.3s-|-%4.0s-|-%+10.3s-|-%+10.3s-|-%4.0s\n",
             "----","----------","----------","----","----------","----------","----");

     for(nneb=0; nneb<5; nneb++)
      {
        if ( MBF->GetEntry(nneb,0) != MHR->GetEntry(nneb,0) )
         printf("** Bawonkatage! %i: (%g,%g)\n", nneb, MBF->GetEntry(nneb,0), MHR->GetEntry(nneb,0) );
      
        printf("%4g | %+10.3e | %+10.3e | %4.0e | %+10.3e | %+10.3e | %4.0e\n",
                MHR->GetEntry(nneb,0),
                MHR->GetEntry(nneb,1),
                MBF->GetEntry(nneb,1),
                RD(MHR->GetEntry(nneb,1), MBF->GetEntry(nneb,1)),
                MHR->GetEntry(nneb,2),
                MBF->GetEntry(nneb,2),
                RD(MHR->GetEntry(nneb,2), MBF->GetEntry(nneb,2)));
      };

   };

}
