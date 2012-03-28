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

#include "libhrutil.h"
#include "libSGJC.h"
#include "libscuff.h"

#define MAXSTR 1000
#define MAXEVAL 10000
#define ABSTOL 1.0e-8
#define RELTOL 1.0e-4

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct GOBFData 
 {
   double *V0, *V1, *V2, *QP, *ZHat;
   double PreFac;
 } GOBFData;

void GetOverlapBFIntegrand(unsigned ndim, const double *x, void *params,
			   unsigned fdim, double *fval)
{
  GOBFData *GOBFD = (GOBFData *)params;
  double u=x[0], v=u*x[1];

  double *V0   = GOBFD->V0;
  double *V1   = GOBFD->V1;
  double *V2   = GOBFD->V2;
  double *QP   = GOBFD->QP;
  double *ZHat = GOBFD->ZHat;
  double PreFac = GOBFD->PreFac;

  double X[3], fA[3], fB[3]; 

  int i;
  for(i=0; i<3; i++)
   { X[i]  = V0[i] + u*(V1[i]-V0[i]) + v*(V2[i]-V1[i]);
     fA[i] = X[i] - V0[i];
     fB[i] = X[i] - QP[i];
   };

  double nxfB[3];
  nxfB[0] = ZHat[1]*fB[2] - ZHat[2]*fB[1];
  nxfB[1] = ZHat[2]*fB[0] - ZHat[0]*fB[2];
  nxfB[2] = ZHat[0]*fB[1] - ZHat[1]*fB[0];

  fval[0] = u*PreFac*(fA[0]*fB[0]   + fA[1]*fB[1]   + fA[2]*fB[2]);
  fval[1] = u*PreFac*(fA[0]*nxfB[0] + fA[1]*nxfB[1] + fA[2]*nxfB[2]);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetOverlapBF(RWGObject *O, int nea, int neb, double OValues[2])
{
  RWGEdge *Ea=O->Edges[nea], *Eb=O->Edges[neb];

  double *QP = O->Vertices + 3*Ea->iQP;
  double *V1 = O->Vertices + 3*Ea->iV1;
  double *V2 = O->Vertices + 3*Ea->iV2;
  double *QM = O->Vertices + 3*Ea->iQM;

  double *QPP = O->Vertices + 3*Eb->iQP;
  double *QMP = O->Vertices + 3*Eb->iQM;

  OValues[0]=OValues[1]=0.0;

  GOBFData MyGOBFD, *GOBFD = &MyGOBFD;
  GOBFD->V1   = V1;
  GOBFD->V2   = V2;

  RWGPanel *PP= O->Panels[Ea->iPPanel];
  double PArea = PP->Area;

  RWGPanel *PM= O->Panels[Ea->iMPanel];
  double MArea = PM->Area;

  double Lower[2]={0.0, 0.0};
  double Upper[2]={1.0, 1.0};
  double I[2], E[2];

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if ( Ea->iPPanel == Eb->iPPanel )
   { 
     GOBFD->V0     = QP; 
     GOBFD->QP     = QPP;
     GOBFD->ZHat   = PP->ZHat;
     GOBFD->PreFac = (Ea->Length * Eb->Length)/(2*PArea);

     adapt_integrate(2, GetOverlapBFIntegrand, (void *)GOBFD, 2, 
		     Lower, Upper, MAXEVAL, ABSTOL, RELTOL, I, E);

     OValues[0] += I[0];
     OValues[1] += I[1];
   };

  if ( Ea->iPPanel == Eb->iMPanel )
   { 
     GOBFD->V0     = QP; 
     GOBFD->QP     = QMP;
     GOBFD->ZHat   = PP->ZHat;
     GOBFD->PreFac = (Ea->Length * Eb->Length)/(2*PArea);

     adapt_integrate(2, GetOverlapBFIntegrand, (void *)GOBFD, 2, 
		     Lower, Upper, MAXEVAL, ABSTOL, RELTOL, I, E);

     OValues[0] -= I[0];
     OValues[1] -= I[1];
   };

  if ( Ea->iMPanel == Eb->iPPanel )
   { 
     GOBFD->V0     = QM; 
     GOBFD->QP     = QPP;
     GOBFD->ZHat   = PM->ZHat;
     GOBFD->PreFac = (Ea->Length * Eb->Length)/(2*MArea);

     adapt_integrate(2, GetOverlapBFIntegrand, (void *)GOBFD, 2, 
		     Lower, Upper, MAXEVAL, ABSTOL, RELTOL, I, E);

     OValues[0] -= I[0];
     OValues[1] -= I[1];
   };

  if ( Ea->iMPanel == Eb->iMPanel )
   { 
     GOBFD->V0     = QM; 
     GOBFD->QP     = QMP;
     GOBFD->ZHat   = PM->ZHat;
     GOBFD->PreFac = (Ea->Length * Eb->Length)/(2*MArea);

     adapt_integrate(2, GetOverlapBFIntegrand, (void *)GOBFD, 2, 
		     Lower, Upper, MAXEVAL, ABSTOL, RELTOL, I, E);

     OValues[0] += I[0];
     OValues[1] += I[1];
   };

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
   { {"geometry",     PA_STRING, (void *)&GeoFileName,   0,     ".scuffgeo file" },
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
  int nea, neb;
  int nnz, NNZ;
  double OValues[2], OHR, OTimesHR;
  HMatrix *MM = new HMatrix(10, 6);
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
      nea = lrand48() % O->NumEdges;

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     for(NNZ=neb=0; NNZ<10 && neb<O->NumEdges; neb++)
      {
        GetOverlapBF(O, nea, neb, OValues);
        OHR=O->GetOverlap(nea, neb, &OTimesHR);
        if (    fabs(OValues[0])>0.0 || fabs(OValues[1])>0.0 
             || fabs(OHR)>0.0 || fabs(OTimesHR)>0.0 )
         { 
           MM->SetEntry(NNZ, 0, (double)neb );
           MM->SetEntry(NNZ, 1, OHR);
           MM->SetEntry(NNZ, 2, OValues[0]);
           MM->SetEntry(NNZ, 3, OTimesHR);
           MM->SetEntry(NNZ, 4, OValues[1]);

           MM->SetEntry(NNZ, 5, O->GetOverlapOld(nea, neb) );

           NNZ++;
         };
      };

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     printf("\n\n");
     printf("%4s | %10s | %10s | %4s | %10s | %10s | %4s\n",
             "neb ","   OHR    ","   OBF    "," RD  ","  OTHR    ","  OTBF    "," RD  ");
     printf("%4s-|-%10s-|-%10s-|-%4s-|-%10s-|-%10s-|-%4s\n",
             "----","----------","----------","-----","----------","----------","-----");

     for(nnz=0; nnz<NNZ; nnz++)
      printf("%4g | %+10.3e | %+10.3e | %5.0e | %+10.3e | %+10.3e | %5.0e\n",
              MM->GetEntryD(nnz,0),
              MM->GetEntryD(nnz,1),
              MM->GetEntryD(nnz,2),
              RD(MM->GetEntryD(nnz,1), MM->GetEntryD(nnz,2)),
              MM->GetEntryD(nnz,3),
              MM->GetEntryD(nnz,4),
              RD(MM->GetEntryD(nnz,3), MM->GetEntryD(nnz,4)));

     printf("\n\n");
     printf("%4s | %10s | %10s | %5s\n",
             "neb ","   OHR    ","  OHROLD  ","  RD ");
     printf("%4s-|-%10s-|-%10s-|-%5s\n",
             "----","----------","----------","-----");

     for(nnz=0; nnz<NNZ; nnz++)
      printf("%4g | %+10.3e | %+10.3e | %5.0e\n",
              MM->GetEntryD(nnz,0), 
              MM->GetEntryD(nnz,1), 
              MM->GetEntryD(nnz,5),
              RD( MM->GetEntryD(nnz,1), MM->GetEntryD(nnz,5) )
            );

   }; // for(;;)

}
