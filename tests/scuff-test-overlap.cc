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
  (void) ndim;
  (void) fdim;

  GOBFData *GOBFD = (GOBFData *)params;
  double u=x[0], v=u*x[1];

  double *V0   = GOBFD->V0;
  double *V1   = GOBFD->V1;
  double *V2   = GOBFD->V2;
  double *QP   = GOBFD->QP;
  double *ZHat = GOBFD->ZHat;
  double PreFac = GOBFD->PreFac;

  double X[3], fA[3], fB[3]; 

  int Mu;
  for(Mu=0; Mu<3; Mu++)
   { X[Mu]  = V0[Mu] + u*(V1[Mu]-V0[Mu]) + v*(V2[Mu]-V1[Mu]);
     fA[Mu] = X[Mu] - V0[Mu];
     fB[Mu] = X[Mu] - QP[Mu];
   };

  double nxfA[3];
  nxfA[0] = ZHat[1]*fA[2] - ZHat[2]*fA[1];
  nxfA[1] = ZHat[2]*fA[0] - ZHat[0]*fA[2];
  nxfA[2] = ZHat[0]*fA[1] - ZHat[1]*fA[0];

  double nxfB[3];
  nxfB[0] = ZHat[1]*fB[2] - ZHat[2]*fB[1];
  nxfB[1] = ZHat[2]*fB[0] - ZHat[0]*fB[2];
  nxfB[2] = ZHat[0]*fB[1] - ZHat[1]*fB[0];

  double DotProd = (fA[0]*fB[0] + fA[1]*fB[1] + fA[2]*fB[2]);

  fval[0] = u*PreFac*DotProd;
  fval[1] = u*PreFac*(fA[0]*nxfB[0] + fA[1]*nxfB[1] + fA[2]*nxfB[2]);

  fval[2] = u*ZHat[0]*PreFac*DotProd;
  fval[3] = u*ZHat[0]*PreFac*4.0;
  fval[4] = u*nxfA[0]*PreFac*2.0;

  fval[5] = u*ZHat[1]*PreFac*DotProd;
  fval[6] = u*ZHat[1]*PreFac*4.0;
  fval[7] = u*nxfA[1]*PreFac*2.0;

  fval[8]  = u*ZHat[2]*PreFac*DotProd;
  fval[9]  = u*ZHat[2]*PreFac*4.0;
  fval[10] = u*nxfA[2]*PreFac*2.0;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetOverlapBF(RWGObject *O, int nea, int neb, double OValues[11])
{
  RWGEdge *Ea=O->Edges[nea], *Eb=O->Edges[neb];

  double *QP = O->Vertices + 3*Ea->iQP;
  double *V1 = O->Vertices + 3*Ea->iV1;
  double *V2 = O->Vertices + 3*Ea->iV2;
  double *QM = O->Vertices + 3*Ea->iQM;

  double *QPP = O->Vertices + 3*Eb->iQP;
  double *QMP = O->Vertices + 3*Eb->iQM;

  GOBFData MyGOBFD, *GOBFD = &MyGOBFD;
  GOBFD->V1   = V1;
  GOBFD->V2   = V2;

  RWGPanel *PP= O->Panels[Ea->iPPanel];
  double PArea = PP->Area;

  RWGPanel *PM= O->Panels[Ea->iMPanel];
  double MArea = PM->Area;

  double Lower[2]={0.0, 0.0};
  double Upper[2]={1.0, 1.0};
  double I[11], E[11];

  int n;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  memset(OValues,0,11*sizeof(double));
  if ( Ea->iPPanel == Eb->iPPanel )
   { 
     GOBFD->V0     = QP; 
     GOBFD->QP     = QPP;
     GOBFD->ZHat   = PP->ZHat;
     GOBFD->PreFac = (Ea->Length * Eb->Length)/(2.0*PArea);

     adapt_integrate(11, GetOverlapBFIntegrand, (void *)GOBFD, 2, 
		     Lower, Upper, MAXEVAL, ABSTOL, RELTOL, I, E);
     for(n=0; n<11; n++)
      OValues[n] += I[n];
   };

  if ( Ea->iPPanel == Eb->iMPanel )
   { 
     GOBFD->V0     = QP; 
     GOBFD->QP     = QMP;
     GOBFD->ZHat   = PP->ZHat;
     GOBFD->PreFac = (Ea->Length * Eb->Length)/(2.0*PArea);

     adapt_integrate(11, GetOverlapBFIntegrand, (void *)GOBFD, 2, 
		     Lower, Upper, MAXEVAL, ABSTOL, RELTOL, I, E);
     for(n=0; n<11; n++)
      OValues[n] -= I[n];

   };

  if ( Ea->iMPanel == Eb->iPPanel )
   { 
     GOBFD->V0     = QM; 
     GOBFD->QP     = QPP;
     GOBFD->ZHat   = PM->ZHat;
     GOBFD->PreFac = (Ea->Length * Eb->Length)/(2.0*MArea);

     adapt_integrate(11, GetOverlapBFIntegrand, (void *)GOBFD, 2, 
		     Lower, Upper, MAXEVAL, ABSTOL, RELTOL, I, E);
     for(n=0; n<11; n++)
      OValues[n] -= I[n];

   };

  if ( Ea->iMPanel == Eb->iMPanel )
   { 
     GOBFD->V0     = QM; 
     GOBFD->QP     = QMP;
     GOBFD->ZHat   = PM->ZHat;
     GOBFD->PreFac = (Ea->Length * Eb->Length)/(2.0*MArea);

     adapt_integrate(11, GetOverlapBFIntegrand, (void *)GOBFD, 2, 
		     Lower, Upper, MAXEVAL, ABSTOL, RELTOL, I, E);
     for(n=0; n<11; n++)
      OValues[n] += I[n];

   };

}

/***************************************************************/
/* Compute the overlap integral between the RWG basis functions*/
/* associated with two edges in an RWG object. (This is the    */
/* older libRWG version of the current GetOverlap function in  */
/* libscuff.)                                                  */
/***************************************************************/
#if 0
double GetOverlapOld(RWGObject *O, int neAlpha, int neBeta)
{ 
  RWGEdge *EAlpha=O->Edges[neAlpha], *EBeta=O->Edges[neBeta];
  double *V1, *V2, *QAlpha, *QBeta;
  double PreFac, Area, Term, Sum;
  int mu;

  V1=O->Vertices + 3*(EAlpha->iV1);
  V2=O->Vertices + 3*(EAlpha->iV2);
  PreFac=EAlpha->Length * EBeta->Length / (24.0);

  Sum=0.0;
  if ( EAlpha->iPPanel == EBeta->iPPanel )
   {  
      QAlpha=O->Vertices + 3*(EAlpha->iQP);
      QBeta =O->Vertices + 3*(EBeta->iQP);
      Area=O->Panels[EAlpha->iPPanel]->Area;

      for(Term=0.0, mu=0; mu<3; mu++)
       Term+= ( V1[mu] - QAlpha[mu] ) * ( V1[mu] + V2[mu] - 2.0*QBeta[mu] )
             +( V2[mu] - QAlpha[mu] ) * ( V2[mu] + QAlpha[mu] - 2.0*QBeta[mu] );

      Sum += PreFac * Term / Area;
   };
  if ( EAlpha->iPPanel == EBeta->iMPanel )
   {  
      QAlpha=O->Vertices + 3*(EAlpha->iQP);
      QBeta =O->Vertices + 3*(EBeta->iQM);
      Area=O->Panels[EAlpha->iPPanel]->Area;

      for(Term=0.0, mu=0; mu<3; mu++)
       Term+= ( V1[mu] - QAlpha[mu] ) * ( V1[mu] + V2[mu] - 2.0*QBeta[mu] )
             +( V2[mu] - QAlpha[mu] ) * ( V2[mu] + QAlpha[mu] - 2.0*QBeta[mu] );

      Sum -= PreFac * Term / Area;
   };
  if ( EAlpha->iMPanel == EBeta->iPPanel )
   {  
      QAlpha=O->Vertices + 3*(EAlpha->iQM);
      QBeta =O->Vertices + 3*(EBeta->iQP);
      Area=O->Panels[EAlpha->iMPanel]->Area;

      for(Term=0.0, mu=0; mu<3; mu++)
       Term+= ( V1[mu] - QAlpha[mu] ) * ( V1[mu] + V2[mu] - 2.0*QBeta[mu] )
             +( V2[mu] - QAlpha[mu] ) * ( V2[mu] + QAlpha[mu] - 2.0*QBeta[mu] );

      Sum -= PreFac * Term / Area;
   };
  if ( EAlpha->iMPanel == EBeta->iMPanel )
   {  
      QAlpha=O->Vertices + 3*(EAlpha->iQM);
      QBeta =O->Vertices + 3*(EBeta->iQM);
      Area=O->Panels[EAlpha->iMPanel]->Area;

      for(Term=0.0, mu=0; mu<3; mu++)
       Term+= ( V1[mu] - QAlpha[mu] ) * ( V1[mu] + V2[mu] - 2.0*QBeta[mu] )
             +( V2[mu] - QAlpha[mu] ) * ( V2[mu] + QAlpha[mu] - 2.0*QBeta[mu] );

      Sum += PreFac * Term / Area;
   };

  return Sum;

}
#endif

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
  int nea, neb, n;
  double OBF[11], OHR[11];
  srand48(time(0));
  for(;;)
   { 
     /*--------------------------------------------------------------*/
     /*- print prompt and get input string --------------------------*/
     /*--------------------------------------------------------------*/
     printf(" options: --nea xx");
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

     double ReferenceArea = O->Panels[O->Edges[nea]->iPPanel]->Area;

     printf(" ** for nea=%i: \n",nea);

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     for(neb=0; neb<O->NumEdges; neb++)
      {
        GetOverlapBF(O, nea, neb, OBF);
        O->GetOverlaps(nea, neb, OHR);
 //       OOld[0]=O->GetOverlapOld(nea, neb, OOld+1);

        if ( fabs(OBF[0])<1.0e-2*ReferenceArea ) 
         continue;

        printf("**neb=%i: \n",neb);
        for(n=0; n<11; n++)
         printf(" %2i:  %+12.4e | %+12.4e | %+.2e\n",n,OBF[n],OHR[n],RD(OBF[n],OHR[n]));
        printf(" \n");
         
      };

   }; // for(;;)

}
