/*
 * GetQDFIPPIs_BruteForce.cc -- compute Q-dependent frequency-independent 
 *                           -- panel-panel integrals by brute force
 * 
 * homer reid -- 11/2005 -- 1/2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include <libSGJC.h>
#include <libPolyFit.h>

#include "libscuff.h"
#include "libscuffInternals.h"

using namespace scuff;

#define ABSTOL 1.0e-12    // absolute tolerance
#define RELTOL 1.0e-8    // relative tolerance

#define NFIPPIS 12

/***************************************************************/
/* data structure used to pass data to integrand routines      */
/***************************************************************/
typedef struct FIPPIBFData
 { 
   double *V0, A[3], B[3], *Q;
   double V0P[3], AP[3], BP[3], QP[3];
 } FIPPIBFData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
static void QDFIPPIBFIntegrand(unsigned ndim, const double *x, void *parms,
                               unsigned nfun, double *fval)
{
  FIPPIBFData *FIPPIBFD=(FIPPIBFData *)parms;
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double u=x[0];
  double v=u*x[1];
  double up=x[2];
  double vp=up*x[3];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double X[3], XmQ[3], XP[3], XPmQP[3], R[3], Scratch[3];
  double r, r2;

  memcpy(X,FIPPIBFD->V0,3*sizeof(double));
  VecPlusEquals(X,u,FIPPIBFD->A);
  VecPlusEquals(X,v,FIPPIBFD->B);
  VecSub(X,FIPPIBFD->Q,XmQ);

  memcpy(XP,FIPPIBFD->V0P,3*sizeof(double));
  VecPlusEquals(XP,up,FIPPIBFD->AP);
  VecPlusEquals(XP,vp,FIPPIBFD->BP);
  VecSub(XP,FIPPIBFD->QP,XPmQP);

  VecSub(X, XP, R);
  r2=R[0]*R[0] + R[1]*R[1] + R[2]*R[2];
  r=sqrt(r2);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double hDot   = u*up*VecDot(XmQ, XPmQP);
  double hNabla = u*up*4.0;
  double hTimes = u*up*VecDot( VecCross(XmQ, XPmQP, Scratch), R ); 
  int nf=0;

  fval[nf++] = hTimes/(r*r*r);

  fval[nf++] = hDot/r;
  fval[nf++] = hNabla/r;
  fval[nf++] = hTimes/r;

  fval[nf++] = hDot;
  fval[nf++] = hNabla;
  fval[nf++] = hTimes;

  fval[nf++] = hDot*r;
  fval[nf++] = hNabla*r;
  fval[nf++] = hTimes*r;

  fval[nf++] = hDot*r2;
  fval[nf++] = hNabla*r2;
} 

/***************************************************************/
/* compute QDFIPPIs using brute-force technique                */
/* (adaptive cubature over both panels).                       */
/***************************************************************/
void ComputeQDFIPPIData_BruteForce(double **Va, double *Qa, 
                                   double **Vb, double *Qb, 
                                   QDFIPPIData *QDFD)
{ 
  /***************************************************************/
  /* setup for call to cubature routine    ***********************/
  /***************************************************************/
  FIPPIBFData MyFIPPIBFData, *FIPPIBFD=&MyFIPPIBFData;
 
  FIPPIBFD->V0 = Va[0];
  VecSub(Va[1], Va[0], FIPPIBFD->A);
  VecSub(Va[2], Va[1], FIPPIBFD->B);
  FIPPIBFD->Q = Qa;

  // note that for Vb[0] and Qb we make copies of the 
  // entries (not just the pointers) because we may need
  // to displace them, see below.
  memcpy(FIPPIBFD->V0P,Vb[0],3*sizeof(double));
  VecSub(Vb[1], Vb[0], FIPPIBFD->AP);
  VecSub(Vb[2], Vb[1], FIPPIBFD->BP);
  memcpy(FIPPIBFD->QP,Qb,3*sizeof(double));
   
  double Lower[4]={0.0, 0.0, 0.0, 0.0};
  double Upper[4]={1.0, 1.0, 1.0, 1.0};

  int nf, fDim=NFIPPIS;
  double Result[fDim], Error[fDim];

  /***************************************************************/
  /* switch off based on whether or not there are any common     */
  /* vertices                                                    */
  /***************************************************************/
  int ncv=AssessPanelPair(Va, Vb);
  if (ncv==0)
   {
     /*--------------------------------------------------------------*/
     /* if there are no common vertices then we can just use naive   */
     /* cubature                                                     */
     /*--------------------------------------------------------------*/
     adapt_integrate(fDim, QDFIPPIBFIntegrand, (void *)FIPPIBFD, 4, Lower, Upper,
                     0, ABSTOL, RELTOL, Result, Error);
   }
  else
   {
     /*--------------------------------------------------------------*/
     /* if there are common vertices then we estimate the panel-panel*/
     /* integrals using a limiting process in which we displace the  */
     /* second of the two panels through a distance Z in the         */
     /* direction of the panel normal and try to fit to Z==0         */
     /*--------------------------------------------------------------*/
     int nz, NZ=10;
     double Z[NZ], F[NZ], FValues[NFIPPIS * NZ];
     PolyFit *PF;

     double Radius, DeltaZ, Centroid[3], BPP[3], ZHat[3];
     Centroid[0] = (Vb[0][0] + Vb[1][0] + Vb[2][0]) / 3.0;
     Centroid[1] = (Vb[0][1] + Vb[1][1] + Vb[2][1]) / 3.0;
     Centroid[2] = (Vb[0][2] + Vb[1][2] + Vb[2][2]) / 3.0;
     VecSub(Vb[2], Vb[0], BPP);
     VecCross(FIPPIBFD->AP, BPP, ZHat);
     VecNormalize(ZHat);
     Radius = VecDistance(Centroid, Vb[0]);
     Radius = fmax(Radius, VecDistance(Centroid, Vb[1]));
     Radius = fmax(Radius, VecDistance(Centroid, Vb[2]));
     DeltaZ = 0.01*Radius;

     for(nz=0; nz<NZ; nz++)
      { 
        Z[nz]=((double)(nz+1))*DeltaZ;
        VecScaleAdd(Vb[0], Z[nz], ZHat, FIPPIBFD->V0P);
        printf("BFing at Z=%g...\n",Z[nz]);

        adapt_integrate(fDim, QDFIPPIBFIntegrand, (void *)FIPPIBFD, 4, Lower, Upper,
                        100000, ABSTOL, RELTOL, Result, Error);

        memcpy(FValues + nz*NFIPPIS, Result, NFIPPIS*sizeof(double));
      };
 
     for(nf=0; nf<NFIPPIS; nf++)
      { for(nz=0; nz<NZ; nz++)
         F[nz] = FValues[ nf + nz*NFIPPIS ];
        PF=new PolyFit(Z, F, NZ, 4);
        Result[nf] = PF->f(0.0);
        delete PF;
      };
     
   }; // if (ncv==0 ... else)

  /***************************************************************/
  /* unpack the results ******************************************/
  /***************************************************************/
  nf=0;
  QDFD->hTimesRM3 = Result[nf++];

  QDFD->hDotRM1   = Result[nf++];
  QDFD->hNablaRM1 = Result[nf++];
  QDFD->hTimesRM1 = Result[nf++];

  QDFD->hDotR0    = Result[nf++];
  QDFD->hNablaR0  = Result[nf++];
  QDFD->hTimesR0  = Result[nf++];

  QDFD->hDotR1    = Result[nf++];
  QDFD->hNablaR1  = Result[nf++];
  QDFD->hTimesR1  = Result[nf++];

  QDFD->hDotR2    = Result[nf++];
  QDFD->hNablaR2  = Result[nf++];

}
