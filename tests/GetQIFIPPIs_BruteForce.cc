/*
 * GetQIFIPPIs_BruteForce.cc -- compute Q-independent frequency-independent
 *                           -- panel-panel integrals by brute force
 *
 * homer reid -- 11/2005   -- 12/2011
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

#define ABSTOL 1.0e-12    // absolute tolerance
#define RELTOL 1.0e-8    // relative tolerance

#define NFIPPIS 33

/***************************************************************/
/* data structure used to pass data to integrand routines      */
/***************************************************************/
typedef struct FIPPIBFData
 { 
   double *V0, A[3], B[3];
   double V0P[3], AP[3], BP[3];
 } FIPPIBFData;

/***************************************************************/
/***************************************************************/
/***************************************************************/
static void FIPPIBFIntegrand(unsigned ndim, const double *x, void *parms,
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
  double X[3], XP[3], R[3], XxXP[3];
  double r, r2, oor, oor3;

  memcpy(X,FIPPIBFD->V0,3*sizeof(double));
  VecPlusEquals(X,u,FIPPIBFD->A);
  VecPlusEquals(X,v,FIPPIBFD->B);

  memcpy(XP,FIPPIBFD->V0P,3*sizeof(double));
  VecPlusEquals(XP,up,FIPPIBFD->AP);
  VecPlusEquals(XP,vp,FIPPIBFD->BP);

  VecSub(X, XP, R);
  VecCross(X, XP, XxXP);
  r=VecNorm(R);
  oor=u*up/r;
  oor3=u*up/(r*r*r);
  r2=u*up*r*r;
  r*=u*up;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  fval[0] = R[0]*oor3;
  fval[1] = R[1]*oor3;
  fval[2] = R[2]*oor3;
  fval[3] = XxXP[0]*oor3;
  fval[4] = XxXP[1]*oor3;
  fval[5] = XxXP[2]*oor3;

  fval[6]  = oor;
  fval[7]  = up*oor;
  fval[8]  = vp*oor;
  fval[9]  = u*oor;
  fval[10] = u*up*oor;
  fval[11] = u*vp*oor;
  fval[12] = v*oor;
  fval[13] = v*up*oor;
  fval[14] = v*vp*oor;

  fval[15] = r;
  fval[16] = up*r;
  fval[17] = vp*r;
  fval[18] = u*r;
  fval[19] = u*up*r;
  fval[20] = u*vp*r;
  fval[21] = v*r;
  fval[22] = v*up*r;
  fval[23] = v*vp*r;

  fval[24] = r2;
  fval[25] = up*r2;
  fval[26] = vp*r2;
  fval[27] = u*r2;
  fval[28] = u*up*r2;
  fval[29] = u*vp*r2;
  fval[30] = v*r2;
  fval[31] = v*up*r2;
  fval[32] = v*vp*r2;

} 

/***************************************************************/
/* compute FIPPIs using brute-force technique                  */
/* (adaptive cubature over both panels).                       */
/***************************************************************/
void ComputeQIFIPPIData_BruteForce(double **Va, double **Vb, QIFIPPIData *QIFD)
{ 
  double rRel;

  /***************************************************************/
  /* setup for call to cubature routine    ***********************/
  /***************************************************************/
  FIPPIBFData MyFIPPIBFData, *FIPPIBFD=&MyFIPPIBFData;
 
  FIPPIBFD->V0 = Va[0];
  VecSub(Va[1], Va[0], FIPPIBFD->A);
  VecSub(Va[2], Va[1], FIPPIBFD->B);

  // note that for Vb[0] we make copies of the 
  // entries (not just the pointer) because we may need
  // to displace them, see below.
  memcpy(FIPPIBFD->V0P,Vb[0],3*sizeof(double));
  VecSub(Vb[1], Vb[0], FIPPIBFD->AP);
  VecSub(Vb[2], Vb[1], FIPPIBFD->BP);
   
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
     adapt_integrate(fDim, FIPPIBFIntegrand, (void *)FIPPIBFD, 4, Lower, Upper,
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

        adapt_integrate(fDim, FIPPIBFIntegrand, (void *)FIPPIBFD, 4, Lower, Upper,
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
  QIFD->xMxpRM3[0]   = Result[nf++];
  QIFD->xMxpRM3[1]   = Result[nf++];
  QIFD->xMxpRM3[2]   = Result[nf++];
  QIFD->xXxpRM3[0]   = Result[nf++];
  QIFD->xXxpRM3[1]   = Result[nf++];
  QIFD->xXxpRM3[2]   = Result[nf++];

  QIFD->uvupvpRM1[0] = Result[nf++];
  QIFD->uvupvpRM1[1] = Result[nf++];
  QIFD->uvupvpRM1[2] = Result[nf++];
  QIFD->uvupvpRM1[3] = Result[nf++];
  QIFD->uvupvpRM1[4] = Result[nf++];
  QIFD->uvupvpRM1[5] = Result[nf++];
  QIFD->uvupvpRM1[6] = Result[nf++];
  QIFD->uvupvpRM1[7] = Result[nf++];
  QIFD->uvupvpRM1[8] = Result[nf++];

  QIFD->uvupvpR1[0] = Result[nf++];
  QIFD->uvupvpR1[1] = Result[nf++];
  QIFD->uvupvpR1[2] = Result[nf++];
  QIFD->uvupvpR1[3] = Result[nf++];
  QIFD->uvupvpR1[4] = Result[nf++];
  QIFD->uvupvpR1[5] = Result[nf++];
  QIFD->uvupvpR1[6] = Result[nf++];
  QIFD->uvupvpR1[7] = Result[nf++];
  QIFD->uvupvpR1[8] = Result[nf++];

  QIFD->uvupvpR2[0] = Result[nf++];
  QIFD->uvupvpR2[1] = Result[nf++];
  QIFD->uvupvpR2[2] = Result[nf++];
  QIFD->uvupvpR2[3] = Result[nf++];
  QIFD->uvupvpR2[4] = Result[nf++];
  QIFD->uvupvpR2[5] = Result[nf++];
  QIFD->uvupvpR2[6] = Result[nf++];
  QIFD->uvupvpR2[7] = Result[nf++];
  QIFD->uvupvpR2[8] = Result[nf++];

}
