/*
 * GetSphericalMoments.cc  -- libscuff class method for computing spherical
 *                         -- multipole moments of induced charge distributions
 *
 * homer reid              -- 2/2010
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <pthread.h>

#include <libhrutil.h>
#include <libSpherical.h>
#include <libTriInt.h>

#include "libscuff.h"

namespace scuff {

#define II cdouble(0,1)

/*******************************************************************/
/* integrand passed to TriInt to compute the contributions to the  */
/* electric and magnetic multipole moments for a single triangle   */
/*                                                                 */ 
/* the integrand has length 4*(lMax+1)*(lMax+1) doubles            */
/*  because there are (lMax+1)*(lMax+1) multipoles with l<=lMax    */
/*  and there is an electric and magnetic moment for each multipole*/
/*  and each moment is a cdouble.                                  */
/*                                                                 */
/* integrand values:                                               */
/*  f[0],  f[1]        == real , imag a^E_{l=0, m=0}               */
/*  f[2],  f[3]        == real , imag a^E_{l=1, m=-1}              */
/*  f[4],  f[5]        == real , imag a^E_{l=1, m=0}               */
/*  f[6],  f[7]        == real , imag a^E_{l=1, m=+1}              */
/*  f[8],  f[9]        == real , imag a^E_{l=2, m=-2}              */
/*   ...                                                           */
/*  f[2*N-2], f[2*N-1] == real , imag a^E_{l=lMax, m=+lMax}        */
/*  f[2*N],   f[2*N+1] == real , imag a^M_{l=0, m=0}               */
/*  f[2*N+2], f[2*N+3] == real , imag a^M_{l=1, m=-1}              */
/*   ...                                                           */
/*  f[4*N-2], f[4*N-1] == real , imag a^M_{l=lMax, m=+lMax}        */
/*                                                                 */
/* (where N=(lMax+1)*(lMax+1))                                     */
/*******************************************************************/
typedef struct SMMIntegrandData
 { 
   double *Q;          // RWG basis function source/sink vertex 
   double PreFac;      // RWG basis function prefactor

   double *X0;         // origin about which we are computing moments

   double Wavevector;  // \sqrt{Eps*Mu} * frequency
   int RealFreq;

   int lMax;           // compute moments for l=0,1, ..., lMax 

 } SMMIData;

static void SphericalMMIntegrand(double *X, void *parms, double *f)
{ 
  SMMIntegrandData *SMMID=(SMMIntegrandData *)parms;
  double fRWG[3], DivfRWG, fRWGSpherical[3];
  double R[3], r, Theta, Phi;
  double w;
  int l, m, nmm, NMM;
  cdouble aE, aM;
  cdouble *zf=(cdouble *)f;

  double Wavevector=SMMID->Wavevector;
  int lMax=SMMID->lMax;
  NMM=(lMax+1)*(lMax+1); // 'number of multipole moments' 

  cdouble Ylm[NMM], dYlmdTheta[NMM];
  double Rl[lMax+2], dRldr[lMax+1]; /* either j_l(kr) or i_l(kr), l=0,...,lMax */
  double kr, ExpFac;

  cdouble iw, ir, PreFac;
  double w2, r2mrdq, dl;

  /* get the value of the RWG basis function at XP */
  VecSub(X,SMMID->Q,fRWG);
  VecScale(fRWG,SMMID->PreFac);

  /* get its divergence */
  DivfRWG=2.0*SMMID->PreFac;

  /* get the spherical coordinates of the evaluation point */
  VecSub(X, SMMID->X0, R);
  CoordinateC2S(R, &r, &Theta, &Phi);
  kr=Wavevector*r;

  /* get the spherical components of the RWG basis function  */
  /* (after this step, fRWGSpherical[0,1,2] = f_r, f_t, f_p) */ 
  VectorC2S(Theta, Phi, fRWG, fRWGSpherical);

  /* get the values of the spherical harmonics and their theta */
  /* derivatives at the evaluation point                       */
  GetYlmDerivArray(lMax, Theta, Phi, Ylm, dYlmdTheta);

  /* get the values of the radial functions at the evaluation point */
#if 0
FIXME to remove gsl bessel functions
  if (SMMID->RealFreq)
   gsl_sf_bessel_jl_steed_array (lMax+1, kr, Rl);
  else
   { gsl_sf_bessel_il_scaled_array(lMax+2, kr, Rl);
     ExpFac=exp(kr);
     for(l=0; l<=lMax; l++)
      Rl[l]*=ExpFac;
   };
#endif

  /* get derivatives of radial functions using  */
  /* d/dx f_l(x) = (l/x)f_l(x) - f_{l+1}(x)     */
  for(l=0; l<=lMax; l++)
   dRldr[l] = ((double)l)*Rl[l]/kr - Rl[l+1];

  /* compute some auxiliary quantities needed to evaluate    */
  /* the integrands                                          */
  iw=SMMID->RealFreq ? II*(SMMID->Wavevector) : -1.0*(SMMID->Wavevector);
  ir=SMMID->RealFreq ? II*r : 1.0*r;
  w2=(SMMID->Wavevector)*(SMMID->Wavevector);
  r2mrdq=r*r-VecDot(X,SMMID->Q);

  /* now loop over all multipoles (l,m) for l=0 ... lMax and */
  /* m=-l to l.                                              */
  /* nmm ('number of multipole moment') is a running index   */
  /* that numbers the multipoles, from nmm=0 for (l=0,m=0)   */
  /* to nmm=NMM-1 for (l=lMax, m=+lMax).                     */
  for(nmm=0, l=0; l<=SMMID->lMax; l++)
   for(m=-l; m<=l; m++, nmm++)
    { 

      dl=(double)l;
      PreFac= 2.0*w2 / sqrt(dl*(dl+1.0));
      if ( SMMID->RealFreq == IMAG_FREQ )
       PreFac*= -1.0;

      /* set aE and aM equal to the integrands for the (l,m)th */
      /* electric and magnetic multipoles                      */
      aE= PreFac * conj(Ylm[nmm]) * ( ( 0.5*iw*r2mrdq - dl/iw)*Rl[l] - (l>0 ? ir*Rl[l-1] : 0.0 ) );
      aM= /*insert contribution here */ 0.0;
   
      /* insert aE and aM into the correct slots in the output vector */
      zf[nmm]=aE;
      zf[NMM + nmm]=aM;
    };

} 

/***************************************************************/
/* get the electric and magnetic spherical multipole moments   */
/* for a single basis function ('1BF') populated with unit     */
/* strength                                                    */
/*                                                             */
/* inputs:                                                     */
/*  ne: basis function index (0...NumEdges-1)                  */
/*  X0: cartesian components of origin about which to compute  */
/*      moments                                                */
/*  lMax: maximum l-value of spherical multipole moments to    */
/*        compute. (moments are computed for l=0,1,...,lMax    */
/*        and all values of m (i.e. m=-l, ... +l).             */
/*  Wavevector: real or imaginary wavevector                   */
/*  RealFreq:   1 or 0 for real or imaginary frequency         */
/*                                                             */
/* outputs:                                                    */
/*                                                             */
/*  aE, aM: must each point to a buffer with space for at      */
/*          least (lMax+1)*(lMax+1) cdoubles. on output, these */
/*          store the electric and magnetic multipole moments  */
/*          of the (electric) current distribution described   */
/*          by basis function #ne populated with unit strength,*/
/*          packed as follows:                                 */
/*           aE[0] : aE_{l=0, m=0}                             */
/*           aE[1] : aE_{l=1, m=-1}                            */
/*           aE[2] : aE_{l=1, m=0}                             */
/*           aE[3] : aE_{l=1, m=+1}                            */
/*           aE[4] : aE_{l=2, m=-2}                            */
/*           ...                                               */
/*           aE[ (lMax+1)*(lMax+1) - 1] : aE_{l=lmax, m=+lmax} */ 
/*          and similarly for the aM.                          */
/***************************************************************/
void RWGObject::Get1BFSphericalMoments(int ne, double *X0, int lMax,
                                       double Wavevector, int RealFreq, 
                                       cdouble *aE, cdouble *aM)
{
  double *QP, *V1, *V2, *QM;
  double PArea, MArea; RWGEdge *E;
  int mu;
  SMMIntegrandData MySMMIData, *SMMID=&MySMMIData;
  int nmm, NMM=(lMax+1)*(lMax+1);
  cdouble PBuf[2*NMM], MBuf[2*NMM];

  /* get edge vertices */
  E=Edges[ne];
  QP=Vertices + 3*(E->iQP);
  V1=Vertices + 3*(E->iV1);
  V2=Vertices + 3*(E->iV2);
  QM=Vertices + 3*(E->iQM);
  PArea=Panels[E->iPPanel]->Area;
  MArea=Panels[E->iMPanel]->Area;

  /* set up data structure passed to SMMIntegrand */
  SMMID->X0=X0;
  SMMID->lMax=lMax;
  SMMID->Wavevector=Wavevector;
  SMMID->RealFreq=RealFreq;

  /* get contribution of positive panel */
  SMMID->Q=QP;
  SMMID->PreFac = E->Length / (2.0*PArea);  
  TriIntFixed(SphericalMMIntegrand, 4*NMM, (void *)SMMID, QP, V1, V2, 25, (double *)PBuf);

  /* contribution of negative panel */
  SMMID->Q=QM;
  SMMID->PreFac = E->Length / (2.0*MArea);
  TriIntFixed(SphericalMMIntegrand, 4*NMM, (void *)SMMID, V1, V2, QM, 25, (double *)MBuf);

  /* pack results into output arrays */
  for(nmm=0; nmm<NMM; nmm++)
   { aE[nmm] = PBuf[nmm ]       - MBuf[NMM];
     aM[nmm] = PBuf[NMM + nmm ] - MBuf[NMM + nmm];
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct ThreadData
 { 
   int nt, nThread;

   RWGGeometry *G;   
   int WhichObject;
   int lMax;
   double *X0;
   double Wavevector;
   int RealFreq;
   HVector *KN;
   cdouble *aE; 
   cdouble *aM;

 } ThreadData;

void *GetSphericalMoments_Thread(void *data)
{ 
  ThreadData *TD=(ThreadData *)data;

  RWGGeometry *G    = TD->G;
  int WhichObject   = TD->WhichObject;
  int lMax          = TD->lMax;
  double *X0        = TD->X0;
  double Wavevector = TD->Wavevector;
  double RealFreq   = TD->RealFreq;
  HVector *KN       = TD->KN;
  cdouble *aE       = TD->aE;
  cdouble *aM       = TD->aM;

  RWGObject *O;
  int nt, ne, mu, Type, Offset;
  cdouble KAlpha, NAlpha;
  double *DKN;
  cdouble *ZKN;
  int IsPEC;

  int nmm, NMM=(lMax+1)*(lMax+1);
  cdouble aE1BF[NMM], aM1BF[NMM];

  O=G->Objects[WhichObject];
  IsPEC=O->MP->IsPEC();
  Offset=G->BFIndexOffset[O->Index];

  if ( RealFreq )
   ZKN=KN->ZV;
  else
   DKN=KN->DV;

  memset(aE, 0, NMM*sizeof(cdouble));
  memset(aM, 0, NMM*sizeof(cdouble));
  
  /*--------------------------------------------------------------*/
  /*- loop over all basis functions on the object                -*/
  /*--------------------------------------------------------------*/
  for(nt=ne=0; ne<O->NumEdges; ne++)
   { 
     /*--------------------------------------------------------------*/
     /*- thread #nt handles edges #nt,lnt+nThread, nt+2*nThread ...  */
     /*--------------------------------------------------------------*/
     nt++;
     if (nt==TD->nThread) nt=0;
     if (nt!=TD->nt) continue;

     /*--------------------------------------------------------------*/
     /* extract basis function weight(s)                             */
     /*--------------------------------------------------------------*/
     if ( IsPEC && RealFreq==1 )
      { KAlpha=ZKN[ Offset + ne ];
        NAlpha=0.0;
      }
     else if ( IsPEC && RealFreq==0 )
      { KAlpha=DKN[ Offset + ne ];
        NAlpha=0.0;
      }
     else if ( !IsPEC && RealFreq==1 )
      { KAlpha=ZKN[ Offset + 2*ne ];
        NAlpha=ZKN[ Offset + 2*ne + 1 ];
      }
     else if ( !IsPEC && RealFreq==0 )
      { KAlpha=DKN[ Offset + 2*ne ];
        NAlpha=DKN[ Offset + 2*ne + 1 ];
      };

     /*--------------------------------------------------------------*/
     /*- get spherical moments of the basis function in question    -*/
     /*- populated with unit strength                               -*/
     /*--------------------------------------------------------------*/
     O->Get1BFSphericalMoments(ne, X0, lMax, Wavevector, RealFreq, aE1BF, aM1BF);
      
     /*--------------------------------------------------------------*/
     /*- add contributions to overall moments -----------------------*/
     /*--------------------------------------------------------------*/
     for(nmm=0; nmm<NMM; nmm++)
      { aE[nmm] += KAlpha * aE1BF[nmm] - NAlpha * aM1BF[nmm];
        aM[nmm] += KAlpha * aM1BF[nmm] - NAlpha * aE1BF[nmm];
      };

   }; // for(nt=ne=0; ne<O->NumEdges; ne++)

  return 0;
 
}

/***************************************************************/
/* get the electric and magnetic spherical multipole moments   */
/* of the current distribution on a single object.             */
/*                                                             */
/* inputs:                                                     */
/*  WhichObject: index of object for which to compute moments  */
/*               (0...NumObjects - 1)                          */
/*  X0: cartesian components of origin about which to compute  */
/*      moments                                                */
/*  lMax: maximum l-value of spherical multipole moments to    */
/*        compute. (moments are computed for l=0,1,...,lMax    */
/*        and all values of m (i.e. m=-l, ... +l).             */
/*  Wavevector: real or imaginary wavevector                   */
/*  RealFreq:   1 or 0 for real or imaginary frequency         */
/*  KN:  solution vector of linear BEM system                  */
/*  nThread: number of threads to use                          */
/*                                                             */
/* outputs:                                                    */
/*                                                             */
/*  aE, aM: must each point to a buffer with space for at      */
/*          least (lMax+1)*(lMax+1) cdoubles. on output, these */
/*          store the electric and magnetic multipole moments  */
/*          of the source distribution on the given object,    */
/*          packed as follows:                                 */
/*                                                             */
/*           aE[0] : aE_{l=0, m=0}                             */
/*           aE[1] : aE_{l=1, m=-1}                            */
/*           aE[2] : aE_{l=1, m=0}                             */
/*           aE[3] : aE_{l=1, m=+1}                            */
/*           aE[4] : aE_{l=2, m=-2}                            */
/*           ...                                               */
/*           aE[ (lMax+1)*(lMax+1) - 1] : aE_{l=lmax, m=+lmax} */ 
/*                                                             */
/*          and similarly for the aM.                          */
/***************************************************************/
void RWGGeometry::GetSphericalMoments(int WhichObject, double *X0, int lMax,
                                      double Frequency, int RealFreq,
                                      HVector *KN, int nThread,
                                      cdouble *aE, cdouble *aM)
{ 
  int nt;
  ThreadData TDS[nThread], *TD;
  pthread_t Threads[nThread];

  int nmm, NMM=(lMax+1)*(lMax+1);
  cdouble aEPartial[nThread*NMM];
  cdouble aMPartial[nThread*NMM];

  cdouble zEps;
  double Eps, Mu;
  double Wavevector;

  if (nThread<=0)
   ErrExit("GetSphericalMoments called with nThread=%i",nThread);

  /***************************************************************/
  /* sanity check on WhichObject *********************************/
  /***************************************************************/
  if (WhichObject<0 || WhichObject > NumObjects)
   { fprintf(stderr,"invalid object %i specified in GetSphericalMoments\n",WhichObject);
     return;
   };

  /***************************************************************/
  /* compute wavevector in external medium                       */
  /***************************************************************/
  ExteriorMP->GetEpsMu(Frequency, RealFreq, &zEps, &Mu); 
  Eps=real(zEps);
  Wavevector=sqrt(Eps*Mu*Frequency);

  /***************************************************************/
  /* fire off threads                                            */
  /*  note: thread #nt writes its contributions to the moments   */
  /*        into aE[ NMM*nt, NMM*nt + 1, ... , NMM*nt + NMM - 1] */
  /*         and aM[ NMM*nt, NMM*nt + 1, ... , NMM*nt + NMM - 1] */
  /*        and afterward we go through and sum them all up.     */
  /***************************************************************/
  for(nt=0; nt<nThread; nt++)
   { 
     TD=&(TDS[nt]);

     TD->nt=nt;
     TD->nThread=nThread;

     TD->G=this;
     TD->WhichObject=WhichObject;
     TD->lMax=lMax;
     TD->X0=X0;
     TD->Wavevector=Wavevector;
     TD->RealFreq=RealFreq;
     TD->KN=KN;
     TD->aE=aE + nt*NMM;
     TD->aM=aM + nt*NMM;
     
     if (nt+1 == nThread)
       GetSphericalMoments_Thread((void *)TD);
     else
       pthread_create( &(Threads[nt]), 0, GetSphericalMoments_Thread, (void *)TD);
   }

  /* wait for threads to complete */
  for(nt=0; nt<nThread-1; nt++)
   pthread_join(Threads[nt],0);

  /* sum contributions from all threads */
  memset(aE,0,NMM*sizeof(cdouble));
  memset(aM,0,NMM*sizeof(cdouble));
  for(nt=0; nt<nThread; nt++)
   for(nmm=0; nmm<NMM; nmm++)
    { aE[nmm] += aEPartial[nt*NMM + nmm];
      aM[nmm] += aMPartial[nt*NMM + nmm];
    };

}

} // namespace scuff
