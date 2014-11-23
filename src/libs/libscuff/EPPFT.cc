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
 * EPPFT.cc    -- 'equivalence-principle power, force, and torque'
 *             -- calculation in scuff-EM
 *
 * homer reid  -- 9/2014
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <libhmat.h>
#include <libhrutil.h>
#include <libTriInt.h>

#include "libscuff.h"
#include "libscuffInternals.h"
#include "TaylorDuffy.h"

#include "cmatheval.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif
#ifdef USE_PTHREAD
#  include <pthread.h>
#endif
#ifdef USE_OPENMP

#endif

#define II cdouble(0.0,1.0)
#define TENTHIRDS (10.0/3.0)
#define NUMPFT 7
#define NUMFT 6

namespace scuff {

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetKNAlpha(RWGGeometry *G, int ns, int ne, HVector *KN,
                cdouble *KAlpha, cdouble *NAlpha)
{
  RWGSurface *S=G->Surfaces[ns];
  int Offset = G->BFIndexOffset[ns];
  if (S->IsPEC)
   { *KAlpha = KN->GetEntry(Offset + ne);
     *NAlpha = 0.0;
   }
  else
   { *KAlpha =       KN->GetEntry(Offset + 2*ne + 0);
     *NAlpha = -ZVAC*KN->GetEntry(Offset + 2*ne + 1);
   };
}

/***************************************************************/
/* Get the surface charge and current densities at a point X,  */
/* where X lies within panel #np on surface #ns.               */
/*                                                             */
/* CCDensities[0]    = SigmaE                                  */
/* CCDensities[1..3] = K                                       */
/* CCDensities[4]    = SigmaM                                  */
/* CCDensities[5..7] = N                                       */
/***************************************************************/
void GetCCDensities(RWGGeometry *G, int ns, int np, HVector *KNVector,
                    cdouble Omega,
                    double X[3], cdouble CCDensities[8])
{
  RWGSurface *S = G->Surfaces[ns];
  RWGPanel *P   = S->Panels[np];

  // loop over panel edges to sum contributions of each edge basis function
  memset(CCDensities, 0, 8*sizeof(cdouble));
  for(int nce=0; nce<3; nce++) // 'number of contributing edges'
   { 
     int ne = P->EI[nce];
     if (ne < 0) continue; // panel edge #nce is an exterior edge

     // get the value of the RWG basis function associated with 
     // panel edge #nce at the observation point
     RWGEdge *E    = S->Edges[ne];
     double *Q     = S->Vertices + 3*(P->VI[nce]);
     double Sign   = ( (np == E->iMPanel) ? -1.0 : 1.0);
     double PreFac = Sign * E->Length / (2.0*P->Area);
     double fRWG[3];
     fRWG[0] = PreFac * (X[0] - Q[0]);
     fRWG[1] = PreFac * (X[1] - Q[1]);
     fRWG[2] = PreFac * (X[2] - Q[2]);
         
     // look up the coefficients of this RWG basis function in the
     // expansions of the electric and magnetic surface currents
     cdouble KAlpha, NAlpha;
     GetKNAlpha(G, ns, ne, KNVector, &KAlpha, &NAlpha);

     // add the contributions of this RWG basis function to
     // the source densities at the panel centroid
     CCDensities[0] += 2.0*KAlpha*PreFac / (II*Omega); // Sigma_E
     CCDensities[1] += KAlpha * fRWG[0];               // K_x
     CCDensities[2] += KAlpha * fRWG[1];               // K_y
     CCDensities[3] += KAlpha * fRWG[2];               // K_z
     CCDensities[4] += 2.0*NAlpha*PreFac / (II*Omega); // Sigma_M
     CCDensities[5] += NAlpha * fRWG[0];               // N_x
     CCDensities[6] += NAlpha * fRWG[1];               // N_y
     CCDensities[7] += NAlpha * fRWG[2];               // N_z

   }; // for(int nce=0; nce<3; nce++) // 'number of contributing edges'

}

/***************************************************************/
/* Get the equivalence-principle force and torque on           */
/* DestSurface using a cubature scheme of degree Order over    */
/* the surface of each panel.                                  */
/***************************************************************/
void RWGGeometry::GetEPFT(int DestSurface, cdouble Omega,
                          HVector *KNVector, IncField *IF, double FT[NUMFT],
                          double **ByEdge, int Order, double Delta)
{
  Log("Computing EPFT for surface %i",DestSurface);

  RWGSurface *S=Surfaces[DestSurface];
  int NE = S->NumEdges;
  int NP = S->NumPanels;

  int NCP; // number of cubature points
  double *TCR=GetTCR(Order, &NCP);

  double XTorque[3]={0.0, 0.0, 0.0}; // torque center
  if (S->OTGT) S->OTGT->Apply(XTorque);
  if (S->GT) S->GT->Apply(XTorque);

  /***************************************************************/
  /* populate a matrix with the coordinates of all cubature      */
  /* points and the surface charge/current densities there.      */
  /*                                                             */
  /* XCCMatrix[0...2] = cubature point x,y,z                     */
  /* XCCMatrix[3    ] = cubature weight                          */
  /* XCCMatrix[4    ] = Sigma_E (electric surface charge)        */
  /* XCCMatrix[5...7] = Kx, Ky, Kz (electric surface current)    */
  /* XCCMatrix[8    ] = Sigma_M (magnetic surface charge)        */
  /* XCCMatrix[9..11] = Nx, Ny, Nz (electric surface current)    */
  /***************************************************************/
  HMatrix *XCCMatrix=new HMatrix(NP*NCP, 12);
  for(int np=0; np<NP; np++)
   { 
     RWGPanel *P = S->Panels[np];
     double Radius = P->Radius;
     double *ZHat  = P->ZHat;
     double  Area  = P->Area;
     double *V1 = S->Vertices + 3*(P->VI[0]);
     double *V2 = S->Vertices + 3*(P->VI[1]);
     double *V3 = S->Vertices + 3*(P->VI[2]);
     double A[3], B[3];
     VecSub(V2, V1, A);
     VecSub(V3, V2, B);

     for(int ncp=0; ncp<NCP; ncp++)
      { 
        double u=TCR[3*ncp+0];
        double v=TCR[3*ncp+1];
        double w=TCR[3*ncp+2] * 2.0*Area;
        u+=v;

        double X[3], XDisplaced[3];
        for(int Mu=0; Mu<3; Mu++)
         { X[Mu]          = V1[Mu] + u*A[Mu] + v*B[Mu];
           XDisplaced[Mu] = X[Mu] + Delta*Radius*ZHat[Mu];
         };

        XCCMatrix->SetEntriesD(np*NCP + ncp,"0:2",XDisplaced);
        XCCMatrix->SetEntry(np*NCP + ncp,3,w);

        cdouble CCDensities[8];
        GetCCDensities(this, DestSurface, np, KNVector, Omega, X, CCDensities);
        XCCMatrix->SetEntries(np*NCP + ncp,"4:end",CCDensities);
      };

    };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  bool UGFVSave=RWGGeometry::UseGetFieldsV2P0;
  RWGGeometry::UseGetFieldsV2P0=true;
  HMatrix *FMatrix=GetFields(IF, KNVector, Omega, XCCMatrix);
  RWGGeometry::UseGetFieldsV2P0=UGFVSave;

  /***************************************************************/
  /*- initialize edge-by-edge contributions to zero -------------*/
  /***************************************************************/
  if (ByEdge)
   for(int n=1; n<NUMPFT; n++)
    if (ByEdge[n]) memset(ByEdge[n],0,NE*sizeof(double));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double FTThisPanel[NUMFT];
  memset(FT, 0, NUMFT*sizeof(double) );
  for(int n=0; n<(NP*NCP); n++)
   {  
     double XmX0[3];
     int np=n/NCP;
     RWGPanel *P=S->Panels[np];
     double Area=P->Area;
     double Radius=P->Radius;
     double *ZHat=P->ZHat;
     XmX0[0] = XCCMatrix->GetEntryD(n,0) - Delta*Radius*ZHat[0] - XTorque[0];
     XmX0[1] = XCCMatrix->GetEntryD(n,1) - Delta*Radius*ZHat[1] - XTorque[1];
     XmX0[2] = XCCMatrix->GetEntryD(n,2) - Delta*Radius*ZHat[2] - XTorque[2];

     double w = XCCMatrix->GetEntryD(n,3);

     cdouble CCDensities[8], EH[6];
     XCCMatrix->GetEntries(n, "4:end" ,CCDensities);

     FMatrix->GetEntries(n, ":" ,EH);

     cdouble SigmaE = CCDensities[0];
     cdouble     *K = CCDensities + 1;
     cdouble SigmaM = CCDensities[4];
     cdouble     *N = CCDensities + 5;
     cdouble     *E = EH + 0;
     cdouble     *H = EH + 3;

     // get the three components of the force per unit area
     double dF[3];
     for(int Mu=0; Mu<3; Mu++)
      { int MP1 = (Mu+1)%3, MP2=(Mu+2)%3;

        cdouble KxH = conj(K[MP1])*H[MP2] - conj(K[MP2])*H[MP1];
        cdouble NxE = conj(N[MP1])*E[MP2] - conj(N[MP2])*E[MP1];

        dF[Mu] = 0.25*w*TENTHIRDS*real(  conj(SigmaE)*E[Mu]
                                        +KxH * ZVAC
                                        +conj(SigmaM)*H[Mu]
                                        -NxE / ZVAC
                                      );
      };

     // add contributions to the total force and torque
     for(int Mu=0; Mu<3; Mu++)
      { int MP1 = (Mu+1)%3, MP2=(Mu+2)%3;
        FT[Mu + 0] += dF[Mu];
        FT[Mu + 3] += XmX0[MP1]*dF[MP2] - XmX0[MP2]*dF[MP1];
      };

     // handle edge-resolved contributions if requested
     if (ByEdge)
      { 
        int ncp=n%NCP;

        if (ncp==0) 
         memset(FTThisPanel, 0, NUMFT*sizeof(double) );

        for(int Mu=0; Mu<3; Mu++)
         { int MP1 = (Mu+1)%3, MP2=(Mu+2)%3;
           FTThisPanel[Mu + 0] += dF[Mu];
           FTThisPanel[Mu + 3] += XmX0[MP1]*dF[MP2] - XmX0[MP2]*dF[MP1];
         };

        if ( ncp == (NCP-1) )
         { for(int nq=0; nq<NUMFT; nq++)
            if (ByEdge[1+nq])
             { for(int nce=0; nce<3; nce++)
                { int ne=P->EI[nce];
                  if (ne<0) continue;
                  ByEdge[1+nq][ne] += FTThisPanel[nq]/(3.0*Area);
                };
             };
         };

      };

   }; // for(int n=0; n<(NP*NCP); n++)

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  delete FMatrix;
  delete XCCMatrix;

}

/***************************************************************/
/* Get the equivalence-principle absorbed power for DestSurface*/
/* using the interior expansion, or the radiated/scattered     */
/* power using the exterior expansion.                         */
/*                                                             */
/* Absorbed=true selects the former option, false the latter.  */
/*                                                             */
/* Either KNVector or SigmaMatrix should be non-null.          */
/*                                                             */
/* If TIn is non-null, it should be the interior portion of the*/
/* BEM matrix for surface DestSurface. (Only used if Absorbed  */
/* ==true).                                                    */
/*                                                             */
/* If TOut is non-null, it should be the exterior portion of   */
/* BEM matrix for surface DestSurface. (Only used if Absorbed  */
/* ==false.)                                                   */
/***************************************************************/
double RWGGeometry::GetEPP(int DestSurface, cdouble Omega,
                           HVector *KNVector, HMatrix *SigmaMatrix,
                           double **ByEdge, bool Absorbed,
                           HMatrix *TIn, HMatrix *TOut)
{
  Log("Computing EPP for surface %i ",DestSurface);

  /*--------------------------------------------------------------*/
  /*- get material parameters of interior or exterior region      */
  /*--------------------------------------------------------------*/
  RWGSurface *S = Surfaces[DestSurface];
  int Offset    = BFIndexOffset[DestSurface];
  int NE        = S->NumEdges;
  int nrOut     = S->RegionIndices[0];
  int nrIn      = S->RegionIndices[1];

  /*--------------------------------------------------------------*/
  /*- initialize edge-by-edge contributions to zero --------------*/
  /*--------------------------------------------------------------*/
  if (ByEdge && ByEdge[0])
   memset(ByEdge[0],0,NE*sizeof(double));

  if (nrIn==-1 || S->IsPEC) 
   return 0.0; // equivalence-principle power vanishes for PEC bodies

  bool Radiated=!Absorbed;

  cdouble k, ZRel;
  if (Absorbed) // use interior expansion
   {
     cdouble EpsIn, MuIn;
     RegionMPs[nrIn]->GetEpsMu(Omega, &EpsIn, &MuIn);
     k = Omega * sqrt(EpsIn*MuIn);
     ZRel = sqrt(MuIn/EpsIn);
   }
  else // use exterior expansion
   { cdouble EpsOut, MuOut;
     RegionMPs[nrOut]->GetEpsMu(Omega, &EpsOut, &MuOut);
     k = Omega * sqrt(EpsOut*MuOut);
     ZRel = sqrt(MuOut/EpsOut);
   };

  cdouble IKZ  = II*k*ZVAC*ZRel;
  cdouble PEE  = +0.5*II*k*ZVAC*ZRel;
  cdouble PEM  = -0.5;
  cdouble PME  = +0.5;
  cdouble PMM  = +0.5*II*k/(ZVAC*ZRel);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  struct GetEEIArgStruct MyArgs, *Args=&MyArgs;
  InitGetEEIArgs(Args);
  Args->Sa = Args->Sb = S;
  Args->k  = k;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int NumThreads=GetNumThreads();
  if ( (Absorbed && TIn) || (Radiated && TOut) )
   { NumThreads=1;
     int NBF = 2*NE;
     if ( TIn && (TIn->NR!=NBF || TIn->NC!=NBF) )
      { Warn("wrong-size TIn matrix passed to GetEPP (ignoring)");
        TIn=0;
      }
     if ( TOut && (TOut->NR!=NBF || TOut->NC!=NBF) )
      { Warn("wrong-size TOut matrix passed to GetEPP (ignoring)");
        TOut=0;
      };
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double P=0.0;
#ifndef USE_OPENMP
  if (LogLevel>=SCUFF_VERBOSE2) Log(" no multithreading...");
#else
  if (LogLevel>=SCUFF_VERBOSE2) Log(" using %i OpenMP threads",NumThreads);
#pragma omp parallel for schedule(dynamic,1),      \
                         collapse(2),              \
                         num_threads(NumThreads),  \
                         reduction(+:P)
#endif
  for(int nea=0; nea<NE; nea++)
    for(int neb=0; neb<NE; neb++)
     { 
       if (neb<nea) continue;

       /*--------------------------------------------------------------*/
       /*- get the matrix elements <b_\alpha | e_\beta> and           -*/
       /*-                         <b_\alpha | h_\beta>               -*/
       /*--------------------------------------------------------------*/
       cdouble be, bh;
       if (Absorbed && TIn)
        { be=TIn->GetEntry(2*nea+0, 2*neb+0)/(II*k*ZRel);
          bh=TIn->GetEntry(2*nea+1, 2*neb+0);
        }
       else if (Radiated && TOut)
        { be=TOut->GetEntry(2*nea+0, 2*neb+0)/(II*k*ZRel);
          bh=TOut->GetEntry(2*nea+1, 2*neb+0);
        }
       else
        { Args->nea = nea;
          Args->neb = neb;
          GetEdgeEdgeInteractions(Args);
          be = Args->GC[0];
          bh = Args->GC[1] * (-1.0*II*k);
        };

       /*--------------------------------------------------------------*/
       /*- extract the surface-current coefficient either from the KN -*/
       /*- vector or the Sigma matrix                                 -*/
       /*--------------------------------------------------------------*/
       cdouble KK, KN, NK, NN;
       if (KNVector)
        { 
          cdouble kAlpha =       KNVector->GetEntry(Offset + 2*nea + 0);
          cdouble nAlpha = -ZVAC*KNVector->GetEntry(Offset + 2*nea + 1);
          cdouble kBeta  =       KNVector->GetEntry(Offset + 2*neb + 0);
          cdouble nBeta  = -ZVAC*KNVector->GetEntry(Offset + 2*neb + 1);

          KK = conj(kAlpha) * kBeta;
          KN = conj(kAlpha) * nBeta;
          NK = conj(nAlpha) * kBeta;
          NN = conj(nAlpha) * nBeta;
        }
       else
        { KK = SigmaMatrix->GetEntry(Offset+2*neb+0, Offset+2*nea+0);
          KN = SigmaMatrix->GetEntry(Offset+2*neb+1, Offset+2*nea+0);
          NK = SigmaMatrix->GetEntry(Offset+2*neb+0, Offset+2*nea+1);
          NN = SigmaMatrix->GetEntry(Offset+2*neb+1, Offset+2*nea+1);
        };

       /*--------------------------------------------------------------*/
       /*- get the contributions of this edge pair to the power and    */
       /*- accumulate contributions to full and by-edge sums           */
       /*--------------------------------------------------------------*/
       double Weight = (nea==neb) ? 1.0 : 2.0;
       double dP = -1.0*Weight*real( (KK*PEE + NN*PMM)*be + (KN*PEM + NK*PME)*bh );

       P += dP;
       if ( ByEdge && ByEdge[0] )
        { 
#pragma omp critical
           ByEdge[0][nea] += dP;
        };

     };  // end of multithreaded loop

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  return P;

}

} // namespace scuff
