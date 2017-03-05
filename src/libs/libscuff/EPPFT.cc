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
void GetEPFT(RWGGeometry *G, int DestSurface, cdouble Omega,
             HVector *KNVector, IncField *IF, double FT[NUMFT],
             double **ByEdge, int Order, double Delta)
{
  Log("Computing EPFT for surface %i",DestSurface);

  RWGSurface *S=G->Surfaces[DestSurface];
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
        GetCCDensities(G, DestSurface, np, KNVector, Omega, X, CCDensities);
        XCCMatrix->SetEntries(np*NCP + ncp,"4:end",CCDensities);
      };

    };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HMatrix *FMatrix=G->GetFields(IF, KNVector, Omega, XCCMatrix);

  /***************************************************************/
  /*- initialize edge-by-edge contributions to zero -------------*/
  /***************************************************************/
  if (ByEdge)
   for(int n=PFT_XFORCE; n<=PFT_ZTORQUE; n++)
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
            if (ByEdge[PFT_XFORCE+nq])
             { for(int nce=0; nce<3; nce++)
                { int ne=P->EI[nce];
                  if (ne<0) continue;
                  ByEdge[PFT_XFORCE+nq][ne] += FTThisPanel[nq]/(3.0*Area);
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
/* Get the absorbed and scattered/radiated power for           */
/* DestSurface using the equivalence-principle method.         */
/*                                                             */
/* Either KNVector *or* DRMatrix should be non-null.           */
/*                                                             */
/* If TInterior and TExterior are non-null, they should be the */
/* interior and exterior contributions to the BEM matrix block */
/* corresponding to the self-interactions of DestSurface.      */
/* (TInterior is not referenced if the object is PEC).         */
/*                                                             */
/* On output,                                                  */
/*  Power[0] = absorbed power                                  */
/*  Power[1] = scattered/radiated power                        */
/***************************************************************/
void GetEPP(RWGGeometry *G, int DestSurface, cdouble Omega,
            HVector *KNVector, HMatrix *DRMatrix, double Power[2],
            double **ByEdge, HMatrix *TInterior, HMatrix *TExterior)
{
  // FIXME the use of symmetry yields erroneous results. something
  //       about the sign and/or phase of the (KN-NK)*ikC terms.
  bool UseSymmetry=false;
  char *s=getenv("SCUFF_EPP_SYMMETRY");
  if ( s && s[0]=='1' )
   UseSymmetry=true;

  Log("Computing EPP for surface %i (%s symmetry)",
       DestSurface,UseSymmetry ? "with" : "without");

  RWGSurface *S = G->Surfaces[DestSurface];
  int Offset    = G->BFIndexOffset[DestSurface];
  int NE        = S->NumEdges;
  int nrOut     = S->RegionIndices[0];
  int nrIn      = S->RegionIndices[1];

  /*--------------------------------------------------------------*/
  /*- get wavenumber and relative wave impedance of interior and  */
  /*- exterior media                                              */
  /*--------------------------------------------------------------*/
  cdouble kOut, ZOutRel, ZOutAbs;
  cdouble EpsOut, MuOut;
  G->RegionMPs[nrOut]->GetEpsMu(Omega, &EpsOut, &MuOut);
  kOut = Omega * sqrt(EpsOut*MuOut);
  ZOutRel = sqrt(MuOut/EpsOut);
  ZOutAbs = ZVAC*ZOutRel;

  cdouble kIn, ZInRel, ZInAbs;
  if (S->IsPEC || nrIn<0)
   { 
     kIn=ZInRel=0.0;
     ZInAbs=1.0;
   }
  else
   { cdouble EpsIn, MuIn;
     G->RegionMPs[nrIn]->GetEpsMu(Omega, &EpsIn, &MuIn);
     kIn    = Omega * sqrt(EpsIn*MuIn);
     ZInRel = sqrt(MuIn/EpsIn);
     ZInAbs = ZVAC*ZInRel;
   };

  /*--------------------------------------------------------------*/
  /*- initialize edge-by-edge contributions to the absorbed and --*/
  /*- scattered power to zero                                   --*/
  /*--------------------------------------------------------------*/
  if (ByEdge)
   { if(ByEdge[PFT_PABS])  memset(ByEdge[PFT_PABS], 0,NE*sizeof(double));
     if(ByEdge[PFT_PSCAT]) memset(ByEdge[PFT_PSCAT],0,NE*sizeof(double));
   };

  /*--------------------------------------------------------------*/
  /*- check if TInterior and TExterior are present and the proper */
  /*- sizes                                                       */
  /*--------------------------------------------------------------*/
  int NBF = S->IsPEC ? NE : 2*NE;
  if ( !(S->IsPEC) && TInterior && (TInterior->NR!=NBF || TInterior->NC!=NBF) )
   { Warn("wrong-size TInterior matrix passed to GetEPP (ignoring)");
     TInterior=0;
   }
  if ( TExterior && (TExterior->NR!=NBF || TExterior->NC!=NBF) )
   { Warn("wrong-size TExterior matrix passed to GetEPP (ignoring)");
     TExterior=0;
   };

  int NumThreads=GetNumThreads();
  if ( TExterior && (TInterior || S->IsPEC) )
   NumThreads=1;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double PAbs=0.0, PScat=0.0;
#ifndef USE_OPENMP
  if (G->LogLevel>=SCUFF_VERBOSE2) Log(" no multithreading...");
#else
  if (G->LogLevel>=SCUFF_VERBOSE2) Log(" using %i OpenMP threads",NumThreads);
#pragma omp parallel for schedule(dynamic,1),      \
                         num_threads(NumThreads),  \
                         reduction(+:PAbs, PScat)
#endif
  for(int neab=0; neab<NE*NE; neab++)
   { 
     int nea=neab / NE;
     int neb=neab % NE;
     if (UseSymmetry && neb<nea) continue;

       /*--------------------------------------------------------------*/
       /*- get the matrix elements                                     */
       /*-  ikG  =  <b_\alpha |  ik*G   | b_\beta>                     */
       /*- mikC  =  <b_\alpha | -ik*C | b_\beta>                       */
       /*--------------------------------------------------------------*/
       cdouble ikGIn, mikCIn;
       if (S->IsPEC)
        { 
          ikGIn=mikCIn=0.0; 
        }
       else if (TInterior)
        { 
           ikGIn=TInterior->GetEntry(2*nea+0, 2*neb+0)/ZInRel;
          mikCIn=TInterior->GetEntry(2*nea+1, 2*neb+0);
        }
       else
        { GetEEIArgStruct MyArgs, *Args=&MyArgs;
          InitGetEEIArgs(Args);
          Args->Sa = Args->Sb = S;
          Args->nea = nea;
          Args->neb = neb;
          Args->k   = kIn;
          GetEdgeEdgeInteractions(Args);
          ikGIn  =  II*kIn*Args->GC[0];
          mikCIn = -II*kIn*Args->GC[1];
        };

       cdouble ikGOut, mikCOut;
       if (TExterior && S->IsPEC )
        {  ikGOut=TExterior->GetEntry(nea, neb)/ZOutRel;
          mikCOut=0.0;
        }
       else if (TExterior && !(S->IsPEC) )
        {  ikGOut=TExterior->GetEntry(2*nea+0, 2*neb+0)/ZOutRel;
          mikCOut=TExterior->GetEntry(2*nea+1, 2*neb+0);
        }
       else
        { GetEEIArgStruct MyArgs, *Args=&MyArgs;
          InitGetEEIArgs(Args);
          Args->Sa = Args->Sb = S;
          Args->nea = nea;
          Args->neb = neb;
          Args->k   = kOut;
          GetEdgeEdgeInteractions(Args);
          ikGOut  =  II*kOut*Args->GC[0];
          mikCOut = -II*kOut*Args->GC[1];
        };

       /*--------------------------------------------------------------*/
       /*- extract the surface-current coefficients either from the KN-*/
       /*- vector or the Rytov matrix                                 -*/
       /*--------------------------------------------------------------*/
       cdouble KK, KN, NK, NN;
       if ( KNVector && S->IsPEC )
        { cdouble kAlpha = KNVector->GetEntry(Offset + nea);
          cdouble kBeta  = KNVector->GetEntry(Offset + neb);
          KK = conj(kAlpha) * kBeta;
          KN = NK = NN  = 0.0;
        }
       else if (KNVector && !(S->IsPEC) )
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
        { KK = DRMatrix->GetEntry(Offset+2*neb+0, Offset+2*nea+0);
          KN = DRMatrix->GetEntry(Offset+2*neb+1, Offset+2*nea+0);
          NK = DRMatrix->GetEntry(Offset+2*neb+0, Offset+2*nea+1);
          NN = DRMatrix->GetEntry(Offset+2*neb+1, Offset+2*nea+1);
        };

       /*--------------------------------------------------------------*/
       /*- get the contributions of this edge pair to the power and    */
       /*- accumulate contributions to full and by-edge sums           */
       /*--------------------------------------------------------------*/
       double Weight = (nea==neb) ? -0.5 : -1.0;
       double dPAbs  = Weight*real( (KK*ZInAbs  + NN/ZInAbs)*real(ikGIn)   + (NK-KN)*imag(mikCIn)  );
       double dPScat = Weight*real( (KK*ZOutAbs + NN/ZOutAbs)*real(ikGOut) + (NK-KN)*imag(mikCOut) );
       if (!UseSymmetry)
        { Weight = -0.5;
          dPAbs  = Weight*real( (KK*ZInAbs  + NN/ZInAbs)*ikGIn   + (NK-KN)*mikCIn  );
          dPScat = Weight*real( (KK*ZOutAbs + NN/ZOutAbs)*ikGOut + (NK-KN)*mikCOut );
        };

       PAbs  += dPAbs;
       PScat += dPScat;
       if ( ByEdge )
        { 
#pragma omp critical
          {
            if (ByEdge[PFT_PABS])  ByEdge[PFT_PABS][nea]  += dPAbs;
            if (ByEdge[PFT_PSCAT]) ByEdge[PFT_PSCAT][nea] += dPScat;
          }
        };

     };  // end of multithreaded loop
    
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  Power[0] = PAbs;
  Power[1] = PScat;

}

} // namespace scuff
