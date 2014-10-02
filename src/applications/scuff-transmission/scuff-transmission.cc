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
 * scuff-transmission.cc -- a command-line tool to compute transmission
 *                       -- through thin films and metamaterial arrays
 *
 * homer reid            -- 9/2012
 */
#include <stdio.h>
#include <stdlib.h>

#include <libhrutil.h>
#include <libscuff.h>
#include <libIncField.h>

using namespace scuff;

#define II cdouble (0.0, 1.0)

#define MAXFREQ 10
#define MAXCACHE 10    // max number of cache files for preload

#define POLARIZATION_TE 0
#define POLARIZATION_TM 1

/*******************************************************************/
/* this is a 9th-order, 17-point cubature rule for the unit square */
/* with corners {(0,0) (1,0) (1,1) (1,0)}.                         */
/* array entries:                                                  */
/*  x_0, y_0, w_0,                                                 */
/*  x_1, y_1, w_1,                                                 */
/*  ...                                                            */
/*  x_16, y_16, w_16                                               */
/* where (x_n, y_n) and w_n are the nth cubature point and weight. */
/*******************************************************************/
double SCR9[]={
  +5.0000000000000000e-01, +5.0000000000000000e-01, +1.3168724279835392e-01,
  +9.8442498318098881e-01, +8.1534005986583447e-01, +2.2219844542549678e-02,
  +9.8442498318098881e-01, +1.8465994013416559e-01, +2.2219844542549678e-02,
  +1.5575016819011134e-02, +8.1534005986583447e-01, +2.2219844542549678e-02,
  +1.5575016819011134e-02, +1.8465994013416559e-01, +2.2219844542549678e-02,
  +8.7513854998945029e-01, +9.6398082297978482e-01, +2.8024900532399120e-02,
  +8.7513854998945029e-01, +3.6019177020215176e-02, +2.8024900532399120e-02,
  +1.2486145001054971e-01, +9.6398082297978482e-01, +2.8024900532399120e-02,
  +1.2486145001054971e-01, +3.6019177020215176e-02, +2.8024900532399120e-02,
  +7.6186791010721466e-01, +7.2666991056782360e-01, +9.9570609815517519e-02,
  +7.6186791010721466e-01, +2.7333008943217640e-01, +9.9570609815517519e-02,
  +2.3813208989278534e-01, +7.2666991056782360e-01, +9.9570609815517519e-02,
  +2.3813208989278534e-01, +2.7333008943217640e-01, +9.9570609815517519e-02,
  +5.3810416409630857e-01, +9.2630786466683113e-01, +6.7262834409945196e-02,
  +5.3810416409630857e-01, +7.3692135333168873e-02, +6.7262834409945196e-02,
  +4.6189583590369143e-01, +9.2630786466683113e-01, +6.7262834409945196e-02,
  +4.6189583590369143e-01, +7.3692135333168873e-02, +6.7262834409945196e-02 
};

/*******************************************************************/
/* get the transmitted and reflected flux by integrating the       */
/* scattered poynting vector over the area of the unit cell.       */
/*                                                                 */
/* more specifically, the transmitted power is the integral of     */
/* the upward-directed poynting vector at ZAbove, while the        */
/* reflected power is the integral of the downward-directed        */
/* poynting vector at ZBelow; the flux is the power divided by     */
/* the area of the unit cell.                                      */
/*                                                                 */
/* Return values: TRFlux[0,1] = transmitted, reflected flux        */
/*******************************************************************/
void GetTRFlux(RWGGeometry *G, IncField *IF, HVector *KN, cdouble Omega, 
               int NQPoints, double *kBloch, 
               double ZAbove, double ZBelow, double *TRFlux)
{
  double *SCR=SCR9;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NCP; // number of cubature points
  if (NQPoints==0)
   NCP=17; 
  else
   NCP=NQPoints*NQPoints;

  /***************************************************************/
  /* on the first invocation we allocate space for the matrices  */
  /* of evaluation points and fields.                            */
  /***************************************************************/
  static HMatrix *XMatrixAbove = 0, *XMatrixBelow = 0;
  static HMatrix *FMatrixAbove = 0, *FMatrixBelow = 0;
  if (XMatrixAbove==0)
   { XMatrixAbove = new HMatrix(NCP, 3 ); 
     XMatrixBelow = new HMatrix(NCP, 3 ); 
     FMatrixAbove = new HMatrix(NCP, 6, LHM_COMPLEX);
     FMatrixBelow = new HMatrix(NCP, 6, LHM_COMPLEX);
   };

  /***************************************************************/ 
  /* fill in coordinates of evaluation points.                   */ 
  /* the first NCP points are for the upper surface; the next    */ 
  /* NCP points are for the lower surface.                       */ 
  /***************************************************************/ 
  double x, y, *LBV[2];
  LBV[0]=G->LBasis[0];
  LBV[1]=G->LBasis[1];
  if (NQPoints==0)
   { for(int ncp=0; ncp<NCP; ncp++)
      { 
        x=SCR[3*ncp+0];
        y=SCR[3*ncp+1];

        XMatrixAbove->SetEntry(ncp, 0, x*LBV[0][0] + y*LBV[1][0]);
        XMatrixAbove->SetEntry(ncp, 1, x*LBV[0][1] + y*LBV[1][1]);
        XMatrixAbove->SetEntry(ncp, 2, ZAbove);

        XMatrixBelow->SetEntry(ncp, 0, x*LBV[0][0] + y*LBV[1][0]);
        XMatrixBelow->SetEntry(ncp, 1, x*LBV[0][1] + y*LBV[1][1]);
        XMatrixBelow->SetEntry(ncp, 2, ZBelow);

      };
   }
  else
   { double Delta = 1.0 / ( (double)NQPoints );
     for(int nqpx=0, ncp=0; nqpx<NQPoints; nqpx++)
      for(int nqpy=0; nqpy<NQPoints; nqpy++, ncp++)
       { 
         x = ((double)nqpx + 0.5)*Delta;
         y = ((double)nqpy + 0.5)*Delta;

         XMatrixAbove->SetEntry(ncp, 0, x*LBV[0][0] + y*LBV[1][0]);
         XMatrixAbove->SetEntry(ncp, 1, x*LBV[0][1] + y*LBV[1][1]);
         XMatrixAbove->SetEntry(ncp, 2, ZAbove);

         XMatrixBelow->SetEntry(ncp, 0, x*LBV[0][0] + y*LBV[1][0]);
         XMatrixBelow->SetEntry(ncp, 1, x*LBV[0][1] + y*LBV[1][1]);
         XMatrixBelow->SetEntry(ncp, 2, ZBelow);
       };
   };

  /***************************************************************/ 
  /* get scattered fields at all cubature points                 */ 
  /***************************************************************/ 
  G->GetFields(IF, KN, Omega, kBloch, XMatrixAbove, FMatrixAbove);
  G->GetFields(0, KN, Omega, kBloch, XMatrixBelow, FMatrixBelow);

  /***************************************************************/
  /* integrate poynting vector over upper and lower surfaces.    */
  /* Note: The jacobian in this cubature is the area of the unit */
  /*       cell, so omitting that factor is equivalent to        */
  /*       dividing the integrated power by the unit-cell area,  */
  /*       which is what we want to do anyway.                   */
  /***************************************************************/
  double w, PTransmitted=0.0, PReflected=0.0;
  cdouble E[3], H[3];
  for(int ncp=0; ncp<NCP; ncp++)
   {
     if (NQPoints==0) 
      w=SCR[3*ncp+2];  // cubature weight
     else
      w=1.0/((double)(NCP));

     E[0]=FMatrixAbove->GetEntry(ncp, 0);
     E[1]=FMatrixAbove->GetEntry(ncp, 1);
     H[0]=FMatrixAbove->GetEntry(ncp, 3);
     H[1]=FMatrixAbove->GetEntry(ncp, 4);
     PTransmitted += 0.5*w*real( E[0]*conj(H[1]) - E[1]*conj(H[0]) );

     E[0]=FMatrixBelow->GetEntry(ncp, 0);
     E[1]=FMatrixBelow->GetEntry(ncp, 1);
     H[0]=FMatrixBelow->GetEntry(ncp, 3);
     H[1]=FMatrixBelow->GetEntry(ncp, 4);
     PReflected -= 0.5*w*real( E[0]*conj(H[1]) - E[1]*conj(H[0]) );

   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  TRFlux[0]=PTransmitted;
  TRFlux[1]=PReflected;

}

/***************************************************************/
/* f1 = \int_0^1 \int_0^u u*e^{-i(uX+vY)} dv du                */
/* f2 = \int_0^1 \int_0^u v*e^{-i(uX+vY)} dv du                */
/***************************************************************/
void f1f2(double X, double Y, cdouble *f1, cdouble *f2)
{
  cdouble ExpIX=exp(II*X);
  cdouble ExpIY=exp(II*Y);
  double X2=X*X; 
  double X3=X2*X;
  double Y2=Y*Y;
  double Y3=Y2*Y;
  double XPY2=(X+Y)*(X+Y);

  if (X==0.0 && Y==0.0)
   { *f1=1.0/3.0;
     *f2=1.0/6.0;
   }
  else if (X==0.0 && Y!=0.0)
   { 
     *f1 = ((II/2.0)*(-2.0 + (2.0 + (2.0*II)*Y)/ExpIY - Y2))/Y3;
     *f2 = -(2.0*II + Y + (-2.0*II + Y)/ExpIY)/Y3;
   }
  else if (Y==0.0 && X!=0.0)
   { *f1= (2.0*II + (-2.0*II + (2.0 + II*X)*X)/ExpIX)/X3;
     *f2= (2.0*II + (-2.0*II + (2.0 + II*X)*X)/ExpIX)/(2.0*X3);
   }
  else if ( fabs(XPY2)<1.0e-20 )
   { *f1 = ((II/2.0)*(-2.0 + (2.0 + (2.0*II)*X)/ExpIX - X2))/X3;
     *f2 = ((-II/2.0)*(-2.0 + 2.0/ExpIX + X*(2.0*II + X)))/X3;
   }
  else
   {
     *f1=((-II)*(-1.0 + (1.0 + II*X)/ExpIX))/(X2*Y) + 
         (II*(-1.0 + (II*(-II + X + Y))/ ExpIX*ExpIY))/(Y*XPY2);

     *f2 = (-(X*(X*(-II + Y) + Y*(-2.0*II + Y))) - II*ExpIY*(-(ExpIX*Y2) + (XPY2)))/
            (ExpIX*ExpIY*X*Y2*XPY2);
   };

}

/***************************************************************/
/* compute the vector-valued integral                          */
/*  \int Exp[-i*(K \cdot X)] b[X] dX                           */
/***************************************************************/
void GetEMiKXRWGIntegral(RWGSurface *S, int ne, double K[3], cdouble Integral[3])
{
  RWGEdge *E    = S->Edges[ne];
  double *QP    = S->Vertices + 3*E->iQP;
  double *V1    = S->Vertices + 3*E->iV1;
  double *V2    = S->Vertices + 3*E->iV2;
  double *QM    = S->Vertices + 3*E->iQM;
  double Length = E->Length;

  double AP[3], AM[3], B[3]; 
  double KQP=0.0, KAP=0.0, KQM=0.0, KAM=0.0, KB=0.0;
  for(int Mu=0; Mu<3; Mu++)
   { AP[Mu] = V1[Mu] - QP[Mu];
     AM[Mu] = V1[Mu] - QM[Mu];
      B[Mu] = V2[Mu] - V1[Mu];
       KQP += K[Mu]*QP[Mu];
       KAP += K[Mu]*AP[Mu];
       KQM += K[Mu]*QM[Mu];
       KAM += K[Mu]*AM[Mu];
        KB += K[Mu]*B[Mu];
   };

  cdouble ExpFac, f1, f2;

  f1f2(KAP, KB, &f1, &f2);
  ExpFac = exp(-II*KQP);
  Integral[0] = Length*ExpFac*( f1*AP[0] + f2*B[0] );
  Integral[1] = Length*ExpFac*( f1*AP[0] + f2*B[0] );
  Integral[2] = Length*ExpFac*( f1*AP[0] + f2*B[0] );

  f1f2(KAM, KB, &f1, &f2);
  ExpFac = exp(-II*KQM);
  Integral[0] -= Length*ExpFac*( f1*AM[0] + f2*B[0] );
  Integral[1] -= Length*ExpFac*( f1*AM[0] + f2*B[0] );
  Integral[2] -= Length*ExpFac*( f1*AM[0] + f2*B[0] );
  
}

/***************************************************************/
/* This routine computes the contributions of currents on a    */
/* single surface to the transmission amplitude.               */
/***************************************************************/
void GetTransmissionAmplitudes(RWGGeometry *G, HVector *KN,
                               int WhichSurface, cdouble EpsPrime,
                               cdouble Omega, double Theta,
                               cdouble *ptTE, cdouble *ptTM)
{
  double nn = real(sqrt(EpsPrime));
  cdouble ZPrime = 1.0/nn;
  double SinThetaPrime = sin(Theta)/nn;
  double CosThetaPrime = sqrt(1.0-SinThetaPrime*SinThetaPrime);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double K[3];
  K[0] = nn*real(Omega)*SinThetaPrime;
  K[1] = 0.0;
  K[2] = nn*real(Omega)*CosThetaPrime;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double EpsTE[3], EpsBarTE[3], EpsTM[3], EpsBarTM[3];
  EpsTE[0]=0.0;     EpsBarTE[0] = -CosThetaPrime;
  EpsTE[1]=1.0;     EpsBarTE[1] = 0.0;
  EpsTE[2]=0.0;     EpsBarTE[2] = +SinThetaPrime;

  EpsTM[0]=+CosThetaPrime;  EpsBarTM[0]=0.0;
  EpsTM[1]=0.0;             EpsBarTM[1]=1.0;
  EpsTM[2]=-SinThetaPrime;  EpsBarTM[2]=0.0;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGSurface *S = G->Surfaces[WhichSurface];
  int BFIndexOffset = G->BFIndexOffset[WhichSurface];

  cdouble tTE=0.0, tTM=0.0;
  cdouble KAlpha, NAlpha=0.0;
  for(int ne=0; ne<S->NumEdges; ne++)
   { 
     if (S->IsPEC)
      { 
        KAlpha=KN->GetEntry( BFIndexOffset + ne );
      }
     else
      { KAlpha=KN->GetEntry( BFIndexOffset + 2*ne + 0 );
        NAlpha=-ZVAC*KN->GetEntry( BFIndexOffset + 2*ne + 1 );
      };

     cdouble EMiQXBAlpha[3]; 
     GetEMiKXRWGIntegral(S, ne, K, EMiQXBAlpha);
     
     tTE += KAlpha*(   EpsTE[0]*EMiQXBAlpha[0]
                     + EpsTE[1]*EMiQXBAlpha[1]
                     + EpsTE[2]*EMiQXBAlpha[2] 
                   )
           +NAlpha*(   EpsBarTE[0]*EMiQXBAlpha[0]
                     + EpsBarTE[1]*EMiQXBAlpha[1]
                     + EpsBarTE[2]*EMiQXBAlpha[2] 
                   );
     
     tTM += KAlpha*(   EpsTM[0]*EMiQXBAlpha[0]
                     + EpsTM[1]*EMiQXBAlpha[1]
                     + EpsTM[2]*EMiQXBAlpha[2] 
                   )
           +NAlpha*(   EpsBarTM[0]*EMiQXBAlpha[0]
                     + EpsBarTM[1]*EMiQXBAlpha[1]
                     + EpsBarTM[2]*EMiQXBAlpha[2] 
                   );

   };

  if (G->LDim!=2)
   ErrExit("%s: %i: internal error",__FILE__,__LINE__);
  double *L1 = G->LBasis[0];
  double *L2 = G->LBasis[1];
  double UnitCellVolume = L1[0]*L2[1]-L1[1]*L2[0];

  tTE *= ZVAC*ZPrime / (2.0*UnitCellVolume*CosThetaPrime);
  tTM *= ZVAC*ZPrime / (2.0*UnitCellVolume*CosThetaPrime);

  if (ptTE) *ptTE = tTE;
  if (ptTM) *ptTM = tTM;
}

/***************************************************************/
/* This routine computes the contributions of currents on ALL  */
/* surfaces bounding the uppermost region to the transmission  */
/* amplitude.                                                  */
/***************************************************************/
void GetTransmissionAmplitudes(RWGGeometry *G, HVector *KN,
                               int UppermostRegionIndex,
                               cdouble Omega, double Theta,
                               cdouble *ptTE, cdouble *ptTM)
{
  double EpsPrime 
   = real( G->RegionMPs[UppermostRegionIndex]->GetEps(Omega) );

  cdouble tTE=0.0, tTM=0.0;
  for(int ns=0; ns<G->NumSurfaces; ns++)
   { 
     double Sign;

     if (G->Surfaces[ns]->RegionIndices[0]==UppermostRegionIndex)
      Sign=1.0;
     else if (G->Surfaces[ns]->RegionIndices[1]==UppermostRegionIndex)
      Sign=-1.0;
     else
      continue;

     cdouble tTEPartial, tTMPartial;
     GetTransmissionAmplitudes(G, KN, ns, EpsPrime, Omega, Theta,
                               &tTEPartial, &tTMPartial);

     tTE += Sign*tTEPartial;
     tTM += Sign*tTMPartial;
   };

  if (ptTE) *ptTE = tTE;
  if (ptTM) *ptTM = tTM;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  EnableAllCPUs();
  InstallHRSignalHandler();

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  char *GeoFileName=0;
  cdouble OmegaVals[MAXFREQ];	int nOmegaVals;
  char *OmegaFile=0;
  double Theta=0.0;
  double ThetaMin=0.0;
  double ThetaMax=0.0;
  int ThetaPoints=25;
  double ZAbove=2.0;
  double ZBelow=-1.0;
  int NQPoints=0;
  char *OutFileName=0;
  char *Cache=0;
  char *ReadCache[MAXCACHE];         int nReadCache;
  char *WriteCache=0;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry",    PA_STRING,  1, 1,       (void *)&GeoFileName,  0,       ".scuffgeo file"},
/**/
     {"Omega",       PA_CDOUBLE, 1, MAXFREQ, (void *)OmegaVals,     &nOmegaVals, "(angular) frequency"},
     {"OmegaFile",   PA_STRING,  1, 1,       (void *)&OmegaFile,    0,       "list of (angular) frequencies"},
/**/
     {"Theta",       PA_DOUBLE,  1, 1,       (void *)&Theta,        0,       "incident angle in degrees"},
     {"ThetaMin",    PA_DOUBLE,  1, 1,       (void *)&ThetaMin,     0,       "minimum incident angle in degrees"}, 
     {"ThetaMax",    PA_DOUBLE,  1, 1,       (void *)&ThetaMax,     0,       "maximum incident angle in degrees"},
     {"ThetaPoints", PA_INT,     1, 1,       (void *)&ThetaPoints,  0,       "number of incident angles"},
/**/
     {"ZAbove",      PA_DOUBLE,  1, 1,       (void *)&ZAbove,       0,       "Z-coordinate of upper integration plane"},
     {"ZBelow",      PA_DOUBLE,  1, 1,       (void *)&ZBelow,       0,       "Z-coordinate of lower integration plane"},
/**/
/**/
     {"NQPoints",    PA_INT,     1, 1,       (void *)&NQPoints,     0,       "number of quadrature points per dimension"},
/**/
     {"OutFile",     PA_STRING,  1, 1,       (void *)&OutFileName,  0,       "output file name"},
/**/
     {"Cache",       PA_STRING,  1, 1,       (void *)&Cache,        0,             "read/write cache"},
     {"ReadCache",   PA_STRING,  1, MAXCACHE,(void *)ReadCache,     &nReadCache,   "read cache"},
     {"WriteCache",  PA_STRING,  1, 1,       (void *)&WriteCache,   0,             "write cache"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFileName==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");
  
  /*******************************************************************/
  /*- create the RWGGeometry                                        -*/
  /*******************************************************************/
  SetLogFileName("scuff-transmission.log");
  RWGGeometry::AssignBasisFunctionsToExteriorEdges=false;
  RWGGeometry *G=new RWGGeometry(GeoFileName);
  if (G->LDim!=2)
   ErrExit("%s: geometry must have two-dimensional lattice periodicity",GeoFileName);
  HMatrix *M   = G->AllocateBEMMatrix();
  HVector *RHS = G->AllocateRHSVector();
  HVector *KN  = G->AllocateRHSVector();

  /*******************************************************************/
  /* determine the index of the uppermost region in the geometry     */
  /*******************************************************************/
  double X[3]={0.0, 0.0, 1.0e6};
  int UpperRegionIndex=G->GetRegionIndex(X);
  printf("Identified uppermost region as region # %i (%s).\n",
          UpperRegionIndex, G->RegionLabels[UpperRegionIndex]);

  /*******************************************************************/
  /* process frequency-related options to construct a list of        */
  /* frequencies at which to run calculations                        */
  /*******************************************************************/
  HVector *OmegaVector=0;
  int nFreq, nOV, NumFreqs=0;
  if (OmegaFile)
   { 
     OmegaVector=new HVector(OmegaFile,LHM_TEXT);
     if (OmegaVector->ErrMsg)
      ErrExit(OmegaVector->ErrMsg);
     NumFreqs=OmegaVector->N;
   };

  // now add any individually specified --Omega options
  if (nOmegaVals>0)
   { 
     NumFreqs += nOmegaVals;
     HVector *OmegaVector0=OmegaVector;
     OmegaVector=new HVector(NumFreqs, LHM_COMPLEX);
     nFreq=0;
     if (OmegaVector0)
      { for(nFreq=0; nFreq<OmegaVector0->N; nFreq++)
         OmegaVector->SetEntry(nFreq, OmegaVector0->GetEntry(nFreq));
        delete OmegaVector0;
      };
     for(nOV=0; nOV<nOmegaVals; nOV++)
      OmegaVector->SetEntry(nFreq+nOV, OmegaVals[nOV]);
   };

  if ( !OmegaVector || OmegaVector->N==0)
   OSUsage(argv[0], OSArray, "you must specify at least one frequency");

  /*******************************************************************/
  /* process incident-angle-related options to construct a list of   */
  /* incident angles at which to run calculations.                   */
  /* Note: The --ThetaMin and --ThetaMax arguments are interpreted   */
  /*       in degrees (i.e. they should be numbers between 0 and 90),*/
  /*       but internally the entries of the ThetaVector vector are  */
  /*       in radians (values between 0 and Pi/2).                   */
  /*******************************************************************/
  HVector *ThetaVector;
  if ( ThetaMax!=0.0 )
   { if ( ThetaMax<ThetaMin )
      OSUsage(argv[0], OSArray, "--ThetaMin must not be greater than --ThetaMax");
     if ( ThetaMax==ThetaMin )  
      ThetaPoints=1;
     ThetaVector=LinSpace(ThetaMin*DEG2RAD, ThetaMax*DEG2RAD, ThetaPoints);
   }
  else 
   { ThetaVector=new HVector(1);
     ThetaVector->SetEntry(0,Theta*DEG2RAD);
   };
  if (ThetaPoints==1)
   Log("Calculating at the single incident angle Theta=%e degrees",ThetaVector->GetEntryD(0)*RAD2DEG);
  else 
   Log("Calculating at %i incident angles in range Theta=(%e,%e) degrees",
       ThetaVector->N, ThetaVector->GetEntryD(0)*RAD2DEG,
       ThetaVector->GetEntryD(ThetaVector->N-1)*RAD2DEG);

  /*******************************************************************/
  /* preload the scuff cache with any cache preload files the user   */
  /* may have specified                                              */
  /*******************************************************************/
  if ( Cache!=0 && WriteCache!=0 )
   ErrExit("--cache and --writecache options are mutually exclusive");
  if (Cache) 
   WriteCache=Cache;
  for (int nrc=0; nrc<nReadCache; nrc++)
   PreloadCache( ReadCache[nrc] );
  if (Cache)
   PreloadCache( Cache );

  /*******************************************************************/
  /*- create the incident field                                      */
  /*******************************************************************/
  cdouble E0[3]={1.0, 0.0, 0.0};
  double nHat[3]={0.0, 0.0, 1.0};
  PlaneWave PW(E0, nHat);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  FILE *f;
  if (OutFileName)
   { f=fopen(OutFileName,"w");
     if (!f) ErrExit("could not open file %s",f);
   }
  else
   { char buffer[1000];
     snprintf(buffer,1000,"%s.transmission",GetFileBase(GeoFileName));
     f=CreateUniqueFile(buffer,1,buffer);
     OutFileName=strdup(buffer);
   };
  fprintf(f,"# data file columns: \n");
  fprintf(f,"# 1:     omega \n");
  fprintf(f,"# 2:     theta (incident angle) (theta=0 --> normal incidence)\n");
  fprintf(f,"# 3:     ktransmitted flux / incident flux (perpendicular polarization)\n");
  fprintf(f,"# 4:     reflected flux   / incident flux (perpendicular polarization)\n");
  fprintf(f,"# 5:     transmitted flux / incident flux (parallel polarization)\n");
  fprintf(f,"# 6:     reflected flux   / incident flux (parallel polarization)\n");
  fprintf(f,"# 7,8:   mag2, phase tPerp\n");
  fprintf(f,"# 9,10:  mag2, phase rPerp\n");
  fprintf(f,"# 11,12: mag2, phase tPar\n");
  fprintf(f,"# 13,14: mag2, phase rPar\n");
  fflush(f);

  cdouble EpsExterior, MuExterior, kExterior;

  /*--------------------------------------------------------------*/
  /*- loop over frequencies and incident angles ------------------*/
  /*--------------------------------------------------------------*/
  double kBloch[2];
  double SinTheta, CosTheta;
  cdouble Omega;
  double FluxTE[2], FluxTM[2], IncFlux;
  cdouble tTETE, tTETM, tTMTE, tTMTM;
  for(int nOmega=0; nOmega<OmegaVector->N; nOmega++)
   for(int nTheta=0; nTheta<ThetaVector->N; nTheta++)
    { 
      Omega = OmegaVector->GetEntry(nOmega);
      Theta = ThetaVector->GetEntryD(nTheta);
      SinTheta=sin(Theta);
      CosTheta=cos(Theta);
      Log("Solving the scattering problem at (Omega,Theta)=(%g,%g)",real(Omega),Theta*RAD2DEG);

      // set bloch wavevector and assemble BEM matrix 
      G->RegionMPs[0]->GetEpsMu(Omega, &EpsExterior, &MuExterior);
      kExterior = csqrt2(EpsExterior*MuExterior)*Omega;
      kBloch[0] = real(kExterior)*SinTheta;
      kBloch[1] = 0.0;
      G->AssembleBEMMatrix(Omega, kBloch, M);
      if (WriteCache)
       { StoreCache( WriteCache );
         WriteCache=0;       
       };
      M->LUFactorize();

      // set plane wave direction 
      nHat[0] = SinTheta;
      nHat[1] = 0.0;
      nHat[2] = CosTheta;
      PW.SetnHat(nHat);

      // solve with E-field perpendicular to plane of incidence  (TE)
      E0[0]=0.0;
      E0[1]=1.0;
      E0[2]=0.0;
      PW.SetE0(E0);
      G->AssembleRHSVector(Omega, kBloch, &PW, RHS);
      KN->Copy(RHS);
      M->LUSolve(KN);
      GetTRFlux(G, &PW, KN, Omega, NQPoints, kBloch, ZAbove, ZBelow, FluxTE);
      GetTransmissionAmplitudes(G, KN, UpperRegionIndex, Omega, Theta, 
                                &tTETE, &tTMTE);

      // solve with E-field parallel to plane of incidence (TM)
      E0[0]=CosTheta;
      E0[1]=0.0;
      E0[2]=-SinTheta;
      PW.SetE0(E0);
      G->AssembleRHSVector(Omega, kBloch, &PW, RHS);
      KN->Copy(RHS);
      M->LUSolve(KN);
      GetTRFlux(G, &PW, KN, Omega, NQPoints, kBloch, ZAbove, ZBelow, FluxTM);
      GetTransmissionAmplitudes(G, KN, UpperRegionIndex, Omega, Theta,
                                &tTETM, &tTMTM);
   
      IncFlux = CosTheta/(2.0*ZVAC);

      fprintf(f,"%s %e ", z2s(Omega), Theta*RAD2DEG);
      fprintf(f,"%e %e ", FluxTE[0]/IncFlux, FluxTE[1]/IncFlux);
      fprintf(f,"%e %e ", FluxTM[0]/IncFlux, FluxTM[1]/IncFlux);
      fprintf(f,"%e %e ", norm(tTETE), arg(tTETE));
      fprintf(f,"%e %e ", norm(tTMTE), arg(tTMTE));
      fprintf(f,"%e %e ", norm(tTETM), arg(tTMTM));
      fprintf(f,"%e %e ", norm(tTMTM), arg(tTMTM));
      fprintf(f,"\n");
      fflush(f);

   }; 
  fclose(f);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  printf("Transmission/reflection data written to %s.\n",OutFileName);
  printf("Thank you for your support.\n");

}
