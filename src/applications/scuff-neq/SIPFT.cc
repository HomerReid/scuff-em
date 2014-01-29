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

/***************************************************************/
/* SIPFT.cc -- SCUFF-EM code for computing                     */
/*          -- surface-integral power, force, and torque       */
/*                                                             */
/* There are two                                               */
/*                                                             */
/*                                                             */
/* Homer Reid -- 1/2014                                        */
/***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#include <libhrutil.h>
#include <libhmat.h>
#include <libscuff.h>
#include <libTriInt.h>
#include <config.h>

#include "CCRules.cc"
 
using namespace scuff;

#define II cdouble(0.0,1.0) 
#define TENTHIRDS 3.33333333333333333333333

// 
#define SIPOWER   0
#define SIXFORCE  1
#define SIYFORCE  2
#define SIZFORCE  3
#define SIXTORQUE 4
#define SIYTORQUE 5
#define SIZTORQUE 6
#define NUMSIPFT  7

/***************************************************************/
/***************************************************************/
/***************************************************************/
void CalcGCFarField(double R[3], cdouble k,
                    cdouble GMuNu[3][3], cdouble CMuNu[3][3])
{ 
  double r2=R[0]*R[0] + R[1]*R[1] + R[2]*R[2];

  if (r2==0.0)
   return;

  double r=sqrt(r2);

  cdouble ExpFac=exp(II*k*r)/(4.0*M_PI*r);

  GMuNu[0][0] =               ExpFac * ( 1.0 - R[0]*R[0]/r2 );
  GMuNu[0][1] = GMuNu[1][0] = ExpFac * (     - R[0]*R[1]/r2 );
  GMuNu[0][2] = GMuNu[2][0] = ExpFac * (     - R[0]*R[2]/r2 );
  GMuNu[1][1] =               ExpFac * ( 1.0 - R[1]*R[1]/r2 );
  GMuNu[1][2] = GMuNu[2][1] = ExpFac * (     - R[1]*R[2]/r2 );
  GMuNu[2][2] =               ExpFac * ( 1.0 - R[2]*R[2]/r2 );

  CMuNu[0][0]=CMuNu[1][1]=CMuNu[2][2]=0.0;

  CMuNu[0][1] = ExpFac * R[2]/r;
  CMuNu[1][2] = ExpFac * R[0]/r;
  CMuNu[2][0] = ExpFac * R[1]/r;

  CMuNu[1][0] = -CMuNu[0][1];
  CMuNu[2][1] = -CMuNu[1][2];
  CMuNu[0][2] = -CMuNu[2][0];

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetFarFieldGCIntegrals(RWGSurface *S, int ne, 
                            double *X, cdouble k, 
                            cdouble GInt[3], cdouble CInt[3])
{
  /***************************************************************/ 
  /***************************************************************/ 
  /***************************************************************/ 
  RWGEdge *E    = S->Edges[ne];
  double *QP    = S->Vertices + 3*E->iQP;
  double *V1    = S->Vertices + 3*E->iV1;
  double *V2    = S->Vertices + 3*E->iV2;
  double *QM    = S->Vertices + 3*E->iQM;
  double Length = E->Length;

  /***************************************************************/ 
  /***************************************************************/ 
  /***************************************************************/
  double AP[3], BP[3], AM[3], BM[3];
  for(int Mu=0; Mu<3; Mu++)
   { AP[Mu] = V1[Mu] - QP[Mu];
     AM[Mu] = V1[Mu] - QM[Mu];
     BP[Mu] = V2[Mu] - QP[Mu];
     BM[Mu] = V2[Mu] - QM[Mu];
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NumPts;
  double *TCR = GetTCR(25, &NumPts);
  GInt[0]=GInt[1]=GInt[2]=CInt[0]=CInt[1]=CInt[2]=0.0;
  for(int np=0, ncp=0; np<NumPts; np++) 
   { 
     double u=TCR[ncp++];
     double v=TCR[ncp++]; 
     double w=Length * TCR[ncp++];

     /***************************************************************/
     /* set X and F=X-Q *********************************************/
     /***************************************************************/
     double XmQP[3], XP[3], XmXP[3], XmQM[3], XM[3], XmXM[3];
     for(int Mu=0; Mu<3; Mu++)
      { XmQP[Mu] = u*AP[Mu] + v*BP[Mu];
          XP[Mu] = XmQP[Mu] + QP[Mu];
        XmXP[Mu] = X[Mu] - XP[Mu];
        XmQM[Mu] = u*AM[Mu] + v*BM[Mu];
          XM[Mu] = XmQM[Mu] + QM[Mu];
        XmXM[Mu] = X[Mu] - XM[Mu];
      };

     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     cdouble G[3][3], C[3][3];
     CalcGCFarField(XmXP, k, G, C);
     GInt[0] += w*( G[0][0]*XmQP[0] + G[0][1]*XmQP[1] + G[0][2]*XmQP[2] );
     GInt[1] += w*( G[1][0]*XmQP[0] + G[1][1]*XmQP[1] + G[1][2]*XmQP[2] );
     GInt[2] += w*( G[2][0]*XmQP[0] + G[2][1]*XmQP[1] + G[2][2]*XmQP[2] );
     CInt[0] += w*( C[0][0]*XmQP[0] + C[0][1]*XmQP[1] + C[0][2]*XmQP[2] );
     CInt[1] += w*( C[1][0]*XmQP[0] + C[1][1]*XmQP[1] + C[1][2]*XmQP[2] );
     CInt[2] += w*( C[2][0]*XmQP[0] + C[2][1]*XmQP[1] + C[2][2]*XmQP[2] );

     CalcGCFarField(XmXM, k, G, C);
     GInt[0] -= w*( G[0][0]*XmQM[0] + G[0][1]*XmQM[1] + G[0][2]*XmQM[2] );
     GInt[1] -= w*( G[1][0]*XmQM[0] + G[1][1]*XmQM[1] + G[1][2]*XmQM[2] );
     GInt[2] -= w*( G[2][0]*XmQM[0] + G[2][1]*XmQM[1] + G[2][2]*XmQM[2] );
     CInt[0] -= w*( C[0][0]*XmQM[0] + C[0][1]*XmQM[1] + C[0][2]*XmQM[2] );
     CInt[1] -= w*( C[1][0]*XmQM[0] + C[1][1]*XmQM[1] + C[1][2]*XmQM[2] );
     CInt[2] -= w*( C[2][0]*XmQM[0] + C[2][1]*XmQM[1] + C[2][2]*XmQM[2] );

   };

}

/***************************************************************/
/* CR matrix stands for 'Cubature rule matrix.'                */
/* A CR matrix is an NCx7 matrix, where the nth row has entries*/
/*  w_n x_n y_n z_n nx_n ny_n nz_n                             */
/* where w is the cubature weight, (x,y,z) is the cubature pt, */
/* and nx_n, ny_n, nz_n is the outward-pointing surface normal.*/
/***************************************************************/
HMatrix *GetCRMatrix(RWGSurface *BS, double R, int NumPoints)
{
  HMatrix *CRMatrix;
  if (BS)
   { CRMatrix = new HMatrix(BS->NumPanels, 7);
     for(int np=0; np<BS->NumPanels; np++)
      { CRMatrix->SetEntry(np, 0, BS->Panels[np]->Area);
        CRMatrix->SetEntry(np, 1, BS->Panels[np]->Centroid[0]);
        CRMatrix->SetEntry(np, 2, BS->Panels[np]->Centroid[1]);
        CRMatrix->SetEntry(np, 3, BS->Panels[np]->Centroid[2]);
        CRMatrix->SetEntry(np, 4, BS->Panels[np]->ZHat[0]);
        CRMatrix->SetEntry(np, 5, BS->Panels[np]->ZHat[1]);
        CRMatrix->SetEntry(np, 6, BS->Panels[np]->ZHat[2]);
      };
   }
  else 
   {
     InitCCRules();
     int NumThetaPoints=NumPoints;
     int NumPhiPoints=2*(NumPoints+1);
     int TotalPoints = NumPhiPoints*(NumThetaPoints-2) + 2;
     CRMatrix = new HMatrix(TotalPoints, 7);
     double dPhi= 2.0*M_PI / NumPhiPoints;
     double *CCRule = CCRules[NumThetaPoints];
     if (CCRule==0) ErrExit("invalid value of NumPoints in GetCRMatrix");
     for(int npTheta=0; npTheta<NumThetaPoints; npTheta++)
      { 
        double CosTheta = CCRule[2*npTheta+0];
        double w = CCRule[2*npTheta+1];
        double Theta = acos(CosTheta);
        double SinTheta = sqrt(1.0-CosTheta*CosTheta);

        // north, south pole
        if ( npTheta==0 || npTheta==(NumThetaPoints-1) )
         { int Index  = (npTheta==0) ? 0 : TotalPoints-1;
           CRMatrix->SetEntry(Index, 0, 2.0*M_PI*R*R*w);
           CRMatrix->SetEntry(Index, 1, 0.0);
           CRMatrix->SetEntry(Index, 2, 0.0);
           CRMatrix->SetEntry(Index, 3, R*CosTheta);
           CRMatrix->SetEntry(Index, 4, 0.0);
           CRMatrix->SetEntry(Index, 5, 0.0);
           CRMatrix->SetEntry(Index, 6, CosTheta);
         }
        else
         { for(int npPhi=0; npPhi<NumPhiPoints; npPhi++)
            { int Index  = 1 + NumPhiPoints*(npTheta-1) + npPhi;
              double Phi = npPhi * dPhi;
              CRMatrix->SetEntry(Index, 0, R*R*w*dPhi);
              CRMatrix->SetEntry(Index, 1, R*SinTheta*cos(Phi));
              CRMatrix->SetEntry(Index, 2, R*SinTheta*sin(Phi));
              CRMatrix->SetEntry(Index, 3, R*CosTheta);
              CRMatrix->SetEntry(Index, 4, SinTheta*cos(Phi));
              CRMatrix->SetEntry(Index, 5, SinTheta*sin(Phi));
              CRMatrix->SetEntry(Index, 6, CosTheta);
            };
         };
      };
   };

  return CRMatrix;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *GetFSVMatrix(RWGSurface *S, HMatrix *CRMatrix,
                      cdouble Omega, cdouble Eps, cdouble Mu)
{ 
  int NE = S->NumEdges;
  int NC = CRMatrix->NR;

  HMatrix *FSVMatrix = new HMatrix(6, 2*NE*NC, LHM_COMPLEX);

  cdouble k = Omega*sqrt(Eps*Mu);
  cdouble iwe = II*Omega*Eps;
  cdouble iwu = II*Omega*Mu;
  
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  bool FarField=false;
  if (getenv("SCUFF_FARFIELD")) 
   { FarField=true;
     Log("Retaining only far-field contributions to FSV matrices.");
   };

  /***************************************************************/
  /* Think of F as a matrix of dimension 6 x N where N is the    */
  /* number of surface-current basis functions (2x number of RWG */
  /* functions) times the number of cubature points on the       */
  /* bounding surface, i.e. C = (2*NE) * NC.                     */
  /* For the pairing of the NCth cubature point with the         */
  /* P-type current in the neth RWG function,                    */
  /* we have the tuple (nc, ne, P), and this is assigned index   */
  /* nc*(2*NE) + 2*ne + P.                                       */
  /***************************************************************/
  Log("Precomputing FSV matrices (%i columns)",CRMatrix->NR);
  int NumThreads;
#ifndef USE_OPENMP
  NumThreads=1;
#else
  NumThreads=GetNumThreads();
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int ne=0; ne<NE; ne++)
   for(int nc=0; nc<NC; nc++)
    { 
      if (nc==0) LogPercent(ne,NE,5);

      // coordinates of ncth cubature point
      double X[3];
      X[0] = CRMatrix->GetEntryD(nc, 1);
      X[1] = CRMatrix->GetEntryD(nc, 2);
      X[2] = CRMatrix->GetEntryD(nc, 3);

      // get fields due to this basis function at this cubature point
      // FSVE = six-vector of fields due to electric current in bf #ne
      // FSVM =                             magnetic 
      cdouble FSVE[6], FSVM[6];
      if (FarField)
       { 
         cdouble GInt[3], CInt[3];
         GetFarFieldGCIntegrals(S, ne, X, k, GInt, CInt);

         FSVE[0] = iwu*ZVAC*GInt[0];
         FSVE[1] = iwu*ZVAC*GInt[1];
         FSVE[2] = iwu*ZVAC*GInt[2];

         FSVE[3] = -II*Omega*CInt[0];
         FSVE[4] = -II*Omega*CInt[1];
         FSVE[5] = -II*Omega*CInt[2];

         FSVM[0] =  II*Omega*CInt[0];
         FSVM[1] =  II*Omega*CInt[1];
         FSVM[2] =  II*Omega*CInt[2];

         FSVM[3] = iwe*GInt[0] / ZVAC;
         FSVM[4] = iwe*GInt[1] / ZVAC;
         FSVM[5] = iwe*GInt[2] / ZVAC;
       }
      else
       { cdouble a[3], Curla[3], Gradp[3];
         S->GetReducedPotentials(ne, X, k, 0, a, Curla, Gradp);
         for(int n=0; n<3; n++)
          { FSVE[0*3 + n] = ZVAC*( iwu*a[n] - Gradp[n]/iwe );
            FSVE[1*3 + n] = Curla[n];
            FSVM[0*3 + n] = -Curla[n];
            FSVM[1*3 + n] = (iwe*a[n] - Gradp[n]/iwu)/ZVAC;
          };
       };
       
      // store field six-vectors as columns in matrices
      int ColumnE = nc*(2*NE) + (2*ne) + 0;
      int ColumnM = nc*(2*NE) + (2*ne) + 1;
      for(int n=0; n<6; n++)
       { FSVMatrix->SetEntry(n, ColumnE, FSVE[n]);
         FSVMatrix->SetEntry(n, ColumnM, FSVM[n]);
       };

    }; // for(int ne=0; ...  for(int nc=0; ...

  return FSVMatrix;

}

static double LeviCivita[3][3][3]=
{ { { 0.0,  0.0,  0.0 }, {  0.0, 0.0, +1.0 }, {  0.0, -1.0, +0.0 }  },
  { { 0.0,  0.0, -1.0 }, {  0.0, 0.0,  0.0 }, { +1.0, +0.0, +0.0 }  },
  { { 0.0, +1.0,  0.0 }, { -1.0, 0.0,  0.0 }, {  0.0,  0.0,  0.0 }  }
};
void GetNMatrices(double nHat[3], double X[3], double XTorque[3],
                  double NMatrix[NUMSIPFT][3][3])
{
  /***************************************************************/
  /* first matrix is the power matrix ****************************/
  /***************************************************************/
  NMatrix[SIPOWER][0][0] = NMatrix[SIPOWER][1][1] = NMatrix[SIPOWER][2][2] = 0.0;
  NMatrix[SIPOWER][1][2] = nHat[0]; NMatrix[SIPOWER][2][1] = -nHat[0];
  NMatrix[SIPOWER][2][0] = nHat[1]; NMatrix[SIPOWER][0][2] = -nHat[1];
  NMatrix[SIPOWER][0][1] = nHat[2]; NMatrix[SIPOWER][1][0] = -nHat[2];

  /***************************************************************/
  /* next three matrices are force matrices                      */
  /***************************************************************/
  NMatrix[SIXFORCE][0][0]                           =  nHat[0]; 
  NMatrix[SIXFORCE][1][1]                           = -nHat[0];
  NMatrix[SIXFORCE][2][2]                           = -nHat[0];
  NMatrix[SIXFORCE][0][1] = NMatrix[SIXFORCE][1][0] =  nHat[1];
  NMatrix[SIXFORCE][0][2] = NMatrix[SIXFORCE][2][0] =  nHat[2];
  NMatrix[SIXFORCE][1][2] = NMatrix[SIXFORCE][2][0] =  0.0;

  NMatrix[SIYFORCE][0][0]                           = -nHat[1]; 
  NMatrix[SIYFORCE][1][1]                           =  nHat[1];
  NMatrix[SIYFORCE][2][2]                           = -nHat[1];
  NMatrix[SIYFORCE][0][1] = NMatrix[SIYFORCE][1][0] =  nHat[0];
  NMatrix[SIYFORCE][1][2] = NMatrix[SIYFORCE][2][1] =  nHat[2];
  NMatrix[SIYFORCE][0][2] = NMatrix[SIYFORCE][2][0] =  0.0;

  NMatrix[SIZFORCE][0][0]                           = -nHat[2]; 
  NMatrix[SIZFORCE][1][1]                           =  nHat[2];
  NMatrix[SIZFORCE][2][2]                           = -nHat[2];
  NMatrix[SIZFORCE][0][2] = NMatrix[SIZFORCE][2][0] =  nHat[0];
  NMatrix[SIZFORCE][1][2] = NMatrix[SIZFORCE][2][1] =  nHat[1];
  NMatrix[SIZFORCE][0][1] = NMatrix[SIZFORCE][1][0] =  0.0;

  /***************************************************************/
  /* last three matrices are the torque matrices                 */
  /***************************************************************/
  double D[3], DxN[3], EpsD[3][3]; 
  D[0] = X[0]-XTorque[0];
  D[1] = X[1]-XTorque[1];
  D[2] = X[2]-XTorque[2];
  DxN[0] = D[1]*nHat[2] - D[2]*nHat[1]; 
  DxN[1] = D[2]*nHat[0] - D[0]*nHat[2];
  DxN[2] = D[0]*nHat[1] - D[1]*nHat[0];

  for(int a=0; a<3; a++)
   for(int b=0; b<3; b++)
    for(int i=0; i<3; i++)
     { NMatrix[SIXTORQUE+i][a][b] = (a==b) ? -DxN[i] : 0.0;
       for(int j=0; j<3; j++)
        NMatrix[SIXTORQUE+i][a][b] 
         +=   LeviCivita[i][j][a]*D[j]*nHat[b]
            + LeviCivita[i][j][b]*D[j]*nHat[a];
     };

}

/***************************************************************/
/* This routine gets the contributions of a single pair of     */
/* RWG basis functions to all SIPFT matrices. Each pair of     */
/* BFs contributes to up to 4 SIPFT matrix entries.            */
/*                                                             */
/* SIPFT = 'surface integral for power / force / torque'       */
/*                                                             */
/* Entries is an array with enough room for 4xNUMSIPFT cdoubles.*/
/* On return,                                                  */
/*  Entries[4*n + 0..3] = (EE, EM, ME, MM) entries of n-th     */
/*                         SIPFT matrix.                       */ 
/***************************************************************/
void GetSIPFTMatrixEntries(RWGSurface *S, int neA, int neB,
                           HMatrix *CRMatrix, HMatrix *FSVMatrix,
                           cdouble Omega, cdouble Eps, cdouble Mu,
                           double *XTorque,
                           bool NeedMatrix[NUMSIPFT], cdouble *SEntries)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble Z = sqrt(Mu/Eps);
  cdouble EpsAbs = TENTHIRDS * Eps / ZVAC;
  cdouble MuAbs  = TENTHIRDS * Mu * ZVAC;

  int NE = S->NumEdges;
  int NC = CRMatrix->NR;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double NMatrix[NUMSIPFT][3][3];
  for(int n=0; n<NUMSIPFT; n++)
   for(int Mu=0; Mu<3; Mu++)
    for(int Nu=0; Nu<3; Nu++)
     NMatrix[n][Mu][Nu] = 0.0;

  /***************************************************************/
  /* if the user didn't specify a torque center, use the origin  */
  /* of coordinates of the meshed surface, possibly translated   */
  /* if the surface has been translated since reading from file  */
  /***************************************************************/
  double DefaultXTorque[3]={0.0, 0.0, 0.0};
  if (XTorque==0)
   { if (S->GT) S->GT->Transform(DefaultXTorque);
     XTorque=DefaultXTorque;
   };
 
  /***************************************************************/
  /* nc loops over cubature points.                              */
  /***************************************************************/
  memset(SEntries, 0, 4*NUMSIPFT*sizeof(cdouble));
  for (int nc=0; nc<NC; nc++)
   {  
     double w, X[3], nHat[3];

     w       = CRMatrix->GetEntryD(nc, 0);
     X[0]    = CRMatrix->GetEntryD(nc, 1);
     X[1]    = CRMatrix->GetEntryD(nc, 2);
     X[2]    = CRMatrix->GetEntryD(nc, 3);
     nHat[0] = CRMatrix->GetEntryD(nc, 4);
     nHat[1] = CRMatrix->GetEntryD(nc, 5);
     nHat[2] = CRMatrix->GetEntryD(nc, 6);

     /***************************************************************/
     /* get 3x3 N matrices at this cubature point *******************/
     /***************************************************************/
     GetNMatrices(nHat, X, XTorque, NMatrix);

     /***************************************************************/
     /* fetch 'F' six-vectors for this cubature point           *****/
     /***************************************************************/
     int ColumnEA = nc*(2*NE) + (2*neA) + 0;
     int ColumnMA = nc*(2*NE) + (2*neA) + 1;
     int ColumnEB = nc*(2*NE) + (2*neB) + 0;
     int ColumnMB = nc*(2*NE) + (2*neB) + 1;
     cdouble FSVEA[6], FSVMA[6], FSVEB[6], FSVMB[6];
     for(int n=0; n<6; n++)
      { FSVEA[n] = FSVMatrix->GetEntry(n, ColumnEA);
        FSVMA[n] = FSVMatrix->GetEntry(n, ColumnMA);
        FSVEB[n] = FSVMatrix->GetEntry(n, ColumnEB);
        FSVMB[n] = FSVMatrix->GetEntry(n, ColumnMB);
      };
 
     /***************************************************************/
     /* entries of power matrix               ***********************/
     /***************************************************************/
     if (NeedMatrix[SIPOWER])
      { for(int m=0; m<3; m++)
         for(int n=0; n<3; n++)
          { SEntries[ 0 ] += w * conj(FSVEA[m]) * NMatrix[SIPOWER][m][n] * FSVEB[3+n];
            SEntries[ 1 ] += w * conj(FSVEA[m]) * NMatrix[SIPOWER][m][n] * FSVMB[3+n];
            SEntries[ 2 ] += w * conj(FSVMB[m]) * NMatrix[SIPOWER][m][n] * FSVEB[3+n];
            SEntries[ 3 ] += w * conj(FSVMB[m]) * NMatrix[SIPOWER][m][n] * FSVMB[3+n];
          }; 
      };

     /***************************************************************/
     /* entries of force and torque matrices ************************/
     /***************************************************************/
     for(int nSIFT=0; nSIFT<NUMSIPFT-1; nSIFT++)
      { int nSIPFT=nSIFT+1;
        if (NeedMatrix[nSIPFT]==false) continue;
        cdouble UpperVMVP[4]={0.0, 0.0, 0.0, 0.0};
        cdouble LowerVMVP[4]={0.0, 0.0, 0.0, 0.0};
        for(int m=0; m<3; m++)
         for(int n=0; n<3; n++)
          { 
            UpperVMVP[ 0 ] += conj(FSVEA[m]) * NMatrix[nSIPFT][m][n] * FSVEB[n];
            UpperVMVP[ 1 ] += conj(FSVEA[m]) * NMatrix[nSIPFT][m][n] * FSVMB[n];
            UpperVMVP[ 2 ] += conj(FSVMA[m]) * NMatrix[nSIPFT][m][n] * FSVEB[n];
            UpperVMVP[ 3 ] += conj(FSVMA[m]) * NMatrix[nSIPFT][m][n] * FSVMB[n];

            LowerVMVP[ 0 ] += conj(FSVEA[3+m]) * NMatrix[nSIPFT][m][n] * FSVEB[3+n];
            LowerVMVP[ 1 ] += conj(FSVEA[3+m]) * NMatrix[nSIPFT][m][n] * FSVMB[3+n];
            LowerVMVP[ 2 ] += conj(FSVMA[3+m]) * NMatrix[nSIPFT][m][n] * FSVEB[3+n];
            LowerVMVP[ 3 ] += conj(FSVMA[3+m]) * NMatrix[nSIPFT][m][n] * FSVMB[3+n];
          };
        SEntries[ 4*nSIPFT + 0 ] += w*( EpsAbs*UpperVMVP[0] + MuAbs*LowerVMVP[0] );
        SEntries[ 4*nSIPFT + 1 ] += w*( EpsAbs*UpperVMVP[1] + MuAbs*LowerVMVP[1] );
        SEntries[ 4*nSIPFT + 2 ] += w*( EpsAbs*UpperVMVP[2] + MuAbs*LowerVMVP[2] );
        SEntries[ 4*nSIPFT + 3 ] += w*( EpsAbs*UpperVMVP[3] + MuAbs*LowerVMVP[3] );
      }; //for(int nSIFT=0; nSIFT<NUMSIPFT-1; nSIFT++)
   }; //for (int nc=0; nc<NC; nc++)

  // put in the factors of 1/4
  for(int nSIPFT=0; nSIPFT<4*NUMSIPFT; nSIPFT++)
   SEntries[nSIPFT] *= 0.25;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetSIPFTMatrices(RWGGeometry *G, int WhichSurface,
                      RWGSurface *BS, int R, int NumPoints,
                      cdouble Omega, bool NeedMatrix[NUMSIPFT],
                      HMatrix *MSIPFT[NUMSIPFT])
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  RWGSurface *S=G->Surfaces[WhichSurface];
  if (BS)
   Log("Computing SIPFT matrix for surface %s, bounding surface %s...",
        S->MeshFileName, BS->MeshFileName);
  else
   Log("Computing SIPFT matrix for surface %s, (R,NPts)=(%e,%i)",
        S->MeshFileName, R, NumPoints);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NE = S->NumEdges;
  int NBF = S->NumBFs;

  if (S->IsPEC) 
   ErrExit("GetSIPFTMatrices() not implemented for PEC bodies");

  cdouble Eps, Mu;
  G->RegionMPs[ S->RegionIndices[0] ] -> GetEpsMu(Omega, &Eps, &Mu);

  double XTorque[3] = {0.0, 0.0, 0.0};
  if (S->GT)
   S->GT->Transform(XTorque);

  /***************************************************************/
  /* (re)allocate SIPFT matrices as necessary *********************/
  /***************************************************************/
  for (int n=0; n<NUMSIPFT; n++)
   { if (NeedMatrix[n]==false) continue;

     HMatrix *M = MSIPFT[n];
     if ( M && M->NR!=NBF || M->NC!=NBF )
      { Warn("incorrect matrix passed to ComputeSMatrix (reallocating)...");
        delete MSIPFT[n];
        M=0;
      };
     if (M==0)
      { M=new HMatrix(NBF, NBF, LHM_COMPLEX);
        MSIPFT[n] = M; 
      };
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HMatrix *CRMatrix  = GetCRMatrix(BS, R, NumPoints);
  HMatrix *FSVMatrix = GetFSVMatrix(S, CRMatrix, Omega, Eps, Mu);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NumThreads;
#ifndef USE_OPENMP
  NumThreads=1;
#else
  NumThreads=GetNumThreads();
  Log("OpenMP multithreading (%i threads)",NumThreads);
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int Alpha=0; Alpha<NE; Alpha++)
   for(int Beta=0; Beta<NE; Beta++)
    { 
      if (Beta==0) LogPercent(Alpha,NE,5);
      cdouble Entries[NUMSIPFT*4];

      GetSIPFTMatrixEntries(S, Alpha, Beta, CRMatrix, FSVMatrix,
                            Omega, Eps, Mu, XTorque, NeedMatrix, Entries);

      for(int n=0; n<NUMSIPFT; n++)
       { if (NeedMatrix[n]==false) continue;
         MSIPFT[n]->SetEntry( 2*Alpha+0, 2*Beta+0, Entries[4*n + 0] );
         MSIPFT[n]->SetEntry( 2*Alpha+0, 2*Beta+1, Entries[4*n + 1] );
         MSIPFT[n]->SetEntry( 2*Alpha+1, 2*Beta+0, Entries[4*n + 2] );
         MSIPFT[n]->SetEntry( 2*Alpha+1, 2*Beta+1, Entries[4*n + 3] );
       };

    };

  delete FSVMatrix;
  delete CRMatrix;
  
}

/***************************************************************/
/* 'hermitian vector/matrix/vector product'                    */
/***************************************************************/
double HVMVP(cdouble V1[3], double M[3][3], cdouble V2[3])
{ 
  double Sum=0.0;

  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    Sum += real( conj(V1[Mu]) * M[Mu][Nu] * V2[Nu]);

  return Sum;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetSIPFT(RWGGeometry *G, IncField *IF, HVector *KN, 
              cdouble Omega, RWGSurface *BS, double R, int NumPoints,
              double SIPFT[NUMSIPFT])
{
  Log("Computing SIPFT (R,NPts)=(%e,%i)", R, NumPoints);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  HMatrix *CRMatrix  = GetCRMatrix(BS, R, NumPoints);

  // we assume that all cubature points lie in the same region 
  // of the scuff geometry, so we use the first point in the rule
  // to determine which region that is and look up its eps/mu
  double X[3];
  X[0]=CRMatrix->GetEntryD(0,1);
  X[1]=CRMatrix->GetEntryD(0,2);
  X[2]=CRMatrix->GetEntryD(0,3);
  int RegionIndex=G->GetRegionIndex(X);
  cdouble EpsRel, MuRel;
  G->RegionMPs[ RegionIndex ] -> GetEpsMu(Omega, &EpsRel, &MuRel);
  double EpsAbs = TENTHIRDS * real(EpsRel) / ZVAC;
  double  MuAbs = TENTHIRDS * real(MuRel) * ZVAC;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  //HMatrix XMatrix = CRMatrix->GetEntries(":","1:3");
  HMatrix *XMatrix = new HMatrix(CRMatrix->NR, 3);
  for(int nr=0; nr<XMatrix->NR; nr++)
   { XMatrix->SetEntry(nr, 0, CRMatrix->GetEntryD(nr,1));
     XMatrix->SetEntry(nr, 1, CRMatrix->GetEntryD(nr,2));
     XMatrix->SetEntry(nr, 2, CRMatrix->GetEntryD(nr,3));
   };
  HMatrix *FMatrix = G->GetFields(IF, KN, Omega, XMatrix);

  double XTorque[3] = {0.0, 0.0, 0.0};

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  memset(SIPFT, 0, NUMSIPFT*sizeof(double));
  for(int nr=0; nr<CRMatrix->NR; nr++)
   { 
     double w, X[3], nHat[3];
           w = CRMatrix->GetEntryD(nr, 0);
        X[0] = CRMatrix->GetEntryD(nr, 1);
        X[1] = CRMatrix->GetEntryD(nr, 2);
        X[2] = CRMatrix->GetEntryD(nr, 3);
     nHat[0] = CRMatrix->GetEntryD(nr, 4);
     nHat[1] = CRMatrix->GetEntryD(nr, 5);
     nHat[2] = CRMatrix->GetEntryD(nr, 6);

     double NMatrix[NUMSIPFT][3][3];
     GetNMatrices(nHat, X, XTorque, NMatrix);

     cdouble E[3], H[3];
     //FMatrix->GetEntries(nr, "0:2", E);
     E[0] = FMatrix->GetEntry(nr, 0);
     E[1] = FMatrix->GetEntry(nr, 1);
     E[2] = FMatrix->GetEntry(nr, 2);
     H[0] = FMatrix->GetEntry(nr, 3);
     H[1] = FMatrix->GetEntry(nr, 4);
     H[2] = FMatrix->GetEntry(nr, 5);

     SIPFT[SIPOWER] += 0.25 * w * (  HVMVP(E, NMatrix[SIPOWER], H)
                                    -HVMVP(H, NMatrix[SIPOWER], E)  
                                  );

     for(int n=SIXFORCE; n<=SIZTORQUE; n++)
      SIPFT[n] += 0.25 * w * ( EpsAbs*HVMVP(E, NMatrix[n], E)
                               +MuAbs*HVMVP(H, NMatrix[n], H)
                             );
   };

  delete FMatrix; 
  delete XMatrix; 
  delete CRMatrix;
  
}
