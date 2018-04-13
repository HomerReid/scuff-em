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
 * TranslationMatrices.cc   -- compute the 'translation matrices' that relate
 *                             scalar and vector helmholtz solutions
 *                             about different origins 
 *
 * more specifically: 
 *
 *  a) consider two points in space xSource, xOrigin (which we think of as the
 *     center of a source distribution and the origin of a spherical 
 *     coordinate system in which we will express the fields of this 
 *     distribution as an expansion in spherical waves).
 *
 *  b) put Xij = xSource - xOrigin.
 *
 *  c) now consider a third point xDest that lies closer to xOrigin than to xSource.
 *
 *  d) let Phi_{lm}(xDest-xSource) be the value at xDest of an outgoing
 *     wave emitted by the source at xSource
 *
 *  e) let Psi_{lm}(xDest-xOrigin) be the value of an INTERIOR helmholtz
 *     solution at xDest as reckoned using a coordinate system with
 *     origin xOrigin.
 *
 *  e) then the translation matrix expresses a single Phi_{lm} as a 
 *     linear combination of Psi_{lm}s: 
 *
 *      Phi_{Alpha}(xDest-xSource)
 *        = \sum_{AlphaP} A_{Alpha,AlphaP} Psi_{AlphaP}(xDest-xOrigin)
 *
 * where Alpha=(lm), AlphaP=(lp,mp) are compound indices, and 
 * where A is the matrix computed by GetTranslationMatrices below.
 *
 * Similarly, B and C are the matrices that relate M-type vector solutions
 * about xp to M- and N-type vector solutions about x: we have 
 * 
 * [M^{ext}(xDest-xSource)]     = [ B   C ]         [ M^{regular}(xDest-xOrigin) ]
 * [N^{ext}(xDest-xSource)]_{A} = [ -C  B ]_{A, AP} [ N^{regular}(xDest-xOrigin) ]_{AP}
 *
 * This calculation follows Wittmann, 'Spherical Wave Operators and the
 * Translation Formulas,' IEEE Transactions on Antennas and Propagation *36* 1078
 * (1988), with the distinction that my N function differs from Wittmann's by 
 * a factor of I. In Wittmann's convention one has N = \curl M / k, whereas 
 * in my convention one has N = \curl M / (-ik).
 *
 * homer reid               -- 4/2005 -- 8/2016
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "libhrutil.h"
#include "libhmat.h"
#include "libSpherical.h"

#define II cdouble(0,1)
#define M1POW(n) ( (n)%2 ? -1.0 : 1.0)
#define IIPOW(n) ( (n)%4 == 0 ? 1 : (n)%4 == 1 ? II : (n)%4 == 2 ? -1 : -1.0*II )

/***************************************************************/
/***************************************************************/
/***************************************************************/
extern "C" {
void drc3jm_(double *L1, double *L2, double *L3, 
             double *M1, double *M2Min, double *M3Max, 
             double *Result, int *NDim, int *ier);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
double ThreeJSymbol(double L1, double L2, double L3, 
                    double M1, double M2, double M3)
{
  (void) M3; // unused 
  double M2Max = L2 > L3-M1     ?  L2 :   L3-M1;
  double M2Min = -L2 > -(L3+M1) ? -L2 : -(L3+M1);
  int Length = (int)(M2Max - M2Min) + 1;
  double *Result = new double[Length];
  int ier;
  drc3jm_(&L1, &L2, &L3, &M1, &M2Min, &M2Max, Result, &Length, &ier);
  double RetVal=0.0;
  for(int im2=M2Min; im2<=M2Max; im2++)
   if ( im2==(int)M2 )
    RetVal = Result[im2-((int)M2Min)];
  delete[] Result;
  return RetVal;
}   

/***************************************************************/
/* Compute the translation matrices that express outgoing      */
/* Helmholtz solutions emanating from a source point xSource   */
/* as linear combinations of regular solutions in a spherical  */
/* coordinate system centered at the origin xOrigin            */
/*                                                             */
/* Scalar waves:                                               */
/*  PsiOut_a(xDest-XSource)                                    */
/*   = \sum A_{ab}(XSource) PsiReg_b(xDest)                    */
/*                                                             */
/*                                                             */
/* inputs:                                                     */
/*                                                             */
/*  Xij = difference vector xSource - xOrigin                  */ 
/*                                                             */
/*    k     = wavenumber that enters into exponent e^{ikr}     */ 
/*                                                             */
/*  LMax    = maximum \Ell-value of translation matrix entry   */
/*            computed                                         */
/*                                                             */
/* Outputs:                                                    */
/*  on return, A->GetEntry(Alpha, AlphaP) is the matrix entry  */
/*  described above, and similarly for B and C.                */
/*                                                             */
/* Note: on entry, A, B, C must point to HMatrices that have   */
/* been preallocated like this:                                */
/*                                                             */
/* int NAlpha=(LMax+1)*(LMax+1);                               */
/*                                                             */
/* A=new HMatrix(NAlpha, NAlpha, LHM_COMCLEX);                 */
/* B=new HMatrix(NAlpha, NAlpha, LHM_COMCLEX);                 */
/* C=new HMatrix(NAlpha, NAlpha, LHM_COMCLEX);                 */
/***************************************************************/
void GetTranslationMatrices(double Xij[3], cdouble k, int LMax,
                            HMatrix *A, HMatrix *B, HMatrix *C)
{
  A->Zero();
  B->Zero();
  C->Zero();
  if ( abs(k)*VecNorm(Xij) < 1.0e-6 )
   { for(int Alpha=0; Alpha<(LMax+1)*(LMax+1); Alpha++)
      { A->SetEntry(Alpha,Alpha,1.0);
        B->SetEntry(Alpha,Alpha,1.0);
      };
     return;
   };

  int LCMax=2*LMax;
  cdouble *R = new cdouble[LCMax+2];
  cdouble *Ylm = new cdouble[(LCMax+1)*(LCMax+1)];

  /***************************************************************/
  /* get all spherical bessel functions and spherical harmonics **/
  /* that we will need for this computation                     **/
  /***************************************************************/
  double r, Theta, Phi;
  CoordinateC2S(Xij, &r, &Theta, &Phi);
  GetRadialFunctions(LCMax, k, r, LS_OUTGOING, R, 0);
  GetYlmArray(LCMax, Theta, Phi, Ylm);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(int Alpha=0, LA=0; LA<=LMax; LA++)
   for(int MA=-LA; MA<=LA; MA++, Alpha++)
    for(int Beta=0, LB=0; LB<=LMax; LB++)
     for(int MB=-LB; MB<=LB; MB++, Beta++)
      { 
        int MC = MA-MB;
        cdouble AA=0.0, BB=0.0;
        for(int LC=abs(LA-LB); LC<=LA+LB; LC++)
         { 
           cdouble Factor= 4.0*M_PI*M1POW(MA)*IIPOW(LA-LB+LC)
                          *sqrt((2*LA+1)*(2*LB+1)*(2*LC+1)/(4.0*M_PI))
                          *ThreeJSymbol(LA,LB,LC,0,0,0)
                          *ThreeJSymbol(LA,LB,LC,-MA,MB,MC)
                          *R[LC]*Ylm[LM2ALPHA(LC,MC)];
           AA+=Factor;
           if (LA>0 && LB>0)
            BB += (double)( LA*(LA+1.0) + LB*(LB+1.0) - LC*(LC+1.0)) * Factor;
         };
        A->SetEntry(Alpha,Beta,AA);
        if (LB>0)
         B->SetEntry(Alpha,Beta,BB/(2.0*sqrt(LA*(LA+1.0)*LB*(LB+1.0))));
      };

  /***************************************************************/
  /* compute the C matrix using the A matrix                     */
  /* Note: My N function = i*Wittman's N function.               */
  /* Thus, my C coefficient = ittman's N function.               */
  /***************************************************************/
  for(int Alpha=1, LA=1; LA<=LMax; LA++)
   for(int MA=-LA; MA<=LA; MA++, Alpha++)
    for(int Beta=1, LB=1; LB<=LMax; LB++)
     for(int MB=-LB; MB<=LB; MB++, Beta++)
      { 
        cdouble CC = Xij[2]*MA*A->GetEntry(Alpha,Beta);

        if (MA<LA) 
         CC += 0.5*(Xij[0] - II*Xij[1])
                  *sqrt((LA-MA)*(LA+MA+1.0))
                  *A->GetEntry(LM2ALPHA(LA,MA+1), Beta);

        if (MA>(-LA)) 
         CC += 0.5*(Xij[0] + II*Xij[1])
                  *sqrt((LA+MA)*(LA-MA+1.0))
                  *A->GetEntry(LM2ALPHA(LA,MA-1), Beta);

        C->SetEntry(Alpha, Beta, -k*CC / sqrt( LA*(LA+1.0)*LB*(LB+1.0) ) );
      };

  delete[] R;
  delete[] Ylm;

}
