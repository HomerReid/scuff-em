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
 *  a) consider two points in space x, xp (which we think of as the 
 *     centers of two source distributions).
 *
 *  b) put Xij = x - xp, the vector pointing from xp to x.
 *
 *  c) now consider a third point xpp that lies closer to xp than to x.
 *
 *  d) let Phi_{lm}(xpp-x) be the value of an EXTERIOR helmholtz 
 *     solution at xpp as reckoned from a coordinate system whose origin
 *     lies at x. 
 *
 *     (for real frequencies, Phi_{lm}(R) = h_l(kr) Y_{lm}(theta, phi);
 *      for imag frequencies, Phi_{lm}(R) = k_l(kr) Y_{lm}(theta, phi))
 *
 *  e) similarly, let Psi_{lm}(xpp-xp) be the value of an INTERIOR helmholtz 
 *     solution at xpp as reckoned from a coordinate system whose origin
 *     lies at xp. 
 *
 *     (for real frequencies, Psi_{lm}(R) = j_l(kr) Y_{lm}(theta, phi);
 *      for imag frequencies, Psi_{lm}(R) = i_l(kr) Y_{lm}(theta, phi))
 *     
 *  e) then the translation matrix expresses Phi_{lm}(xpp-x) in terms
 *     of Psi_{lm}(xpp-xp), as follows:
 *
 *      Phi_{Alpha}(xpp-x) 
 *        = \sum_{AlphaP} A_{Alpha,AlphaP} Psi_{AlphaP}(xpp-xp) 
 *
 * where Alpha=(lm), AlphaP=(lp,mp) are compound indices, and 
 * where A is the matrix computed by GetTranslationMatrices below.
 *
 * Similarly, B and C are the matrices that relate M-type vector solutions
 * about xp to M- and N-type vector solutions about x: we have 
 * 
 * [M^{exterior}(xpp-x)]     = [ B   C ]         [ M^{interior}(xpp-xp) ]
 * [N^{exterior}(xpp-x)]_{A} = [ -C  B ]_{A, AP} [ N^{interior}(xpp-xp) ]_{AP}
 *
 * This calculation follows Chew, Fields and Waves in Inhomogeneous Media,
 * Appendix D, with the difference that my vector-N function is equal to 
 * i times Chews's vector-N function, and thus my 'C' matrix coefficient 
 * is equal to -i times Chew's 'C' matrix coefficient. (My convention agrees 
 * with Jackson's.)
 *
 * homer reid               -- 4/2005 -- 4/2011
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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
/* Compute the translation matrices that relate scalar and     */
/* vector Helmholtz solutions about different origins, as      */
/* described in detail above.                                  */
/*                                                             */
/* inputs:                                                     */
/*                                                             */
/*  Xij = difference vector                                    */ 
/*                                                             */
/*    k     = wavenumber that enters into exponent e^{ikr}     */ 
/*                                                             */
/*     lMax = maximum l-value of translation matrix entry      */
/*            computed                                         */
/*                                                             */
/* Outputs:                                                    */
/*  on return, A->GetEntry(Alpha, AlphaP) is the matrix entry  */
/*  described above, and similarly for B and C.                */
/*                                                             */
/* Note: on entry, A, B, C must point to HMatrices that have   */
/* been preallocated like this:                                */
/*                                                             */
/* int nAlpha=(lMax+1)*(lMax+1);                               */
/*                                                             */
/* A=new HMatrix(nAlpha, nAlpha, LHM_COMPLEX);                 */
/* B=new HMatrix(nAlpha, nAlpha, LHM_COMPLEX);                 */
/* C=new HMatrix(nAlpha, nAlpha, LHM_COMPLEX);                 */
/***************************************************************/
void GetTranslationMatrices(double Xij[3], cdouble k,
                            int lMax, HMatrix *A, HMatrix *B, HMatrix *C)
{
  int l, m, Alpha;
  int lP, mP, AlphaP;
  int lPP, mPP, AlphaPP;
  double r, Theta, Phi, Sign;
  cdouble ik, Factor, AA, BB, CC, F1, F2, F3;

  int lPPMax=2*lMax;
  cdouble *R = new cdouble[lPPMax+2];
  cdouble *Ylm = new cdouble[(lPPMax+1)*(lPPMax+1)];

  /***************************************************************/
  /* get all spherical bessel functions and spherical harmonics **/
  /* that we will need for this computation                     **/
  /***************************************************************/
  CoordinateC2S(Xij, &r, &Theta, &Phi);
  GetRadialFunctions(lPPMax, k, r, LS_OUTGOING, R, 0);
  GetYlmArray(lPPMax, Theta, Phi, Ylm);

  ik = II*k;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  A->Zero();
  B->Zero();
  C->Zero();
  for(Alpha=l=0; l<=lMax; l++)
   for(m=-l; m<=l; m++, Alpha++)
    for(AlphaP=lP=0; lP<=lMax; lP++)
     for(mP=-lP; mP<=lP; mP++, AlphaP++)
      { 
        mPP=m-mP;
        for(AA=BB=0.0, lPP=abs(l-lP); lPP<=l+lP; lPP++)
         { AlphaPP=LM2ALPHA(lPP, mPP);
 
           cdouble IIPow = IIPOW(lP + lPP - l);

           Factor= 4.0*M_PI*IIPow*M1POW(m)
                   *sqrt((2*l+1)*(2*lP+1)*(2*lPP+1)/(4.0*M_PI))
                   *ThreeJSymbol(l,lP,lPP,0,0,0)*ThreeJSymbol(l,lP,lPP,-m,mP,mPP)
                   *R[lPP]*Ylm[AlphaPP];
           AA+=Factor;
           if (l>0 && lP>0)
            BB += (double)((l+1) + lP*(lP+1) - lPP*(lPP+1)) * Factor;
         };
        A->SetEntry(Alpha,AlphaP,AA);
        if (lP>0) B->SetEntry(Alpha,AlphaP,BB/(2.0*lP*(lP+1)));
      };

  /***************************************************************/
  /* compute the C matrix using the A matrix *********************/
  /***************************************************************/
  cdouble krC = k*r*cos(Theta); // in Chew these have a prefactor of II 
  cdouble krS = k*r*sin(Theta); // (not here bc different convention)
  cdouble epiP  = exp(+II*Phi); 
  cdouble emiP  = exp(-II*Phi);
  for(Alpha=l=1; l<=lMax; l++)
   for(m=-l; m<=l; m++, Alpha++)
    for(AlphaP=lP=1; lP<=lMax; lP++)
     for(mP=-lP; mP<=lP; mP++, AlphaP++)
      { 
        cdouble CC = krC*(2.0*mP)*A->GetEntry(Alpha,AlphaP);

        if (mP<lP) 
         CC+= krS*epiP*sqrt((lP-mP)*(lP+mP+1) )*A->GetEntry(Alpha, LM2ALPHA(lP,mP+1));

        if (mP>(-lP)) 
         CC+= krS*emiP*sqrt((lP+mP)*(lP-mP+1) )*A->GetEntry(Alpha, LM2ALPHA(lP,mP-1));

        C->SetEntry(Alpha, AlphaP, CC / ((double)(2*lP*(lP+1))));
      };

  delete[] R;
  delete[] Ylm;

}
