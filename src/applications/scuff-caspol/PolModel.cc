/*
 * PolModel.cc   -- implementation of a simple class for describing
 *               -- the frequency-dependent polarizability of atoms
 *               -- and molecules
 *
 * homer reid    -- 1/2012
 *
 */

#include "scuff-caspol.h"
#include <libhmat.h>

/***************************************************************/
/* PolModel class constructor: construct an instance of        */
/* PolModel from a data file.                                  */
/***************************************************************/
PolModel::PolModel(const char *PolFileName)
{

   // the current implementation reads in a table of polarizability
   // data from a file and constructs an interpolator,  
   // which is used in the GetPolarizability() routine below 
   // to compute the polarizability at arbitrary frequencies

   // read in the data file as an HMatrix
   HMatrix *PolDataMatrix=new HMatrix("PolData.dat",LHM_TEXT);
   if (PolDataMatrix->ErrMsg)
    ErrExit(PolDataMatrix->ErrMsg);

   // extract the 1st and 7th columns as the X and Y data
   // (the 7th column is for rubidium)
   int n, N=PolDataMatrix->NR;
   double XValues[N], YValues[N];
   for(n=0; n<N; n++)
    { XValues[n]=PolDataMatrix->GetEntryD(n,0);
      YValues[n]=PolDataMatrix->GetEntryD(n,6);
    };

   // initialize the interpolator
   PolInterp = new Interp1D(XValues, YValues, N, 1);

   delete PolDataMatrix;
  
}


/***************************************************************/
/* routine to compute polarizability tensor at a given         */
/* frequency.                                                  */
/*                                                             */
/* on input, Xi is the imaginary frequency in units of         */
/* 3e14 rad/sec.                                               */
/*                                                             */
/* Alpha points to a vector of size 9 that the routine must    */
/* fill in with the components of the polarizability tensor,   */
/* as follows:                                                 */
/*                                                             */
/*  Alpha[0] = Alpha_{xx}                                      */
/*  Alpha[1] = Alpha_{xy}                                      */
/*  Alpha[2] = Alpha_{xz}                                      */
/*  Alpha[3] = Alpha_{yx}                                      */
/*  ...                                                        */
/*                                                             */
/* etc. (in general: Alpha[ i + 3*j ] = Alpha_{ij} )           */
/***************************************************************/
void PolModel::GetPolarizability(double Xi, double *Alpha) 
 {
  
   // the constant here convertes Xi from my units, 
   // in which '1' == 3e14 rad/sec, to atomic units, 
   // in which '1' == 2.598e+17 rad/sec
   double AlphaDiag;
   PolInterp->Evaluate(0.00115493*Xi , &AlphaDiag);

   // note we only have to fill in the nonzero entries of Alpha
   Alpha[0]=Alpha[4]=Alpha[9]=AlphaDiag;

 }
