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
 * PFTOptions.h -- interface for passing non-default options
 *              -- to SCUFF-EM's algorithms for computing 
 *              -- power, force, and torque
 *
 * homer reid   -- 1/2015
 */
#ifndef PFTOPTIONS_H
#define PFTOPTIONS_H

#include "libhmat.h"
#include "libhrutil.h"
#include "libIncField.h"

namespace scuff {

/***************************************************************/
/* number of power, force, and torque quantities, plus indices */
/* into arrays of them                                         */
/***************************************************************/
#define PFT_PABS    0  // absorbed power
#define PFT_PSCAT   1  // scattered power
#define PFT_XFORCE  2  // force components
#define PFT_YFORCE  3  // 
#define PFT_ZFORCE  4  // 
#define PFT_XTORQUE 5  // torque components
#define PFT_YTORQUE 6  // 
#define PFT_ZTORQUE 7  // 
#define NUMPFT      8

// for NEQ calculations, radiated power goes in the PSCAT slot
#define PFT_PRAD    1

#define NUMSRFLUX   12 // number of spatially-resolved flux quantities

// 10/3 is the value of the conversion factor that
// converts forces from internal SCUFF units of force
// into the units in which SCUFF reports forces
// (nanonewtons).
// internal SCUFF force units:
//    1 volt * 1 amp / (1 micron * omega_0 )
//  = 1 watt / 3e8 m
//  = (10/3) nanoNewtons
// (where omega_0 = 3e14 rad/sec = c/1 micron)
#define TENTHIRDS 3.33333333333333333333333

/***************************************************************/
/* values for the PFTMethod field of PFTOptions ****************/
/***************************************************************/
#define SCUFF_PFT_OVERLAP       0   // overlap method
#define SCUFF_PFT_DSI           1   // displaced-surface-integral method
#define SCUFF_PFT_EP            2   // equivalence-principle method

// these two options mean: compute force and torque using 
// overlap / DSI, but compute power (both absorbed and scattered) using EP
#define SCUFF_PFT_EPOVERLAP     3
#define SCUFF_PFT_EPDSI         4

#define SCUFF_PFT_EMT           5   // energy/momentum transfer method
#define SCUFF_PFT_EMT_EXTERIOR  5   //  (exterior is default)
#define SCUFF_PFT_EMT_INTERIOR  6
#define SCUFF_PFT_MOMENTS       7   // dipole-moment method
#define SCUFF_PFT_NUMMETHODS    8

#define SCUFF_PFT_DEFAULT       SCUFF_PFT_EMT_EXTERIOR

// the following are various methods for computing
// the 4-dimensional integrals needed to get scattering
// contributions to EMTPFT
#define SCUFF_EMTPFTI_EHDERIVATIVES1 0
#define SCUFF_EMTPFTI_EHDERIVATIVES2 1
#define SCUFF_EMTPFTI_EHVALUES1      2
#define SCUFF_EMTPFTI_EHVALUES2      3
#define SCUFF_EMTPFTI_NUMMETHODS     4
#define SCUFF_EMTPFTI_DEFAULT SCUFF_EMTPFTI_EHDERIVATIVES2

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct PFTOptions 
 {
   // general options
   int PFTMethod;
   char *FluxFileName;
   double *kBloch;
   HMatrix *DRMatrix;
   IncField *IF;
   HVector *RHSVector;

   // options affecting DSI PFT computation
   char *DSIMesh;
   double DSIRadius;
   int DSIPoints;
   bool DSIFarField;
   bool NeedQuantity[NUMPFT];

   // options affecting EMT PFT computation
   bool Itemize;
   bool Interior;
   int EMTPFTIMethod;
   HMatrix *TInterior, *TExterior;

   bool GetRegionPFTs;

 } PFTOptions;

/***************************************************************/
/* routine for initializing a PFTOptions structure to default  */
/* values; creates and returns a new default structure if      */
/* called with Options=NULL or with no argument                */
/***************************************************************/
PFTOptions *InitPFTOptions(PFTOptions *Options=0);

/***************************************************************/
/***************************************************************/
/***************************************************************/
class RWGGeometry;
HMatrix *GetSRFluxTrace(RWGGeometry *G, HMatrix *XMatrix, cdouble Omega,
                   HMatrix *DRMatrix, HMatrix *FMatrix=0);

void GetKNBilinears(HVector *KNVector, HMatrix *DRMatrix,
                    bool IsPECA, int KNIndexA,
                    bool IsPECB, int KNIndexB,
                    cdouble Bilinears[4]);

} // namespace scuff 

#endif // #ifndef PFTOPTIONS_H
