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

namespace scuff {

/***************************************************************/
/* number of power, force, and torque quantities, plus indices */
/* into arrays of them                                         */
/***************************************************************/
#define NUMPFT        8
#define SCUFF_PABS    0  // absorbed power
#define SCUFF_PSCAT   1  // scattered power
#define SCUFF_XFORCE  2  // force components
#define SCUFF_YFORCE  3  // 
#define SCUFF_ZFORCE  4  // 
#define SCUFF_XTORQUE 5  // torque components
#define SCUFF_YTORQUE 6  // 
#define SCUFF_ZTORQUE 7  // 

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

#define SCUFF_PFT_DEFAULT       SCUFF_PFT_EPOVERLAP

/***************************************************************/
/***************************************************************/
/***************************************************************/
typedef struct PFTOptions 
 {
   // general options
   int PFTMethod;
   char *FluxFileName;
   HMatrix *RytovMatrix;
   double *kBloch;

   // options affecting overlap PFT computation
   HVector *RHSVector;

   // options affecting DSI PFT computation
   char *DSIMesh;
   double DSIRadius;
   int DSIPoints;
   bool DSIFarField;
   bool NeedQuantity[NUMPFT];

   // options affecting EP power computation
   HMatrix *TInterior, *TExterior;

   // options affecting EP force/torque computation
   int EPFTOrder;
   double EPFTDelta;

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
HMatrix *GetSRFlux(RWGGeometry *G, HMatrix *XMatrix, cdouble Omega,
                   HVector *KNVector, HMatrix *RytovMatrix,
                   HMatrix *FMatrix=0, bool FarField=false);

} // namespace scuff 

#endif // #ifndef PFTOPTIONS_H
