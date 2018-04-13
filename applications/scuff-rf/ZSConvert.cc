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
 * ZSConvert.cc -- routines for converting a matrix of multiport 
 *              -- S-parameters into Z-parameters and viceversa
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "libscuff.h"
#include <libhrutil.h>
#include <libhmat.h>

using namespace scuff;

/***************************************************************/
/* a function that converts the Z-matrix of a multiport network*/
/* into the S-matrix of that network                           */
/*                                                             */
/* formula: S = (Z-ZC)*(Z+ZC)^{-1}                             */
/*                                                             */
/*  where ZC = characteristic impedance * identity matrix      */
/***************************************************************/
void ZToS(HMatrix *ZMatrix, HMatrix *SMatrix, double ZCharacteristic)
{
  HMatrix *ZpZC=new HMatrix(ZMatrix->NR, ZMatrix->NC, LHM_COMPLEX);  
  HMatrix *ZmZC=new HMatrix(ZMatrix->NR, ZMatrix->NC, LHM_COMPLEX);  

  ZpZC->Copy(ZMatrix);
  ZmZC->Copy(ZMatrix);
  int nr;
  for(nr=0; nr<ZpZC->NR; nr++)
   { ZpZC->AddEntry(nr, nr, ZCharacteristic);
     ZmZC->AddEntry(nr, nr, -ZCharacteristic);
   };
  ZpZC->LUFactorize();
  ZpZC->LUInvert();
  ZmZC->Multiply(ZpZC, SMatrix);

  delete ZpZC;
  delete ZmZC;
    
}

void ZToS(HMatrix *ZMatrix, HMatrix *SMatrix)
{ ZToS(ZMatrix, SMatrix, 50.0); }

/***************************************************************/
/* a function that converts the S-matrix of a multiport network*/
/* into the Z-matrix of that network                           */
/*                                                             */
/* formula: Z = ZCharacteristic * (1-S)^{-1}(1+S)              */
/*                                                             */
/*  where S  = S matrix                                        */
/*        1  = identity matrix                                 */
/***************************************************************/
void SToZ(HMatrix *SMatrix, HMatrix *ZMatrix, double ZCharacteristic)
{
  // 'OmS = 'one minus S'
  HMatrix *OmS=new HMatrix(SMatrix->NR, SMatrix->NC, LHM_COMPLEX);  
  HMatrix *OpS=new HMatrix(SMatrix->NR, SMatrix->NC, LHM_COMPLEX);  

  OpS->Copy(SMatrix);
  OmS->Copy(SMatrix);
  OmS->Scale(-1.0);
  int nr;
  for(nr=0; nr<SMatrix->NR; nr++)
   { 
     OpS->AddEntry(nr, nr, 1.0);
     OmS->AddEntry(nr, nr, 1.0);
   };
  OmS->LUFactorize();
  OmS->LUInvert();
  OmS->Multiply(OpS, ZMatrix);

  ZMatrix->Scale(ZCharacteristic);

  delete OpS;
  delete OmS;

}

void SToZ(HMatrix *SMatrix, HMatrix *ZMatrix)
{ SToZ(SMatrix, ZMatrix, 50.0); }
