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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include <libhrutil.h>
#include "libhmat.h"

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool Equal(SMatrix *S1, SMatrix *S2)
{
  if (    (S1->NR!=S2->NR) 
       || (S1->NC!=S2->NC) 
       || (S1->nnz!=S2->nnz) 
       || (S1->RealComplex!=S2->RealComplex) 
     )
   { printf("Basic statistics disagree!\n");
     printf(" S1=(%i,%i,%i,%i)\n",S1->NR,S1->NC,S1->nnz,S1->RealComplex);
     printf(" S2=(%i,%i,%i,%i)\n",S2->NR,S2->NC,S2->nnz,S2->RealComplex);
     return false;
   }

  for(int nr=0; nr<S1->NR; nr++)
   if( S1->RowStart[nr] != S2->RowStart[nr])
    { printf("RowStart[%i]: (S1,S2) = (%i,%i)\n",nr,S1->RowStart[nr],S2->RowStart[nr]);
      return false;
    }

  for(int nnz=0; nnz<S1->nnz; nnz++)
   { 
     if( S1->ColIndices[nnz] != S2->ColIndices[nnz])
      { printf("ColIndices[%i]: (S1,S2) = (%i,%i)\n",nnz,S1->ColIndices[nnz],S2->ColIndices[nnz]);
        return false;
      };
   };

  if ( S1->RealComplex==LHM_REAL )
   { for(int nnz=0; nnz<S1->nnz; nnz++)
      if( !EqualFloat(S1->DM[nnz],S2->DM[nnz]))
       { printf("DM[%i]: (S1,S2) = (%e,%e)\n",nnz,S1->DM[nnz],S2->DM[nnz]);
         return false;
       }
   }
  else
   { for(int nnz=0; nnz<S1->nnz; nnz++)
      if( !EqualFloat(S1->ZM[nnz],S2->ZM[nnz]))
       { printf("ZM[%i]: (S1,S2) = (%s,%s)\n",nnz,CD2S(S1->DM[nnz]),CD2S(S2->DM[nnz]));
         return false;
       };
   };

  printf("Matrices agree!\n");
  return true;
}

/***************************************************************/ 
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{ 
  /*--------------------------------------------------------------*/
  /*- process command-line arguments -----------------------------*/
  /*--------------------------------------------------------------*/
  char *GeoFile;
  int NR=10;
  int NC=12;
  double pnnz=0.2; // probability of nonzero entry
  bool Complex=false;
  /* name, type, #args, max_instances, storage, count, description*/
  OptStruct OSArray[]=
   { {"NR",       PA_INT,     1, 1, (void *)&NR,       0,  "NR"},
     {"NC",       PA_INT,     1, 1, (void *)&NC,       0,  "NC"},
     {"pnnz",     PA_DOUBLE,  1, 1, (void *)&pnnz,     0,  "pnnz"},
     {"Complex",  PA_BOOL,    0, 1, (void *)&Complex,  0,  "Complex"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  SMatrix *S1=new SMatrix(NR, NC, Complex ? LHM_COMPLEX : LHM_REAL);
  S1->BeginAssembly( (int)(ceil(pnnz*NR*NC)) );
  for(int nr=0; nr<NR; nr++)
   for(int nc=0; nc<NC; nc++)
    if( drand48() < pnnz )
     S1->SetEntry(nr,nc,drand48());
  S1->EndAssembly();
  
  void *pHC=HMatrix::OpenHDF5Context("tSMatrix.hdf5");
  S1->ExportToHDF5(pHC, "S1");
  HMatrix::CloseHDF5Context(pHC);

  SMatrix *S2=new SMatrix("tSMatrix.hdf5","S1");
  if (S2->ErrMsg)
   ErrExit(S2->ErrMsg);

  return Equal(S1, S2);

}
