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

#include "libhrutil.h"
#include "libhmat.h"

#define DIM 10
#define II cdouble(0.0,1.0)

int main()
{ 
  HMatrix *M1, *M2;
  HVector *V1, *V2;
  double Norm, DNorm;
  int nr, nc;

  srand48(time(0));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  M1=new HMatrix(DIM,DIM,LHM_COMPLEX,LHM_SYMMETRIC);
  for(nr=0; nr<DIM; nr++) 
   for(nc=nr; nc<DIM; nc++) 
    M1->SetEntry(nr,nc,drand48() + II*drand48());

  V1=new HVector(DIM,LHM_COMPLEX);
  for(nr=0; nr<DIM; nr++) 
   V1->SetEntry(nr,drand48() + II*drand48());

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  void *hPC=HMatrix::OpenHDF5Context("MV.hdf5");
  M1->ExportToHDF5(hPC,"M1");
  V1->ExportToHDF5(hPC,"V1");
  HMatrix::CloseHDF5Context(hPC);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  M2=new HMatrix("MV.hdf5",LHM_HDF5, "M1");
  if (M2->ErrMsg)
   printf("Error 1: %s \n",M2->ErrMsg);

  V2=new HVector("MV.hdf5",LHM_HDF5, "V2");
  if (V2->ErrMsg)
   printf("This should say error: %s \n",V2->ErrMsg);

  V2=new HVector("MV.hdf5",LHM_HDF5,"V1");

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  printf("Checking matrix entries...\n");
  for(nr=0; nr<DIM; nr++)
   for(nc=0; nc<DIM; nc++)
    if ( norm( M1->GetEntry(nr,nc) - M2->GetEntry(nr,nc)) > 1.0e-10 )
     printf("whoops! M1(%i,%i) = %s != %s = M2(%i,%i)\n",
             nr,nc,CD2S(M1->GetEntry(nr,nc)),CD2S(M2->GetEntry(nr,nc)),nr,nc);

  printf("\nChecking vector entries...\n");
  for(nr=0; nr<DIM; nr++)
   if ( norm( V1->GetEntry(nr) - V2->GetEntry(nr)) > 1.0e-10 )
    printf("whoops! V1(%i) = %s != %s = V2(%i)\n",
             nr,CD2S(V1->GetEntry(nr)),CD2S(V2->GetEntry(nr)),nr);

}
