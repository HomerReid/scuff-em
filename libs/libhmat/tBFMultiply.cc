
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

#define NRA 2
#define NCA 5

#define NRB 5
#define NCB 3

#define NRC 2
#define NCC 3

#if defined(_WIN32)
#  define srand48 srand
#  define drand48 my_drand48
static double my_drand48(void) {
  return rand() * 1.0 / RAND_MAX;
}
#endif

int main()
{ 
  int m, n, k;
  HMatrix *A, *B, *C, *CBF;
  cdouble X;

  A=new HMatrix(NRA,NCA,LHM_COMPLEX);
  B=new HMatrix(NRB,NCB,LHM_COMPLEX);
  C=new HMatrix(NRC,NCC,LHM_COMPLEX);
  CBF=new HMatrix(NRC,NCC,LHM_COMPLEX);

  srand48(time(0));
  for(m=0; m<NRA; m++)
   for(n=0; n<NCA; n++)
     A->SetEntry(m,n, cdouble(5.0*(drand48()-2.5), 5.0*(drand48()-2.5) ));
  for(m=0; m<NRB; m++)
   for(n=0; n<NCB; n++)
     B->SetEntry(m,n, cdouble(5.0*(drand48()-2.5), 5.0*(drand48()-2.5) ));
 
  A->Multiply(B,C);

  if ( NCA!=NRB || NRA!=NRC || NCB!=NCC )
   ErrExit("badcallatage");

  for(m=0; m<NRC; m++)
   for(n=0; n<NCC; n++)
    { X=0.0;
      for(k=0; k<NCA; k++)
       X+= A->GetEntry(m,k) * B->GetEntry(k,n);
      CBF->SetEntry(m,n,X);
    };

  printf("\n\n* libhmat: \n*\n");
  SetDefaultCD2SFormat("(%+4.2e,%+4.2e)");
  for(m=0; m<NRC; m++)
   for(n=0; n<NCC; n++)
    printf("%s%c",CD2S(C->GetEntry(m,n)), n==(NCC-1) ? '\n' : ' ');

  printf("\n\n* BF: \n*\n");
  SetDefaultCD2SFormat("(%+4.2e,%+4.2e)");
  for(m=0; m<NRC; m++)
   for(n=0; n<NCC; n++)
    printf("%s%c",CD2S(CBF->GetEntry(m,n)), n==(NCC-1) ? '\n' : ' ');


}
