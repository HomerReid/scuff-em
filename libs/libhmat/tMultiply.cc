
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

#if defined(_WIN32)
#  define srand48 srand
#  define drand48 my_drand48
static double my_drand48(void) {
  return rand() * 1.0 / RAND_MAX;
}
#endif

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{ 
  /*--------------------------------------------------------------*/
  /*- process command-line arguments -----------------------------*/
  /*--------------------------------------------------------------*/
  int N1=100;
  int N2=200;
  ArgStruct ASArray[]=
   { {"N1",         PA_INT,  (void *)&N1,        "100", "dimension 1"},
     {"N2",         PA_INT,  (void *)&N2,        "200", "dimension 2"},
     {0,0,0,0,0}
   };
  ProcessArguments(argc, argv, ASArray);

  HMatrix *A=new HMatrix(N2,N1,LHM_COMPLEX);
  HMatrix *B=new HMatrix(N2,N1,LHM_COMPLEX);

  srand48(time(0));
  for(int m=0; m<A->NR; m++)
   for(int n=0; n<A->NC; n++)
     A->SetEntry(m,n, cdouble(5.0*(drand48()-2.5), 5.0*(drand48()-2.5) ));
  for(int m=0; m<B->NR; m++)
   for(int n=0; n<B->NC; n++)
     B->SetEntry(m,n, cdouble(5.0*(drand48()-2.5), 5.0*(drand48()-2.5) ));
 
  /***************************************************************/
  /* test multiplication with the first factor adjointed *********/
  /***************************************************************/
  HMatrix *C1=new HMatrix(N1,N1,LHM_COMPLEX);
  Tic();
  A->Multiply(B,C1,"--transA C");
  printf("--transA C: %e s\n",Toc());

  HMatrix *C2=new HMatrix(N1,N1,LHM_COMPLEX);
  Tic();
  A->Adjoint();
  A->Multiply(B,C2);
  A->Adjoint(); // undo the adjointing
  printf("--notrans:  %e s\n",Toc());

  double NormC=0.0, NormError=0.0;
  for(int nr=0; nr<C1->NR; nr++)
   for(int nc=0; nc<C1->NC; nc++)
    { NormC     += norm(C1->GetEntry(nr,nc));
      NormError += norm(C1->GetEntry(nr,nc) - C2->GetEntry(nr,nc));
    };
  printf("rel norm error=%e\n",sqrt(NormError/NormC));

  delete C1; 
  delete C2; 
 
  /***************************************************************/
  /* test multiplication with the second factor adjointed ********/
  /***************************************************************/
  C1=new HMatrix(N2,N2,LHM_COMPLEX);
  Tic();
  A->Multiply(B,C1,"--transB C");
  printf("--transB C: %e s\n",Toc());

  C2=new HMatrix(N2,N2,LHM_COMPLEX);
  Tic();
  B->Adjoint();
  A->Multiply(B,C2);
  B->Adjoint(); // undo the adjointing
  printf("--notrans:  %e s\n",Toc());

  NormC=0.0, NormError=0.0;
  for(int nr=0; nr<C1->NR; nr++)
   for(int nc=0; nc<C1->NC; nc++)
    { NormC     += norm(C1->GetEntry(nr,nc));
      NormError += norm(C1->GetEntry(nr,nc) - C2->GetEntry(nr,nc));
    };
  printf("rel norm error=%e\n",sqrt(NormError/NormC));

  delete C1; 
  delete C2; 
}
