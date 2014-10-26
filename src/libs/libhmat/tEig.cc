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
 * tEig.cc 
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <libhrutil.h>
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
  /*- process options  -------------------------------------------*/
  /*--------------------------------------------------------------*/
  int N=1000;
  bool Complex=false;
  bool NS=false;
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { {"N",       PA_INT,     1, 1, (void *)&N,       0, "dimension "},
     {"Complex", PA_BOOL,    0, 1, (void *)&Complex, 0, "complex-valued matrix"},
     {"NS",      PA_BOOL,    0, 1, (void *)&NS,      0, "non-symmetric matrix"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  printf("Creating a random %ix%i matrix ... \n",N,N);
  HMatrix *M=new HMatrix(N, N, Complex ? LHM_COMPLEX : LHM_REAL);

  cdouble II = Complex ? cdouble(0.0,1.0) : cdouble(0.0,0.0);
  srand48(time(0));

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (NS==true)
   { for(int m=0; m<N; m++)
      for(int n=0; n<N; n++)
       M->SetEntry(m, n, drand48() + II*drand48());
   }
  else
   { for(int m=0; m<N; m++)
      { M->SetEntry(m, m, drand48());
        for(int n=m+1; n<N; n++)
         { cdouble Entry = drand48() + II*drand48();
           M->SetEntry(m, n, Entry);
           M->SetEntry(n, m, conj(Entry));
         };
      };
   };

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  SetDefaultCD2SFormat("{%+10.3e , %+10.3e}\n");
  printf("(1,1)   = %s \n",CD2S(M->GetEntry(0,0)));
  printf("(%i,1)  = %s \n",N,CD2S(M->GetEntry(N-1,0)));
  printf("(1,%i)  = %s \n",N,CD2S(M->GetEntry(0,N-1)));
  printf("(%i,%i) = %s \n",N,N,CD2S(M->GetEntry(N-1,N-1)));
  M->ExportToHDF5("tEig.hdf5","M");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double Elapsed;
  Tic();
  HVector *Lambda = NS ? M->NSEig() : M->Eig();
  Elapsed=Toc();
  printf("...%.3f s\n",Elapsed);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int n=0; n<N; n++)
   printf("%i %s \n",n,CD2S(Lambda->GetEntry(n)));

}
