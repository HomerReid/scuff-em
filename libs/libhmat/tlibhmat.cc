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
#define II cdouble (0.0,1.0)

int main()
{ 
  HMatrix *MNorm, *MHerm, *MCSym;
  HVector *B, *XNormal, *XPacked;
  int nr, nc, n, info;
  double Norm, DNorm;
  cdouble Z;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  void *pCC = HMatrix::OpenMATLABContext("tlibhmat");

  MNorm=new HMatrix(DIM,DIM,LHM_COMPLEX);
  MHerm=new HMatrix(DIM,DIM,LHM_COMPLEX,LHM_HERMITIAN);
  MCSym=new HMatrix(DIM,DIM,LHM_COMPLEX,LHM_SYMMETRIC);

  B=new HVector(DIM,LHM_COMPLEX);
  XNormal=new HVector(DIM,LHM_COMPLEX);
  XPacked=new HVector(DIM,LHM_COMPLEX);

  srand48(time(0));

  /* create a random complex vector */
  for(n=0; n<DIM; n++)
   B->SetEntry(n, 2.0*(drand48()-0.5) + 2.0*II*(drand48()-0.5) );

  /*--------------------------------------------------------------*/
  /*- phase 1: test lu-factorization of packed hermitian matrix --*/
  /*--------------------------------------------------------------*/
 
  /* create a random hermitian matrix, stored using both */
  /* normal storage and packed storage                   */
  for(nr=0; nr<DIM; nr++)
   for(nc=nr; nc<DIM; nc++)
    { 
      Z = 2.0*(drand48()-0.5) + 2.0*II*(drand48()-0.5);

      if ( nr==nc )
       imag(Z)=0.0;

       MNorm->SetEntry(nr,nc,Z);
       MNorm->SetEntry(nc,nr,conj(Z));

       MHerm->SetEntry(nr,nc,Z);  
    };
  MNorm->ExportToMATLAB(pCC,"M1_Normal");
  MHerm->ExportToMATLAB(pCC,"M1_Herm");

  /* solve the system using the normal and packed        */
  /* matrices and compare the result                     */
  info=MNorm->LUFactorize();
  printf("LU(Normal 1) : info = %i \n",info); 
  memcpy(XNormal->ZV,B->ZV,DIM*sizeof(double));
  MNorm->LUSolve(XNormal);

  info=MHerm->LUFactorize();
  printf("LU(Herm   1) : info = %i \n",info); 
  memcpy(XPacked->ZV,B->ZV,DIM*sizeof(double));
  MHerm->LUSolve(XPacked);

  for(Norm=DNorm=0.0, n=0; n<DIM; n++)
   { Norm += norm( XNormal->ZV[n] );
     DNorm += norm( XNormal->ZV[n] - XPacked->ZV[n] );
   };
  Norm=sqrt(Norm) / DIM;
  DNorm=sqrt(DNorm) / DIM;
  printf("Hermitian:  D = %7.3e (RD=%7.3e)\n",DNorm,DNorm/Norm);
 

  /*--------------------------------------------------------------*/
  /*- phase 2: test lu-factorization of complex-symmetric matrix  */
  /*--------------------------------------------------------------*/
 
  /* create a random complex-symmetric matrix, stored using both */
  /* normal storage and packed storage                           */
  for(nr=0; nr<DIM; nr++)
   for(nc=nr; nc<DIM; nc++)
    { 
       Z=2.0*(drand48()-0.5) + 2.0*II*(drand48()-0.5);
 
       MNorm->SetEntry(nr,nc,Z);
       MNorm->SetEntry(nc,nr,Z);

       MCSym->SetEntry(nr,nc,Z);  
    };
  MNorm->ExportToMATLAB(pCC,"M2_Normal");
  MCSym->ExportToMATLAB(pCC,"M2_CSym");

  /* solve the system using the normal and packed        */
  /* matrices and compare the result                     */
  info=MNorm->LUFactorize();
  printf("LU(Normal 2) : info = %i \n",info); 
  memcpy(XNormal->ZV,B->ZV,DIM*sizeof(double));
  MNorm->LUSolve(XNormal);

  info=MCSym->LUFactorize();
  printf("LU(CSym   2) : info = %i \n",info); 
  memcpy(XPacked->ZV,B->ZV,DIM*sizeof(double));
  MCSym->LUSolve(XPacked);

  for(Norm=DNorm=0.0, n=0; n<DIM; n++)
   { Norm += norm( XNormal->ZV[n] );
     DNorm += norm( XNormal->ZV[n] - XPacked->ZV[n] );
   };
  Norm=sqrt(Norm) / DIM;
  DNorm=sqrt(DNorm) / DIM;
  printf("CSymmetric: D = %7.3e (RD=%7.3e)\n",DNorm,DNorm/Norm);

  /*--------------------------------------------------------------*/
  /*- phase 3: test cholesky factorization of packed hermitian matrix */
  /*--------------------------------------------------------------*/
  /* create a random hermitian matrix, stored using both */
  /* normal storage and packed storage                   */
  for(nr=0; nr<DIM; nr++)
   for(nc=nr; nc<DIM; nc++)
    { 
      Z = 2.0*(drand48()-0.5) + 2.0*II*(drand48()-0.5);

      if ( nr==nc )
       { imag(Z)=0.0;
         real(Z)+=2.0*((double)DIM); /* to make it positive-definite */
       };

       MNorm->SetEntry(nr,nc,Z);
       MNorm->SetEntry(nc,nr,conj(Z));

       MHerm->SetEntry(nr,nc,Z);  
    };
  MNorm->ExportToMATLAB(pCC,"M3_Normal");
  MHerm->ExportToMATLAB(pCC,"M3_Herm");

  /* solve the system using the normal and packed        */
  /* matrices and compare the result                     */
  info=MNorm->CholFactorize();
  printf("Chol(Normal 3) : info = %i \n",info); 
  memcpy(XNormal->ZV,B->ZV,DIM*sizeof(double));
  MNorm->CholSolve(XNormal);

  info=MHerm->CholFactorize();
  printf("Chol(Herm   3) : info = %i \n",info); 
  memcpy(XPacked->ZV,B->ZV,DIM*sizeof(double));
  MHerm->CholSolve(XPacked);

  for(Norm=DNorm=0.0, n=0; n<DIM; n++)
   { Norm += norm( XNormal->ZV[n] );
     DNorm += norm( XNormal->ZV[n] - XPacked->ZV[n] );
   };
  Norm=sqrt(Norm) / DIM;
  DNorm=sqrt(DNorm) / DIM;
  printf("Hermitian:  D = %7.3e (RD=%7.3e)\n",DNorm,DNorm/Norm);
 
  /*--------------------------------------------------------------*/
  /*- phase 2: test cholesky-factorization of complex-symmetric  -*/
  /*-          (which must fail)                                  */
  /*--------------------------------------------------------------*/
  /* create a random complex-symmetric matrix, stored using both */
  /* normal storage and packed storage                           */
  for(nr=0; nr<DIM; nr++)
   for(nc=nr; nc<DIM; nc++)
    { 
       Z=2.0*(drand48()-0.5) + 2.0*II*(drand48()-0.5);

       if (nr==nc)
        real(Z)+=2.0*((double)DIM); /* to make it positive-definite */
 
       MNorm->SetEntry(nr,nc,Z);
       MNorm->SetEntry(nc,nr,Z);

       MCSym->SetEntry(nr,nc,Z);  
    };
  MNorm->ExportToMATLAB(pCC,"M4_Normal");
  MCSym->ExportToMATLAB(pCC,"M4_CSym");

  /* solve the system using the normal and packed        */
  /* matrices and compare the result                     */
  info=MNorm->CholFactorize();
  printf("Chol(Normal 4) : info = %i \n",info); 
  MNorm->ExportToMATLAB(pCC,"M4_Normal_Chol");
  memcpy(XNormal->ZV,B->ZV,DIM*sizeof(double));
  MNorm->CholSolve(XNormal);

  info=MCSym->CholFactorize();
  printf("Chol(CSym   4) : info = %i \n",info); 
  memcpy(XPacked->ZV,B->ZV,DIM*sizeof(double));
  MCSym->CholSolve(XPacked);

  for(Norm=DNorm=0.0, n=0; n<DIM; n++)
   { Norm += norm( XNormal->ZV[n] );
     DNorm += norm( XNormal->ZV[n] - XPacked->ZV[n] );
   };
  Norm=sqrt(Norm) / DIM;
  DNorm=sqrt(DNorm) / DIM;
  printf("CSymmetric: D = %7.3e (RD=%7.3e)\n",DNorm,DNorm/Norm);

  HMatrix::CloseMATLABContext(pCC);
  printf("\nThank you for your support.\n");
}
