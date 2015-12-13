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
 * GetLDOS.cc -- routine to compute the LDOS at a single 
 *               (Omega, kBloch) point
 *
 * homer reid -- 3/2015
 */
#include <libhrutil.h>
#include <libhmat.h>
#include "libscuff.h"
#include "scuff-ldos.h"

#define ABSTOL 1.0e-10

using namespace scuff;

/***************************************************************/
/* routine to compute the LDOS at a single (Omega, kBloch)     */
/* point (but typically multiple spatial evaluation points)    */
/***************************************************************/
void GetLDOS(void *pData, cdouble Omega, double *kBloch, 
             double *Result)
{
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  SLDData *Data        = (SLDData *)pData;
  RWGGeometry *G       = Data->G;
  HMatrix *M           = Data->M;
  HMatrix *XMatrix     = Data->XMatrix;
  HMatrix *GMatrix     = Data->GMatrix;
  void **ABMBCache     = Data->ABMBCache;
  MatProp *HalfSpaceMP = Data->HalfSpaceMP;

  int LDim = G->LDim;
  double *LBasis[2]={0,0};
  if (LDim>=1) LBasis[0]=G->LBasis[0];
  if (LDim>=2) LBasis[1]=G->LBasis[1];
  double LBV[2][2]={{0.0,0.0},{0.0,0.0}};
  if (LDim>=1)
   { LBV[0][0]=LBasis[0][0]; LBV[0][1]=LBasis[0][1]; }
  if (LDim>=2)
   { LBV[1][0]=LBasis[1][0]; LBV[1][1]=LBasis[1][1]; }

  switch(LDim)
   { case 0: Log("Computing LDOS at Omega=%s",z2s(Omega)); break; 
     case 1: Log("Computing LDOS at (Omega,kx)=(%s,%e)",z2s(Omega),kBloch[0]);
     case 2: Log("Computing LDOS at (Omega,kx,ky)=(%s,%e),",z2s(Omega),kBloch[0],kBloch[1]);
   };

  /*--------------------------------------------------------------*/
  /*- assemble the BEM matrix at this frequency and Bloch vector, */
  /*- then get DGFs at all evaluation points                      */
  /*--------------------------------------------------------------*/
  if (HalfSpaceMP==0)
   { int NS = G->NumSurfaces;
     for(int ns=0, nb=0; ns<NS; ns++)
      for(int nsp=ns; nsp<NS; nsp++, nb++)
       { 
         int RowOffset = G->BFIndexOffset[ns];
         int ColOffset = G->BFIndexOffset[nsp];
         G->AssembleBEMMatrixBlock(ns, nsp, Omega, kBloch,
                                   M, 0, RowOffset, ColOffset,
                                   ABMBCache[nb], false);

         if (nsp>ns)
          G->AssembleBEMMatrixBlock(nsp, ns, Omega, kBloch,
                                    M, 0, ColOffset, RowOffset,
                                    ABMBCache[nb], true);
       };
     M->LUFactorize();

     G->GetDyadicGFs(Omega, kBloch, XMatrix, M, GMatrix);

   };
                               
  /*--------------------------------------------------------------*/
  /*- get LDOS at all evaluation points.                          */
  /*- Note: The LDOS is defined as                                */
  /*-  \Rho = (\omega / \pi c^2) * Im Tr G                        */
  /*- where G is the dyadic GF. OTOH, the vacuum LDOS is          */
  /*-  \Rho_0 = \omega^3 / (\pi^2 c^3).                           */
  /*- Thus the LDOS enhancement is                                */
  /*-  \Rho/Rho_0 = Im Tr G / ( \pi \omega^2 c )                  */
  /*- which explains the prefactor below.                         */
  /*--------------------------------------------------------------*/
  FILE *DataFile=0;
  if (Data->ByKFileName)
   DataFile=fopen(Data->ByKFileName, "a");
  double PreFac = 1.0 / (M_PI * real(Omega) * real(Omega) );
  for(int nx=0; nx<XMatrix->NR; nx++)
   { 
     double X[3];
     XMatrix->GetEntriesD(nx, ":", X);
  
     /***************************************************************/
     /* get the DGFs at this evaluation point                       */
     /***************************************************************/
     cdouble GE[3][3], GM[3][3];
     if ( HalfSpaceMP && HalfSpaceMP->IsPEC())
      GetGroundPlaneDGFs(X, Omega, kBloch, LDim, LBasis, GE, GM);
     else if (HalfSpaceMP)
      GetHalfSpaceDGFs(Omega, kBloch, X[2], LBV, HalfSpaceMP,
                       Data->RelTol, ABSTOL, Data->MaxEvals, GE, GM);
     else
      for(int i=0; i<3; i++)
       for(int j=0; j<3; j++)
        { GE[i][j] = GMatrix->GetEntry(nx, 0 + 3*i+j);
          GM[i][j] = GMatrix->GetEntry(nx, 9 + 3*i+j);
        };

     /***************************************************************/
     /* figure out what to return                                   */
     /***************************************************************/
     if (Data->LDOSOnly)
      {
        Result[2*nx+0] = PreFac * imag(GE[0][0] + GE[1][1] + GE[2][2]);
        Result[2*nx+1] = PreFac * imag(GM[0][0] + GM[1][1] + GM[2][2]);
      }
     else
      { int NFun=20, nf=0;
        Result[NFun*nx + nf++] = PreFac * imag(GE[0][0] + GE[1][1] + GE[2][2]);
        Result[NFun*nx + nf++] = PreFac * imag(GM[0][0] + GM[1][1] + GM[2][2]);
        for(int Mu=0; Mu<3; Mu++)
         for(int Nu=Mu; Nu<3; Nu++)
          { Result[NFun*nx + nf++] = real(GE[Mu][Nu]);
            Result[NFun*nx + nf++] = imag(GE[Mu][Nu]);
          };
        for(int Mu=0; Mu<3; Mu++)
         { int Nu=(Mu+1)%3;
           Result[NFun*nx + nf++] = real(GM[Mu][Nu]);
           Result[NFun*nx + nf++] = imag(GM[Mu][Nu]);
         };
      };

     /***************************************************************/
     /* write output to data file if we have one ********************/
     /***************************************************************/
     if (DataFile)
      { 
        fprintf(DataFile,"%e %e %e %s ",X[0],X[1],X[2],z2s(Omega));
        if (G->LDim>=1) fprintf(DataFile,"%e ",kBloch[0]);
        if (G->LDim>=2) fprintf(DataFile,"%e ",kBloch[1]);
        fprintf(DataFile,"%e %e ",Result[2*nx+0],Result[2*nx+1]);

        fprintf(DataFile,"%s %s %s %s %s %s ",
                         CD2S(GE[0][0]),CD2S(GE[0][1]),CD2S(GE[0][2]),
                         CD2S(GE[1][1]),CD2S(GE[1][2]),CD2S(GE[2][2]));
        fprintf(DataFile,"%s %s %s %s %s %s ",
                         CD2S(GM[0][0]),CD2S(GM[0][1]),CD2S(GM[0][2]),
                         CD2S(GM[1][1]),CD2S(GM[1][2]),CD2S(GM[2][2]));

        fprintf(DataFile,"\n");
        fflush(DataFile);
      };

   }; // for(int nr=0, nr<XMatrix->NR; nr++)

  if (DataFile) fclose(DataFile);

}
