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
 * CasimirIntegrand.cc    -- evaluate the casimir energy, force, and/or torque
 *                        -- integrand at a single Xi point or a single (Xi,kBloch)
 *                        -- point.
 *
 * homer reid  -- 10/2008 -- 6/2010
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdarg.h>
#include <unistd.h>
#include <ctype.h>

#include "scuff-cas3D.h"
#include "libscuffInternals.h"

extern "C" {
int dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
int zgetrf_(int *m, int *n, cdouble *a, int *lda, int *ipiv, int *info);
}

using namespace scuff;

#define II cdouble(0.0,1.0)

/***************************************************************/
/* compute \log \det \{ M^{-1} MInfinity \}                    */
/***************************************************************/
double GetLNDetMInvMInf(SC3Data *SC3D)
{ 
  HMatrix *M              = SC3D->M;
  int N                   = SC3D->N;

  double LNDet=0.0;
  if (SC3D->NewEnergyMethod==false)
   {  
     /*--------------------------------------------------------------*/
     /*- calculation method 1  --------------------------------------*/
     /*--------------------------------------------------------------*/
     HVector *MInfLUDiagonal = SC3D->MInfLUDiagonal;

     for(int n=0; n<N; n++)
      LNDet+=log( abs( MInfLUDiagonal->GetEntryD(n) / M->GetEntry(n,n) ) );
   }
  else
   {
     /*--------------------------------------------------------------*/
     /*- calculation method 2 ---------------------------------------*/
     /*--------------------------------------------------------------*/
     Log("Computing Casimir energy by new method...");
     RWGGeometry *G = SC3D->G;
     HMatrix *MM1MInf = SC3D->MM1MInf;
     MM1MInf->Zero();
     for(int ns=0; ns<G->NumSurfaces; ns++)
      MM1MInf->InsertBlock(SC3D->TBlocks[ns], G->BFIndexOffset[ns], G->BFIndexOffset[ns]);
     M->LUSolve(MM1MInf); 
     MM1MInf->LUFactorize();

     for(int n=0; n<N; n++)
      LNDet+=log( abs( MM1MInf->GetEntry(n,n) ) );
   }

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  // paraphrasing the physicists of the 1930s, 'just because
  // something is infinite doesn't mean that it's zero.' and yet...
  if (!IsFinite(LNDet))
   LNDet=0.0;
  return -LNDet/(2.0*M_PI);

} 

/***************************************************************/
/* compute \trace \{ M^{-1} dMdAlpha\},                        */
/* where Alpha=x, y, z                                         */
/***************************************************************/
double GetTraceMInvdM(SC3Data *SC3D, char XYZT)
{ 
  /***************************************************************/
  /* unpack fields from workspace structure **********************/
  /***************************************************************/
  RWGGeometry *G     = SC3D->G;
  HMatrix *M         = SC3D->M;
  HMatrix *dM        = SC3D->dM;
  HMatrix **dUBlocks = SC3D->dUBlocks;

  int Mu=0;
  switch(XYZT)
   {  case 'X': Mu   = 0; break;
      case 'Y': Mu   = 1; break;
      case 'Z': Mu   = 2; break;
      case '1': Mu   = 3; break; 
      case '2': Mu   = 4; break; 
      case '3': Mu   = 5; break; 
   };

  if ( '1'<=XYZT && XYZT<='3' )
   Log("  Computing torque about axis #%c...",XYZT);
  else
   Log("  Computing %cForce...",XYZT);

  /***************************************************************/
  /* stamp derivative blocks into dM matrix, compute M^{-1} dM,  */
  /* then sum the diagonals of the upper matrix block            */
  /***************************************************************/
  dM->Zero();
  for(int ns=1; ns<G->NumSurfaces; ns++)
   dM->InsertBlockAdjoint(dUBlocks[ 6*(ns-1) + Mu ], G->BFIndexOffset[ns], 0);

  M->LUSolve(dM);

  double Trace=0.0;
  for(int n=0; n<dM->NC; n++)
   Trace+=dM->GetEntryD(n,n);
  Trace*=2.0;

  // paraphrasing the physicists of the 1930s, 'just because
  // something is infinite doesn't mean that it's zero.' and yet...
  if (!IsFinite(Trace))
   Trace=0.0;

  return -Trace/(2.0*M_PI);
} 


/***************************************************************/
/* stamp T and U blocks into the BEM matrix, then LU-factorize.*/
/***************************************************************/
void Factorize(SC3Data *SC3D)
{ 
  RWGGeometry *G = SC3D->G;
  HMatrix *M     = SC3D->M;

  /***************************************************************/
  /* stamp blocks into M matrix                                  */
  /***************************************************************/
  /* T blocks */
  int Offset;
  for(int ns=0; ns<G->NumSurfaces; ns++)
   { 
     Offset=G->BFIndexOffset[ns];
     M->InsertBlock(SC3D->TBlocks[ns], Offset, Offset);
   };

  /* U blocks */
  int RowOffset, ColOffset;
  for(int nb=0, ns=0; ns<G->NumSurfaces; ns++)
   for(int nsp=ns+1; nsp<G->NumSurfaces; nsp++, nb++)
    { 
      RowOffset=G->BFIndexOffset[ns];
      ColOffset=G->BFIndexOffset[nsp];
      M->InsertBlock(SC3D->UBlocks[nb], RowOffset, ColOffset);
      M->InsertBlockAdjoint(SC3D->UBlocks[nb], ColOffset, RowOffset);
    };

  /***************************************************************/
  /* LU factorize                                                */
  /***************************************************************/
  M->LUFactorize();

} 

/***************************************************************/
/* evaluate the casimir energy, force, and/or torque integrand */
/* at a single Xi point, or a single (Xi,kBloch) point for PBC */
/* geometries, possibly under multiple spatial transformations.*/
/*                                                             */
/* the output vector 'EFT' stands for 'energy, force, and      */
/* torque.' the output quantities are packed into this vector  */
/* in the following order:                                     */
/*                                                             */
/*  energy  integrand (spatial transformation 1)               */
/*  xforce  integrand (spatial transformation 1)               */
/*  yforce  integrand (spatial transformation 1)               */
/*  zforce  integrand (spatial transformation 1)               */
/*  torque1 integrand (spatial transformation 1)               */
/*  torque2 integrand (spatial transformation 1)               */
/*  torque3 integrand (spatial transformation 1)               */
/*  energy  integrand (spatial transformation 2)               */
/* ...                                                         */
/*                                                             */
/* where only requested quantities are present, i.e. if the    */
/* user requested only energy and yforce then we would have    */
/*                                                             */
/*  EFT[0] = energy integrand (spatial transformation 1)       */
/*  EFT[1] = yforce integrand (spatial transformation 1)       */
/*  EFT[2] = energy integrand (spatial transformation 2)       */
/*  EFT[3] = yforce integrand (spatial transformation 2)       */
/*  ...                                                        */
/*  EFT[2*N-1] = yforce integrand (spatial transformation N)   */
/***************************************************************/
void GetCasimirIntegrand(SC3Data *SC3D, double Xi, 
                         double *kBloch, double *EFT)
{ 
  /***************************************************************/
  /* attempt to bypass the calculation by reading data from a    */
  /* cache file                                                  */
  /***************************************************************/
  if ( CacheRead(SC3D, Xi, kBloch, EFT) )
   return;

  RWGGeometry *G = SC3D->G;
  bool PBC = (G->LDim > 0);

  if (G->LDim==0)
   Log("Computing Casimir integrand at Xi=%g\n",Xi);
  else if (G->LDim==1)
   Log("Computing Casimir integrand at (Xi,kx)=(%g,%g)\n",Xi,kBloch[0]);
  else // (G->LDim==2)
   Log("Computing Casimir integrand at (Xi,kx,ky)=(%g,%g,%g)\n",Xi,kBloch[0],kBloch[1]);

  /***************************************************************/
  /* SurfaceNeverMoved[ns] is initialized true and remains true  */
  /* as long as surface #ns is not displaced or rotated by any   */
  /* geometrical transformation, but switches to false as soon as*/
  /* that surface is moved, and then remains false for the rest  */
  /* of the calculations done at this frequency                  */
  /***************************************************************/
  static bool *SurfaceNeverMoved=0;
  if (SurfaceNeverMoved==0)
   SurfaceNeverMoved=(bool *)mallocEC(G->NumSurfaces*sizeof(bool));
  for(int ns=0; ns<G->NumSurfaces; ns++)
   SurfaceNeverMoved[ns]=true;

  /***************************************************************/
  /* assemble T matrices                                         */
  /***************************************************************/
  cdouble Omega = cdouble(0.0, Xi);
  for(int ns=0; ns<G->NumSurfaces; ns++)
   { 
     /* skip if this surface is identical to a previous surface */
     int nsp = G->Mate[ns];
     if ( nsp != -1 )
      { Log("Surface %i is identical to object %i (skipping)",ns,nsp);
        continue;
      };

     Log("Assembling T%i at Xi=%e...",ns+1,Xi);
     void *Accelerator = PBC ? SC3D->TAccelerators[ns] : 0;
     G->AssembleBEMMatrixBlock(ns, ns, Omega, kBloch, SC3D->TBlocks[ns], 0,
                               0, 0, Accelerator, false);

   }; // for(ns=0; ns<G->NumSurfaces; ns++)

  /***************************************************************/
  /* if an energy calculation was requested, compute and save    */
  /* the diagonals of the LU factorization of the T blocks       */
  /* (which collectively constitute the diagonal of the LU       */
  /* factorization of the M_{\infinity} matrix).                 */
  /***************************************************************/
  if ( SC3D->WhichQuantities & QUANTITY_ENERGY )
   {
     HMatrix *M=SC3D->M;
     HVector *V=SC3D->MInfLUDiagonal;
     for(int ns=0; ns<G->NumSurfaces; ns++)
      { 
        int Offset=G->BFIndexOffset[ns];
        int NBF=G->Surfaces[ns]->NumBFs;

        int nsp=G->Mate[ns];
        if ( nsp != -1 )
         { 
           /* if this object has a mate, just copy the diagonals of the */
           /* mate's LU factor                                          */
           int MateOffset=G->BFIndexOffset[nsp];
           memcpy(V->DV+Offset,V->DV+MateOffset,NBF*sizeof(double));
         }
        else
         { 
           /* LU factorize the matrix (non-PEC) case or cholesky-factorize */
           /* the negative of the matrix (all-PEC case).                   */
           /* KIND OF HACKY: we use the data buffer inside M as temporary  */
           /* storage for the content of T; this means that we have to     */
           /* call the lapack routines directly instead of using the nice  */
           /* wrappers provided by libhmat.                                */
           Log("LU-factorizing T%i at Xi=%g...",ns+1,Xi);
           int info;
           if (PBC)
            { for(int nbf=0; nbf<NBF; nbf++)
               for(int nbfp=0; nbfp<NBF; nbfp++)
                M->ZM[nbf + nbfp*NBF] = SC3D->TBlocks[ns]->GetEntry(nbf,nbfp);
              zgetrf_(&NBF, &NBF, M->ZM, &NBF, SC3D->ipiv, &info);
              for(int nbf=0; nbf<NBF; nbf++)
               V->SetEntry(Offset+nbf, abs(M->ZM[nbf+nbf*NBF]) );
            }
           else
            { for(int nbf=0; nbf<NBF; nbf++)
               for(int nbfp=0; nbfp<NBF; nbfp++)
                M->DM[nbf + nbfp*NBF] = SC3D->TBlocks[ns]->GetEntryD(nbf,nbfp);
              dgetrf_(&NBF, &NBF, M->DM, &NBF, SC3D->ipiv, &info);
              for(int nbf=0; nbf<NBF; nbf++)
               V->SetEntry(Offset+nbf, M->DM[nbf+nbf*NBF]);
            };
           if (info!=0)
            Log("...FAILED with info=%i (N=%i)",info,NBF);
         };
      };
   };
     
  /***************************************************************/
  /* for each line in the TransFile, apply the specified         */
  /* transformation, then calculate all quantities requested.    */
  /***************************************************************/
  FILE *ByXiKFile=0;
  if (SC3D->ByXiKFileName)
   ByXiKFile=fopen(SC3D->ByXiKFileName,"a");
  for(int ntnq=0, nt=0; nt<SC3D->NumTransformations; nt++)
   { 
     char *Tag=SC3D->GTCList[nt]->Tag;

     /******************************************************************/
     /* skip if all quantities are already converged at this transform */
     /******************************************************************/
     bool AllConverged=true;
     for(int nq=0; AllConverged && nq<SC3D->NumQuantities; nq++)
      { 
        bool ThisConverged;

        if (SC3D->XiConverged[ntnq+nq]) 
         ThisConverged=true;
        else if (PBC && SC3D->BZConverged[ntnq+nq])
         ThisConverged=true;
        else
         ThisConverged=false;

        if (ThisConverged==false) 
         AllConverged=false;

      };
     if (AllConverged)
      { Log("All quantities already converged at Tag %s",Tag);

        for(int nq=0; nq<SC3D->NumQuantities; nq++)
         EFT[ntnq++]=0.0;

        if (ByXiKFile)
         { 
           if (G->LDim==1)
            fprintf(ByXiKFile,"%s %6e %6e ",Tag,Xi,kBloch[0]);
           else // (LDim==2)
            fprintf(ByXiKFile,"%s %6e %6e %6e ",Tag,Xi,kBloch[0],kBloch[1]);
           for(int nq=0; nq<SC3D->NumQuantities; nq++)
            fprintf(ByXiKFile,"%8e ",0.0);
           fprintf(ByXiKFile,"\n");
           fflush(ByXiKFile);
         };

        continue;
      };

     /******************************************************************/
     /* apply the geometrical transform                                */
     /******************************************************************/
     Log("Applying transform %s...",Tag);
     G->Transform( SC3D->GTCList[nt] );
     for(int ns=0; ns<G->NumSurfaces; ns++)
      if (G->SurfaceMoved[ns]) SurfaceNeverMoved[ns]=false;

     /***************************************************************/
     /* assemble U_{a,b} blocks and dUdXYZT_{0,b} blocks            */
     /***************************************************************/
     for(int nb=0, ns=0; ns<G->NumSurfaces; ns++)
      for(int nsp=ns+1; nsp<G->NumSurfaces; nsp++, nb++)
       { 
         /* if we already computed the interaction between objects ns  */
         /* and nsp once at this frequency, and if neither object has  */
         /* moved, then we do not need to recompute the interaction    */
         if ( nt>0 && SurfaceNeverMoved[ns] && SurfaceNeverMoved[nsp] )
          continue;

         Log(" Assembling U(%i,%i)",ns,nsp);
         void *Accelerator = PBC ? SC3D->UAccelerators[nt][nb] : 0;
         if (ns==0)
          G->AssembleBEMMatrixBlock(ns, nsp, Omega, kBloch,
                                    SC3D->UBlocks[nb], SC3D->dUBlocks + 6*nb,
                                    0, 0, Accelerator, false, 
                                    SC3D->NumTorqueAxes, SC3D->dUBlocks + 6*nb + 3, SC3D->GammaMatrix);
         else
          G->AssembleBEMMatrixBlock(ns, nsp, Omega, kBloch, SC3D->UBlocks[nb], 0,
                                    0, 0, Accelerator, false);
       };

     /***************************************************************/
     /* factorize the M matrix and compute casimir quantities       */
     /***************************************************************/
     Factorize(SC3D);
     if ( SC3D->WhichQuantities & QUANTITY_ENERGY )
      EFT[ntnq++]=GetLNDetMInvMInf(SC3D);
     if ( SC3D->WhichQuantities & QUANTITY_XFORCE )
      EFT[ntnq++]=GetTraceMInvdM(SC3D,'X');
     if ( SC3D->WhichQuantities & QUANTITY_YFORCE )
      EFT[ntnq++]=GetTraceMInvdM(SC3D,'Y');
     if ( SC3D->WhichQuantities & QUANTITY_ZFORCE )
      EFT[ntnq++]=GetTraceMInvdM(SC3D,'Z');
     if ( SC3D->WhichQuantities & QUANTITY_TORQUE1 )
      EFT[ntnq++]=GetTraceMInvdM(SC3D,'1');
     if ( SC3D->WhichQuantities & QUANTITY_TORQUE2 )
      EFT[ntnq++]=GetTraceMInvdM(SC3D,'2');
     if ( SC3D->WhichQuantities & QUANTITY_TORQUE3 )
      EFT[ntnq++]=GetTraceMInvdM(SC3D,'3');

     /******************************************************************/
     /* write results to .byXi file (for non-periodic geometries) or   */
     /* to .byXiK file (for periodic geometries).                      */
     /* Note that, for periodic geometries, data are also written to   */
     /* the .byXi file, but this happens at one level higher up in the */
     /* calling hierarchy, in the GetXiIntegrand() routine.            */
     /******************************************************************/
     FILE *f=0;
     if (G->LDim==0)
      { f = fopen(SC3D->ByXiFileName, "a");
        fprintf(f,"%s %.6e ",Tag,Xi);
      }
     else if (G->LDim==1)
      { f = fopen(SC3D->ByXiKFileName, "a");
        fprintf(f,"%s %.6e %.6e ",Tag,Xi,kBloch[0]);
      }
     else if (G->LDim==2)
      { f = fopen(SC3D->ByXiKFileName, "a");
        fprintf(f,"%s %.6e %.6e %.6e ",Tag,Xi,kBloch[0],kBloch[1]);
      }
     else
      ErrExit("%s:%i: internal error",__FILE__,__LINE__);

     for(int nq=SC3D->NumQuantities; nq>0; nq--)
      fprintf(f,"%.8e ",EFT[ntnq-nq]);
     fprintf(f,"\n");
     fclose(f);

     /******************************************************************/
     /* undo the geometrical transform                                 */
     /******************************************************************/
     G->UnTransform();

   }; // for(ntnq=nt=0; nt<SC3D->NumTransformations; nt++)

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (SC3D->WriteCache)
   { StoreCache(SC3D->WriteCache); 
     SC3D->WriteCache=0;
   };

}
