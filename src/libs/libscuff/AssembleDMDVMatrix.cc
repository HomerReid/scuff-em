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
 * AssembleDMDVMatrix.cc
 * 
 * homer reid   -- 05/2011
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

#include <libhmat.h>
#include <libhrutil.h>

#include "libscuff.h"

namespace scuff {

#define II cdouble(0,1)

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGGeometry::GetMEVertexDerivative(void *pLFW, 
                                        RWGObject *O, int ne, 
                                        RWGObject *OP, int nep,
                                        int VertexIndex, int Mu, 
                                        double AverageRadius,
                                        double Frequency, int RealFreq,
                                        cdouble *dmdv)
{

  /***************************************************************/
  /* for the purposes of this calculation we need to mask the ****/
  /* fact that O might have a table of static panel-panel      ***/
  /* integrals                                                 ***/
  /***************************************************************/
  StaticPPIDataTable *SPPIDTableSave;
  SPPIDTableSave=O->SPPIDTable;
  O->SPPIDTable=0;
  
  /***************************************************************/
  /* use a four-point finite-difference scheme to estimate the   */
  /* derivatives of the matrix element w.r.t. displacement of    */
  /* this vertex                                                 */
  /***************************************************************/
  double Delta;
  static HMatrix *MM2D=0, *MM1D, *MP1D, *MP2D;

  if (MM2D==0)
   { MM2D=new HMatrix(2,2,LHM_COMPLEX);
     MM1D=new HMatrix(2,2,LHM_COMPLEX);
     MP1D=new HMatrix(2,2,LHM_COMPLEX);
     MP2D=new HMatrix(2,2,LHM_COMPLEX);
   };

  Delta=1.0e-2*AverageRadius;

  O->Vertices[3*VertexIndex + Mu] -= 2.0*Delta;
  StampMatrixElements(pLFW, O, ne, 0, OP, nep, 0, Frequency, RealFreq, 
                      0, 0, MM2D, 0, 0, 0, 0, 0, 0);

  O->Vertices[3*VertexIndex + Mu] += Delta;
  StampMatrixElements(pLFW, O, ne, 0, OP, nep, 0, Frequency, RealFreq, 
                      0, 0, MM1D, 0, 0, 0, 0, 0, 0);
  O->Vertices[3*VertexIndex + Mu] += 2.0*Delta;
  StampMatrixElements(pLFW, O, ne, 0, OP, nep, 0, Frequency, RealFreq, 
                      0, 0, MP1D, 0, 0, 0, 0, 0, 0);
  O->Vertices[3*VertexIndex + Mu] += Delta;
  StampMatrixElements(pLFW, O, ne, 0, OP, nep, 0, Frequency, RealFreq, 
                      0, 0, MP2D, 0, 0, 0, 0, 0, 0);
  O->Vertices[3*VertexIndex + Mu] -= 2.0*Delta;
 
  dmdv[0] = (  -1.0*MP2D->GetEntry(0,0) 
               +8.0*MP1D->GetEntry(0,0)
               -8.0*MM1D->GetEntry(0,0)
               +1.0*MM2D->GetEntry(0,0) ) / (12.0);

  dmdv[1] = (  -1.0*MP2D->GetEntry(0,1) 
               +8.0*MP1D->GetEntry(0,1)
               -8.0*MM1D->GetEntry(0,1)
               +1.0*MM2D->GetEntry(0,1) ) / (12.0);

  dmdv[2] = (  -1.0*MP2D->GetEntry(1,1) 
               +8.0*MP1D->GetEntry(1,1)
               -8.0*MM1D->GetEntry(1,1)
               +1.0*MM2D->GetEntry(1,1) ) / (12.0);

  /***************************************************************/
  /* unmask static panel-panel-integral table inside O ***********/
  /***************************************************************/
  O->SPPIDTable=SPPIDTableSave;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void RWGGeometry::AssembleDMDVMatrix(int ObjectIndex, int VertexIndex, int Mu, 
                                     double Frequency, int RealFreq,
                                     HMatrix *DMDV, int nThread)
{ 
  int ne, nep, nop, Offset, OffsetP;
  int Formulation;
  RWGObject *O, *OP;
  RWGEdge *E, *EP;
  cdouble dmdv[3];
  static void *pLFW=0;

  if (nThread <= 0) nThread = GetNumThreads();

  if (pLFW==0)
   pLFW=CreateLFWorkspace();

  /***************************************************************/
  /* update the cached values of Epsilon and Mu inside all       */
  /* objects and inside (this)                                   */
  /***************************************************************/
  int no;
  MP->GetEpsMu(Frequency, RealFreq, &EpsThisFreq, &MuThisFreq);
  for(no=0; no<NumObjects; no++)
   Objects[no]->MP->GetEpsMu(Frequency, RealFreq, 
                             &(Objects[no]->EpsThisFreq), 
                             &(Objects[no]->MuThisFreq));

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  O=Objects[ObjectIndex];
  Offset=BFIndexOffset[ObjectIndex];
  DMDV->Zero();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  double AverageRadius;
  int np, NumTouchingPanels;
  RWGPanel *P;
  AverageRadius=0.0;
  NumTouchingPanels=0;
  for(np=0, P=O->Panels[0]; np<O->NumPanels; P=O->Panels[++np])
   if (P->VI[0]==VertexIndex || P->VI[1]==VertexIndex || P->VI[2]==VertexIndex )
    { AverageRadius+=P->Radius;
      NumTouchingPanels++;
    };
  AverageRadius/=((double)NumTouchingPanels);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  for(ne=0, E=O->Edges[ne]; ne<O->NumEdges; E=O->Edges[++ne])
   { 
     if ( (E->iQP != VertexIndex) && (E->iQM != VertexIndex ) &&
          (E->iV1 != VertexIndex) && (E->iV2 != VertexIndex ) )
      continue;

     for(nop=0, OP=Objects[nop]; nop<NumObjects; OP=Objects[++nop])
      { 
        OffsetP=BFIndexOffset[nop];

        for(nep=(OP==O ? ne : 0); nep<O->NumEdges; nep++)
         { 
           GetMEVertexDerivative(pLFW, O, ne, OP, nep, VertexIndex, Mu, AverageRadius, 
                                 Frequency, RealFreq, dmdv);

           if (AllPEC)
            { 
              DMDV->SetEntry(Offset+ne, OffsetP+nep, dmdv[0] );
            }
           else
            { DMDV->SetEntry(Offset+2*ne,   OffsetP+2*nep,   dmdv[0] );
              DMDV->SetEntry(Offset+2*ne,   OffsetP+2*nep+1, dmdv[1] );
              DMDV->SetEntry(Offset+2*ne+1, OffsetP+2*nep,   dmdv[1] );
              DMDV->SetEntry(Offset+2*ne+1, OffsetP+2*nep+1, dmdv[2] );
            };
         };

      };

   };

} 

/***************************************************************/
/* Allocate a new HMatrix to store the dMdV matrix.    *********/
/***************************************************************/
HMatrix *RWGGeometry::AllocateDMDVMatrix(int RealFreq)
{ 
  HMatrix *DMDV;

  if (RealFreq)
   DMDV=new HMatrix(TotalBFs,TotalBFs,LHM_COMPLEX,LHM_SYMMETRIC);
  else
   DMDV=new HMatrix(TotalBFs,TotalBFs,LHM_REAL,LHM_SYMMETRIC);

  return DMDV;

} 

} // namespace scuff
