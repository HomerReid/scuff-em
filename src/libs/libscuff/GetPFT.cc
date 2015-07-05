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
 * GetPFT.cc   -- master (switchboard) routines for SCUFF algorithms
 *             -- for computing power, force, and torque
 *
 * homer reid  -- 1/2015
 */

#include "libscuff.h"
#include "PFTOptions.h"

namespace scuff { 

/***************************************************************/
/* function prototypes for the various PFT algorithms.         */
/*                                                             */
/* Note: In an effort to improve modularity and readability,   */
/*       the specific PFT algorithms are implemented as        */
/*       standalone (non-class-method) functions.              */
/*       Only the master GetPFT() routine is a class method in */
/*       RWGGeometry; it is just a switchboard routine that    */
/*       hands off to the various non-class-method functions   */
/*       to do the computation.                                */
/***************************************************************/

// PFT by overlap method
void GetOPFT(RWGGeometry *G, int SurfaceIndex, cdouble Omega,
             HVector *KNVector, HVector *RHS, HMatrix *RytovMatrix,
             double PFT[NUMPFT], double **ByEdge=0);

// PFT by displaced-surface-integral method
void GetDSIPFT(RWGGeometry *G, cdouble Omega, double *kBloch,
               HVector *KN, IncField *IF, double PFT[NUMPFT],
               char *BSMesh, double R, int NumPoints,
               bool UseCCQ, bool FarField,
               char *PlotFileName, GTransformation *GT);

void GetDSIPFTTrace(RWGGeometry *G, cdouble Omega,
                    HMatrix *RytovMatrix,
                    double PFT[NUMPFT], bool NeedQuantity[NUMPFT],
                    char *BSMesh, double R, int NumPoints,
                    bool UseCCQ, bool FarField,
                    char *PlotFileName, GTransformation *GT);

// absorbed and scattered/radiated power by equivalence-principle method
void GetEPP(RWGGeometry *G, int SurfaceIndex, cdouble Omega,
            HVector *KNVector, HMatrix *RytovMatrix, double Power[2],
            double **ByEdge=0, HMatrix *TInterior=0, HMatrix *TExterior=0);

// force/torque by equivalence-principle method
void GetEPFT(RWGGeometry *G, int SurfaceIndex, cdouble Omega,
             HVector *KNVector, IncField *IF, double FT[6],
             double **ByEdge=0, int Order=1, double Delta=1.0e-5);

/***************************************************************/
/* Get the full GTransformation that takes an object/surface   */
/* from its native configuration (as described in its .msh     */
/* file) to its current configuration in a scuff-em calculation.*/
/* This is a composition of 0, 1, or 2 GTransformations        */
/* depending on whether (a) the object was DISPLACED/ROTATED   */
/* in the .scuffgeo file, and (b) the object has been          */
/* transformed since it was read in from that file.            */
/***************************************************************/
GTransformation *GetFullSurfaceTransformation(RWGGeometry *G,
                                              int ns,
                                              bool *CreatedGT)
{
  /*--------------------------------------------------------------*/
  /*- If the surface in question has been transformed since we   -*/
  /*- read it in from the meshfile (including any "one-time"     -*/
  /*- transformation specified in the .scuffgeo file) we need    -*/
  /*- to transform the cubature rule accordingly.                -*/
  /*--------------------------------------------------------------*/
  RWGSurface *S=G->Surfaces[ns];
  GTransformation *GT=0;
  *CreatedGT=false;
  if ( (S->OTGT!=0) && (S->GT==0) ) 
   GT=S->OTGT;
  else if ( (S->OTGT==0) && (S->GT!=0) ) 
   GT=S->GT;
  else if ( (S->OTGT!=0) && (S->GT!=0) )
   { *CreatedGT=true;
     GT=new GTransformation(S->GT);
     GT->Transform(S->OTGT);
   };

  return GT;
}

/***************************************************************/
/* Get the power, force, and torque on surface #SurfaceIndex.  */
/***************************************************************/
void RWGGeometry::GetPFT(int SurfaceIndex, HVector *KN,
                         cdouble Omega, double PFT[NUMPFT], PFTOptions *Options)
{
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  PFTOptions DefaultOptions;
  if (Options==0)
   { Options=&DefaultOptions;
     InitPFTOptions(Options);
   };
  int PFTMethod      = Options->PFTMethod;
  char *FluxFileName = Options->FluxFileName;
  RWGSurface *S      = Surfaces[SurfaceIndex];
  IncField *IF       = Options->IF;

  /***************************************************************/
  /* allocate arrays for the edge-by-edge contributions if the   */
  /* user requested flux plots                                   */
  /***************************************************************/
  double **ByEdge=0;
  if (FluxFileName && PFTMethod!=SCUFF_PFT_DSI )
   { 
     int NE = S->NumEdges;
     ByEdge=(double **)mallocEC(NUMPFT*sizeof(double *));
     ByEdge[0]=(double *)mallocEC(NUMPFT*NE*sizeof(double));
     for(int nq=1; nq<NUMPFT; nq++)
      ByEdge[nq] = ByEdge[nq-1] + NE;
   }

  /***************************************************************/
  /* hand off to the individual PFT algorithms to do the         */
  /* computation                                                 */
  /***************************************************************/
  if (     PFTMethod==SCUFF_PFT_OVERLAP
        || PFTMethod==SCUFF_PFT_EPOVERLAP
     )
   { 
     HVector *RHSVector   = Options->RHSVector;
     HMatrix *RytovMatrix = Options->RytovMatrix;
     GetOPFT(this, SurfaceIndex, Omega, KN,
             RHSVector, RytovMatrix, PFT, ByEdge);
   }
  else if (     PFTMethod==SCUFF_PFT_DSI
             || PFTMethod==SCUFF_PFT_EPDSI
          )
   { 
     bool CreatedGT;
     GTransformation *GT
      =GetFullSurfaceTransformation(this, SurfaceIndex, &CreatedGT);

     char *DSIMesh      = Options->DSIMesh;
     double DSIRadius   = Options->DSIRadius;
     int DSIPoints      = Options->DSIPoints;
     bool DSIFarField   = Options->DSIFarField;
     double *kBloch     = Options->kBloch;
     HMatrix *RytovMatrix = Options->RytovMatrix;
     bool *NeedQuantity = Options->NeedQuantity;

     if (RytovMatrix==0)
      GetDSIPFT(this, Omega, kBloch, KN, IF, PFT,
                DSIMesh, DSIRadius, DSIPoints,
                false, DSIFarField, FluxFileName, GT);
     else 
      GetDSIPFTTrace(this, Omega, RytovMatrix,
                     PFT, NeedQuantity,
                     DSIMesh, DSIRadius, DSIPoints,
                     false, DSIFarField, FluxFileName, GT);

     if (CreatedGT) delete(GT);
   }
  else if (PFTMethod==SCUFF_PFT_EP)
   { 
     // EP force, torque
     int EPFTOrder=Options->EPFTOrder;
     double EPFTDelta=Options->EPFTDelta;
     GetEPFT(this, SurfaceIndex, Omega, KN, IF, PFT + 2,
             ByEdge, EPFTOrder, EPFTDelta);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( PFTMethod==SCUFF_PFT_EP             ||
       PFTMethod==SCUFF_PFT_EPOVERLAP      ||
       PFTMethod==SCUFF_PFT_EPDSI
     )
   {
     double Power[2];
     HMatrix *TInterior =  Options->TInterior;
     HMatrix *TExterior =  Options->TExterior;
     HMatrix *RytovMatrix = Options->RytovMatrix;
     GetEPP(this, SurfaceIndex, Omega, KN, RytovMatrix, Power,
            ByEdge, TInterior, TExterior);

     // replace absorbed and scattered power with EP calculations
     PFT[0] = Power[0];
     PFT[1] = Power[1];
   };


  /***************************************************************/
  /* produce flux plots if that was requested ********************/
  /***************************************************************/
  if (ByEdge)
   { 
     static const char *PFTNames[NUMPFT]
      ={"PAbs","PScat","FX","FY","FZ","TX","TY","TZ"};

     for(int nq=0; nq<NUMPFT; nq++)
      { char Tag[20];
        snprintf(Tag,20,"%s(%s)",PFTNames[nq],z2s(Omega));
        S->PlotScalarDensity(ByEdge[nq],true,FluxFileName,Tag);
      };

     free(ByEdge[0]);
     free(ByEdge);

   };
 
}

/***************************************************************/
/* routine for initializing a PFTOptions structure to default  */
/* values; creates and returns a new default structure if      */
/* called with Options=NULL or with no argument                */
/***************************************************************/
PFTOptions *InitPFTOptions(PFTOptions *Options)
{
  if (Options==0)
   Options = (PFTOptions *)mallocEC( sizeof(PFTOptions) );

  // general options
  Options->PFTMethod = SCUFF_PFT_DEFAULT;
  Options->FluxFileName=0;
  Options->RytovMatrix=0;
  Options->IF=0;
  Options->kBloch=0;

  // options affecting overlap PFT computation
  Options->RHSVector = 0;

  // options affecting DSI PFT computation
  Options->DSIMesh=0;
  Options->DSIRadius=10.0;
  Options->DSIPoints=302;
  Options->DSIFarField=false;
  for(int nq=0; nq<NUMPFT; nq++) 
   Options->NeedQuantity[nq]=true;
 
  // options affecting EP power computation
  Options->TInterior=0;
  Options->TExterior=0;

  // options affecting EP force / torque computation
  Options->EPFTOrder=1;
  Options->EPFTDelta=1.0e-5;

  Options->GetRegionPFTs=false;

  return Options;
}

} // namespace scuff
