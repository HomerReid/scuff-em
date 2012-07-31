
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
 * GetSphericalMoments.cc  -- libscuff class method for computing the 
 *                         -- spherical multipole moments induced by 
 *                         -- the incident field on the scattering objects
 *
 * homer reid              -- 3/2007  -- 7/2012
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <libhrutil.h>
#include <libSpherical.h>
#include <libTriInt.h>

#include "libscuff.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

using namespace scuff;

#define II cdouble(0,1)

/***************************************************************/
/* cubature integrand for the GetMNProjections routine         */
/***************************************************************/
typedef struct GMNPData
 {
   cdouble k; 
   int lMax;
   double *Q;
   double RWGPreFac;
   cdouble *MArray;
   cdouble *NArray;
  
 } GMNPData;

void GMNPIntegrand(double *X, void *parms, double *f)
{
  /*--------------------------------------------------------------*/
  /*- extract parameters from data structure. --------------------*/
  /*- note: the MArray and NArray fields must point to caller-    */
  /*- allocated buffers with enough room for (lMax+1)*(lMax+1)    */
  /*- cdoubles each.                                              */
  /*--------------------------------------------------------------*/
  GMNPData *Data    = (GMNPData *)parms;
  cdouble k         = Data->k;
  int lMax          = Data->lMax;
  double *Q         = Data->Q;
  double RWGPreFac  = Data->RWGPreFac;
  cdouble *MArray   = Data->MArray;
  cdouble *NArray   = Data->NArray;

  /*--------------------------------------------------------------*/
  /*- convert the eval point to spherical coordinates, then get   */
  /*- the M, N vector helmholtz solutions at this eval point      */
  /*--------------------------------------------------------------*/
  double r, Theta, Phi;
  CoordinateC2S(X, &r, &Theta, &Phi);
  GetMNlmArray(lMax, k, r, Theta, Phi, LS_REGULAR, MArray, NArray);

  /*--------------------------------------------------------------*/
  /*- get the vector-valued RWG basis function at this eval point */
  /*--------------------------------------------------------------*/
  double fRWG[3];
  fRWG[0] = RWGPreFac*(X[0] - Q[0]); 
  fRWG[1] = RWGPreFac*(X[1] - Q[1]); 
  fRWG[2] = RWGPreFac*(X[2] - Q[2]); 

  /*--------------------------------------------------------------*/
  /*- and convert it to spherical components                      */
  /*--------------------------------------------------------------*/
  double fRWGS[3];
  VectorC2S(Theta, Phi, fRWG, fRWGS);

  /*--------------------------------------------------------------*/
  /*- loop over all elements in the M, N array to compute the     */
  /*- elements of the integrand vector                            */
  /*--------------------------------------------------------------*/
  cdouble MDotF, NDotF;
  int nf=0;
  int NumLMs=(lMax+1)*(lMax+1);
  for(int nLM=0; nLM<NumLMs; nLM++)
   { 
     MDotF = MArray[3*nLM + 0]*fRWGS[0] + MArray[3*nLM + 1]*fRWGS[1] + MArray[3*nLM + 2]*fRWGS[2];
     NDotF = NArray[3*nLM + 0]*fRWGS[0] + NArray[3*nLM + 1]*fRWGS[1] + NArray[3*nLM + 2]*fRWGS[2];

     f[nf++] = real( MDotF );
     f[nf++] = imag( MDotF );
     f[nf++] = real( NDotF );
     f[nf++] = imag( NDotF );

   };

}

/***************************************************************/
/* get the projections of a single RWG basis function onto the */
/* M and N spherical waves, up to a maximum l-value of lMax.   */
/* Workspace is a caller-allocated workspace with enough room  */
/* to store at least 4*(lMax+1)*(lMax+1) cdoubles.             */
/***************************************************************/
void GetMNProjections(RWGObject *O, int ne, cdouble k, int lMax,
                      cdouble *Workspace, 
                      cdouble *MProjection, cdouble *NProjection)
{
  /* extract edge vertices */
  RWGEdge *E=O->Edges[ne];
  double *QP = O->Vertices + 3*(E->iQP);
  double *V1 = O->Vertices + 3*(E->iV1);
  double *V2 = O->Vertices + 3*(E->iV2);
  double *QM = O->Vertices + 3*(E->iQM);
  double PArea = O->Panels[E->iPPanel]->Area;
  double MArea = O->Panels[E->iMPanel]->Area;

  /* set up data structure passed to GMNPIntegrand.         */
  /* note: the MArray and NArray fields are scratch buffers */
  /* used internally by the integrand routine to store      */
  /* temporary data; we pass for these the MProjection and  */
  /* NProjection output buffers since they happen to have   */ 
  /* the correct size.                                      */ 
  GMNPData MyData, *Data=&MyData;
  Data->k=k;
  Data->lMax=lMax;
  Data->MArray = MProjection;
  Data->NArray = NProjection;

  int NumLMs = (lMax+1)*(lMax+1);
  cdouble *IP = Workspace + 0;
  cdouble *IM = Workspace + 2*NumLMs;

  /* contribution of positive panel */
  Data->Q=QP;
  Data->RWGPreFac = E->Length / (2.0*PArea);
  TriIntFixed(GMNPIntegrand, 4*NumLMs, (void *)Data, QP, V1, V2, 25, (double *)IP);

  /* contribution of negative panel */
  Data->Q=QM;
  Data->RWGPreFac = E->Length / (2.0*MArea);
  TriIntFixed(GMNPIntegrand, 4*NumLMs, (void *)Data, V1, V2, QM, 25, (double *)IM);

  for(int nLM=0; nLM<NumLMs; nLM++)
   { MProjection[nLM] = IP[2*nLM+0] - IM[2*nLM+0];
     NProjection[nLM] = IP[2*nLM+1] - IM[2*nLM+1];
   };

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
HVector *GetSphericalMoments(RWGObject *O, cdouble k, int lMax,
                             HVector *KN, int BFIndexOffset, 
                             HVector *AVector)
{ 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  int NumLMs = (lMax+1)*(lMax+1);
  int NumMoments = 2*NumLMs; // a^E and a^M moments for each l,m
  if ( AVector->N != NumMoments )
   { Warn("wrong-size AVector passed to GetSphericalMoments (reallocating...)");
     AVector=0;
   };
  if ( AVector==0 )
   AVector=new HVector(NumMoments, LHM_COMPLEX);
  AVector->Zero();

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble *Workspace   = new cdouble[4*NumLMs];
  cdouble *MProjection = new cdouble[NumLMs];
  cdouble *NProjection = new cdouble[NumLMs];
 
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  cdouble KAlpha, NAlpha;
  int nbf, ne, nLM;
  for(nbf=BFIndexOffset, ne=0; ne<O->NumEdges; ne++)
    { 
      /***************************************************************/
      /* get the projections of the basis function onto the M and N  */
      /* spherical waves                                             */
      /***************************************************************/
      GetMNProjections(O, ne, k, lMax, Workspace, MProjection, NProjection);

      /***************************************************************/
      /* add the contributions of the electric and magnetic currents */
      /* described by this basis function to the a^{M,E} moments     */
      /***************************************************************/
      KAlpha=KN->GetEntry( nbf++ );
      if ( O->MP->Type==MP_PEC )
       NAlpha = 0.0;
      else 
       NAlpha = -ZVAC*KN->GetEntry( nbf++ );

      for(nLM=0; nLM<NumLMs; nLM++)
       { 
         // contributions to a^M moment 
         AVector->AddEntry( 2*nLM + 0, -k*k*( ZVAC*MProjection[nLM]*KAlpha
                                                  -NProjection[nLM]*NAlpha));

         // contributions to a^E moment
         AVector->AddEntry( 2*nLM + 1, -k*k*( ZVAC*NProjection[nLM]*KAlpha
                                                  +MProjection[nLM]*NAlpha));

       };
 
    };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  delete[] Workspace;
  delete[] MProjection;
  delete[] NProjection;
  return AVector;

}
