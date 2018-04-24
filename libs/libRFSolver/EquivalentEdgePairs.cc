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
 * EquivalentEdgePairs.cc
 */

#include "EquivalentEdgePairs.h"

namespace scuff {

/***************************************************************/
/* Compute an orthonormal basis of vectors {xHat, yHat, zHat}  */
/* for the positive panel of an RWG edge.                      */
/* QVVQ[0,1,2,3] = QP, V1, V2, QM / centroid                   */
/* P = positive panel                                          */
/***************************************************************/
bool AdjustXYZAxes=true;
void GetXYZSystem(double *QVVQ[4], double *Centroid, RWGPanel *P, double xyzHat[3][3])
{
  xyzHat[2][0] = P->ZHat[0];
  xyzHat[2][1] = P->ZHat[1];
  xyzHat[2][2] = P->ZHat[2];
  double LengthScale = VecNormalize( VecSub(QVVQ[1], Centroid, xyzHat[0]) );
  VecCross(xyzHat[2], xyzHat[0], xyzHat[1]);

  // ensure that the Y axis points into the positive panel, i.e. towards QP
  if (AdjustXYZAxes)
   { double EP[3];
     VecScaleAdd(Centroid, 0.1*LengthScale, xyzHat[1], EP);
     if ( VecDistance(EP,QVVQ[0]) > VecDistance(Centroid,QVVQ[0]) )
      { VecScale(xyzHat[1], -1.0); 
        VecScale(xyzHat[2], -1.0); 
      }
   }
}

/***************************************************************/
/* compute a vector associated with an RWG edge that helps to  */
/* determine when two edges are equivalent.                    */
/* This vector is proportional to the dipole moment of the RWG */
/* basis function, which is why I call it 'P.'                 */
/***************************************************************/
void GetPVector(RWGSurface *S, RWGEdge *E, double *V[4], double P[3])
{ 
  V[0] = S->Vertices + 3*E->iQP;
  V[1] = S->Vertices + 3*E->iV1;
  V[2] = S->Vertices + 3*E->iV2;
  V[3] = (E->iQM==-1) ? E->Centroid : S->Vertices + 3*E->iQM;
  VecSub(V[0], V[3], P);
  VecScale(P, E->Length);
}

// possible results of a call to TestEdgeEquivalence: either the
// edges are equivalent (result==0) or not (result>=1),
// with the index of the result indicating which criterion tipped
// us off to detecting the inequivalence (lower=faster)
static const char *EEResultNames[]=
{ "Equivalent", "Length",    "Radius",   "Full/Half", "Dipole",
  "QPbPrime",   "V1bPrime",  "V2bPrime", "QMbPrime"
};
#define NUMEERESULTS (sizeof(EEResultNames)/sizeof(EEResultNames[0]))

/***************************************************************/
/* Given two (full or half) RWG edges {Ea, Eb}, determine      */
/* whether or not they are equivalent, i.e Ea is the result    */
/* of applying a rigid geometric transformation (displacement+ */
/* rotation) to Eb.                                            */
/* If the edges are equivalent, return 0 and set Gb2a = the    */
/* transform, i.e. Ea = Gb2a(Eb).                              */
/* if the edges are inequivalent, return a non-zero integer    */
/* code indicating how we determined this.                     */
/* If PPFileName is non-null, write visualization data to      */
/* file (for debugging purposes).                              */
/***************************************************************/
#define EERELTOL 1.0e-6
int TestEdgeEquivalence(RWGSurface *Sa, int nea, RWGSurface *Sb, int neb,
                        GTransformation &Gb2a, double *Distance)
{ 
  RWGEdge *Ea = Sa->GetEdgeByIndex(nea);
  RWGEdge *Eb = Sb->GetEdgeByIndex(neb);

  double LengthTol = EERELTOL*Ea->Length;

  // stage-1 check: gross statistics must match
  if ( fabs(Ea->Length - Eb->Length) > LengthTol ) return 1;
  if ( fabs(Ea->Radius - Eb->Radius) > LengthTol ) return 2;
  if ( (Ea->iMPanel==-1) != (Eb->iMPanel==-1) )    return 3;

  // stage-2 check: magnitudes of dipole moments must match
  double *Va[4], Pa[3], *Vb[4], Pb[3];
  GetPVector(Sa, Ea, Va, Pa);
  GetPVector(Sb, Eb, Vb, Pb);
  double Na = VecNorm(Pa), Nb = VecNorm(Pb);
  if ( fabs(Na-Nb) > EERELTOL*0.5*(Na+Nb) ) return 4;

  // stage-3 check: Compute the G-transformation that (a) displaces Eb
  // until its centroid agrees with that of Ea, (b) rotates Eb so that
  // the X,Y,Z coordinate system of its positive panel is aligned with
  // that of Ea. Then confirm that all 4 vertices of the transformed Eb 
  // coincide with those of Ea.
  double aBasis[3][3], bBasis[3][3];
  GetXYZSystem(Va, Ea->Centroid, Sa->Panels[Ea->iPPanel], aBasis);
  GetXYZSystem(Vb, Eb->Centroid, Sb->Panels[Eb->iPPanel], bBasis);

  // to transform Eb into Ea, first translate through the centroid-centroid displacement
  // vector, then rotate the xyzb system into the xyza system
  Gb2a.Reset();
  VecSub(Ea->Centroid, Eb->Centroid, Gb2a.DX);
  GTransformation GT2;
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    GT2.M[Mu][Nu] = aBasis[0][Mu]*bBasis[0][Nu] + aBasis[1][Mu]*bBasis[1][Nu] + aBasis[2][Mu]*bBasis[2][Nu];
  Gb2a.Transform(GT2);
 
  // pause to make visualization plots if requested
#if 0
  if (PPFileName)
   { double *aBasis[3]; aBasis[0] = xaHat; aBasis[1] = yaHat; aBasis[2] = zaHat;
     double *bBasis[3]; bBasis[0] = xbHat; bBasis[1] = ybHat; bBasis[2] = zbHat;
     PlotEdgeEquivalence(PPFileName, Gb2a, Ea, Va, aBasis, Eb, Vb, bBasis);
   }
#endif
   
  for(int nv=0; nv<(Ea->iMPanel==-1 ? 3 : 4); nv++)
   { double VbPrime[3];
     Gb2a.Apply(Vb[nv], VbPrime);
     if ( VecDistance(VbPrime, Va[nv]) > LengthTol )
      return 5 + nv;
   }
  *Distance = VecDistance(Ea->Centroid, Eb->Centroid);
  return 0; // edges are equivalent
}

int TestEdgeEquivalence(RWGGeometry *G, int nfea, int nfeb, GTransformation &Gb2a, double *Distance)
{ int nea, nsa; RWGSurface *Sa = G->ResolveEdge(nfea, &nsa, &nea);
  int neb, nsb; RWGSurface *Sb = G->ResolveEdge(nfeb, &nsb, &neb);
  return TestEdgeEquivalence(Sa, nea, Sb, neb, Gb2a, Distance);
}

/******************************************************************/
/* Given a single (full or half) RWG edge in a geometry, identify */
/* all edges equivalent to that edge.                             */
/* On entry, nfeParent is the full edge index of the parent edge. */
/* Returns a d-sorted list of tuples {d, nfeChild, GT}, where     */
/*  nfeChild = (full) edge index of child edge                    */
/*  GTs      = transform that operates on child edge to yield     */
/*             parent edge                                        */
/*  d        = child-parent centroid-centroid distance, quantized */
/*             in units of 1/10 the minimal panel radius          */
/******************************************************************/
#define LOOP_OVER_CHILD_EDGES(List,it) for (ChildEdgeList::iterator it=List.begin(); it!=List.end(); it++)

#define LOOP_OVER_CHILD_EDGES_WITH_DISTANCE(List,it,qd)                                     \
 std::pair <ChildEdgeList::iterator, ChildEdgeList::iterator> range = List.equal_range(qd); \
 for (ChildEdgeList::iterator it=range.first; it!=range.second; it++)

typedef struct ChildEdgeData
 { int nfeChild;
   GTransformation GT;
 } ChildEdgeData;

typedef std::multimap<size_t, ChildEdgeData> ChildEdgeList;

ChildEdgeList FindChildEdges(RWGGeometry *G, int nfeParent, double DistanceQuantum, int Results[NUMEERESULTS])
{ 
  ChildEdgeList ChildEdges;
  for(int nfeChild=nfeParent+1; nfeChild<G->TotalEdges; nfeChild++)
   { GTransformation GT;
     double Distance;
     int Status=TestEdgeEquivalence(G, nfeParent, nfeChild, GT, &Distance);
     Results[Status]++;
     size_t QDistance = (size_t)round(Distance/DistanceQuantum);
     if (Status==0)
      { ChildEdgeData CEData;
        CEData.nfeChild = nfeChild;
        CEData.GT = GT;
	ChildEdges.insert(std::pair<size_t, ChildEdgeData>(QDistance,CEData));
      }
   }
  return ChildEdges;
}

/****************************************************************************/
/* constructor helper routine to add a single equivalent edge pair **********/
/****************************************************************************/
void EquivalentEdgePairTable::AddEquivalentEdgePair(int ParentPair, int ChildPair)
{ 
  if (ParentPairArray)
   { if (ParentPairArray[ChildPair]!=-1) return;
     ParentPairArray[ChildPair]=ParentPair;
     ChildPairListArray[ParentPair].push_back(ChildPair);
   }
  else
   { if ( ParentPairMap.find(ChildPair) != ParentPairMap.end() )
      return;
     ParentPairMap[ChildPair]=ParentPair;
     ChildPairListMap[ParentPair].push_back(ChildPair);
   }
}

void EquivalentEdgePairTable::AddEquivalentEdgePair(int nfeaParent, int nfebParent,
                                                    int nfeaChild, int nfebChild)
{ AddEquivalentEdgePair(nfeaParent*G->TotalEdges + nfebParent, nfeaChild*G->TotalEdges + nfebChild); }

/******************************************************************/
/******************************************************************/
/******************************************************************/
double GetMinPanelRadius(RWGGeometry *G)
{
  double MPR=HUGE_VAL;
  for(int ns=0; ns<G->NumSurfaces; ns++)
   for(int np=0; np<G->Surfaces[ns]->NumPanels; np++)
    MPR = fmin(MPR, G->Surfaces[ns]->Panels[np]->Radius);
  return MPR;
}

/******************************************************************/
/******************************************************************/
/******************************************************************/
EquivalentEdgePairTable::EquivalentEdgePairTable(RWGGeometry *_G): G(_G)
{
  int NFE = G->TotalEdges, NFEPairs=NFE*(NFE+1)/2;
  Log("Detecting equivalences among %i edges (%i pairs)...",NFE,NFEPairs);

  ChildPairListArray=0;
  ParentPairArray=0;
  bool UseArrays = !CheckEnv("SCUFF_USE_EEPTMAPS");
  if (UseArrays)
   { ChildPairListArray = new iVec[NFE*NFE];
     ParentPairArray    = new int[NFE*NFE];
     for(int n=0; n<NFE*NFE; n++)
      ParentPairArray[n]=-1;
   };
  
  /*----------------------------------------------------------------*/
  /* loops in this routine are single-threaded by default           */
  /*----------------------------------------------------------------*/
  int NumThreads=1;
  CheckEnv("SCUFF_EEP_THREADS", &NumThreads);
 
  /*----------------------------------------------------------------*/
  /* first pass to identify equivalent edges. If edge #nfe is       */
  /* equivalent to edge #nfePrime>nfe, we say nfePrime is a child   */
  /* of nfe (and nfe is a parent of nfePrime). After this loop,     */
  /* we have e.g.                                                   */
  /*                                                                */
  /* Children[13][7] = edge index of the 7th child of edge #13      */
  /*                                                                */
  /*      GTs[13][7] = GTransform that operates on the 7th          */
  /*                   child of edge #13 to yield edge #13          */
  /*                                                                */
  /* Note that the full edge index of a child is always strictly    */
  /* greater than the edge index of any of its parents.             */
  /*----------------------------------------------------------------*/
  Log("  Step 1: fetching lists of equivalent edges...");
  ChildEdgeList *Children= new ChildEdgeList[NFE];
  int EEResults[NUMEERESULTS]; // histogram of edge-equivalence test results
  memset(EEResults, 0, NUMEERESULTS*sizeof(int)); 
  int NumParentEdges=0;
  double DistanceQuantum = 0.1 * GetMinPanelRadius(G);
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
  for(int nfe=0; nfe<NFE; nfe++)
   { Children[nfe] = FindChildEdges(G, nfe, DistanceQuantum, EEResults);
     if (Children[nfe].size() > 0) NumParentEdges +=1;
   }
  Log("  %i/%i edges are parents (%e %%)",NumParentEdges,NFE,100.0*((double)NumParentEdges)/((double)NFE));

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  int TotalResults=0;
  Log(" results of equivalent-edge detection:");
  for(size_t nr=0; nr<NUMEERESULTS; nr++)
   if (EEResults[nr]) 
    { TotalResults+=EEResults[nr];
      Log("  %15s: %6i",EEResultNames[nr],EEResults[nr]);
    }
  Log("-------------------------------------------");
  Log("  %15s: %6i (should be %i)","Total",TotalResults,NFE);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /*----------------------------------------------------------------*/
  /* second pass to identify equivalent off-diagonal edge pairs.    */
  /* If nfeaPrime is a child of nfea, while nfebPrime is a child of */
  /* nfeb, then (nfeaPrime,nfebPrime) is equivalent to (nfea,nfeb)  */
  /* if the GTransform that takes nfeaPrime to nfea is identical    */
  /* to the GTransform that takes nfebPrime to nfeb.                */
  /*----------------------------------------------------------------*/
  Log("  Step 2: Identifying equivalent edge pairs...");

  // handle diagonal pairs first separately
  for(int nfeaParent=0; nfeaParent<NFE; nfeaParent++)
   if ( !HasParent(nfeaParent*NFE + nfeaParent) )
    LOOP_OVER_CHILD_EDGES(Children[nfeaParent],it)
     AddEquivalentEdgePair(nfeaParent, nfeaParent,it->second.nfeChild,it->second.nfeChild);

  // off-diagonal pairs
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
  for(int nfeaParent=0; nfeaParent<NFE; nfeaParent++)
   for(int nfebParent=nfeaParent+1; nfebParent<NFE; nfebParent++)
    {  
      if (nfebParent==nfeaParent+1) LogPercent(nfeaParent, NFE, 10);
      int ParentPair = nfeaParent*NFE + nfebParent;
      if (HasParent(ParentPair)) continue;
      LOOP_OVER_CHILD_EDGES(Children[nfeaParent],ita)
       { LOOP_OVER_CHILD_EDGES_WITH_DISTANCE(Children[nfebParent],itb,ita->first)
          { int nfeaChild = ita->second.nfeChild, nfebChild = itb->second.nfeChild;
            int ChildPair = nfeaChild*NFE + nfebChild;
            if (HasParent(ChildPair)) continue;
            if ( ita->second.GT.IsIdentical(itb->second.GT) )
             { 
#pragma omp critical
               AddEquivalentEdgePair(ParentPair, ChildPair);
             }
          }
       }
    }

  /*----------------------------------------------------------------*/
  /* report some statistics on equivalent edge pairs to log file    */
  /*----------------------------------------------------------------*/
  int NumParentPairs = CountParentPairs();
  int NumChildPairs  = CountChildPairs();
  Log(" Of %u total edge-edge pairs: ",NFEPairs);
  Log("    %u are children (savings of %.1f %%)",NumChildPairs,100.0*((double)NumChildPairs)/((double)NFEPairs));
  Log("    %u are parents (%.1f %%)",NumParentPairs, 100.0*((double)NumParentPairs) / ((double)NFEPairs));
  Log("    %u are unicorns",NFEPairs - NumParentPairs - NumChildPairs);
}


/***************************************************************/
/* destructor **************************************************/
/***************************************************************/
EquivalentEdgePairTable::~EquivalentEdgePairTable()
{
 if (ChildPairListArray) delete[] ChildPairListArray;
 if (ParentPairArray) delete[] ParentPairArray;
}

/****************************************************************************/
/* API routines              ************************************************/
/****************************************************************************/
bool EquivalentEdgePairTable::HasParent(int ChildPair)
{ return ParentPairArray ? (ParentPairArray[ChildPair]!=-1)
                         : (ParentPairMap.find(ChildPair) != ParentPairMap.end());
} 

bool EquivalentEdgePairTable::HasParent(int nfeaChild, int nfebChild)
{ return HasParent(nfeaChild*G->TotalEdges + nfebChild); }

bool EquivalentEdgePairTable::HasParent(int nsa, int neaChild, int nsb, int nebChild)
{ return HasParent( G->UnResolveEdge(nsa, neaChild), G->UnResolveEdge(nsb, nebChild) ); } 

iVec EquivalentEdgePairTable::GetChildren(int ParentPair)
{ return ChildPairListArray ? ChildPairListArray[ParentPair] : ChildPairListMap[ParentPair]; } 

iVec EquivalentEdgePairTable::GetChildren(int nfeaParent, int nfebParent)
{ return GetChildren(nfeaParent*G->TotalEdges + nfebParent); }

iVec EquivalentEdgePairTable::GetChildren(int nsa, int neaParent, int nsb, int nebParent)
{ return GetChildren( G->UnResolveEdge(nsa, neaParent), G->UnResolveEdge(nsb, nebParent) ); }

int EquivalentEdgePairTable::CountParentPairs()
{ 
  if (!ChildPairListArray) return ChildPairListMap.size();

  int NumParentPairs=0, NFE=G->TotalEdges;
  for(int n=0; n<NFE*NFE; n++)
   if (ChildPairListArray[n].size()>0) 
    NumParentPairs++;
  return NumParentPairs;
}

int EquivalentEdgePairTable::CountChildPairs()
{ 
   if (!ParentPairArray) return ParentPairMap.size();
   int NumChildPairs[2]={0,0}, NFE=G->TotalEdges;
   for(int n=0; n<NFE*NFE; n++)
    if (ParentPairArray[n]!=-1)
     NumChildPairs[0]++;
   for(int n=0; n<NFE*NFE; n++)
    NumChildPairs[1]+=(ChildPairListArray[n].size());
  if (NumChildPairs[0]!=NumChildPairs[1])
   Log("Whoops! NumChildPairs={%i,%i} by different methods",NumChildPairs[0],NumChildPairs[1]); 
  else 
   Log("NumChildPairs={%i,%i} by different methods",NumChildPairs[0],NumChildPairs[1]); 
  return NumChildPairs[0]; 
} 

} // namespace scuff
