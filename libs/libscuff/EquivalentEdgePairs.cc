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
 * EquivalentEdgePairs.cc -- automatic detection of "equivalent edge pairs," i.e.
 *                        -- pairs of RWG basis functions with identical SIE matrix
 *                        -- elements
 */

#include "EquivalentEdgePairs.h"

#define HAVE_TR1

#ifdef HAVE_CXX11
  #include <unordered_map>
#elif defined(HAVE_TR1)
  #include <tr1/unordered_map>
#endif
#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

namespace scuff {

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- The first portion of this file contains routines for       -*/
/*- analyzing equivalences between edges on single surfaces.   -*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/***************************************************************/
/* given a (full or half) RWG edge E, compute the standard     */
/* rotation that operates on E to put it in standard position. */
/* Standard position means the edge lies on the $x$ axis with  */
/* its centroid at the origin and its positive panel lying in  */
/* the XY plane with Q1.y > 0.                                 */
/***************************************************************/
GTransformation GetStandardRotation(RWGSurface *S, int ne, bool Invert=false)
{
  RWGEdge *E = S->Edges[ne];
  double *QP = S->Vertices + 3*E->iQP;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
if (Invert && E->iQM!=-1) QP = S->Vertices + 3*E->iQM;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  double *V1 = S->Vertices + 3*E->iV1;
  double *X0 = E->Centroid;
  
  double MX0[3];
  MX0[0]=-X0[0];
  MX0[1]=-X0[1];
  MX0[2]=-X0[2];
  GTransformation GTA(MX0);
  double VScratch1[3], VScratch2[3], nHat[3];
  VecNormalize( VecCross( VecSub(V1,X0,VScratch1), VecSub(QP, V1, VScratch2), nHat));
  double Phi   = atan2(nHat[1], nHat[0]);
  double Theta = atan2(sqrt(nHat[0]*nHat[0] + nHat[1]*nHat[1]), nHat[2]);

  //double XHat[3]={1.0, 0.0, 0.0};
  double YHat[3]={0.0, 1.0, 0.0};
  double ZHat[3]={0.0, 0.0, 1.0};
  GTransformation GTB(ZHat, -Phi   * RAD2DEG);
  GTransformation GTC(YHat, -Theta * RAD2DEG);
  GTransformation GT = GTC + GTB + GTA;
  double V1Prime[3];
  GT.Apply(V1, V1Prime);
  double Xi = atan2(V1Prime[1], V1Prime[0]);
  GTransformation GTD(ZHat, -Xi * RAD2DEG);
  return GTD + GT;
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

/***************************************************************/
/* 9-integer "signature" of a G-transformation, useful for     */
/* comparing transformations.                                  */
/* Actually only 6 numbers of needed to specify a G-transform  */
/* (displacement and Euler angles) but I kept getting tripped  */
/* up by the signs of the Euler angles in degenerate cases so  */
/* we'll just be a little wasteful.                            */
/***************************************************************/
#define GTSIGLEN 9
typedef struct { int Signature[GTSIGLEN]; } GTSignature;

#define ANGLEQUANTUM (1.0/(1000.0*2.0*M_PI))
void GetGTSignature(GTransformation GT, double DistanceQuantum, GTSignature *GTSig)
{ 
/*
  // euler angles, defined by M=R_z(ThetaZ) * R_y(ThetaY) * R_x(Theta_x)
  double ThetaX=0.0, ThetaY=0.0, ThetaZ=0.0;
  if ( EqualFloat( GT.M[0][0], 1.0 ) )
   ThetaX = atan2( GT.M[2][1], GT.M[1][1] );
  else if ( EqualFloat( GT.M[1][1], 1.0 ) )
   ThetaY = -atan2( GT.M[2][0], GT.M[0][0] );
  else if ( EqualFloat( GT.M[2][2], 1.0 ) )
   ThetaZ = atan2( GT.M[0][1], GT.M[0][0] );
  else
   { ThetaX = atan2(-1.0*GT.M[1][2], GT.M[1][1]);
     ThetaY = -asin(GT.M[2][0]);
     ThetaZ = atan2( GT.M[1][0], GT.M[0][0]);
   }
*/
  double AbsQuantum=1.0e-4;
  GTSig->Signature[0] = (int)lround(GT.DX[0]/DistanceQuantum);
  GTSig->Signature[1] = (int)lround(GT.DX[1]/DistanceQuantum);
  GTSig->Signature[2] = (int)lround(GT.DX[2]/DistanceQuantum);
  GTSig->Signature[3] = (int)lround(GT.M[0][0]/AbsQuantum);
  GTSig->Signature[4] = (int)lround(GT.M[0][1]/AbsQuantum);
  GTSig->Signature[5] = (int)lround(GT.M[0][2]/AbsQuantum);
  GTSig->Signature[6] = (int)lround(GT.M[1][1]/AbsQuantum);
  GTSig->Signature[7] = (int)lround(GT.M[1][2]/AbsQuantum);
  GTSig->Signature[8] = (int)lround(GT.M[2][2]/AbsQuantum);
}


/***************************************************************/
/* Given two (full or half) RWG edges {Ea, Eb}, determine      */
/* whether or not they are equivalent, i.e Ea is the result    */
/* of applying a rigid geometric transformation (displacement+ */
/* rotation) to Eb.                                            */
/* If the edges are equivalent, return 0 and set Gb2a = the    */
/* transform, i.e. Ea = Gb2a(Eb).                              */
/* if the edges are inequivalent, return a non-zero size_teger */
/* code indicating how we determined this.                     */
/* If PPFileName is non-null, write visualization data to      */
/* file (for debugging purposes).                              */
// possible results of a call to TestEdgeEquivalence: either   */
// edges are equivalent as is or after a sign flip, or else    */
// one of several possible criteria identify them as ineqvlnt. */
/***************************************************************/
#define EQUIV_PLUS        0
#define EQUIV_MINUS       1
#define INEQUIV_LENGTH    2
#define INEQUIV_RADIUS    3
#define INEQUIV_FULLHALF  4
#define INEQUIV_DIPOLE    5
#define INEQUIV_V12       6
#define INEQUIV_QPM       7
static const char *EEResultNames[]=
{ "Equivalent+", "Equivalent-", 
  "length", "radius", "full/half", "dipole", "V12", "QPM" };
#define NUMEERESULTS (sizeof(EEResultNames)/sizeof(EEResultNames[0]))
int TestEdgeEquivalence(RWGSurface *Sa, size_t nea, RWGSurface *Sb, size_t neb,
                        double DistanceQuantum, GTSignature *Gb2aSig, bool bInvert=false)
{ 
  RWGEdge *Ea = Sa->GetEdgeByIndex(nea);
  RWGEdge *Eb = Sb->GetEdgeByIndex(neb);

#define EERELTOL 1.0e-6
  double LengthTol = EERELTOL*Ea->Length;

  // stage-1 check: gross statistics must match
  if ( fabs(Ea->Length - Eb->Length) > LengthTol ) return INEQUIV_LENGTH;
  if ( fabs(Ea->Radius - Eb->Radius) > LengthTol ) return INEQUIV_RADIUS;
  if ( (Ea->iMPanel==-1) != (Eb->iMPanel==-1) )    return INEQUIV_FULLHALF;

  // stage-2 check: magnitudes of dipole moments must match
  double *Va[4], Pa[3], *Vb[4], Pb[3];
  GetPVector(Sa, Ea, Va, Pa);
  GetPVector(Sb, Eb, Vb, Pb);
  double Na = VecNorm(Pa), Nb = VecNorm(Pb);
  if ( fabs(Na-Nb) > EERELTOL*0.5*(Na+Nb) ) return INEQUIV_DIPOLE;

  // stage-3 check: Compute the "standard transformations" GTa, GTb
  // for edges Ea, Eb, then check that Ea = Ta^{-1} Tb(Eb).
   GTransformation GStandard_b = GetStandardRotation(Sb, neb, bInvert);
   GTransformation GStandard_a = GetStandardRotation(Sa, nea);
   GTransformation Gb2a = -GStandard_a + GStandard_b;
   GetGTSignature(Gb2a,DistanceQuantum,Gb2aSig);

  // FIXME this leaves some equivalences undetected....
  // make sure {V1,V2} = {V1,V2}, possibly with twist
  double VbPrime[4][3];
  Gb2a.Apply(Vb[1], VbPrime[1]);
  Gb2a.Apply(Vb[2], VbPrime[2]);
  bool V12Direct = (    VecDistance(VbPrime[1], Va[1])<=LengthTol
                    &&  VecDistance(VbPrime[2], Va[2])<=LengthTol
                   );
  bool V12Twist  = (    VecDistance(VbPrime[1], Va[2])<=LengthTol
                    &&  VecDistance(VbPrime[2], Va[1])<=LengthTol
                   );
  if (!V12Direct && !V12Twist) return INEQUIV_V12;
 
  // make sure {QP,QM} = {QP,QM}, possibly with twist
  bool QPMDirect=true, QPMTwist=false;
  Gb2a.Apply(Vb[0], VbPrime[0]);
  if (Ea->iMPanel==-1)
   QPMDirect = ( VecDistance(VbPrime[0], Va[0]) < LengthTol );
  else
   { Gb2a.Apply(Vb[3], VbPrime[3]);
     QPMDirect = (     VecDistance(VbPrime[0], Va[0])<=LengthTol
                   &&  VecDistance(VbPrime[3], Va[3])<=LengthTol
                 );
     QPMTwist  = (     VecDistance(VbPrime[0], Va[3])<=LengthTol
                   &&  VecDistance(VbPrime[3], Va[0])<=LengthTol
                 );
   }

  if (!QPMDirect && !QPMTwist) return INEQUIV_QPM;

  return QPMTwist ? EQUIV_MINUS : EQUIV_PLUS;
}

/******************************************************************/
/* Given a single (full or half) RWG edge (the "parent" edge) in  */
/* a geometry, construct a list of all "child" edges, i.e.        */
/* edges equivalent to the parent with strictly greater index.    */
/*                                                                */
/* On entry, neParent is the edge index of the parent edge.       */
/* Returns a list of tuples {neChild, Sign, GTSig}, where         */
/*  neChild  = index of child edge                                */
/*  Sign:    = edges are equivalent with a sign flip              */
/*  GTSig    = transform that operates on child edge to yield     */
/*             parent edge                                        */
/******************************************************************/
typedef struct ChildEdgeData
 { int neChild;
   int Sign;
   GTSignature GTSig;
 } ChildEdgeData;

typedef vector<ChildEdgeData> ChildEdgeList;

ChildEdgeList GetChildEdgeList(RWGSurface *S, int neParent, double DistanceQuantum, size_t Results[NUMEERESULTS])
{ 
  ChildEdgeList ChildEdges;
  ChildEdgeData CEData;
  for(CEData.neChild=neParent+1; CEData.neChild<S->NumEdges; CEData.neChild++)
   { 
     int Status=TestEdgeEquivalence(S, neParent, S, CEData.neChild, DistanceQuantum, &(CEData.GTSig));
     if (Status==EQUIV_PLUS || Status==EQUIV_MINUS)
      { CEData.Sign = (Status==EQUIV_MINUS ? -1 : 1);
        ChildEdges.push_back(CEData);
      }
     else
      { Status=TestEdgeEquivalence(S, neParent, S, CEData.neChild, DistanceQuantum, &(CEData.GTSig), true);
        if (Status==EQUIV_PLUS)
         { CEData.Sign = -1;
           ChildEdges.push_back(CEData);
         }
      }
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
#pragma omp critical
     Results[Status]++;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
   }
  return ChildEdges;
}

ChildEdgeList *GetChildren(RWGGeometry *G, int ns, double DistanceQuantum, size_t EEResults[NUMEERESULTS])
{
  RWGSurface *S = G->Surfaces[ns];
  int NE = S->NumEdges;
  ChildEdgeList *Children = new ChildEdgeList[NE];
  int NumParents=0;
  for(int ne=0; ne<NE; ne++)
   { Children[ne] = GetChildEdgeList(S, ne, DistanceQuantum, EEResults);
     if (Children[ne].size() > 0) NumParents +=1;
   }
  if (G->LogLevel>=SCUFF_VERBOSE2)
   Log(" Surface %s: %i/%i edges are parents (%.0g %%)",S->Label, NumParents,NE,100.0*((double)NumParents)/((double)NE));
  return Children;
}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- Now we move into the portion of the file that looks at     -*/
/*- relationships between *pairs* of edges.                    -*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/******************************************************************/
/* General-purpose helper routines that convert back and forth    */
/* between {edge, edge} pairs and integers indexing the list of   */
/* all possible edge pairs.                                       */
/* Note that the two indices that define an edge pair are always  */
/* *full* edge indices, i.e. indices that specify both a surface  */
/* within an RWG geometry as well as an edge on that surface.     */
/******************************************************************/
int EquivalentEdgePairTable::GetEdgePairIndex(int neParent, int neChild)
{
  return neParent*NERadix + neChild;
}

int iabs(int x) { return (x<0 ? -x : x); }

void EquivalentEdgePairTable::ResolveEdgePairIndex(int nPair, int *neParent, int *neChild)
{
  int Sign = (nPair < 0 ? -1 : 1);
  nPair = iabs(nPair);
 *neParent =         nPair / NERadix;
  *neChild = Sign * (nPair % NERadix);
}

/******************************************************************************/
/* If GTSig is the signature of a GTransformation GT, then                    */
/* ParentChildTable[GTSig] is a list of GT-related {parent, child} pairs,     */
/* i.e. for each {parent, child} pair in the list we have parent = GT(child). */
/******************************************************************************/
long JenkinsHash(const char *key, size_t len); // in FIBBICache.cc 

struct GTSigHash 
 { long operator() (const GTSignature &GTSig) const 
    { return JenkinsHash( (const char *)GTSig.Signature, sizeof(GTSignature)); } 
 }; 

typedef struct 
 { bool operator()(const GTSignature &GTSig1, const GTSignature &GTSig2) const  
     { return !memcmp( (const void *)GTSig1.Signature, 
                       (const void *)GTSig2.Signature, 
                       sizeof(GTSignature)); 
     } 
 } GTSigCmp; 

#ifdef HAVE_CXX11 
 typedef unordered_map<GTSignature, iVec *, GTSigHash, GTSigCmp> ParentChildTable; 
#elif defined(HAVE_TR1) 
 typedef tr1::unordered_map<GTSignature, iVec *, GTSigHash, GTSigCmp> ParentChildTable;
#else 
 typedef map<GTSignature, iVec *, GTSigCmp> ParentChildTable;
#endif

void AddParentChildPair(ParentChildTable &PCTable, GTSignature GTSig, int PCPair)
{ 
  ParentChildTable::iterator pctIterator=PCTable.find(GTSig);
  if (pctIterator!=PCTable.end())
   pctIterator->second->push_back(PCPair);
  else
   { iVec *PCPairList = new iVec(1, PCPair);
     PCTable.insert( std::pair<GTSignature, iVec *>(GTSig, PCPairList));
   }
}

ParentChildTable CreateParentChildTable(RWGGeometry *G, int ns, int NERadix,
                                        ChildEdgeList *Children)
{
  // construct list of {parent,child=G(parent)} pairs ordered by G
  ParentChildTable PCTable;
  int NE = G->Surfaces[ns]->NumEdges;
  for(int neParent=0; neParent<NE; neParent++)
   for(size_t nc=0; nc<Children[neParent].size(); nc++)
    { ChildEdgeData CEData = Children[neParent][nc];
      AddParentChildPair(PCTable, CEData.GTSig, CEData.Sign*(neParent*NERadix + CEData.neChild));
    }
  return PCTable;
}

iVec *GetParentChildPairs(ParentChildTable &PCTable, GTSignature GTSig)
{ 
  ParentChildTable::iterator pctIterator=PCTable.find(GTSig);
  return (pctIterator==PCTable.end() ? 0 : pctIterator->second);
}

void DestroyParentChildTable(ParentChildTable &PCTable)
{
  for(ParentChildTable::iterator it=PCTable.begin(); it!=PCTable.end(); it++)
   delete it->second;
}

/****************************************************************************/
/* constructor helper routine to add a single equivalent edge pair **********/
/****************************************************************************/
void EquivalentEdgePairTable::AddEquivalentEdgePair(int ParentPair, int ChildPair,
                                                    bool SignFlip)
{
  if (ChildPair <= ParentPair) return;

#pragma omp critical
 { if (!HasEquivalentPairFlag[ParentPair] && !HasEquivalentPairFlag[ChildPair])
    { HasEquivalentPairFlag[ChildPair]=true;
      if (SignFlip) ChildPair*=-1;
      ChildPairListMap::iterator it=CPLMap.find(ParentPair);
      if (it!=CPLMap.end())
       it->second.push_back(ChildPair);
      else
       { iVec NewChildPairList(1,ChildPair);
         CPLMap.insert(std::pair<int, iVec>(ParentPair,NewChildPairList));
       }
    }
 }
}

void EquivalentEdgePairTable::AddEquivalentEdgePair(int neaParent, int nebParent,
                                                    int neaChild, int nebChild,
                                                    bool SignFlip)
{ 
  if( SameSurface && ( (nebParent<neaParent) || (nebChild < neaChild) ) ) return;

  AddEquivalentEdgePair( GetEdgePairIndex(neaParent, nebParent),
                         GetEdgePairIndex(neaChild, nebChild),
                         SignFlip
                       );
}

void EquivalentEdgePairTable::AddEquivalentEdgePair(int neaParent, int neaChild,
                                                    int nebPair)
{ 
  int aSign = ( (neaChild < 0) ? -1 : 1);
  int bSign = ( (nebPair  < 0) ? -1 : 1);
  int nebParent, nebChild;
  ResolveEdgePairIndex( iabs(nebPair), &nebParent, &nebChild);
  AddEquivalentEdgePair( neaParent, nebParent, iabs(neaChild), nebChild,
                         aSign!=bSign );
}

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

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/* EquivalentEdgePairTable class constructor: Construct a table */
/* of equivalent edge pairs for two surfaces in an RWG geometry.*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
EquivalentEdgePairTable::EquivalentEdgePairTable(RWGGeometry *_G, int nsa, int nsb)
 : G(_G)
{
  int NEA=G->Surfaces[nsa]->NumEdges, NEB=G->Surfaces[nsb]->NumEdges;

  // set some internal class data fields
  DistanceQuantum = 0.1 * GetMinPanelRadius(G);
  NERadix = (NEA > NEB ? NEA : NEB);
  SameSurface = (nsa==nsb);

  /*----------------------------------------------------------------*/
  /* loops in this routine are single-threaded by default           */
  /*----------------------------------------------------------------*/
  int NumThreads=1;
  CheckEnv("SCUFF_EEP_THREADS", &NumThreads);

  /*----------------------------------------------------------------*/
  /* Step 1: identify pairs of equivalent edges within surfaces.    */
  /* If edge ne is equivalent to nePrime>ne, we say nePrime is a    */
  /* child of ne (and ne is a parent of nePrime). After this loop,  */
  /* we have e.g.                                                   */
  /*                                                                */
  /* Children[13][7] = edge index of the 7th child of edge #13      */
  /*                                                                */
  /*      GTs[13][7] = "signature" of the GTransform that transforms*/
  /*                   the 7th child of edge #13 into #13           */
  /*                                                                */
  /* Note that the edge index of a child is always strictly greater */
  /* than the edge index of any of its parents.                     */
  /*                                                                */
  /* For the second surface we further build a "Parent-Child table" */
  /* (bPCTable).  This is a table of {Parent,Child} edge pairs      */
  /* sorted by geometrical transform:                               */
  /*  bPCTable[GT] = list of all GT-related parent-child pairs on   */
  /*                 surface #nsb.                                  */
  /* A parent-child pair is "GT-related" if GT(Child) = Parent.     */
  /*----------------------------------------------------------------*/
  if (G->LogLevel>=SCUFF_VERBOSE2)
   { Log("Detecting edge-pair equivalences(%i,%i)...",nsa,nsb);
     Log("  Step 1: fetching lists of equivalent edges...");
   }
  size_t EEResults[NUMEERESULTS]; // histogram of edge-equivalence test results
  memset(EEResults, 0, NUMEERESULTS*sizeof(size_t)); 
  ChildEdgeList *aChildren = GetChildren(G, nsa, DistanceQuantum, EEResults);
  ChildEdgeList *bChildren = (nsa==nsb) ? aChildren : GetChildren(G, nsb, DistanceQuantum, EEResults);
  ParentChildTable bPCTable = CreateParentChildTable(G, nsb, NERadix, bChildren);

  if (G->LogLevel >= SCUFF_VERBOSE2)
   {
     Log(" results of equivalent-edge detection:");
     int TotalResults=EEResults[EQUIV_PLUS] + EEResults[EQUIV_MINUS];
     Log("  Equivalent            : %6i",EEResults[EQUIV_PLUS]);
     Log("  Equivalent with flip  : %6i",EEResults[EQUIV_MINUS]);
     for(size_t nr=2; nr<NUMEERESULTS; nr++)
      if (EEResults[nr]) 
       { TotalResults+=EEResults[nr];
         Log("Inequivalent(%10s): %6i",EEResultNames[nr],EEResults[nr]);
       }
     Log("-------------------------------------------");
     int TotalPairsTested = NEA*(NEA-1)/2  + (nsa==nsb ? 0 : NEB*(NEB-1)/2);
     Log("  %15s: %6i (should be %i)","Total",TotalResults,TotalPairsTested);
   }

  /*----------------------------------------------------------------*/
  /* second pass to identify equivalent off-diagonal edge pairs.    */
  /* If neaChild is a child of neaParent while nebChild is a child */
  /* of nebParent, then (neaChild,nebChild) is equivalent to        */
  /* (neaParent, nebParent) iff the GTransform that transforms      */
  /* neaChild into neaParent is identical to the GTransform that    */
  /* transforms nebChild into nebParent.                            */
  /*----------------------------------------------------------------*/
  if (G->LogLevel>=SCUFF_VERBOSE2) Log("  Step 2: Identifying equivalent edge pairs...");
  HasEquivalentPairFlag.resize(NERadix*NERadix,false);

  // handle diagonal pairs first separately
  if (nsa==nsb)
   for(int neaParent=0; neaParent<NEA; neaParent++)
    for(size_t nc=0; nc<aChildren[neaParent].size(); nc++)
     AddEquivalentEdgePair(neaParent, neaParent,
                           aChildren[neaParent][nc].neChild,
                           aChildren[neaParent][nc].neChild);

  /*----------------------------------------------------------------*/
  /*- For all {Parent,Child} pairs on surface #nsa, look for pairs  */
  /*- on surface #nsb that are related by the same G-transformation.*/
  /*----------------------------------------------------------------*/
  for(int neaParent=0; neaParent<NEA; neaParent++)
   for(size_t nc=0; nc < aChildren[neaParent].size(); nc++)
    { int neaChild  = aChildren[neaParent][nc].Sign * aChildren[neaParent][nc].neChild;
      iVec *bPairs = GetParentChildPairs(bPCTable, aChildren[neaParent][nc].GTSig);
      if (bPairs==0) continue;
      for(size_t n=0; n<bPairs->size(); n++)
       AddEquivalentEdgePair(neaParent, neaChild, (*bPairs)[n]);
    }

  /*----------------------------------------------------------------*/
  /*- Construct the "irreducible list" of inequivalent matrix       */
  /*- elements, i.e. pairs of edges  (nea, neb) for which there is  */
  /*- no previous equivalent pair (neaPrime, nebPrime) and for which*/
  /*- we will thus need to compute the matrix element.              */
  /*- Each entry in the IEPList is stored with a table of           */
  /*- equivalent edge pairs that come after, so that once we        */
  /*- compute the matrix element we can stamp in place everywhere   */
  /*- else it needs to. The length of IEPList is bounded above by   */
  /*- the total number of edge pairs (nea, neb); the extent to      */
  /*- which IEPList is *shorter* than this upper bound is the       */
  /*- extent of the cost savings we achieve in matrix assembly.     */
  /*----------------------------------------------------------------*/
  IEPList = new IrreducibleEdgePairList;
  for(int neaParent=0; neaParent<NEA; neaParent++)
   for(int nebParent=(nsa==nsb ? neaParent : 0); nebParent<NEB; nebParent++)
    { 
      int ParentPair = GetEdgePairIndex(neaParent, nebParent);
      if (HasEquivalentPairFlag[ParentPair]) continue;
      IrreducibleEdgePair IEP;
      IEP.ParentPair = ParentPair;
      IEP.ChildPairs = GetEquivalentEdgePairs(ParentPair);
      IEPList->push_back(IEP);
    }

  DestroyParentChildTable(bPCTable);
 
  /*----------------------------------------------------------------*/
  /* report some statistics on equivalent edge pairs to log file    */
  /*----------------------------------------------------------------*/
  int NumParentPairs = CountParentPairs();
  int NumChildPairs  = CountChildPairs();
  int NEPairs = (nsa==nsb ? NEA*(NEA+1)/2 : NEA*NEB);
  Log(" Of %u total edge-edge pairs on surfaces (%i,%i) (%s,%s):",NEPairs,nsa,nsb,G->Surfaces[nsa]->Label,G->Surfaces[nsb]->Label);
  Log("    %u are children (savings of %.1f %%)",NumChildPairs,100.0*((double)NumChildPairs)/((double)NEPairs));
  if (G->LogLevel>SCUFF_VERBOSE2)
   { Log("    %u are parents (%.1f %%)",NumParentPairs, 100.0*((double)NumParentPairs) / ((double)NEPairs));
     Log("    %u are unicorns (should be %u)",IEPList->size(), NEPairs - NumParentPairs - NumChildPairs);
   }
}

/***************************************************************/
/* destructor **************************************************/
/***************************************************************/
EquivalentEdgePairTable::~EquivalentEdgePairTable()
{
}

/****************************************************************************/
/* API routines              ************************************************/
/****************************************************************************/
bool EquivalentEdgePairTable::HasEquivalentEdgePair(int ChildPair)
{ return HasEquivalentPairFlag[ChildPair]; } 

bool EquivalentEdgePairTable::HasEquivalentEdgePair(int neaChild, int nebChild)
{ return HasEquivalentPairFlag[GetEdgePairIndex(neaChild, nebChild)]; }

iVec *EquivalentEdgePairTable::GetEquivalentEdgePairs(int ParentPair)
{ 
  ChildPairListMap::iterator it=CPLMap.find(ParentPair);
  return (it==CPLMap.end() ? 0 : &(it->second));
}

iVec *EquivalentEdgePairTable::GetEquivalentEdgePairs(int neaParent, int nebParent)
{ return GetEquivalentEdgePairs( GetEdgePairIndex(neaParent, nebParent) ); }

int EquivalentEdgePairTable::CountParentPairs()
{ return CPLMap.size(); }

int EquivalentEdgePairTable::CountChildPairs()
{ 
   // first way of counting
   int NumChildPairs[2]={0,0};
   for(int nePair=0; nePair<NERadix*NERadix; nePair++)
    if (HasEquivalentPairFlag[nePair])
     NumChildPairs[0]++;

   // second way of counting
   for(ChildPairListMap::iterator it=CPLMap.begin(); it!=CPLMap.end(); it++)
    NumChildPairs[1]+=it->second.size();

  if (G->LogLevel>=SCUFF_VERBOSE2)
   { if (NumChildPairs[0]!=NumChildPairs[1])
      Log("Whoops! NumChildPairs={%i,%i} by different methods",NumChildPairs[0],NumChildPairs[1]); 
     else 
      Log("NumChildPairs={%i,%i} by different methods",NumChildPairs[0],NumChildPairs[1]); 
   }
  return NumChildPairs[0];
}

IrreducibleEdgePairList *EquivalentEdgePairTable::GetIrreducibleEdgePairList()
{ return IEPList; }

/***************************************************************/
/***************************************************************/
/***************************************************************/
void EquivalentEdgePairTable::Export(char *FileName)
{
  char FileNameBuffer[100];
  if (!FileName)
   { FileName = FileNameBuffer;
     snprintf(FileNameBuffer,100,"%s.EEPTable",GetFileBase(G->GeoFileName));
   }
  FILE *f = fopen(FileName,"w");
  if (!f) 
   { Warn("could not open file %s (skipping edge-pair table export)",FileName); 
     return;
   }

  for(int ParentPair=0; ParentPair<NERadix*NERadix; ParentPair++)
   { iVec *ChildPairs = GetEquivalentEdgePairs(ParentPair);
     if (ChildPairs ==0 || ChildPairs->size()==0) continue;
     int neaParent, nebParent;
     ResolveEdgePairIndex(ParentPair, &neaParent, &nebParent);
     fprintf(f,"{%i,%i} ",neaParent, nebParent);
     for(size_t n=0; n<ChildPairs->size(); n++)
      { int neaChild, nebChild;
        ResolveEdgePairIndex( (*ChildPairs)[n], &neaChild, &nebChild);
        fprintf(f,"{%i,%i} ",neaChild,nebChild);
      }
     fprintf(f,"\n");
   }
  fclose(f);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
#if 0
void EquivalentEdgePairTable::Import(char *FileName)
{
}
#endif

} // namespace scuff

#if 0
  bool WriteEEPLogs = CheckEnv("SCUFF_WRITE_EEPLOGS");
  if (WriteEEPLogs)
   { FILE *f=vfopen("/tmp/%s.ChildEdgeLists","w",GetFileBase(G->GeoFileName));
     for(int nfe=0; nfe<NFE; nfe++)
      { fprintf(f,"Parent edge %i: %lu children\n",nfe,Children[nfe].size());
        for(size_t nc=0; nc<Children[nfe].size(); nc++)
         { fprintf(f," %6i {",Children[nfe][nc].Sign*Children[nfe][nc].nfeChild);
           for(int n=0; n<GTSIGLEN; n++) fprintf(f,"%i%s",Children[nfe][nc].GTSig.Signature[n],n==(GTSIGLEN-1) ? "}\n" : ",");
         }
      }
     fclose(f);
   }
#endif

#if 0
  if (WriteEEPLogs)
   { FILE *f=vfopen("/tmp/%s.ParentChildPairLists","w",GetFileBase(G->GeoFileName));
     for(ParentChildTable::iterator pctIterator=PCTable.begin(); pctIterator!=PCTable.end(); pctIterator++)
      { GTSignature Sig = pctIterator->first;
        iVec *ParentChildPairs = pctIterator->second;
        size_t NCP = ParentChildPairs ? ParentChildPairs->size() : 0;
        fprintf(f,"%lu {Parent,Child} pairs with GT {",NCP); 
        for(int n=0; n<GTSIGLEN; n++) fprintf(f,"%i%s",Sig.Signature[n],n==(GTSIGLEN-1) ? "}:\n" : ",");
        for(size_t n=0; n<NCP; n++)
         { int PCPair = (*ParentChildPairs)[n];
           int nfeParent, nfeChild; ResolveEdgePairIndex(PCPair, &nfeParent, &nfeChild);
           fprintf(f,"%i(%i,%i) ",PCPair,nfeParent,nfeChild);
         }
        fprintf(f,"\n");
      }
     fclose(f);
   }
#endif
