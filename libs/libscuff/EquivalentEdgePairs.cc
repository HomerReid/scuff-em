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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

namespace scuff {

long JenkinsHash(const char *key, size_t len); // in FIBBICache.cc 

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/* The first portion of this file contains routines for         */
/* analyzing equivalences between pairs of RWG functions, i.e.  */
/* relationships of the form (neParent<-->neChild).             */
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/***************************************************************/
/* given a (full or half) RWG edge E, compute the standard     */
/* rotation that operates on E to put it in standard position. */
/* Standard position means the edge lies on the $x$ axis with  */
/* its centroid at the origin and its positive panel lying in  */
/* the XY plane with Q1.x > 0.                                 */
/***************************************************************/
GTransformation GetStandardRotation(RWGSurface *S, int ne, bool FlipQPM=false, bool FlipV12=false)
{
  RWGEdge *E = S->Edges[ne];
  double *QP = S->Vertices + 3*E->iQP;
  if (FlipQPM && E->iQM!=-1) QP=S->Vertices + 3*E->iQM;
  double *V1 = S->Vertices + 3*E->iV1;
  double *V2 = S->Vertices + 3*E->iV2;
  double *X0 = E->Centroid;

  if (FlipV12)
   { double *Temp=V1; V1=V2; V2=Temp; }
  if (VecDistance(QP,V2) < VecDistance(QP,V1))
   { double *Temp=V1; V1=V2; V2=Temp; }
  
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
void GetPVector(RWGSurface *S, RWGEdge *E, double *V[4], double P[3], bool FlipQPM=false, bool FlipV12=false)
{ 
  V[0] = S->Vertices + 3*E->iQP;
  V[1] = S->Vertices + 3*E->iV1;
  V[2] = S->Vertices + 3*E->iV2;
  V[3] = (E->iQM==-1) ? E->Centroid : S->Vertices + 3*E->iQM;
  if (FlipV12)
   { double *Temp=V[1]; V[1]=V[2]; V[2]=Temp; }
  if (FlipQPM && E->iQM!=-1)
   { double *Temp=V[0]; V[0]=V[3]; V[3]=Temp; }
  VecSub(V[0], V[3], P);
  VecScale(P, E->Length);
}

#define EDGESIGLEN 3
typedef struct { int Signature[EDGESIGLEN]; } EdgeSignature;
void GetEdgeSignature(RWGSurface *S, int ne, double DistanceQuantum, EdgeSignature &ES)
 { 
   RWGEdge *E = S->Edges[ne];
   double *QP = S->Vertices + 3*E->iQP;
   double *V1 = S->Vertices + 3*E->iV1;
   double *V2 = S->Vertices + 3*E->iV2;

   double Length    = E->Length;
   double Perimeter = VecDistance(QP, V1) + VecDistance(QP, V2);
   double Area      = S->Panels[E->iPPanel]->Area;
   if ( E->iQM!=-1 )
    { double *QM = S->Vertices + 3*E->iQM;
      Perimeter += VecDistance(QM, V1) + VecDistance(QM, V2);
      Area += S->Panels[E->iMPanel]->Area;
    }

   ES.Signature[0] = round( Length     / DistanceQuantum );
   ES.Signature[1] = round( Perimeter  / DistanceQuantum );
   ES.Signature[2] = round( sqrt(Area) / DistanceQuantum );
}

struct EdgeSigHash 
 { long operator() (const EdgeSignature &EdgeSig) const 
    { return scuff::JenkinsHash( (const char *)EdgeSig.Signature, sizeof(EdgeSignature)); }
 };

typedef struct
 { bool operator()(const EdgeSignature &EdgeSig1, const EdgeSignature &EdgeSig2) const  
     { return !memcmp( (const void *)EdgeSig1.Signature, 
                       (const void *)EdgeSig2.Signature, 
                       sizeof(EdgeSignature)); 
     } 
 } EdgeSigCmp; 

#ifdef HAVE_CXX11 
 typedef unordered_map<EdgeSignature, iVec, EdgeSigHash, EdgeSigCmp> SimilarEdgeTable;
#elif defined(HAVE_TR1) 
 typedef tr1::unordered_map<EdgeSignature, iVec, EdgeSigHash, EdgeSigCmp> SimilarEdgeTable;
#else 
 typedef map<EdgeSignature, iVec, EdgeSigCmp> SimilarEdgeTable;
#endif

void AddSimilarEdge(SimilarEdgeTable &SETable, EdgeSignature EdgeSig, int ne)
{ 
  SimilarEdgeTable::iterator it = SETable.find(EdgeSig);
  if (it!=SETable.end())
   it->second.push_back(ne);
  else
   SETable.insert( std::pair<EdgeSignature, iVec>(EdgeSig, iVec(1,ne)));
}

iVec *GetSimilarEdges(SimilarEdgeTable &SETable, EdgeSignature EdgeSig)
{
  SimilarEdgeTable::iterator it=SETable.find(EdgeSig);
  return (it==SETable.end() ? 0 : &(it->second));
}

iVec *GetSimilarEdges(SimilarEdgeTable &SETable, RWGSurface *S, int ne, double DistanceQuantum)
{ EdgeSignature ES;
  GetEdgeSignature(S, ne, DistanceQuantum, ES);
  return GetSimilarEdges(SETable, ES);
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
#define EQUIV             0
#define INEQUIV_LENGTH    1
#define INEQUIV_RADIUS    2
#define INEQUIV_FULLHALF  3
#define INEQUIV_DIPOLE    4
#define INEQUIV_V12       5
#define INEQUIV_QPM       6
static const char *EEResultNames[]=
{ "Equivalent", "length", "radius", "full/half", "dipole", "V12", "QPM" };
#define NUMEERESULTS (sizeof(EEResultNames)/sizeof(EEResultNames[0]))
int TestEdgeEquivalence(RWGSurface *Sa, size_t nea, RWGSurface *Sb, size_t neb,
                        double DistanceQuantum, GTSignature *Gb2aSig, bool FlipQPM=false, bool FlipV12=false)
{ 
  RWGEdge *Ea = Sa->GetEdgeByIndex(nea);
  RWGEdge *Eb = Sb->GetEdgeByIndex(neb);

  if (FlipQPM && Eb->iMPanel==-1) return INEQUIV_FULLHALF;

#define EERELTOL 1.0e-6
  double LengthTol = EERELTOL*Ea->Length;

  // stage-1 check: gross statistics must match
  if ( fabs(Ea->Length - Eb->Length) > LengthTol ) return INEQUIV_LENGTH;
  if ( fabs(Ea->Radius - Eb->Radius) > LengthTol ) return INEQUIV_RADIUS;
  if ( (Ea->iMPanel==-1) != (Eb->iMPanel==-1) )    return INEQUIV_FULLHALF;

  // stage-2 check: magnitudes of dipole moments must match
  double *Va[4], Pa[3], *Vb[4], Pb[3];
  GetPVector(Sa, Ea, Va, Pa);
  GetPVector(Sb, Eb, Vb, Pb, FlipQPM, FlipV12);
  double Na = VecNorm(Pa), Nb = VecNorm(Pb);
  if ( fabs(Na-Nb) > EERELTOL*0.5*(Na+Nb) ) return INEQUIV_DIPOLE;

  // stage-3 check: Compute the "standard transformations" GTa, GTb
  // for edges Ea, Eb, then check that Ea = Ta^{-1} Tb(Eb).
   GTransformation GStandard_b = GetStandardRotation(Sb, neb, FlipQPM, FlipV12);
   GTransformation GStandard_a = GetStandardRotation(Sa, nea);
   GTransformation Gb2a = -GStandard_a + GStandard_b;
   GetGTSignature(Gb2a,DistanceQuantum,Gb2aSig);

  // make sure {V1a,V2a} = {V1bPrime,V2bPrime} as sets
  double VbPrime[4][3];
  Gb2a.Apply(Vb[1], VbPrime[1]);
  Gb2a.Apply(Vb[2], VbPrime[2]);
printf("VA1={%+g,%+g,%+g} , VA2={%+g,%+g,%+g}\n",Va[1][0],Va[1][1],Va[1][2],Va[2][0],Va[2][1],Va[2][2]);
printf("VB1={%+g,%+g,%+g} , VB2={%+g,%+g,%+g}\n",VbPrime[1][0],VbPrime[1][1],VbPrime[1][2],VbPrime[2][0],VbPrime[2][1],VbPrime[2][2]);
  bool V12Direct = (    VecDistance(VbPrime[1], Va[1])<=LengthTol
                    &&  VecDistance(VbPrime[2], Va[2])<=LengthTol
                   );
  bool V12Twist  = (    VecDistance(VbPrime[1], Va[2])<=LengthTol
                    &&  VecDistance(VbPrime[2], Va[1])<=LengthTol
                   );
  if (!V12Direct && !V12Twist) return INEQUIV_V12;
 
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
Gb2a.Apply(Vb[0], VbPrime[0]);
Gb2a.Apply(Vb[3], VbPrime[3]);
printf("VAP={%+g,%+g,%+g} , VAM={%+g,%+g,%+g}\n",Va[0][0],Va[0][1],Va[0][2],Va[3][0],Va[3][1],Va[3][2]);
printf("VBP={%+g,%+g,%+g} , VBM={%+g,%+g,%+g}\n",VbPrime[0][0],VbPrime[0][1],VbPrime[0][2],VbPrime[3][0],VbPrime[3][1],VbPrime[3][2]);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  // make sure QPa==QPb, QMa==QMb
  Gb2a.Apply(Vb[0], VbPrime[0]);
  if ( VecDistance(VbPrime[0], Va[0]) > LengthTol ) return INEQUIV_QPM;
  if (Ea->iMPanel!=-1)
   { Gb2a.Apply(Vb[3], VbPrime[3]);
     if ( VecDistance(VbPrime[3], Va[3]) > LengthTol ) return INEQUIV_QPM;
   }

  return EQUIV;
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

ChildEdgeList GetChildEdgeList(RWGSurface *S, int neParent, double DistanceQuantum, SimilarEdgeTable &SETable, size_t Results[NUMEERESULTS])
{ 
  ChildEdgeList ChildEdges;
  ChildEdgeData CEData;
  if (TestEdgeEquivalence(S, neParent, S, neParent, DistanceQuantum, &(CEData.GTSig), true))
   { CEData.neChild=neParent;
     CEData.Sign=-1;
     ChildEdges.push_back(CEData);
   }
  int NE = S->NumEdges;

  iVec *neChildCandidates=GetSimilarEdges(SETable, S, neParent, DistanceQuantum);
  int ncMax = neChildCandidates ? neChildCandidates->size() : NE;
  for(int nc=0; nc<ncMax; nc++)
   { 
     CEData.neChild = neChildCandidates ? (*neChildCandidates)[nc] : nc;
     if (CEData.neChild==neParent) continue;

     int Status1=TestEdgeEquivalence(S, neParent, S, CEData.neChild, DistanceQuantum, &(CEData.GTSig), false);
     if (Status1==EQUIV)
      { CEData.Sign=1;
        ChildEdges.push_back(CEData);
      }
     int Status2=TestEdgeEquivalence(S, neParent, S, CEData.neChild, DistanceQuantum, &(CEData.GTSig), true);
     if (Status2==EQUIV)
      { CEData.Sign=-1;
        ChildEdges.push_back(CEData);
      }

     int Status3=TestEdgeEquivalence(S, neParent, S, CEData.neChild, DistanceQuantum, &(CEData.GTSig), false, true);
     if (Status3==EQUIV)
      { CEData.Sign=1;
        ChildEdges.push_back(CEData);
      }
     int Status4=TestEdgeEquivalence(S, neParent, S, CEData.neChild, DistanceQuantum, &(CEData.GTSig), true, true);
     if (Status4==EQUIV)
      { CEData.Sign=-1;
        ChildEdges.push_back(CEData);
      }

#pragma omp critical
{
     Results[Status1]++;
     Results[Status2]++;
     Results[Status3]++;
     Results[Status4]++;
}
   }
  return ChildEdges;
}

// this routine is used only for testing
ChildEdgeList GetChildEdgeList(RWGSurface *S, int neParent, double DistanceQuantum, size_t Results[NUMEERESULTS])
{ 
  SimilarEdgeTable SETable;
  for(int ne=0; ne<S->NumEdges; ne++)
   { EdgeSignature ES;
     GetEdgeSignature(S, ne, DistanceQuantum, ES);
     AddSimilarEdge(SETable, ES, ne);
   }
  return GetChildEdgeList(S, neParent, DistanceQuantum, SETable, Results);
}

ChildEdgeList *GetChildren(RWGGeometry *G, int ns, double DistanceQuantum)
{
  RWGSurface *S = G->Surfaces[ns];
  int NE = S->NumEdges;

  // initial order-N step to construct edge signature map
  // SETable[ES][0, 1, ...] = indices of edges with edge signature ES
  SimilarEdgeTable SETable;
  for(int ne=0; ne<NE; ne++)
   { EdgeSignature ES;
     GetEdgeSignature(S, ne, DistanceQuantum, ES);
     AddSimilarEdge(SETable, ES, ne);
   }

  int NumThreads=1;
  CheckEnv("SCUFF_EEP_THREADS", &NumThreads);
  size_t EEResults[NUMEERESULTS]; // histogram of edge-equivalence test results
  memset(EEResults, 0, NUMEERESULTS*sizeof(size_t)); 
  ChildEdgeList *Children = new ChildEdgeList[NE];
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int ne=0; ne<NE; ne++)
   Children[ne] = GetChildEdgeList(S, ne, DistanceQuantum, SETable, EEResults);

  if (G->LogLevel>=SCUFF_VERBOSE2)
   { 
     int NumParents=0, NumChildren[2]={0,0};
     for(int ne=0; ne<NE; ne++)
      { NumParents += (Children[ne].size() > 0 ? 1 : 0);
        for(size_t nc=0; nc<Children[ne].size(); nc++)
         NumChildren[Children[ne][nc].Sign==1 ? 0 : 1]++;
      }

     Log(" Surface %s: %i/%i edges are parents (%.0g %%)",S->Label, NumParents,NE,100.0*((double)NumParents)/((double)NE));
     Log(" results of equivalent-edge detection:");
     Log("  Equivalent            : %6i (%i,%i)",EEResults[EQUIV],NumChildren[0],NumChildren[1]);
     for(size_t nr=EQUIV+1; nr<NUMEERESULTS; nr++)
      if (EEResults[nr]) 
       Log("Inequivalent(%10s): %6i",EEResultNames[nr],EEResults[nr]);
   }

  return Children;
}

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/* The remainder of the file contains routines for analyzing    */
/* relationships between *pairs* of edge pairs, i.e. relations  */
/* of the form (neaParent,nebParent) <--> (neaChild,nebChild).  */
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/******************************************************************/
/* General-purpose helper routines that convert back and forth    */
/* between pairs of edges indices and indices of edge pairs.      */
/******************************************************************/
EdgePair EquivalentEdgePairTable::GetEdgePair(int neParent, int neChild)
{ return neParent*NERadix + neChild; }

int iabs(int x) { return (x<0 ? -x : x); }

bool EquivalentEdgePairTable::ResolveEdgePair(EdgePair Pair, int *neParent, int *neChild)
{
  bool SignFlip = Pair<0;
  if (SignFlip) Pair*=-1;
  *neParent = Pair / NERadix;
  *neChild =  Pair % NERadix;
  return SignFlip;
}

/******************************************************************************/
/* If GTSig is the signature of a GTransformation GT, then                    */
/* ParentChildTable[GTSig] is a list of GT-related {parent, child} pairs,     */
/* i.e. for each {parent, child} pair in the list we have parent = GT(child). */
/******************************************************************************/
struct GTSigHash 
 { long operator() (const GTSignature &GTSig) const 
    { return scuff::JenkinsHash( (const char *)GTSig.Signature, sizeof(GTSignature)); } 
 };

typedef struct 
 { bool operator()(const GTSignature &GTSig1, const GTSignature &GTSig2) const  
     { return !memcmp( (const void *)GTSig1.Signature, 
                       (const void *)GTSig2.Signature, 
                       sizeof(GTSignature)); 
     } 
 } GTSigCmp; 

#ifdef HAVE_CXX11 
 typedef unordered_map<GTSignature, EdgePairList, GTSigHash, GTSigCmp> ParentChildTable; 
#elif defined(HAVE_TR1) 
 typedef tr1::unordered_map<GTSignature, EdgePairList, GTSigHash, GTSigCmp> ParentChildTable;
#else 
 typedef map<GTSignature, EdgePairList, GTSigCmp> ParentChildTable;
#endif

void AddParentChildPair(ParentChildTable &PCTable, GTSignature GTSig, EdgePair PCPair)
{ 
  ParentChildTable::iterator pctIterator=PCTable.find(GTSig);
  if (pctIterator!=PCTable.end())
   pctIterator->second.push_back(PCPair);
  else
   PCTable.insert( std::pair<GTSignature, EdgePairList>(GTSig, EdgePairList(1,PCPair)));
}

EdgePairList *GetParentChildPairs(ParentChildTable &PCTable, GTSignature GTSig)
{ 
  ParentChildTable::iterator pctIterator=PCTable.find(GTSig);
  return (pctIterator==PCTable.end() ? 0 : &(pctIterator->second));
}

/****************************************************************************/
/* constructor helper routine to add a single equivalent edge pair **********/
/****************************************************************************/
void EquivalentEdgePairTable::AddEquivalentEdgePair(EdgePair ParentPair, EdgePair ChildPair,
                                                    bool SignFlip)
{
  if (ChildPair == ParentPair) return;

#pragma omp critical
 { if (!IsReduced[ParentPair] && !IsReduced[ChildPair])
    { IsReduced[ChildPair]=true;
      if (SignFlip) ChildPair*=-1;
      IrreducibleEdgePairMap::iterator it=IEPMap.find(ParentPair);
      if (it!=IEPMap.end())
       it->second.push_back(ChildPair);
      else
       IEPMap.insert(std::pair<EdgePair, EdgePairList>(ParentPair,EdgePairList(1,ChildPair)));
    }
 }
}

void EquivalentEdgePairTable::AddEquivalentEdgePair(int neaParent, int nebParent,
                                                    int neaChild, int nebChild,
                                                    bool SignFlip)
{ 
  if(nsa==nsb)
   { if (nebParent<neaParent)
      { int temp=nebParent; nebParent=neaParent; neaParent=temp; }
     if (nebChild<neaChild)
      { int temp=nebChild; nebChild=neaChild; neaChild=temp; }
   }

  AddEquivalentEdgePair( GetEdgePair(neaParent, nebParent),
                         GetEdgePair(neaChild, nebChild),
                         SignFlip
                       );
}

void EquivalentEdgePairTable::AddEquivalentEdgePair(int neaParent, int neaChild,
                                                    EdgePair nebPair)
{ 
  bool aFlipped = (neaChild < 0);
  int nebParent, nebChild;
  bool bFlipped = ResolveEdgePair( iabs(nebPair), &nebParent, &nebChild);
  AddEquivalentEdgePair( neaParent, nebParent, iabs(neaChild), nebChild,
                         aFlipped!=bFlipped );
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

/******************************************************************/
/* EquivalentEdgePairTable class constructor: Construct a table   */
/* of equivalent edge pairs for two surfaces in an RWG geometry.  */
/******************************************************************/
EquivalentEdgePairTable::EquivalentEdgePairTable(RWGGeometry *_G, int _nsa, int _nsb, char *EEPTFileName)
 : G(_G), nsa(_nsa), nsb(_nsb)
{
  int NEA=G->Surfaces[nsa]->NumEdges, NEB=G->Surfaces[nsb]->NumEdges;

  // set some internal class data fields
  NERadix = (NEA > NEB ? NEA : NEB);

  // try to import table from file
  if (EEPTFileName)
   { char *ErrMsg=Import(EEPTFileName);
     if (ErrMsg) 
      Warn(ErrMsg);
     else
      { Log("successfully read EEPTable from file %s",EEPTFileName);
        return;
      }
   }

  /*----------------------------------------------------------------*/
  /* loops in this routine are single-threaded by default           */
  /*----------------------------------------------------------------*/
  int NumThreads=1;
  CheckEnv("SCUFF_EEP_THREADS", &NumThreads);

  /*----------------------------------------------------------------*/
  /* Step 1: For each of the two surfaces, identify pairs           */
  /*         of equivalent edges within the surface.                */
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
  /*----------------------------------------------------------------*/
  double DistanceQuantum = 0.1 * GetMinPanelRadius(G);
  ChildEdgeList *aChildren  = GetChildren(G, nsa, DistanceQuantum);
  ChildEdgeList *bChildren  = (nsa==nsb) ? aChildren : GetChildren(G, nsb, DistanceQuantum);

  /*----------------------------------------------------------------*/
  /* Step 2: To accelerate Step 4 below, we pause to build a        */
  /*         "Parent-Child table" for surface *nsb. This is a table */
  /*         of equivalent {Parent,Child} edge pairs sorted by      */
  /*         G-transformation:                                      */
  /*         bPCTable[GT] = list of all GT-related parent-child     */
  /*         edge pairs on surface #nsb.                            */
  /*         A pair is "GT-related" if GT(Child) = Parent.          */
  /*----------------------------------------------------------------*/
  ParentChildTable bPCTable;
  for(int nebParent=0; nebParent<NEB; nebParent++)
   for(size_t nc=0; nc<bChildren[nebParent].size(); nc++)
    AddParentChildPair(bPCTable, bChildren[nebParent][nc].GTSig,
                       bChildren[nebParent][nc].Sign * GetEdgePair(nebParent,bChildren[nebParent][nc].neChild));

  /*----------------------------------------------------------------*/
  /* second pass to identify equivalent off-diagonal edge pairs.    */
  /* If neaChild is a child of neaParent while nebChild is a child  */
  /* of nebParent, then (neaChild,nebChild) is equivalent to        */
  /* (neaParent, nebParent) iff the GTransform that transforms      */
  /* neaChild into neaParent is identical to the GTransform that    */
  /* transforms nebChild into nebParent.                            */
  /*----------------------------------------------------------------*/
  IsReduced.resize(NERadix*NERadix,false);

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
  if (G->LogLevel>=SCUFF_VERBOSE2) 
   Log("  Step 2: Identifying equivalent edge pairs (%i threads)",NumThreads);
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int neaParent=0; neaParent<NEA; neaParent++)
   for(size_t nc=0; nc<aChildren[neaParent].size(); nc++)
    { int neaChild  = aChildren[neaParent][nc].Sign * aChildren[neaParent][nc].neChild;
      EdgePairList *bPairs = GetParentChildPairs(bPCTable, aChildren[neaParent][nc].GTSig);
      if (bPairs==0) continue;
      for(size_t n=0; n<bPairs->size(); n++)
       AddEquivalentEdgePair(neaParent, neaChild, (*bPairs)[n]);
    }

  /*----------------------------------------------------------------*/
  /* report some statistics on equivalent edge pairs to log file    */
  /*----------------------------------------------------------------*/
  int NumParentPairs = IEPMap.size();
  int NumChildPairs  = 0, NumChildPairs2 = 0;
  for(int nePair=0; nePair<NERadix*NERadix; nePair++)
   if (IsReduced[nePair]) NumChildPairs++;
  for(IrreducibleEdgePairMap::iterator it=IEPMap.begin(); it!=IEPMap.end(); it++)
   NumChildPairs2+=it->second.size();
 
  int NEPairs = (nsa==nsb ? NEA*(NEA+1)/2 : NEA*NEB);
  Log(" Of %u total edge-edge pairs on surfaces (%i,%i) (%s,%s):",NEPairs,nsa,nsb,G->Surfaces[nsa]->Label,G->Surfaces[nsb]->Label);
  Log("    %u are children (savings of %.1f %%)",NumChildPairs,100.0*((double)NumChildPairs)/((double)NEPairs));
  if (NumChildPairs2!=NumChildPairs)
   Log(" ** warning: child pair counts disagree (%i,%i)",NumChildPairs,NumChildPairs2);
  if (G->LogLevel>SCUFF_VERBOSE2)
   { Log("    %u are parents (%.1f %%)",NumParentPairs, 100.0*((double)NumParentPairs) / ((double)NEPairs));
//     Log("    %u are unicorns (should be %u)",IEPMap->size(), NEPairs - NumParentPairs - NumChildPairs);
   }
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void EquivalentEdgePairTable::Export(char *FileName)
{
  char FileNameBuffer[100];
  if (!FileName)
   { if (G==0) return;
     FileName = FileNameBuffer;
     snprintf(FileNameBuffer,100,"%s.EEPTable",GetFileBase(G->GeoFileName));
   }
  FILE *f = fopen(FileName,"w");
  if (!f) 
   { Warn("could not open file %s (skipping edge-pair table export)",FileName); 
     return;
   }

  for(IrreducibleEdgePairMap::iterator it=IEPMap.begin(); it!=IEPMap.end(); it++)
   { int neaParent, nebParent;
     ResolveEdgePair(it->first, &neaParent, &nebParent);
     fprintf(f,"{%i,%i}",neaParent,nebParent);
     for(size_t nREP=0; nREP<it->second.size(); nREP++)
      { int neaChild, nebChild;
        bool SignFlip = ResolveEdgePair(it->second[nREP], &neaChild, &nebChild);
        fprintf(f," %c{%i,%i}",SignFlip ? '-' : '+',neaChild,nebChild);
      }
     fprintf(f,"\n");
   }

  fclose(f);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
// read a string of the form {n1,n2} from f.
// return codes: 
// 0: fail because end of file
// 1: fail for other reason
// 2: success, at end of file
// 3: success, at end of line
// 4: success, neither of the above
int ReadIndexPair(FILE *f, int *n1, int *n2)
{
  if ( 2 != fscanf(f," {%i,%i}",n1,n2) )
   return feof(f) ? 0 : 1;
  char c=' ';
  while ( isspace(c) && c!='\n' )
   c=fgetc(f);
  if (feof(f))
   return 2;
  if (c=='\n')
   return 3;
  ungetc(c,f);
  return 4;
}

char *EquivalentEdgePairTable::Import(char *FileName)
{ 
  FILE *f=fopen(FileName,"r");
  if (!f) return vstrdup("could not open equivalent edge-pair file %s",FileName);
  int LineNum=1;
  bool AtTopOfLine=true;
  int neaParent, nebParent;
  char *ErrMsg=0;
  IsReduced.resize(NERadix*NERadix,false);
  while( !feof(f) )
   { int nea, neb;
     int Status=ReadIndexPair(f, &nea, &neb);
     if (Status==1)
      ErrMsg=vstrdup("%s:%i: syntax error",FileName,LineNum);
     if (Status<=1)
      break;
     if (AtTopOfLine)
      neaParent=nea, nebParent=neb;
     else
      AddEquivalentEdgePair(neaParent, nebParent, nea, neb);
     if (Status==2) break;
     AtTopOfLine=(Status==3);
     if (AtTopOfLine) LineNum++;
   }
  fclose(f);
  return ErrMsg;
}

} // namespace scuff
