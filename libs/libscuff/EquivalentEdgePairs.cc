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
      { VecScale(xyzHat[0], -1.0); 
        VecScale(xyzHat[1], -1.0); 
      }
   }
}

/***************************************************************/
/* given a (full or half) RWG edge E, compute the standard     */
/* rotation that operates on E to put it into standard position*/
/***************************************************************/
GTransformation GetStandardRotation(RWGSurface *S, int ne)
{
  RWGEdge *E = S->Edges[ne];
  double *QP = S->Vertices + 3*E->iQP;
  double *V1 = S->Vertices + 3*E->iV1;
  //double *V2 = S->Vertices + 3*E->iV2;
  double *X0 = E->Centroid;
  
  double MX0[3];
  MX0[0]=-X0[0];
  MX0[1]=-X0[1];
  MX0[2]=-X0[2];
  GTransformation GTA(MX0);
  double VScratch1[3], VScratch2[3], nHat[3];
  VecCross( VecSub(V1,X0,VScratch1), VecSub(QP, V1, VScratch2), nHat);
  VecNormalize(nHat);
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
/* comparing transformations                                   */
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

// possible results of a call to TestEdgeEquivalence: either the
// edges are equivalent (result==0) or not (result>=1),
// with the index of the result indicating which criterion tipped
// us off to detecting the inequivalence (lower=faster)
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

/***************************************************************/
/* Given two (full or half) RWG edges {Ea, Eb}, determine      */
/* whether or not they are equivalent, i.e Ea is the result    */
/* of applying a rigid geometric transformation (displacement+ */
/* rotation) to Eb.                                            */
/* If the edges are equivalent, return 0 and set Gb2a = the    */
/* transform, i.e. Ea = Gb2a(Eb).                              */
/* if the edges are inequivalent, return a non-zero size_teger    */
/* code indicating how we determined this.                     */
/* If PPFileName is non-null, write visualization data to      */
/* file (for debugging purposes).                              */
/***************************************************************/
#define EERELTOL 1.0e-6
bool OldEEPMethod=true;
int TestEdgeEquivalence(RWGSurface *Sa, size_t nea, RWGSurface *Sb, size_t neb,
                        double DistanceQuantum, GTSignature *Gb2aSig)
{ 
  RWGEdge *Ea = Sa->GetEdgeByIndex(nea);
  RWGEdge *Eb = Sb->GetEdgeByIndex(neb);

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

  // stage-3 check: Compute the G-transformation that (a) displaces Eb
  // until its centroid agrees with that of Ea, (b) rotates Eb so that
  // the X,Y,Z coordinate system of its positive panel is aligned with
  // that of Ea. Then confirm that all 4 vertices of the transformed Eb 
  // coincide with those of Ea.
  double aBasis[3][3], bBasis[3][3];
  GetXYZSystem(Va, Ea->Centroid, Sa->Panels[Ea->iPPanel], aBasis);
  GetXYZSystem(Vb, Eb->Centroid, Sb->Panels[Eb->iPPanel], bBasis);

  // to transform Eb into Ea:
  //  (1) translate through -1.0 * Eb->Centroid (so Eb is centered at the origin)
  //  (2) Rotate the xyz system of Eb into that of Ea
  //  (3) translate through +1.0 * Ea->Centroid
static bool Init=true;
if (Init)
 { Init=false;
   OldEEPMethod = !CheckEnv("SCUFF_NEW_EEPMETHOD");
 }

GTransformation Gb2a;
if (OldEEPMethod)
{
  GTransformation GT1(Eb->Centroid);
  GTransformation GT2;
  for(int Mu=0; Mu<3; Mu++)
   for(int Nu=0; Nu<3; Nu++)
    GT2.M[Mu][Nu] = aBasis[0][Mu]*bBasis[0][Nu] + aBasis[1][Mu]*bBasis[1][Nu] + aBasis[2][Mu]*bBasis[2][Nu];
  GTransformation GT3(Ea->Centroid);
  Gb2a = GT3 + GT2 - GT1;
  //GetGTSignature(Gb2a,DistanceQuantum,Gb2aSig);

  VecSub(Ea->Centroid, Eb->Centroid, GT2.DX);
  GetGTSignature(GT2,DistanceQuantum,Gb2aSig);
}
else
 { GTransformation GStandard_b = GetStandardRotation(Sb, neb);
   GTransformation GStandard_a = GetStandardRotation(Sa, nea);
   Gb2a = -GStandard_a + GStandard_b;
   GetGTSignature(Gb2a,DistanceQuantum,Gb2aSig);
 }

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

int TestEdgeEquivalence(RWGGeometry *G, int nfea, int nfeb, double DistanceQuantum, GTSignature *Gb2aSig)
{ int nea, nsa; RWGSurface *Sa = G->ResolveEdge(nfea, &nsa, &nea);
  int neb, nsb; RWGSurface *Sb = G->ResolveEdge(nfeb, &nsb, &neb);
  return TestEdgeEquivalence(Sa, nea, Sb, neb, DistanceQuantum, Gb2aSig);
}

/******************************************************************/
/* Given a single (full or half) RWG edge (the "parent" edge) in  */
/* a geometry, construct a list of all "child" edges, i.e.        */
/* edges equivalent to the parent and lying after it in the global*/
/* list of edges.                                                 */
/*                                                                */
/* On entry, nfeParent is the full edge index of the parent edge. */
/* Returns a d-sorted list of tuples {d, nfeChild, GTSig}, where  */
/*  nfeChild = (full) edge index of child edge                    */
/*  GTSig    = transform that operates on child edge to yield     */
/*             parent edge                                        */
/*  d        = child-parent centroid-centroid distance, quantized */
/*             in units of 1/10 the minimal panel radius          */
/******************************************************************/
typedef struct ChildEdgeData
 { int nfeChild;
   int Sign;
   GTSignature GTSig;
 } ChildEdgeData;

typedef vector<ChildEdgeData> ChildEdgeList;

ChildEdgeList GetChildEdgeList(RWGGeometry *G, int nfeParent, double DistanceQuantum, int Results[NUMEERESULTS])
{ 
  ChildEdgeList ChildEdges;
  ChildEdgeData CEData;
  for(CEData.nfeChild=nfeParent+1; CEData.nfeChild<G->TotalEdges; CEData.nfeChild++)
   { int Status=TestEdgeEquivalence(G, nfeParent, CEData.nfeChild, DistanceQuantum, &(CEData.GTSig));
     if (Status==EQUIV_PLUS || Status==EQUIV_MINUS)
      { CEData.Sign = (Status==EQUIV_MINUS ? -1 : 1);
        ChildEdges.push_back(CEData);
      }
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
#pragma omp critical
     Results[Status]++;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
   }
  return ChildEdges;
}

/******************************************************************/
/* General-purpose helper routines that convert back and forth    */
/* between {edge, edge} pairs and integers indexing the list of   */
/* all possible edge pairs.                                       */
/* Note that the two indices that define an edge pair are always  */
/* *full* edge indices, i.e. indices that specify both a surface  */
/* within an RWG geometry as well as an edge on that surface.     */
/******************************************************************/
int EquivalentEdgePairTable::GetEdgePairIndex(int nfeParent, int nfeChild)
{
  return nfeParent*NFE + nfeChild;
}

int iabs(int x) { return (x<0 ? -x : x); }

void EquivalentEdgePairTable::ResolveEdgePairIndex(int nPair, int *nfeParent, int *nfeChild)
{
  int Sign = (nPair < 0 ? -1 : 1);
  nPair = iabs(nPair);
 *nfeParent =         nPair / NFE;
  *nfeChild = Sign * (nPair % NFE);
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
 { if (!HasParentFlag[ParentPair] && !HasParentFlag[ChildPair])
    { HasParentFlag[ChildPair]=true;
      if (SignFlip) ChildPair*=-1;
      if (ChildPairListArray)
       ChildPairListArray[ParentPair].push_back(ChildPair);
      else
       { ChildPairListMap::iterator it=CPLMap.find(ParentPair);
         if (it!=CPLMap.end())
          it->second.push_back(ChildPair);
         else
          { iVec NewChildPairList(1,ChildPair);
            CPLMap.insert(std::pair<int, iVec>(ParentPair,NewChildPairList));
          }
       }
    }
 }
}

void EquivalentEdgePairTable::AddEquivalentEdgePair(int nfeaParent, int nfebParent,
                                                    int nfeaChild, int nfebChild,
                                                    bool SignFlip)
{ 
  if ( (nfebParent<nfeaParent) || (nfebChild < nfeaChild) )
   return;

  AddEquivalentEdgePair( GetEdgePairIndex(nfeaParent, nfebParent),
                         GetEdgePairIndex(nfeaChild, nfebChild),
                         SignFlip
                       );
}

void EquivalentEdgePairTable::AddEquivalentEdgePair(int nfeaParent, int nfeaChild,
                                                    int nfebPair)
{ 
  int aSign = ( (nfeaChild < 0) ? -1 : 1);
  int bSign = ( (nfebPair  < 0) ? -1 : 1);
  int nfebParent, nfebChild;
  ResolveEdgePairIndex( iabs(nfebPair), &nfebParent, &nfebChild);
  AddEquivalentEdgePair( nfeaParent, nfebParent, iabs(nfeaChild), nfebChild,
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

/******************************************************************/
/******************************************************************/
/******************************************************************/
EquivalentEdgePairTable::EquivalentEdgePairTable(RWGGeometry *_G)
 : G(_G)
{
  // set some internal class data fields
  DistanceQuantum = 0.1 * GetMinPanelRadius(G);
  NFE = G->TotalEdges;

  int NFEPairs=NFE*(NFE+1)/2;
  Log("Detecting equivalences among %i edges (%i pairs)...",NFE,NFEPairs);

  ChildPairListArray=0;
  bool UseArrays = CheckEnv("SCUFF_USE_EEPARRAYS");
  if (UseArrays)
   ChildPairListArray = new iVec[NFE*NFE];
  
  /*----------------------------------------------------------------*/
  /* loops in this routine are single-threaded by default           */
  /*----------------------------------------------------------------*/
  int NumThreads=1;
  CheckEnv("SCUFF_EEP_THREADS", &NumThreads);

  bool WriteEEPLogs = CheckEnv("SCUFF_WRITE_EEPLOGS");
 
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
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
  for(int nfe=0; nfe<NFE; nfe++)
   { Children[nfe] = GetChildEdgeList(G, nfe, DistanceQuantum, EEResults);
     if (Children[nfe].size() > 0) NumParentEdges +=1;
   }
  Log("  %i/%i edges are parents (%e %%)",NumParentEdges,NFE,100.0*((double)NumParentEdges)/((double)NFE));

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
  if (G->LogLevel >= SCUFF_VERBOSE2)
   { Log(" results of equivalent-edge detection:");
     int TotalResults=EEResults[EQUIV_PLUS] + EEResults[EQUIV_MINUS];
     Log("  Equivalent plus       : %6i",EEResults[EQUIV_PLUS]);
     Log("  Equivalent minus      : %6i",EEResults[EQUIV_MINUS]);
     for(size_t nr=2; nr<NUMEERESULTS; nr++)
      if (EEResults[nr]) 
       { TotalResults+=EEResults[nr];
         Log("Inequivalent(%10s): %6i",EEResultNames[nr],EEResults[nr]);
       }
     Log("-------------------------------------------");
     Log("  %15s: %6i (should be %i)","Total",TotalResults,NFEPairs);
   }
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

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

  /*----------------------------------------------------------------*/
  /* second pass to identify equivalent off-diagonal edge pairs.    */
  /* If nfeaPrime is a child of nfea, while nfebPrime is a child of */
  /* nfeb, then (nfeaPrime,nfebPrime) is equivalent to (nfea,nfeb)  */
  /* if the GTransform that takes nfeaPrime to nfea is identical    */
  /* to the GTransform that takes nfebPrime to nfeb.                */
  /*----------------------------------------------------------------*/
  Log("  Step 2: Identifying equivalent edge pairs...");
  HasParentFlag.resize(NFE*NFE,false);

  // handle diagonal pairs first separately
  for(int nfeaParent=0; nfeaParent<NFE; nfeaParent++)
   for(size_t nc=0; nc<Children[nfeaParent].size(); nc++)
    AddEquivalentEdgePair(nfeaParent, nfeaParent,
                          Children[nfeaParent][nc].nfeChild,
                          Children[nfeaParent][nc].nfeChild);

  // construct list of {parent,child=G(parent)} pairs ordered by G
  ParentChildTable PCTable;
  for(int nfeParent=0; nfeParent<NFE; nfeParent++)
   for(size_t nc=0; nc<Children[nfeParent].size(); nc++)
    AddParentChildPair(PCTable, Children[nfeParent][nc].GTSig,
                       Children[nfeParent][nc].Sign
                        * GetEdgePairIndex(nfeParent, Children[nfeParent][nc].nfeChild)
                      );

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

  for(int nfeaParent=0; nfeaParent<NFE; nfeaParent++)
   for(size_t nc=0; nc<Children[nfeaParent].size(); nc++)
    { int nfeaChild  = Children[nfeaParent][nc].Sign * Children[nfeaParent][nc].nfeChild;
      iVec *nfebPairs = GetParentChildPairs(PCTable, Children[nfeaParent][nc].GTSig);
      if (nfebPairs==0) continue;
      for(size_t n=0; n<nfebPairs->size(); n++)
       AddEquivalentEdgePair(nfeaParent, nfeaChild, (*nfebPairs)[n]);
    }

  DestroyParentChildTable(PCTable);
 
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
}

/****************************************************************************/
/* API routines              ************************************************/
/****************************************************************************/
bool EquivalentEdgePairTable::HasParent(int ChildPair)
{ return HasParentFlag[ChildPair]; } 

bool EquivalentEdgePairTable::HasParent(int nfeaChild, int nfebChild)
{ return HasParentFlag[GetEdgePairIndex(nfeaChild, nfebChild)]; }

bool EquivalentEdgePairTable::HasParent(int nsa, int neaChild, int nsb, int nebChild)
{ return HasParentFlag[ GetEdgePairIndex( G->UnResolveEdge(nsa, neaChild),
                                          G->UnResolveEdge(nsb, nebChild) 
                                        )
                      ];
}

iVec *EquivalentEdgePairTable::GetChildren(int ParentPair)
{ 
  if (ChildPairListArray) 
   return &(ChildPairListArray[ParentPair]);

  ChildPairListMap::iterator it=CPLMap.find(ParentPair);
  return (it==CPLMap.end() ? 0 : &(it->second));
}

iVec *EquivalentEdgePairTable::GetChildren(int nfeaParent, int nfebParent)
{ return GetChildren(nfeaParent*NFE + nfebParent); }

iVec *EquivalentEdgePairTable::GetChildren(int nsa, int neaParent, int nsb, int nebParent)
{ return GetChildren( G->UnResolveEdge(nsa, neaParent), G->UnResolveEdge(nsb, nebParent) ); }

int EquivalentEdgePairTable::CountParentPairs()
{ 
  if (!ChildPairListArray) return CPLMap.size();

  int NumParentPairs=0;
  for(int n=0; n<NFE*NFE; n++)
   if (ChildPairListArray[n].size()>0) 
    NumParentPairs++;
  return NumParentPairs;
}

int EquivalentEdgePairTable::CountChildPairs()
{ 
   // first way of counting
   int NumChildPairs[2]={0,0};
   for(int n=0; n<NFE*NFE; n++)
    if (HasParentFlag[n])
     NumChildPairs[0]++;

   // second way of counting
   if (ChildPairListArray)
    for(int n=0; n<NFE*NFE; n++)
     NumChildPairs[1]+=(ChildPairListArray[n].size()>0 ? 1 : 0);
   else
    for(ChildPairListMap::iterator it=CPLMap.begin(); it!=CPLMap.end(); it++)
     NumChildPairs[1]+=it->second.size();

  if (NumChildPairs[0]!=NumChildPairs[1])
   Log("Whoops! NumChildPairs={%i,%i} by different methods",NumChildPairs[0],NumChildPairs[1]); 
  else 
   Log("NumChildPairs={%i,%i} by different methods",NumChildPairs[0],NumChildPairs[1]); 
  return NumChildPairs[0]; 
} 

/***************************************************************/
/***************************************************************/
/***************************************************************/
IrreducibleEdgePairList *EquivalentEdgePairTable::GetIrreducibleEdgePairList(int nsa, int nsb)
{
  IrreducibleEdgePairList *IEPList = new IrreducibleEdgePairList;
  RWGSurface *Sa = G->Surfaces[nsa], *Sb = G->Surfaces[nsb];
  int NEA = Sa->NumEdges, NEB = Sb->NumEdges;
   
  for(int neaParent=0; neaParent<NEA; neaParent++)
   for(int nebParent=(nsa==nsb ? neaParent : 0); nebParent<NEB; nebParent++)
    { 
      int nfeaParent = G->UnResolveEdge(nsa, neaParent);
      int nfebParent = G->UnResolveEdge(nsb, nebParent);
      int ParentPair = GetEdgePairIndex(nfeaParent, nfebParent);

      if (HasParentFlag[ParentPair]) continue;
      IrreducibleEdgePair IEP;
      IEP.ParentPair = ParentPair;
      IEP.ChildPairs = GetChildren(ParentPair);
      IEPList->push_back(IEP);
    }
  return IEPList;
}

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
   { Warn("could not open file %s (skipping edge-pair table export)",FileName); return;}
  for(int ParentPair=0; ParentPair<NFE*NFE; ParentPair++)
   { iVec *Children = GetChildren(ParentPair);
     if (Children==0 || Children->size()==0) continue;
     int neaParent, nebParent;
     ResolveEdgePairIndex(ParentPair, &neaParent, &nebParent);
     fprintf(f,"{%i,%i} ",neaParent, nebParent);
     for(size_t n=0; n<Children->size(); n++)
      { int neaChild, nebChild;
        ResolveEdgePairIndex( (*Children)[n], &neaChild, &nebChild);
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
