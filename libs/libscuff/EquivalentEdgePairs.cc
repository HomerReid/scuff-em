/*
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
 * Homer Reid `           -- 4/2018
 */

#include "EquivalentEdgePairs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <map>

namespace scuff {

#define ISWAP(a,b) { int temp=a; a=b; b=temp; }
#define VSWAP(a,b) { double *temp=a; a=b; b=temp; }

long JenkinsHash(const char *key, size_t len); // in FIBBICache.cc 

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/* Part 1 of this file contains routines for analyzing          */
/* equivalences between pairs of RWG functions, i.e.            */
/* relationships of the form (neParent<-->neChild).             */
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/

/***************************************************************/
/* Part 1A: Routines for handling similar edges.               */
/*  "Similarity" of edges is a coarse version of "equivalence".*/
/*  Two edges are similar if they have identical "signatures", */
/*  where the signature of an edge is a geometric statistic    */
/*  that is fast to compute and which yields a useful key for  */
/*  organizing tables of similar edges. Similarity is a        */
/*  necessary but not sufficient condition for equivalence.    */
/*  A SignedEdgeTable is a list of edges sorted by signature.  */
/***************************************************************/
#define EDGESIGLEN 3
typedef struct { float Signature[EDGESIGLEN]; } EdgeSignature;
EdgeSignature GetEdgeSignature(RWGSurface *S, int ne)
 { 
   RWGEdge *E = S->Edges[ne];
   double *QP = S->Vertices + 3*E->iQP;
   double *V1 = S->Vertices + 3*E->iV1;
   double *V2 = S->Vertices + 3*E->iV2;
   double *QM = E->iQM==-1 ? E->Centroid : S->Vertices + 3*E->iQM;
   
   double QmQ[3];   VecSub(QP, QM, QmQ);
   double VmV[3];   VecSub(V1, V2, VmV);
   double Cross[3]; VecCross(QmQ, VmV, Cross);
   
   EdgeSignature ES;
   ES.Signature[0] = (float) VecNorm(QmQ);
   ES.Signature[1] = (float) VecNorm(VmV);
   ES.Signature[2] = (float) VecNorm(Cross);
   return ES;
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
 typedef unordered_map<EdgeSignature, iVec, EdgeSigHash, EdgeSigCmp> SignedEdgeTable;
#elif defined(HAVE_TR1) 
 typedef tr1::unordered_map<EdgeSignature, iVec, EdgeSigHash, EdgeSigCmp> SignedEdgeTable;
#else 
 typedef map<EdgeSignature, iVec, EdgeSigCmp> SignedEdgeTable;
#endif

void AddSignedEdge(SignedEdgeTable &SETable, int ne, EdgeSignature EdgeSig)
{ 
  SignedEdgeTable::iterator it = SETable.find(EdgeSig);
  if (it!=SETable.end())
   it->second.push_back(ne);
  else
   SETable.insert( std::pair<EdgeSignature, iVec>(EdgeSig, iVec(1,ne)));
}

iVec *GetSimilarEdges(SignedEdgeTable &SETable, EdgeSignature EdgeSig)
{
  SignedEdgeTable::iterator it=SETable.find(EdgeSig);
  return (it==SETable.end() ? 0 : &(it->second));
}

iVec *GetSimilarEdges(SignedEdgeTable &SETable, RWGSurface *S, int ne)
{ return GetSimilarEdges(SETable, GetEdgeSignature(S,ne));
}

/***************************************************************/
/* Part 1B: Routines for detecting equivalent edges.           */
/***************************************************************/

/*--------------------------------------------------------------*/
/* Given a (full or half) RWG edge E, return the rigid          */
/* transformation (rotation+displacement) that operates on E to */
/* put it in "standard position." Standard position means the   */
/* edge lies on the $x$ axis with centroid at the origin and    */
/* positive panel lying in the XY plane with Q1.x > 0.          */
/*--------------------------------------------------------------*/
GTransformation GetStandardTransformation(RWGSurface *S, int ne, bool FlipQPM=false, bool FlipV12=false)
{
  RWGEdge *E = S->Edges[ne];
  double *QP = S->Vertices + 3*E->iQP;
  if (FlipQPM && E->iQM!=-1) QP=S->Vertices + 3*E->iQM;
  double *V1 = S->Vertices + 3*E->iV1;
  double *V2 = S->Vertices + 3*E->iV2;
  double *X0 = E->Centroid;

  if (FlipV12)
   { VSWAP(V1,V2) }
  else if (VecDistance(QP,V2) < VecDistance(QP,V1))
   { VSWAP(V1,V2) };
  
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

/*--------------------------------------------------------------*/
/* TestEdgeMatch applies a GTransformation to edge #ne2 and asks*/
/* if the transformed edge is equivalent to edge #ne1, possibly */
/* with a sign flip.                                            */
/*--------------------------------------------------------------*/
bool Nearby(double *X, double *Y, double Tolerance)
{ return VecDistance(X,Y) <= Tolerance; }

#define EQUIV_PLUS       0
#define EQUIV_MINUS      1
#define INEQUIV_FULLHALF 2
#define INEQUIV_V12      3
#define INEQUIV_QPM      4
const char *EdgeMatchResults[]={"Equiv+", "Equiv-", "FullHalf", "V12", "QPM"};

int TestEdgeMatch(RWGSurface *S1, int ne1, RWGSurface *S2, int ne2, GTransformation T, double *Distance=0)
{
  RWGEdge *E1=S1->Edges[ne1], *E2=S2->Edges[ne2];
  if (Distance) *Distance = VecDistance(E1->Centroid, E2->Centroid);

  if ( (E1->iQM==-1) != (E2->iQM==-1)) return INEQUIV_FULLHALF;

  double *V1[4], *V2[4], V2P[4][3];
  V1[0] = S1->Vertices + 3*E1->iQP;
  V1[1] = S1->Vertices + 3*E1->iV1;
  V1[2] = S1->Vertices + 3*E1->iV2;
  V1[3] = (E1->iQM==-1 ? 0 : S1->Vertices + 3*E1->iQM);
  V2[0] = S2->Vertices + 3*E2->iQP;
  V2[1] = S2->Vertices + 3*E2->iV1;
  V2[2] = S2->Vertices + 3*E2->iV2;
  V2[3] = (E2->iQM==-1 ? 0 : S2->Vertices + 3*E2->iQM);

  for(int nv=0; nv<4; nv++)
   if (V2[nv])
    T.Apply(V2[nv], V2P[nv]);

  double Tolerance = 1.0e-4*fmin(E1->Radius, E2->Radius);

  bool V12Direct = ( Nearby(V1[1], V2P[1], Tolerance) && Nearby(V1[2], V2P[2], Tolerance) );
  bool V12Twist  = ( Nearby(V1[1], V2P[2], Tolerance) && Nearby(V1[2], V2P[1], Tolerance) );
  if (!V12Direct && !V12Twist) return INEQUIV_V12;

  bool QPMatch = Nearby(V1[0], V2P[0], Tolerance);
  bool QMMatch = (E1->iQM==-1 ? true : Nearby(V1[3], V2P[3], Tolerance));

  if (QPMatch && QMMatch) return EQUIV_PLUS;

  if ( E1->iQM!=-1 && Nearby(V1[0],V2P[3],Tolerance) && Nearby(V1[3],V2P[0],Tolerance) )
   return EQUIV_MINUS;

  return INEQUIV_QPM;
}

int TestEdgeMatch(RWGGeometry *G, int ns1, int ne1, int ns2, int ne2, GTransformation T, double *Distance=0)
{ return TestEdgeMatch(G->Surfaces[ns1], ne1, G->Surfaces[ns2],ne2,T,Distance); }

int TestEdgePairPair(RWGGeometry *G, int nsa, int nsb,
                     int neaParent, int neaChild, int nebParent, int nebChild)
{
  RWGSurface *Sa = G->Surfaces[nsa], *Sb = G->Surfaces[nsb];

  GTransformation Ta =
   -GetStandardTransformation(Sa, neaParent) + GetStandardTransformation(Sa,neaChild);
  int Status=TestEdgeMatch(Sb, nebParent, Sb, nebChild, Ta);
  if (Status==EQUIV_PLUS || Status==EQUIV_MINUS)
   return Status;

  GTransformation Tb =
   -GetStandardTransformation(Sb, nebParent) + GetStandardTransformation(Sb,nebChild);
  Status=TestEdgeMatch(Sa, neaParent, Sa, neaChild, Tb);
  return Status;
}

/***************************************************************/
/* Part 1C: GTSignatures.                                      */
/***************************************************************/
#define GTSIGLEN 6
typedef struct { float Signature[GTSIGLEN]; } GTSignature;

GTSignature GetGTSignature(GTransformation GT)
{
  // vertices of standard triangle
  const static double V0[4][3]=
   { 0.1,   1.0, 0.0,
     0.9,   0.0, 0.0,
    -0.9,   0.0, 0.0,
    -0.05, -0.9, 0.0,
   };

  double V[4][3];
  for(int nv=0; nv<4; nv++) 
   GT.Apply(V0[nv], V[nv]);

  double QmQ[3], X0[3];
  VecSub(V[0], V[3], QmQ);
  X0[0] = 0.5*(V[1][0]+V[2][0]);
  X0[1] = 0.5*(V[1][1]+V[2][1]);
  X0[2] = 0.5*(V[1][2]+V[2][2]);

  GTSignature GTSig;
  double Quantum=1.0e-6;
  GTSig.Signature[0] = (float)(Quantum*round(X0[0]/Quantum));
  GTSig.Signature[1] = (float)(Quantum*round(X0[1]/Quantum));
  GTSig.Signature[2] = (float)(Quantum*round(X0[2]/Quantum));
  GTSig.Signature[3] = (float)(Quantum*round(QmQ[0]/Quantum));
  GTSig.Signature[4] = (float)(Quantum*round(QmQ[1]/Quantum));
  GTSig.Signature[5] = (float)(Quantum*round(QmQ[2]/Quantum));
  return GTSig;
}

struct GTSigHash 
 { long operator() (const GTSignature &GTSig) const 
    { return JenkinsHash( (const char *)GTSig.Signature, sizeof(GTSignature)); } 
 }; 

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
int GTSigLen=GTSIGLEN;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
typedef struct 
 { bool operator()(const GTSignature &GTSig1, const GTSignature &GTSig2) const  
     { return !memcmp( (const void *)GTSig1.Signature, 
                       (const void *)GTSig2.Signature, 
                       GTSigLen*sizeof(float)
         //              sizeof(GTSignature)
                     ); 
     }
 } GTSigCmp; 

/***************************************************************/
/* Part 1D: Data structures for pairs of equivalent edges.     */
/*  An "EdgePairData" describes a single pair of equivalent    */
/*  edges. An "EdgePairList" is an unsorted list of            */
/*  EdgePairData structures.                                   */
/***************************************************************/
typedef struct EdgePairData
 { int neParent, neChild;
   bool Flipped;
   GTransformation T;
   GTSignature TSig;
   EdgePairData(int _neParent, int _neChild, bool _Flipped, GTransformation _T, GTSignature _TSig):
    neParent(_neParent), neChild(_neChild), Flipped(_Flipped), T(_T), TSig(_TSig) {}
 } EdgePairData;

typedef vector< vector<EdgePairData> > EdgePairTable;

EdgePairTable CreateEdgePairTable(RWGGeometry *G, int ns)
{
  RWGSurface *S = G->Surfaces[ns];
  int NE = S->NumEdges;

  // initial order-N step to construct a table of edges organized by edge signature
  // SETable[ES][0, 1, ...] = indices of edges with edge signature ES
  SignedEdgeTable SETable;
  for(int ne=0; ne<NE; ne++)
   AddSignedEdge(SETable, ne, GetEdgeSignature(S,ne));

  int NumThreads=1;
  CheckEnv("SCUFF_EEP_THREADS", &NumThreads);
  EdgePairTable Children(NE);
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1), num_threads(NumThreads)
#endif
  for(int neParent=0; neParent<NE; neParent++)
   { GTransformation InvTParent = -GetStandardTransformation(S, neParent);
     iVec *neChildCandidates = GetSimilarEdges(SETable, S, neParent);
     for(size_t n=0; neChildCandidates && n<neChildCandidates->size(); n++)
      { int neChild = (*neChildCandidates)[n];
        //if (neChild<neParent) continue;
        GTransformation T = InvTParent + GetStandardTransformation(S, neChild);
        int Result = TestEdgeMatch(S, neParent, S, neChild, &T);
        if (Result==EQUIV_PLUS || Result==EQUIV_MINUS)
         Children[neParent].push_back(EdgePairData(neParent, neChild, (Result==EQUIV_MINUS), T, GetGTSignature(T)));
      }
   }

  if (G->LogLevel>=SCUFF_VERBOSE2)
   { 
     int NumParents=0, NumChildren[2]={0,0};
     for(int ne=0; ne<NE; ne++)
      { NumParents += (Children[ne].size() > 0 ? 1 : 0);
        for(size_t nc=0; nc<Children[ne].size(); nc++)
         NumChildren[Children[ne][nc].Flipped ? 1 : 0]++;
      }
     int TotalChildren=NumChildren[0] + NumChildren[1];

     Log(" Surface %s: %i/%i edges are parents (%.0g %%)",S->Label, NumParents,NE,100.0*((double)NumParents)/((double)NE));
     Log("             %i/%i edges are children (%.0g %%) {%i,%i}", TotalChildren,NE,100.0*((double)TotalChildren)/((double)NE),NumChildren[0],NumChildren[1]);
/*
     for(size_t nr=EQUIV+1; nr<NUMEERESULTS; nr++)
      if (EEResults[nr]) 
       Log("Inequivalent(%10s): %6i",EEResultNames[nr],EEResults[nr]);
*/
   }

  return Children;
}


#if defined(HAVE_CXX11) 
  typedef unordered_multimap<GTSignature, EdgePairData *, GTSigHash, GTSigCmp> EdgePairMap;
#elif defined(HAVE_TR1) 
  typedef tr1::unordered_multimap<GTSignature, EdgePairData *, GTSigHash, GTSigCmp> EdgePairMap;
#else
  typedef multimap<GTSignature, EdgePairData *, GTSigHash, GTSigCmp> EdgePairMap;
#endif

typedef EdgePairMap::iterator EPMIterator;

EdgePairMap CreateEdgePairMap(EdgePairTable &EPTable)
{ 
  EdgePairMap EPMap;
  for(size_t ne=0; ne<EPTable.size(); ne++)
   for(size_t nc=0; nc<EPTable[ne].size(); nc++)
    EPMap.insert( pair<GTSignature, EdgePairData *>(EPTable[ne][nc].TSig, &(EPTable[ne][nc])));
  return EPMap;
}

/******************************************************************/
/* Part 2: of the file contains routines for analyzing          */
/* relationships between *pairs* of edge pairs, i.e. relations  */
/* of the form (neaParent,nebParent) <--> (neaChild,nebChild).  */
/******************************************************************/

/*--------------------------------------------------------------*/
/* General-purpose helper routines that convert back and forth  */
/* between pairs of edges indices and indices of edge pairs.    */
/*--------------------------------------------------------------*/
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

/****************************************************************************/
/* constructor helper routine to add a single equivalent edge pair **********/
/****************************************************************************/
bool EquivalentEdgePairTable::AddReducedPair(EdgePair ParentPair, EdgePair ChildPair, bool SignFlip)
{
  if (ChildPair == ParentPair) return false;
  if (ChildPair < ParentPair) ISWAP(ChildPair, ParentPair);

//#pragma omp critical
 { if (!IsReduced[ParentPair] && !IsReduced[ChildPair])
    { IsReduced[ChildPair]=true;
      if (SignFlip) ChildPair*=-1;
      ReducedEdgePairMap::iterator it=REPMap.find(ParentPair);
      if (it!=REPMap.end())
       it->second.push_back(ChildPair);
      else
       REPMap.insert(std::pair<EdgePair, EdgePairList>(ParentPair,EdgePairList(1,ChildPair)));
      return true;
    }
 }
 return false;
}

bool EquivalentEdgePairTable::AddReducedPair(int neaParent, int nebParent,
                                             int neaChild, int nebChild,
                                             bool SignFlip)
{ 
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
bool Special=false;
if (neaParent==0 && neaChild==0 && nebParent==1 && nebChild==8)
 { Special=true;
   printf("**Bawonkatage %i %i %i %i!\n",neaParent,neaChild,nebParent,nebChild);
 }
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  if(nsa==nsb)
   { if (nebParent<neaParent) ISWAP(nebParent, neaParent);
     if (nebChild<neaChild) ISWAP(nebChild, neaChild);
   }

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
 if(Special)
  printf("**Helloatage %i %i %i %i!\n",neaParent,neaChild,nebParent,nebChild);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  if (nsa==nsb && (neaParent==nebChild && neaChild==nebParent)) return false;

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
 if(Special)
  printf("**Goofnatage %i %i %i\n",
GetEdgePair(neaParent,nebParent),
GetEdgePair(neaChild,nebChild),SignFlip ? -1 : 1 );
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  return AddReducedPair( GetEdgePair(neaParent, nebParent), GetEdgePair(neaChild, nebChild), SignFlip);
}

bool EquivalentEdgePairTable::AddReducedPair(int neaParent, int neaChild, EdgePair nebPair)
{ 
  bool aFlipped = (neaChild < 0);
  int nebParent, nebChild;
  bool bFlipped = ResolveEdgePair( nebPair, &nebParent, &nebChild);
  return AddReducedPair( neaParent, nebParent, iabs(neaChild), nebChild, aFlipped!=bFlipped );
}

bool EquivalentEdgePairTable::EvaluatePairPair(EdgePairData *aPair, EdgePairData *bPair)
{
  if (    IsReduced[GetEdgePair(aPair->neParent,bPair->neParent)]
       || IsReduced[GetEdgePair(aPair->neChild, bPair->neChild)]
     ) return false;

  int Result=0; // 1,-1,0 for positive match, negative match, non-match

  int aSign=TestEdgeMatch(G, nsa, aPair->neParent, nsb, aPair->neChild, bPair->T);
  if (aSign==EQUIV_PLUS || aSign==EQUIV_MINUS)
   { bool aFlipped = (aSign==EQUIV_MINUS);
     Result = (aFlipped == bPair->Flipped) ? 1 : -1;
   }
  else
   { int bSign=TestEdgeMatch(G, nsb, bPair->neParent, nsb, bPair->neChild, aPair->T);
     if (bSign==EQUIV_PLUS || bSign==EQUIV_MINUS)
      { bool bFlipped = (bSign==EQUIV_MINUS);
        Result = (bFlipped == aPair->Flipped) ? 1 : -1;
      }
   }
if (aPair->neParent==0 && aPair->neChild==0 && bPair->neParent==1 && bPair->neChild==8)
 printf("aSign=%i, Result=%i\n",aSign,Result);
  if (Result==0) return false;
  return AddReducedPair(aPair->neParent, bPair->neParent, aPair->neChild, bPair->neChild, Result==-1);
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
  CheckEnv("GTSIGLEN",&GTSigLen);

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
  //double DistanceQuantum = 0.1 * GetMinPanelRadius(G);
  EdgePairTable aPairTable = CreateEdgePairTable(G, nsa);
  EdgePairTable bPairTable;
  if (nsa!=nsb) 
   bPairTable = CreateEdgePairTable(G,nsb);
  EdgePairMap bPairMap = CreateEdgePairMap( nsa==nsb ? aPairTable : bPairTable );
  
  /*----------------------------------------------------------------*/
  /* Step 2: To accelerate Step 4 below, we pause to build a        */
  /*         "Parent-Child table" for surface *nsb. This is a table */
  /*         of equivalent {Parent,Child} edge pairs sorted by      */
  /*         G-transformation:                                      */
  /*         bPCTable[GT] = list of all GT-related parent-child     */
  /*         edge pairs on surface #nsb.                            */
  /*         A pair is "GT-related" if GT(Child) = Parent.          */
  /*----------------------------------------------------------------*/

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
    for(size_t nc=0; nc<aPairTable[neaParent].size(); nc++)
     AddReducedPair(neaParent, neaParent,
                           aPairTable[neaParent][nc].neChild,
                           aPairTable[neaParent][nc].neChild);

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
bool BFMethod=CheckEnv("SCUFF_BF_EEPS");
FILE **LogFiles = new FILE*[NumThreads];
for(int nt=0; nt<NumThreads; nt++)
 LogFiles[nt]=vfopen("/tmp/P%i.log","w",nt);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

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
   for(size_t naPair=0; naPair<aPairTable[neaParent].size(); naPair++)
    { EdgePairData *aPair = &(aPairTable[neaParent][naPair]);
      pair<EPMIterator, EPMIterator> range = bPairMap.equal_range(aPair->TSig);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
int Hits=0, Candidates=0;
if (BFMethod) { range.first=bPairMap.begin(); range.second=bPairMap.end(); }
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
      for(EPMIterator it=range.first; it!=range.second; it++, Candidates++)
       //EvaluatePairPair(aPair, it->second);
       if (EvaluatePairPair(aPair, it->second)) Hits++;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
FILE *f=LogFiles[GetThreadNum()];
fprintf(f,"{%5i,%5i} (",neaParent, aPair->neChild);
for(int n=0; n<GTSigLen; n++)
 fprintf(f,"%g ",aPair->TSig.Signature[n]);
fprintf(f,") %i/%i hits/candidates\n",Hits,Candidates);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
    }

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
for(int nt=0; nt<NumThreads; nt++) fclose(LogFiles[nt]);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

  /*----------------------------------------------------------------*/
  /*----------------------------------------------------------------*/
  /*----------------------------------------------------------------*/
  tr1::unordered_map<int, int> Consolidators;
  for(ReducedEdgePairMap::iterator it1=REPMap.begin(); it1!=REPMap.end(); it1++)
   for(ReducedEdgePairMap::iterator it2=it1; it2!=REPMap.end(); it2++)
    { if (it2==it1) continue;
      int ParentPair1 = it1->first, ParentPair2 = it2->first;
      if (      Consolidators.find(ParentPair1) != Consolidators.end()
           ||   Consolidators.find(ParentPair2) != Consolidators.end()
         ) continue;
      int nea1, neb1; ResolveEdgePair(ParentPair1, &nea1, &neb1);
      int nea2, neb2; ResolveEdgePair(ParentPair2, &nea2, &neb2);
      int Result=TestEdgePairPair(G, nsa, nsb, nea1, nea2, neb1, neb2);
      if (Result==EQUIV_PLUS)
       Consolidators.insert( std::pair<int,int>(ParentPair2, ParentPair1));
      else if (Result==EQUIV_MINUS)
       Consolidators.insert( std::pair<int,int>(ParentPair2, -ParentPair1));
    }

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
{
Export("/tmp/Before.EEPTable");
  /*----------------------------------------------------------------*/
  /* report some statistics on equivalent edge pairs to log file    */
  /*----------------------------------------------------------------*/
  int NumParentPairs = REPMap.size();
  int NumChildPairs  = 0, NumChildPairs2 = 0;
  for(int nePair=0; nePair<NERadix*NERadix; nePair++)
   if (IsReduced[nePair]) NumChildPairs++;
  for(ReducedEdgePairMap::iterator it=REPMap.begin(); it!=REPMap.end(); it++)
   NumChildPairs2+=it->second.size();
 
  int NEPairs = (nsa==nsb ? NEA*(NEA+1)/2 : NEA*NEB);
  Log(" Of %u total edge-edge pairs on surfaces (%i,%i) (%s,%s):",NEPairs,nsa,nsb,G->Surfaces[nsa]->Label,G->Surfaces[nsb]->Label);
  Log("    %u are children (savings of %.1f %%)",NumChildPairs,100.0*((double)NumChildPairs)/((double)NEPairs));
  if (NumChildPairs2!=NumChildPairs)
   Log(" ** warning: child pair counts disagree (%i,%i)",NumChildPairs,NumChildPairs2);
  if (G->LogLevel>SCUFF_VERBOSE2)
   { Log("    %u are parents (%.1f %%)",NumParentPairs, 100.0*((double)NumParentPairs) / ((double)NEPairs));
//     Log("    %u are unicorns (should be %u)",REPMap->size(), NEPairs - NumParentPairs - NumChildPairs);
   }
}
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

   for( tr1::unordered_map<int,int>::iterator it=Consolidators.begin(); it!=Consolidators.end(); it++)
    { 
      int BranchPair = it->first, RootPair = it->second;
      int neaRoot, nebRoot; bool Flipped=ResolveEdgePair(RootPair, &neaRoot, &nebRoot);
      int neaBranch, nebBranch; ResolveEdgePair(BranchPair, &neaBranch, &nebBranch);
      printf("Consolidate {%i,%i}  <--> %c{%i,%i}\n",neaRoot,nebRoot,Flipped ? '-' : '+',neaBranch,nebBranch);

      if (RootPair<0) RootPair*=-1;
      ReducedEdgePairMap::iterator itRoot = REPMap.find(RootPair), itBranch = REPMap.find(BranchPair);
      if (itRoot==REPMap.end() || itBranch==REPMap.end())
       { Warn("%s:%i: internal error (%i,%i)",RootPair,BranchPair);
         continue;
       }
      int Sign = (Flipped ? -1 : 1);
      itRoot->second.push_back( Sign * BranchPair );
      for(size_t n=0; n<itBranch->second.size(); n++) 
       itRoot->second.push_back( Sign * itBranch->second[n] );
      REPMap.erase(BranchPair);
    }

Export("/tmp/After.EEPTable");

  /*----------------------------------------------------------------*/
  /* report some statistics on equivalent edge pairs to log file    */
  /*----------------------------------------------------------------*/
  int NumParentPairs = REPMap.size();
  int NumChildPairs  = 0, NumChildPairs2 = 0;
  for(int nePair=0; nePair<NERadix*NERadix; nePair++)
   if (IsReduced[nePair]) NumChildPairs++;
  for(ReducedEdgePairMap::iterator it=REPMap.begin(); it!=REPMap.end(); it++)
   NumChildPairs2+=it->second.size();
 
  int NEPairs = (nsa==nsb ? NEA*(NEA+1)/2 : NEA*NEB);
  Log(" Of %u total edge-edge pairs on surfaces (%i,%i) (%s,%s):",NEPairs,nsa,nsb,G->Surfaces[nsa]->Label,G->Surfaces[nsb]->Label);
  Log("    %u are children (savings of %.1f %%)",NumChildPairs,100.0*((double)NumChildPairs)/((double)NEPairs));
  if (NumChildPairs2!=NumChildPairs)
   Log(" ** warning: child pair counts disagree (%i,%i)",NumChildPairs,NumChildPairs2);
  if (G->LogLevel>SCUFF_VERBOSE2)
   { Log("    %u are parents (%.1f %%)",NumParentPairs, 100.0*((double)NumParentPairs) / ((double)NEPairs));
//     Log("    %u are unicorns (should be %u)",REPMap->size(), NEPairs - NumParentPairs - NumChildPairs);
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

  for(ReducedEdgePairMap::iterator it=REPMap.begin(); it!=REPMap.end(); it++)
   { int neaParent, nebParent;
     ResolveEdgePair(it->first, &neaParent, &nebParent);
     int Negatives=0;
     for(size_t nREP=0; nREP<it->second.size(); nREP++)
      if (it->second[nREP]<0) 
       Negatives++;
     int Positives = it->second.size() - Negatives;
     fprintf(f,"{%i,%i} %i %i ",neaParent,nebParent,Positives,Negatives);
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
int ReadIndexPair(FILE *f, bool AtTopOfLine, int *n1, int *n2, bool *SignFlip)
{
  bool Success;
  char SignChar;
  int nPos, nNeg;
  if( AtTopOfLine )
   Success = ( fscanf(f," {%i,%i} %i %i ",n1,n2,&nPos,&nNeg)==4 );
  else
   Success = ( fscanf(f," %c{%i,%i}",&SignChar,n1,n2)==3 );
  if (!Success)
   return feof(f) ? 0 : 1;

  if (SignChar=='-')
   *SignFlip=true;
  else if (SignChar=='+')
   *SignFlip=false;
  else 
   return 1;

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
     bool SignFlip;
     int Status=ReadIndexPair(f, AtTopOfLine, &nea, &neb, &SignFlip);
     if (Status==1)
      ErrMsg=vstrdup("%s:%i: syntax error",FileName,LineNum);
     if (Status<=1)
      break;
     if (AtTopOfLine)
      neaParent=nea, nebParent=neb;
     else
      AddReducedPair(neaParent, nebParent, nea, neb, SignFlip);
     if (Status==2) break;
     AtTopOfLine=(Status==3);
     if (AtTopOfLine) LineNum++;
   }
  fclose(f);
  return ErrMsg;
}

} // namespace scuff
