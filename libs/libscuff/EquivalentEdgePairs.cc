/*
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

/******************************************************************************/
/* EquivalentEdgePairs.cc -- automatic detection of "equivalent edge pairs,"  */
/*                        -- i.e. pairs of RWG basis functions with identical */
/*                        -- SIE matrix elements                              */
/* Homer Reid             -- 4/2018                                           */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <set>
#include <map>

#ifdef HAVE_CONFIG_H
  #include "config.h"
#endif

#include "EquivalentEdgePairs.h"
#include "PanelCubature.h"

#if defined(HAVE_TR1)
  #include <tr1/unordered_map>
#elif defined(HAVE_CXX11)
  #include <unordered_map>
#endif

namespace scuff {

#define UNFLIPPED        0
#define FLIPPED          1
#define NUM_ORIENTATIONS 2

#define ISWAP(a,b) { int temp=a; a=b; b=temp; }
#define VSWAP(a,b) { double *temp=a; a=b; b=temp; }

#define NSME_G1   0
#define NSME_G2   1
#define NSME_IKC  2
#define NUM_NSMES 3
typedef struct { float Data[NUM_NSMES]; } EdgePairSignature;

#define EEP_GSIGN 0
#define EEP_CSIGN 1
#define EEP_SIGNS 2
#define EEP_SIGN_PATTERNS (1<<EEP_SIGNS)

static double EEPQuantum=1.0e-10;
static double EEPRelTol=1.0e-6;

/*****************************************************************/
/*****************************************************************/
/* Part 1: general-purpose utility routines                      */
/*****************************************************************/
/*****************************************************************/
float Quantize(float d, float Unit=0.0)
{ if (Unit==0.0) Unit=EEPQuantum;
  return ( fabs(d) < 0.5*Unit ? 0.0 : Unit*round(d/Unit) ); 
}

float sgn(float x) { return (x<0.0 ? -1.0 : 1.0); }

char PMSign(bool SignFlip) { return SignFlip ? '-' : '+'; }

/*--------------------------------------------------------------*/
/*- Hash function and equality test for general data types      */
/*--------------------------------------------------------------*/
long JenkinsHash(const char *key, size_t len); // in FIBBICache.cc
template<typename KeyType>
 struct TypeHash
  { long operator() (const KeyType &K) const
     { return scuff::JenkinsHash( (const char *)&K, sizeof(KeyType)); }
  };

template<typename KeyType>
 struct TypeEq
  { bool operator() (const KeyType &K1, const KeyType &K2) const
     { return 0==memcmp( (const void *)&K1, (const void *)&K2, sizeof(KeyType) ); }
  };

/*--------------------------------------------------------------*/
/*- Edge pairs -------------------------------------------------*/
/*--------------------------------------------------------------*/
class EdgePair
 { 
public:
   EdgePair(int _nea=0, int _neb=0); 
   
   void Resolve(int *pnea, int *pneb) const;
   char *Str(char *s=0) const;
   bool LessThan(const EdgePair EP) const;
   bool Equal(const EdgePair EP) const;
//private:
   int nea, neb;
   SignPattern Signs;
 };

void GetRelativeSigns(const SignPattern *SP1, const SignPattern *SP2, SignPattern *SPRel)
{ for(int n=0; n<NUMKERNELS; n++) SPRel->Flipped[n] = SP1->Flipped[n] ^ SP2->Flipped[n];
}

SignPattern RelativeSignPattern(const SignPattern *SP1, const SignPattern *SP2)
{ SignPattern RelSigns;
  GetRelativeSigns(SP1, SP2, &RelSigns);
  return RelSigns;
}

SignPattern RelativeSignPattern(const EdgePair *EP1, const EdgePair *EP2)
{ return RelativeSignPattern( &(EP1->Signs), &(EP2->Signs) ); }

void DefaultSignPattern(SignPattern *Signs)
{ for(int n=0; n<NUMKERNELS; n++) Signs->Flipped[n]=false; }

EdgePair::EdgePair(int _nea, int _neb): nea(_nea), neb(_neb)
{ DefaultSignPattern(&Signs); }
  
void EdgePair::Resolve(int *pnea, int *pneb) const
{ *pnea=nea; *pneb=neb; }

char *EdgePair::Str(char *s) const
{ static char sBuffer[20]; 
  if (s==0) s=sBuffer;
  snprintf(s,20,"%c%c{%i,%i}",PMSign(Signs.Flipped[0]),PMSign(Signs.Flipped[1]),nea,neb);
  return s;
}

bool EdgePair::LessThan(const EdgePair EP) const
{ return (nea!=EP.nea ? nea<EP.nea : neb<EP.neb); }

bool EdgePair::Equal(const EdgePair EP) const
{ return (nea==EP.nea && neb==EP.neb); }

struct EdgePairCmp
  { bool operator() (const EdgePair &EP1, const EdgePair &EP2) const
     { return EP1.LessThan(EP2); }
  };

struct EdgePairEq
  { bool operator() (const EdgePair &EP1, const EdgePair &EP2) const
     { return EP1.Equal(EP2); }
  };

struct EdgePairHash
  { long operator() (const EdgePair &EP) const
     { int Key[2];
       Key[0]=EP.nea;
       Key[1]=EP.neb;
       return JenkinsHash( (const char *)Key, sizeof(Key) );
     }
  };


bool GetCanonicalPairIndices(int neaParent, int neaChild, int nebParent, int nebChild,
                             EdgePair *ParentPair, EdgePair *ChildPair)
{
  //if (neaParent > neaChild)  { ISWAP(neaParent, neaChild); }
  //if (nebParent > nebChild)  { ISWAP(nebParent, nebChild); }
  bool Modified=false;
  if (neaParent > nebParent) { ISWAP(neaParent, nebParent); Modified=true; }
  if (neaChild > nebChild)   { ISWAP(neaChild, nebChild);   Modified=true; }

  if ( neaParent>neaChild || (neaParent==neaChild && nebParent<nebChild) )
   { Modified=true;
     EdgePair *Temp=ParentPair; ParentPair=ChildPair; ChildPair=Temp;
   }

  ParentPair->nea = neaParent; ParentPair->neb = nebParent;
  ChildPair->nea = neaChild; ChildPair->neb = nebChild;
  return Modified;
}

/*****************************************************************/
/* Part 2: Routines for detecting and collecting similar edges.  */
/*  "Similarity" of edges is a coarse version of "equivalence".  */
/*  Two edges are similar if they have identical "signatures",   */
/*  where the signature of an edge is a geometric statistic that */
/*  is fast to compute and which yields a useful key for         */
/*  organizing tables of similar edges. Similarity is a          */
/*  necessary but not sufficient condition for equivalence.      */
/*****************************************************************/
#define EDGESIGLEN 3
typedef struct { float Data[3]; } EdgeSignature;
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
   ES.Data[0] = (float) VecNorm(QmQ);
   ES.Data[1] = (float) VecNorm(VmV);
   ES.Data[2] = (float) VecNorm(Cross);
   return ES;
}

bool Less(EdgeSignature ES1, EdgeSignature ES2)
{ for(int n=0; n<EDGESIGLEN; n++) 
   if (ES1.Data[n] < ES2.Data[n]) 
    return true;
  return false;
}

// Given two edges {Ea,Eb}, this routine computes an inexpensive hash value H(Ea,Eb)
// that is guaranteed to have the same value for all equivalent pairs, i.e. if
// {Ea,Eb} <--> {EaPrime,EbPrime} then H(Ea,Eb) = H(EaPrime,EbPrime).
long GetEdgePairHash(RWGSurface *Sa, int nea, RWGSurface *Sb, int neb)
{
 // EdgePairData = (greater signature) + (lesser signature) + (centroid-centroid distance)
  float EdgePairData[2*EDGESIGLEN + 1];

  EdgeSignature ES[2];
  ES[0] = GetEdgeSignature(Sa, nea);
  ES[1] = GetEdgeSignature(Sb, neb);
  int iLesser = ( Less(ES[0], ES[1]) ? 0 : 1 ), iGreater = 1-iLesser;
  for(int n=0; n<EDGESIGLEN; n++)
   { EdgePairData[0*EDGESIGLEN + n] = ES[iGreater].Data[n];
     EdgePairData[1*EDGESIGLEN + n] = ES[iLesser].Data[n];
   }
  EdgePairData[2*EDGESIGLEN + 0] 
   = Quantize(VecDistance(Sa->Edges[nea]->Centroid, Sb->Edges[neb]->Centroid), 1.0e-6);
  long HashVal=scuff::JenkinsHash( (const char *)EdgePairData, sizeof(EdgePairData) ); return HashVal<0 ? -HashVal : HashVal;
  return HashVal<0 ? -HashVal : HashVal;
}

/*****************************************************************/
/* Part 3: Routines for characterizing equivalent edge pairs.    */
/*****************************************************************/
void NSMEIntegrand(double x[3], double b[3], double Divb,
                   double xp[3], double bp[3], double Divbp,
                   void *UserData, double Weight, double *Integral)
{ (void) UserData;
  double R[3]; VecSub(x, xp, R); double r=VecNorm(R);
  double bxbp[3]; VecCross(b,bp,bxbp);
  Integral[NSME_G1]  += Weight*r*VecDot(b,bp);
  Integral[NSME_G2]  += Weight*r*Divb*Divbp;
  Integral[NSME_IKC] += Weight*r*VecDot(R,bxbp);
}

void GetNSMEs(RWGSurface *Sa, int nea, RWGSurface *Sb, int neb, int Order, double NSMEs[NUM_NSMES])
{  GetBFBFCubature2(Sa, nea, Sb, neb, NSMEIntegrand, 0, NUM_NSMES, Order, NSMEs); }

EdgePairSignature GetEdgePairSignature(RWGSurface *Sa, int nea, RWGSurface *Sb, int neb,
                                       SignPattern *Signs=0, int Order=4)
{ 
  double NSMEs[NUM_NSMES];
  GetNSMEs(Sa, nea, Sb, neb, Order, NSMEs);
  EdgePairSignature EPSig;
  double Scale = fmax( fmax( fabs(NSMEs[0]), fabs(NSMEs[1]) ), fabs(NSMEs[2]) );
  for(int n=0; n<NUM_NSMES; n++)
   EPSig.Data[n] = Quantize(NSMEs[n], EEPRelTol*Scale);
  if (Signs)
   { Signs->Flipped[GKERNEL]   = (NSMEs[NSME_G1]   < 0.0);
     Signs->Flipped[IKCKERNEL] = (NSMEs[NSME_IKC]  < 0.0);
   }
  return EPSig;
}

char *ToStr(EdgePairSignature EPSig, char *s=0)
{ static char sBuffer[100];
  if (s==0) s=sBuffer;
  snprintf(s,100,"{%g %g %g}",(double)(EPSig.Data[0]),(double)(EPSig.Data[1]),(double)(EPSig.Data[2]));
  return s;
}

char *EdgePairSignatureStr(RWGSurface *Sa, int nea, RWGSurface *Sb, int neb, char *s=0)
{ return ToStr(GetEdgePairSignature(Sa, nea, Sb, neb), s );
}

/************************************************************************/
/* Part 4: Tables for storing information on sets of equivalent edge pairs.*/
/************************************************************************/

/*********************************************************************************/
/* An EdgePair is a pair of edge indices.                                        */
/* An EdgePairSet is a set of EdgePairs, all equivalent to each other.           */
/*********************************************************************************/
typedef set<EdgePair, EdgePairCmp> EdgePairSet;
typedef map<EdgePair, EdgePairSet, EdgePairCmp> ChildPairMap;

struct EPSigHash
  { long operator() (const EdgePairSignature &EPSig) const
     { EdgePairSignature AbsEPSig;
       for(int n=0; n<NUM_NSMES; n++) AbsEPSig.Data[n] = (float)fabs(EPSig.Data[n]);
       return scuff::JenkinsHash( (const char *)&AbsEPSig, sizeof(EdgePairSignature)); 
     }
  };

struct EPSigEqual
{ bool operator() (const EdgePairSignature &EPSig1, const EdgePairSignature &EPSig2) const
   { for(int n=0; n<NUM_NSMES; n++)
      { float f1=fabs(EPSig1.Data[n]), f2=fabs(EPSig2.Data[n]);
        if ( fabs(f1-f2) > EEPQuantum*f1 )
         return false;
      }
     return true;
   }
};

struct EPSigCmp
{ bool operator() (const EdgePairSignature &EPSig1, const EdgePairSignature &EPSig2) const
   { for(int n=0; n<NUM_NSMES; n++)
      { float f1=fabs(EPSig1.Data[n]), f2=fabs(EPSig2.Data[n]);
        if ( fabs(f1-f2) > EEPQuantum*f1 )
         return true;
      }
     return false;
   }
};

#if defined(HAVE_TR1)
  typedef tr1::unordered_map <EdgePair, EdgePair, EdgePairHash, EdgePairEq > ParentPairMap;
  typedef tr1::unordered_map <EdgePairSignature, EdgePair, EPSigHash, EPSigEqual> EPSigMap;
#elif defined(HAVE_CXX11)
  typedef unordered_map <EdgePair, EdgePair, EdgePairHash, EdgePairEq> ParentPairMap;
  typedef unordered_map <EdgePairSignature, EdgePair, EPSigHash, EPSigEqual> EPSigMap;
#else
  typedef map<EdgePair, EdgePair, EdgePairCmp> ParentPairMap;
  typedef map<EdgePairSignature, EdgePair , EPSigCmp> EPSigMap;
#endif

typedef struct EquivalentEdgePairSubTable
 { RWGSurface *Sa, *Sb;
   ChildPairMap   Children;
   ParentPairMap  Parents;
   EPSigMap       EPRepresentatives;
 } EquivalentEdgePairSubTable;

void ExportEEPSubTable(EquivalentEdgePairSubTable *Table, const char *FileName)
{
  FILE *f = (!strcmp(FileName,"stdout") ? stdout : fopen(FileName,"w"));
  if (!f) 
   { Warn("could not open file %s (skipping edge-pair table export)",FileName); 
     return;
   }
  fprintf(f,"%s %i \n",Table->Sa->MeshFileName,Table->Sa->NumEdges);
  fprintf(f,"%s %i \n",Table->Sb->MeshFileName,Table->Sb->NumEdges);
  for(ChildPairMap::iterator it=Table->Children.begin(); it!=Table->Children.end(); it++)
   { EdgePair ParentPair  = it->first;
     EdgePairSet Children = it->second;
     iVec ChildrenBySign(EEP_SIGN_PATTERNS,0);
     for(EdgePairSet::iterator p=Children.begin(); p!=Children.end(); p++)
      { int SignPatternIndex = 2*(p->Signs.Flipped[0] ? 1 : 0) + (p->Signs.Flipped[1] ? 1 : 0);
        ChildrenBySign[SignPatternIndex]++;
      }
     fprintf(f,"#_%s_[%lu:%i,%i,%i,%i]\n",ParentPair.Str()+2,Children.size(),ChildrenBySign[0],ChildrenBySign[1],ChildrenBySign[2],ChildrenBySign[3]);
     fprintf(f,"%s ",ParentPair.Str()+2);
     for(EdgePairSet::iterator p=Children.begin(); p!=Children.end(); p++)
       fprintf(f,"%s ",p->Str());
     fprintf(f,"\n");
   }

  if (f!=stdout) fclose(f); 
}

bool ParentInTable(EquivalentEdgePairSubTable *Table, EdgePair ParentPair, EdgePairSet **Children=0)
{ 
  ChildPairMap::iterator it = Table->Children.find(ParentPair);
  if (it==Table->Children.end()) return false;
  if (Children) *Children=&(it->second);
  return true;
}

bool ChildInTable(EquivalentEdgePairSubTable *Table, EdgePair ChildPair,
                  int *neaParent=0, int *nebParent=0, SignPattern *Signs=0)
{
  ParentPairMap::iterator it=Table->Parents.find(ChildPair);
  if (it==Table->Parents.end()) return false;
  if (neaParent && nebParent) it->second.Resolve(neaParent, nebParent);
  if (Signs) GetRelativeSigns( &(it->first.Signs), &(it->second.Signs), Signs);
  return true;
}

bool PairInTable(EquivalentEdgePairSubTable *Table, EdgePair Pair,
                 EdgePairSet **ChildPairs=0, EdgePair *ParentPair=0)
{
  if (ChildPairs) *ChildPairs=0;

  ChildPairMap::iterator it1 = Table->Children.find(Pair);
  if (it1!=Table->Children.end())
   { if (ChildPairs) *ChildPairs = &(it1->second);
     return true;
   }
  ParentPairMap::iterator it2 = Table->Parents.find(Pair);
  if (it2!=Table->Parents.end())
   { if (ParentPair) *ParentPair = it2->second;
     return true;
   }
  return false;
}

void AddChildPair(EquivalentEdgePairSubTable *Table, EdgePair ParentPair, EdgePair ChildPair, EdgePairSet *ChildPairs=0)
{ 
  if (ChildPairs==0) ChildPairs = &(Table->Children[ParentPair]);
  ChildPairs->insert(ChildPair);
  Table->Parents[ChildPair]=ParentPair;
}

void ReParent(EquivalentEdgePairSubTable *Table, EdgePair ExistingParentPair, EdgePair NewParentPair)
{
  EdgePairSet *OldChildPairs =&(Table->Children[ExistingParentPair]);
  EdgePairSet NewChildPairs = EdgePairSet( *OldChildPairs ); 
  for(EdgePairSet::iterator it=NewChildPairs.begin(); it!=NewChildPairs.end(); it++)
   Table->Parents[*it]=NewParentPair;
  Table->Children.erase(ExistingParentPair);
  Table->Children[NewParentPair] = NewChildPairs;
  AddChildPair(Table, NewParentPair, ExistingParentPair);

  EdgePairSignature EPSig = GetEdgePairSignature(Table->Sa, NewParentPair.nea, Table->Sb, NewParentPair.neb);
  Table->EPRepresentatives[EPSig]=NewParentPair;
}

bool AddEquivalentPair(EquivalentEdgePairSubTable *Table, EdgePair ParentPair, EdgePair ChildPair)
{
  if ( ParentPair.Equal(ChildPair) ) return false;

  EdgePairSet *ChildrenOfParentInTable;
  EdgePair ParentOfParentInTable;
  bool ParentInTable=PairInTable(Table, ParentPair, &ChildrenOfParentInTable, &ParentOfParentInTable);
  bool ParentInTableAsParent = ParentInTable && (ChildrenOfParentInTable != 0);
  bool ParentInTableAsChild  = ParentInTable && !ParentInTableAsParent;

  EdgePairSet *ChildrenOfChildInTable;
  EdgePair ParentOfChildInTable;
  bool ChildInTable=PairInTable(Table, ChildPair, &ChildrenOfChildInTable, &ParentOfChildInTable);
  bool ChildInTableAsParent = ChildInTable && (ChildrenOfChildInTable != 0);

  if (ParentInTable && ChildInTable)
   return false;
  else if (!ParentInTable && !ChildInTable)
   AddChildPair(Table, ParentPair, ChildPair);
  else if (ParentInTableAsParent)
   AddChildPair(Table, ParentPair, ChildPair, ChildrenOfParentInTable);
  else if (ParentInTableAsChild)
   AddChildPair(Table, ParentOfParentInTable, ChildPair);
  else // Child in table, either as parent or child
   { EdgePair ParentPairInTable = ChildInTableAsParent ? ChildPair : ParentOfChildInTable;
     if ( ParentPairInTable.LessThan(ParentPair) )
      AddChildPair(Table, ParentPairInTable, ParentPair);
     else
      ReParent(Table, ParentPairInTable, ParentPair);
   }
  return true;
}

bool AddEquivalentPair(EquivalentEdgePairSubTable *Table,
                       int nea1, int nea2, int neb1, int neb2,
                       bool *RearrangedIndices)
{
  EdgePair ParentPair, ChildPair;
  bool Swapped=GetCanonicalPairIndices(nea1, nea2, neb1, neb2, &ParentPair, &ChildPair);
  if (RearrangedIndices) *RearrangedIndices=Swapped;
  return AddEquivalentPair(Table, ParentPair, ChildPair);
}

void MergeEEPSubTables(EquivalentEdgePairSubTable *Tree,
                       EquivalentEdgePairSubTable *Branch)
{
  for(ChildPairMap::iterator it=Branch->Children.begin(); it!=Branch->Children.end(); it++)
   { 
     EdgePair BranchParentPair     = it->first;
     EdgePairSet *BranchChildPairs = &(it->second);

     int neaBranch, nebBranch; BranchParentPair.Resolve(&neaBranch, &nebBranch);
     EdgePairSignature BranchEPSig = GetEdgePairSignature(Tree->Sa, neaBranch, Tree->Sb, nebBranch);
     EdgePair TreeParentPair = BranchParentPair;
     EPSigMap::iterator TreeRepresentative = Tree->EPRepresentatives.find(BranchEPSig);
     if (TreeRepresentative!=Tree->EPRepresentatives.end())
      { 
        TreeParentPair = TreeRepresentative->second;
        if (BranchParentPair.LessThan(TreeParentPair))
         { ReParent(Tree, TreeParentPair, BranchParentPair);
           TreeParentPair=BranchParentPair;
         }
      }
     else
      Tree->EPRepresentatives[BranchEPSig]=BranchParentPair;
     
     for(EdgePairSet::iterator ChildPair=BranchChildPairs->begin(); ChildPair!=BranchChildPairs->end(); ChildPair++)
      AddEquivalentPair(Tree, TreeParentPair, *ChildPair);

     if ( !TreeParentPair.Equal(BranchParentPair) )
      AddEquivalentPair(Tree, TreeParentPair, BranchParentPair);
   }
  delete Branch;
}

/******************************************************************/
/* EquivalentEdgePairTable class constructor: Construct a table   */
/* of equivalent edge pairs for two surfaces in an RWG geometry.  */
/******************************************************************/
EquivalentEdgePairTable::EquivalentEdgePairTable(RWGGeometry *_G, int _nsa, int _nsb, char *EEPTFileName)
 : G(_G), nsa(_nsa), nsb(_nsb)
{
  CheckEnv("SCUFF_EEP_QUANTUM",&EEPQuantum);
  CheckEnv("SCUFF_EEP_RELTOL", &EEPRelTol);

  RWGSurface *Sa = G->Surfaces[nsa], *Sb=G->Surfaces[nsb];
  int NEA=Sa->NumEdges, NEB=Sb->NumEdges;
  //NERadix = (NEA > NEB ? NEA : NEB);
  //double DistanceQuantum=0.1*GetMinPanelRadius(G);

  // try to import table from file
#if 1
(void) EEPTFileName;
#else
  if (EEPTFileName)
   { char *ErrMsg=Import(EEPTFileName);
     if (ErrMsg) 
      Warn(ErrMsg);
     else
      { Log("successfully read EEPTable from file %s",EEPTFileName);
        return;
      }
   }
#endif

  int NumThreads = GetNumThreads();
  vector<EquivalentEdgePairSubTable *> BranchTables(NumThreads);
  for(int nt=0; nt<NumThreads; nt++)
   { BranchTables[nt] = new EquivalentEdgePairSubTable;
     BranchTables[nt]->Sa = G->Surfaces[nsa];
     BranchTables[nt]->Sb = G->Surfaces[nsb];
   }

  Log("Identifying equivalent pairs (%i threads...)",NumThreads);
#ifdef USE_OPENMP
#pragma omp parallel for num_threads(NumThreads)
#endif
  for(int nt=0; nt<NumThreads; nt++)
   for(int nea=0; nea<Sa->NumEdges; nea++)
    for(int neb=((nsa==nsb) ? nea: 0); neb<Sb->NumEdges; neb++)
     {
       int nTask = GetEdgePairHash(Sa, nea, Sb, neb) % NumThreads;
       if (nt != nTask) continue;

       EdgePair Pair(nea, neb);
       if (ChildInTable(BranchTables[nt],Pair)) continue;

       EdgePairSignature EPSig           = GetEdgePairSignature(Sa, nea, Sb, neb, &(Pair.Signs));
       EPSigMap::iterator Representative = BranchTables[nt]->EPRepresentatives.find(EPSig);
       if (Representative==BranchTables[nt]->EPRepresentatives.end())
        { BranchTables[nt]->EPRepresentatives[EPSig]=Pair;
          BranchTables[nt]->Children[Pair] = EdgePairSet(); // empty
        }
       else
        { EdgePair TablePair = Representative->second;
          AddEquivalentPair(BranchTables[nt], TablePair, Pair);
        }
     }

  for(int nt=0; nt<NumThreads; nt++)
   { int NumParentPairs      = BranchTables[nt]->Children.size();
     int NumChildPairs       = BranchTables[nt]->Parents.size();
     Log("Thread %i: %i/%i parents/children",nt,NumParentPairs,NumChildPairs);
     if (nt==0) continue;
     MergeEEPSubTables(BranchTables[0], BranchTables[nt]);
   }
   
 // FIXME
 // Log("Before pruning: %lu parents",BranchTables[0]->Children.size());
 // for(ChildPairMap::iterator it=BranchTables[0]->Children.begin(); it!=BranchTables[0]->Children.end(); it++)
 // if ( it->second.size() == 0)
 //   BranchTables[0]->Children.erase(it);
 // Log("After pruning: %lu parents",BranchTables[0]->Children.size());

  int NEPairs = (nsa==nsb ? NEA*(NEA+1)/2 : NEA*NEB);
  int NumParentPairs     = BranchTables[0]->Children.size();
  int NumChildPairs      = BranchTables[0]->Parents.size();
  Log(" Of %i total edge-edge pairs on surfaces (%i,%i) (%s,%s):",NEPairs,nsa,nsb,Sa->Label,Sb->Label);
  Log("    %i are children (savings of %.0f %%)",NumChildPairs,100.0*((double)NumChildPairs)/((double)NEPairs));
  Log("    %i are parents (%.1f %%)",NumParentPairs, 100.0*((double)NumParentPairs) / ((double)NEPairs));
  MasterTable = (void *)BranchTables[0];
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool EquivalentEdgePairTable::HasParent(int neaChild, int nebChild, int *neaParent, int *nebParent, SignPattern *Signs)
{ 
  EquivalentEdgePairSubTable *Table=(EquivalentEdgePairSubTable *)MasterTable;
  bool Status=ChildInTable(Table, EdgePair(neaChild,nebChild), neaParent, nebParent, Signs);
  return Status;
}

ParentPairList EquivalentEdgePairTable::GetParents()
{ 
  EquivalentEdgePairSubTable *Table=(EquivalentEdgePairSubTable *)MasterTable;
  ParentPairList ParentPairs;
  for(ChildPairMap::iterator it=Table->Children.begin(); it!=Table->Children.end(); it++)
   ParentPairs.push_back( ParentPairData(it->first.nea, it->first.neb) );
  return ParentPairs;
}

ChildPairList EquivalentEdgePairTable::GetChildren(int neaParent, int nebParent, bool IncludeParent)
{ 
  EquivalentEdgePairSubTable *Table=(EquivalentEdgePairSubTable *)MasterTable;
  ChildPairList ChildPairs;
  if (IncludeParent) ChildPairs.push_back(ChildPairData(neaParent, nebParent));
  ChildPairMap::iterator it=Table->Children.find( EdgePair(neaParent, nebParent) );
  if (it!=Table->Children.end())
   { const EdgePair *ParentEdgePair = &(it->first);
     EdgePairSet *ChildEdgePairs    = &(it->second);
     for(EdgePairSet::iterator ChildEdgePair=ChildEdgePairs->begin(); ChildEdgePair!=ChildEdgePairs->end(); ChildEdgePair++)
      ChildPairs.push_back( ChildPairData(ChildEdgePair->nea, ChildEdgePair->neb, RelativeSignPattern( &(*ChildEdgePair), ParentEdgePair) ) );
   }
  return ChildPairs;
}

int EquivalentEdgePairTable::NumParents()
{ EquivalentEdgePairSubTable *Table = (EquivalentEdgePairSubTable *)MasterTable;
  return Table->Children.size();
}

int EquivalentEdgePairTable::NumChildren()
{ EquivalentEdgePairSubTable *Table = (EquivalentEdgePairSubTable *)MasterTable;
  return Table->Parents.size();
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
char *GetStandardEEPTFilePath(RWGSurface *Sa, RWGSurface *Sb)
{
  static char Path[1000];
  char PWD[2]=".";

  char *Dir = Sa->MeshFileDir;
  if (Dir==0) Dir=PWD;
  if ( !strcmp(Sa->MeshFileName, Sb->MeshFileName) )
   snprintf(Path,1000,"%s/%s.EEPTable",Dir,GetFileBase(Sa->MeshFileName));
  else
   snprintf(Path,100,"%s/%s_%s.EEPTable",Dir,GetFileBase(Sa->MeshFileName),GetFileBase(Sb->MeshFileName));
  return Path;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void EquivalentEdgePairTable::Export(const char *FileName)
{ if (FileName==0) FileName=GetStandardEEPTFilePath(G->Surfaces[nsa], G->Surfaces[nsb]);
  ExportEEPSubTable( (EquivalentEdgePairSubTable *)MasterTable, FileName); 
  Log("Exported EEPTable(%i,%i) to file %s.\n",nsa,nsb,FileName);
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
// read a string of the form {n1,n2} or %c%c{n1,n2} from f.
// return codes:
// 0: fail because end of file
// 1: fail for other reason
// 2: success, at end of file
// 3: success, at end of line
// 4: success, neither of the above
#if 0
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

/***************************************************************/
/***************************************************************/
/***************************************************************/
bool EquivalentEdgePairTable::Import(const char *FileName)
{ 
  /*--------------------------------------------------------------*/
  /* try to open file and check header ---------------------------*/
  /*--------------------------------------------------------------*/
  if (FileName==0) FileName=GetStandardEEPTFilePath(Sa, Sb);
  Log("Trying to import EEPTable from %s...",FileName);
  FILE *f=fopen(FileName,"r");
  if (f==0) return false;
  char Line[1000]; 
  for(int ns=0; ns<2; ns++)
   { char LineShouldBe[1000];
     RWGSurface *S = (ns==0) ? Sa : Sb;
     snprintf(LineShouldBe,1000,"%s %i\n",S->MeshFileName,S->NumEdges);
     if ( !fgets(Line,1000,f) || strcmp(Line, LineShouldBe) )
      { Log("failed: line %i is %s, should be %s",ns+1,Line,LineShouldBe);
        fclose(f); 
        return false;
      }
   }
 
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int LineNum=3;
  bool AtTopOfLine=true;
  int neaParent, nebParent;
  char *ErrMsg=0;
  while( !feof(f) )
   { int nea, neb;
     SignPattern Signs[2];
     int Status=ReadIndexPair(f, AtTopOfLine, &nea, &neb, &SignPattern);
     if (Status==1)
      ErrMsg=vstrdup("%s:%i: syntax error",FileName,LineNum);
     if (Status<=1)
      break;
     if (AtTopOfLine)
      neaParent=nea, nebParent=neb;
     else
      AddEquivalentPair(EPSMap, neaParent, nebParent, nea, neb, SignFlip);
      //AddEquivalentPair(FileEPS, neaParent, nebParent, nea, neb, SignFlip);
     if (Status==2) break;
      Warn("%s:%i: internal error A",__FILE__,__LINE__);
     AtTopOfLine=(Status==3);
     if (AtTopOfLine) LineNum++;
   }
  fclose(f);
  //MergeEPSMaps(&EPSMap, &FileEPS);
  return ErrMsg;
}
#endif

} // namespace scuff
