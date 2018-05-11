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
 * EquivalentEdgePairs.h 
 */

#ifndef EQUIVALENT_EDGE_PAIRS_H
#define EQUIVALENT_EDGE_PAIRS_H

#include <libscuff.h>

#ifdef HAVE_CXX11
  #include <unordered_map>
#elif defined(HAVE_TR1)
  #include <tr1/unordered_map>
#endif
#include <map>

typedef vector<bool> bVec;

typedef struct IrreducibleEdgePair
 { int ParentPair;
   iVec *ChildPairs;
 } IrreducibleEdgePair;

typedef vector<IrreducibleEdgePair> IrreducibleEdgePairList;

namespace scuff {

/****************************************************************************/
/* an EquivalentEdgePairTable is a table of equivalent (edge,edge) pairs    */
/****************************************************************************/
class EquivalentEdgePairTable
{ 
public:
   EquivalentEdgePairTable(RWGGeometry *G, int nsa, int nsb);
  ~EquivalentEdgePairTable();

   IrreducibleEdgePairList *GetIrreducibleEdgePairList();

   bool HasEquivalentEdgePair(int ChildPair);
   bool HasEquivalentEdgePair(int neaChild, int nebChild);

   iVec *GetEquivalentEdgePairs(int ParentPair);
   iVec *GetEquivalentEdgePairs(int neaParent, int nebParent);

   int CountParentPairs();
   int CountChildPairs();

   void Export(char *FileName=0);

//private:
// private methods 

   int GetEdgePairIndex(int neParent, int neChild);
   void ResolveEdgePairIndex(int nPair, int *neParent, int *neChild);
  
   void AddEquivalentEdgePair(int ParentPair, int ChildPair, bool SignFlip=false);
   void AddEquivalentEdgePair(int neaParent, int nebParent, int neaChild, int nebChild, bool SignFlip=false);
   void AddEquivalentEdgePair(int neaParent, int neaChild, int nebPair);

// private data fields

   RWGGeometry *G;
   int NERadix;
   bool SameSurface;
   double DistanceQuantum;
   IrreducibleEdgePairList *IEPList;
   bVec HasEquivalentPairFlag;

/* CPLMap[ParentPair] = list of edge pairs equivalent to ParentPair */


#if defined(HAVE_CXX11) 
  typedef std::unordered_map<int, iVec> ChildPairListMap;
#elif defined(HAVE_TR1) 
  typedef std::tr1::unordered_map<int, iVec> ChildPairListMap;
#else
  typedef std::map<int, iVec> ChildPairListMap;
#endif 
  ChildPairListMap CPLMap;

};

} // namespace scuff

#endif // #ifndef EQUIVALENT_EDGE_PAIRS_H
