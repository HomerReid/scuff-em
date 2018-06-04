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
 * EquivalentEdges.h 
 */

#ifndef EQUIVALENT_EDGE_PAIRS_H
#define EQUIVALENT_EDGE_PAIRS_H

#include <libscuff.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#if defined(HAVE_TR1)
  #include <tr1/unordered_map>
#elif defined(HAVE_CXX11)
  #include <unordered_map>
#endif
#include <map>

namespace scuff {

/*********************************************************************************/
/* An EdgePair is a pair of edge indices.                                        */
/* An EdgePairList is a list of EdgePairs.                                       */
/*********************************************************************************/
typedef int EdgePair;
typedef vector<EdgePair> EdgePairList;

typedef map<int, bool> EquivalentPairSet;
typedef map<int, EquivalentPairSet* > EquivalentPairSetMap;

/****************************************************************************/
/* an EquivalentEdgePairTable is a table of matrix-element redundancies for */
/* a given pair of RWG surfaces.                                            */
/****************************************************************************/
struct EdgePairData; // defined in EquivalentEdges.cc

class EquivalentEdgePairTable
{ 
public:
   EquivalentEdgePairTable(RWGGeometry *G, int nsa, int nsb, char *EEPTFile=0);
   bool HasParent(int nea, int neb);
   EquivalentPairSet GetChildren(int neaParent, int nebParent);

   void Export(char *EEPTFile=0);
   char *Import(char *EEPTFile);

//private:
// private methods 

   // conversion between pairs of edge indices and indices of edge pairs
   EdgePair GetEdgePair(int neParent, int neChild);
   bool ResolveEdgePair(EdgePair nPair, int *neParent, int *neChild);
  
   // helper methods for constructing the table
   int TestEdgePairPair(int ParentPair, int ChildPair);
   bool EvaluatePairPair(EquivalentPairSetMap &EPSetMap, struct EdgePairData *aPair, struct EdgePairData *bPair);
   void AddEquivalentPair(EquivalentPairSetMap &EPSetMap,
                          int neaParent, int neaChild, int nebParent, int nebChild,
                          bool Flipped=false);
   void MergeEPSMaps(EquivalentPairSetMap *Tree, EquivalentPairSetMap *Branch);

// private data fields

   RWGGeometry *G;
   int nsa, nsb;
   int NERadix;

   bVec IsReduced;

   /* EPSetMap[ParentPair] = set of edge pairs equivalent to ParentPair */
   EquivalentPairSetMap EPSMap;
};


} // namespace scuff

#endif // #ifndef EQUIVALENT_EDGE_PAIRS_H
