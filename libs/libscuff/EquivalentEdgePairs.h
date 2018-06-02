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
#else
  #include <map>
#endif

namespace scuff {

/*********************************************************************************/
/* An EdgePair is a pair of edge indices.                                        */
/* An EdgePairList is a list of EdgePairs.                                       */
/* A ReducedEdgePairMap is a set of (Key=EdgePair, Value=EdgePairList) pairs,    */
/* where the key is a single pair of RWG edges (the parent pair) and the value   */
/* is a list of edge pairs (the children) where each child pair has matrix       */
/* elements equal to those of the parent pair.                                   */
/*********************************************************************************/
typedef int EdgePair;
typedef vector<EdgePair> EdgePairList;

#if defined(HAVE_CXX11) 
  typedef std::unordered_map<EdgePair, EdgePairList> ReducedEdgePairMap;
#elif defined(HAVE_TR1) 
  typedef std::tr1::unordered_map<EdgePair, EdgePairList> ReducedEdgePairMap;
#else
  typedef std::map<EdgePair, EdgePairList> ReducedEdgePairMap;
#endif

/****************************************************************************/
/* an EquivalentEdgePairTable is a table of matrix-element redundancies for */
/* a given pair of RWG surfaces.                                            */
/****************************************************************************/
struct EdgePairData; // defined in EquivalentEdges.cc

class EquivalentEdgePairTable
{ 
public:
   EquivalentEdgePairTable(RWGGeometry *G, int nsa, int nsb, char *EEPTFile=0);

   void Export(char *EEPTFile=0);
   char *Import(char *EEPTFile);

//private:
// private methods 

   // conversion between pairs of edge indices and indices of edge pairs
   EdgePair GetEdgePair(int neParent, int neChild);
   bool ResolveEdgePair(EdgePair nPair, int *neParent, int *neChild);
  
   // helper methods for constructing the table
   bool AddReducedPair(EdgePair ParentPair, EdgePair ChildPair, bool SignFlip=false);
   bool AddReducedPair(int neaParent, int nebParent, int neaChild, int nebChild, bool SignFlip=false);
   bool AddReducedPair(int neaParent, int neaChild, EdgePair nebPair);
   bool EvaluatePairPair(struct EdgePairData *aPair, struct EdgePairData *bPair);

// private data fields

   RWGGeometry *G;
   int nsa, nsb;
   int NERadix;

   bVec IsReduced;

   /* REPMap[ParentPair] = list of edge pairs equivalent to ParentPair */
   ReducedEdgePairMap REPMap;
};

} // namespace scuff

#endif // #ifndef EQUIVALENT_EDGE_PAIRS_H
