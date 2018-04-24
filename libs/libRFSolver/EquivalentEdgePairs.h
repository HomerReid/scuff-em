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

namespace scuff {

/****************************************************************************/
/* an EquivalentEdgePairTable is a table of equivalent (edge,edge) pairs    */
/****************************************************************************/
class EquivalentEdgePairTable
{ 
public:
   EquivalentEdgePairTable(RWGGeometry *G);
  ~EquivalentEdgePairTable();

   bool HasParent(int ChildPair);
   bool HasParent(int nfeaChild, int nfebChild);
   bool HasParent(int nsa, int neaChild, int nsb, int nsbChild);

   iVec GetChildren(int ParentPair);
   iVec GetChildren(int nfeaParent, int nfebParent);
   iVec GetChildren(int nsa, int neaParent, int nsb, int nebParent);

   int CountParentPairs();
   int CountChildPairs();

//private:
  void AddEquivalentEdgePair(int ParentPair, int ChildPair);
  void AddEquivalentEdgePair(int nfeaParent, int nfebParent, int nfeaChild, int nfebChild);

  RWGGeometry *G;

/* ChildPairLists[ParentPair] = list of edge pairs equivalent to ParentPair */
/* ParentPair[ChildPair]      = parent edge pair for ChildPair              */

  iVec *ChildPairListArray;
  int *ParentPairArray;

#if defined(HAVE_CXX11) 
  std::unordered_map<int, iVec> ChildPairListMap;
  std::unordered_map<int, int>  ParentPairMap;
#elif defined(HAVE_TR1) 
  std::tr1::unordered_map<int, iVec> ChildPairListMap;
  std::tr1::unordered_map<int, int>  ParentPairMap;
#else
  std::map<int, iVec> ChildPairListMap;
  std::map<int, int>  ParentPairMap;
#endif 

};

} // namespace scuff

#endif // #ifndef EQUIVALENT_EDGE_PAIRS_H
