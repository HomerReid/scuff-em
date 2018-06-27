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
 * EquivalentEdgePairs.h -- definitions for EquivalentEdgePairs module
 */

#ifndef EQUIVALENT_EDGE_PAIRS_H
#define EQUIVALENT_EDGE_PAIRS_H

#ifdef HAVE_CONFIG_H
  #include <config.h>
#endif
#include <vector>

#include "libscuff.h"

namespace scuff {

typedef struct ParentPairData
 { int nea, neb;
   ParentPairData(int _nea, int _neb): nea(_nea), neb(_neb) {}
 } ParentPairData;

typedef vector<ParentPairData> ParentPairList;

#define GKERNEL 0
#define CKERNEL 1
#define NUMKERNELS 2
typedef struct { bool Flipped[NUMKERNELS]; } SignPattern;

typedef struct ChildPairData
 { int nea, neb;
   double GCSign[NUMKERNELS];
   ChildPairData(int _nea, int _neb, SignPattern Signs): nea(_nea), neb(_neb)
    { GCSign[0] = Signs.Flipped[0] ? -1.0 : 1.0, GCSign[1] = Signs.Flipped[1] ? -1.0 : 1.0; }
   ChildPairData(int _nea, int _neb): nea(_nea), neb(_neb) { GCSign[0]=GCSign[1]=1.0; }
 } ChildPairData;
typedef vector<ChildPairData> ChildPairList;

class EquivalentEdgePairTable
 {
public:
    EquivalentEdgePairTable(RWGGeometry *G, int nsa, int nsb, char *EEPTFile=0);
    void Export(const char *EEPTFile=0);

    bool HasParent(int neaChild, int nebChild, int *neaParent=0, int *nebParent=0, SignPattern *Signs=0);
    ParentPairList GetParents();
    ChildPairList GetChildren(int neaParent, int nebParent, bool IncludeParent=false);

// private data fields
// private:
   RWGGeometry *G;
   int nsa, nsb;
   void *MasterTable;
 };

} // namespace scuff 
#endif // #ifndef EQUIVALENT_EDGE_PAIRS_H
