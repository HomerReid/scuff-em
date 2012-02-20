/*
 * AssessPanelPair.cc 
 * 
 * homer reid -- 11/2005 -- 1/2012
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>

#include "libscuff.h"

// the 'common vertex threshold:' two vertices are considered to be 
// the same if their distance is less than CVTHRESHOLD*the panel radius 
#define CVTHRESHOLD 1.0e-6

/***************************************************************/
/* look for vertices that are shared between panels.           */
/* NOTE: in an earlier incarnation of this code, i looked for  */
/*       common vertices by simple integer comparisons         */
/*       (comparing indices within a table of vertices), but   */
/*       i specifically DON'T want to do that here for several */
/*       reasons. ultimately it would be nice to avoid doing   */
/*       9 separate comparisons here, although in practice it  */
/*       won't matter much since most cases will be caught by  */
/*       the rRel>2 check above.                               */
/***************************************************************/
int AssessPanelPair(double **Va, double **Vb, double Radius)
{
  int ncv=0;
  int CVIa[3], CVIb[3];
  int ia, ib;
  double ThresholdDistance2=CVTHRESHOLD*CVTHRESHOLD*Radius*Radius;
  for(ia=0; ia<3; ia++)
   for(ib=0; ib<3; ib++)
    { if ( VecDistance2(Va[ia],Vb[ib]) < ThresholdDistance2 )
       { CVIa[ncv]=ia;    
         CVIb[ncv]=ib;
         ncv++;
       };
    };
  
  /***************************************************************/
  /* if there were any common vertices, reorganize the Va and Vb */
  /* arrays so that common vertices appear first                 */
  /***************************************************************/
  if (ncv==0)
   return 0; // vertices are already in acceptable order

  double *OVa[3]={Va[0], Va[1], Va[2]}; // 'original vertices, a' 
  double *OVb[3]={Vb[0], Vb[1], Vb[2]}; // 'original vertices, b'

  if (ncv==1)
   { Va[0] = OVa[ CVIa[0]       ];
     Va[1] = OVa[ (CVIa[0]+1)%3 ];
     Va[2] = OVa[ (CVIa[0]+2)%3 ];
     Vb[0] = OVb[ CVIb[0]       ];
     Vb[1] = OVb[ (CVIb[0]+1)%3 ];
     Vb[2] = OVb[ (CVIb[0]+2)%3 ];
   }
  else if (ncv==2)
   { Va[0] = OVa[  CVIa[0] ];
     Va[1] = OVa[  CVIa[1] ];
     Va[2] = OVa[  3-CVIa[0]-CVIa[1] ];
     Vb[0] = OVb[  CVIb[0] ];
     Vb[1] = OVb[  CVIb[1] ];
     Vb[2] = OVb[  3-CVIb[0]-CVIb[1] ];
   }
  else if (ncv==3)
   { Va[0] = OVa[ CVIa[0] ];
     Va[1] = OVa[ CVIa[1] ];
     Va[2] = OVa[ CVIa[2] ];
     Vb[0] = OVb[ CVIb[0] ];
     Vb[1] = OVb[ CVIb[1] ];
     Vb[2] = OVb[ CVIb[2] ];
   };

  return ncv;
 
}

/***************************************************************/
/* alternate entry point for the previous routine for the case */
/* in which the maximum panel radius is not known a priori     */
/***************************************************************/
int AssessPanelPair(double **Va, double **Vb)
{ 
  int Mu;
  double CentroidA[3], CentroidB[3];

  for(Mu=0; Mu<3; Mu++)
   { CentroidA[Mu]=(Va[0][Mu] + Va[1][Mu] + Va[2][Mu])/3.0;
     CentroidB[Mu]=(Vb[0][Mu] + Vb[1][Mu] + Vb[2][Mu])/3.0;
   };

  double rMax;
  rMax=VecDistance(Va[0],CentroidA);
  rMax=fmax(rMax, VecDistance(Va[1],CentroidA));
  rMax=fmax(rMax, VecDistance(Va[2],CentroidA));
  rMax=fmax(rMax, VecDistance(Vb[0],CentroidB));
  rMax=fmax(rMax, VecDistance(Vb[1],CentroidB));
  rMax=fmax(rMax, VecDistance(Vb[2],CentroidB));

  return AssessPanelPair(Va, Vb, rMax);

}

/***************************************************************/
/* this routine gathers some information on the pair of panels */
/* (Oa, npa) -- (Ob,npb).                                      */
/*                                                             */
/* on return from this routine,                                */
/*  a) the return value is the # of common vertices (0,1,2,3)  */
/*  b) *rRel is set to the 'relative distance'                 */
/*  c) Va[0..2] and Vb[0..2] are pointers to the vertices of   */
/*     the two panels                                          */
/*  d) if there are any common vertices, then the ordering of  */
/*     the Va and Vb arrays is such that any common vertices   */
/*     come first; for example, if there are 2 common vertices */
/*     then Va[0] = Vb[0] and Va[1] = Vb[1].                   */
/*     (the pointers aren't necessarily equal, but the         */
/*      three-vectors to which they point are equal.)          */
/***************************************************************/
int AssessPanelPair(RWGObject *Oa, int npa, RWGObject *Ob, int npb,
                    double *rRel, double **Va, double **Vb)
{
  RWGPanel *Pa=Oa->Panels[npa];
  RWGPanel *Pb=Ob->Panels[npb];

  Va[0] = Oa->Vertices + 3*Pa->VI[0];
  Va[1] = Oa->Vertices + 3*Pa->VI[1];
  Va[2] = Oa->Vertices + 3*Pa->VI[2];

  Vb[0] = Ob->Vertices + 3*Pb->VI[0];
  Vb[1] = Ob->Vertices + 3*Pb->VI[1];
  Vb[2] = Ob->Vertices + 3*Pb->VI[2];

  double DC, rRel2, rMax=fmax(Pa->Radius, Pb->Radius);

  DC=(Pa->Centroid[0]-Pb->Centroid[0]); rRel2=DC*DC;
  DC=(Pa->Centroid[1]-Pb->Centroid[1]); rRel2+=DC*DC;
  DC=(Pa->Centroid[2]-Pb->Centroid[2]); rRel2+=DC*DC;
  *rRel=sqrt(rRel2) / rMax;
  if ( *rRel > 2.0 ) // there can be no common vertices in this case 
   return 0;

  return AssessPanelPair(Va, Vb, rMax);

}

/***************************************************************/
/* alternate entry point in which the caller wants only a      */
/* subset of the full information returned by AssessPanelPair  */
/***************************************************************/
int AssessPanelPair(RWGObject *Oa, int npa, 
                    RWGObject *Ob, int npb, 
                    double *rRel)
{ 
  double *Va[3], *Vb[3];
  return AssessPanelPair(Oa, npa, Ob, npb, rRel, Va, Vb);
} 

/***************************************************************/
/* this routine simply counts the number of common vertices    */
/* without reordering the vertices                             */
/***************************************************************/
int NumCommonVertices(RWGObject *Oa, int npa, 
                      RWGObject *Ob, int npb)
{ 
  double rRel;
  double *Va[3], *Vb[3];
  return AssessPanelPair(Oa, npa, Ob, npb, &rRel, Va, Vb);
} 

/***************************************************************/
/* 'vertex less than.' returns 1 if V1<=V2, 0 otherwise.       */
/* vertices are sorted using a fairly obvious sorting scheme.  */
/***************************************************************/
static int VLT(double *V1, double *V2)
{
  double DV;

  DV=V1[0]-V2[0];
  if ( fabs(DV) > 1.0e-6*fabs(V1[0]) )
   return DV<0.0 ? 1 : 0;

  DV=V1[1]-V2[1];
  if ( fabs(DV) > 1.0e-6*fabs(V1[1]) )
   return DV<0.0 ? 1 : 0;
 
  DV=V1[2]-V2[2];
  if ( fabs(DV) > 1.0e-6*fabs(V1[2]) )
   return DV<0.0 ? 1 : 0;

  return 1;
}

/***************************************************************/
/* given two sets of panel vertices, this routine puts them    */
/* in a canonical ordering, based on the vertex-ordering       */
/* algorithm of the above routine. the canonical ordering is   */
/* used to construct a search key into a FIPPI table.          */
/*                                                             */
/* inputs: Va[0..2], Vb[0..2]: panel vertices                  */
/*         ncv:                number of common vertices       */
/*                                                             */
/* outputs: OVa[0..2], OVb[0..2]: canonically ordered          */
/*                                panel vertices               */
/*                                                             */
/* return value = 1 if the panels were swapped to achieve      */
/*                  the canonical ordering                     */
/*                                                             */
/*                0 if the panels were not swapped.            */
/*                                                             */
/* important: this routine is intended to be used as an        */
/* optional follow-up to a call to AssessPanelPair. in         */
/* particular, if there are any common vertices (ncv>0), the   */
/* Va and Vb arrays are assumed to be in the order returned    */
/* by AssessPanelPair. (note that AssessPanelPair ensures that */
/* the vertex sets are in acceptable order for passage to      */
/* TaylorMaster, but does not take the further step of         */
/* putting the vertices into the CANONICAL order.)             */
/*                                                             */
/* my 'canonical ordering' for panel vertices obeys the        */
/* following properties.                                       */
/*                                                             */
/*  (a) if the panels have no common vertices, then the        */
/*      vertices of each panel are sorted in ascending order   */
/*      and the smallest vertex of OVa is less than OVb, i.e.  */
/*      we have                                                */ 
/*        OVa[0] < OVa[1] < OVa[2]                             */ 
/*        OVb[0] < OVb[1] < OVb[2]                             */
/*      and                                                    */
/*        OVa[0] < OVb[0].                                     */
/*                                                             */
/*  (b) if the panels have any common vertices, then those     */
/*      common vertices appear first in both lists; within     */
/*      each panel, the subset of vertices that are common     */
/*      with the other panel are sorted in ascending order,    */
/*      the subset of vertices that are not common are         */
/*      separately sorted in ascending order, and the first    */
/*      noncommon vertex in OVa is less than the first         */
/*      noncommon vertex in OVb.                               */
/*                                                             */
/*      more specifically, what this means is:                 */
/*                                                             */
/*       (b1) if there are three common vertices, we have      */
/*                                                             */
/*             OVa[0] < OVa[1] < OVa[2]                        */
/*                                                             */
/*            and                                              */
/*                                                             */
/*             OVb[i] = OVa[i]  for i=0,l,2.                   */
/*                                                             */
/*       (b2) if there are two common vertices, we have        */
/*                                                             */
/*             OVa[0] < OVa[1],                                */
/*                                                             */
/*             OVb[i] = OVa[i]  for i=0,1                      */
/*                                                             */
/*             and OVa[2] < OVb[2].                            */
/*                                                             */
/*       (b3) if there is one common vertex, we have           */
/*                                                             */
/*             OVa[0] = OVb[0],                                */
/*             OVa[1] < OVa[2],                                */
/*             OVb[1] < OVb[2],                                */
/*            and                                              */
/*             OVa[1] < OVb[1].                                */
/***************************************************************/
int CanonicallyOrderVertices(double **Va, double **Vb, int ncv,
                             double **OVa, double **OVb)
{
  int iMina, iMeda, iMaxa;
  int iMinb, iMedb, iMaxb;
  double *TV;

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if ( ncv == 3 )
   { 
     if ( VLT(Va[0], Va[1]) )
      { iMina=0; iMaxa=1; }
     else
      { iMina=1; iMaxa=0; }

     if ( VLT(Va[2], Va[iMina]) )
      { iMeda=iMina; iMina=2; }
     else if ( VLT(Va[2], Va[iMaxa]) )
      { iMeda=2; }
     else
      { iMeda=iMaxa; iMaxa=2; }

     OVa[0]=Va[iMina]; OVa[1]=Va[iMeda]; OVa[2]=Va[iMaxa];
     OVb[0]=Vb[iMina]; OVb[1]=Vb[iMeda]; OVb[2]=Vb[iMaxa];
     return 0;
   }
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  else if ( ncv == 2 )
   {
     int a01Min;

     a01Min = VLT(Va[0], Va[1]) ? 0 : 1;

     if ( VLT( Va[2], Vb[2] ) )
      { 
        OVa[0]=Va[a01Min]; OVa[1]=Va[1-a01Min]; OVa[2]=Va[2];
        OVb[0]=Vb[a01Min]; OVb[1]=Vb[1-a01Min]; OVb[2]=Vb[2];
        return 0;
      }
     else
      { 
        OVa[0]=Vb[a01Min]; OVa[1]=Vb[1-a01Min]; OVa[2]=Vb[2];
        OVb[0]=Va[a01Min]; OVb[1]=Va[1-a01Min]; OVb[2]=Va[2];
        return 1;
      };
   }
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  else if ( ncv == 1 )
   {
     int a12Min, b12Min; 

     a12Min = VLT(Va[1], Va[2]) ? 1 : 2;
     b12Min = VLT(Vb[1], Vb[2]) ? 1 : 2;

     if ( VLT( Va[a12Min], Vb[b12Min] ) )
      { 
        OVa[0]=Va[0]; OVa[1]=Va[a12Min]; OVa[2]=Va[3-a12Min];
        OVb[0]=Vb[0]; OVb[1]=Vb[b12Min]; OVb[2]=Vb[3-b12Min];
        return 0;
      }
     else
      { OVa[0]=Vb[0]; OVa[1]=Vb[b12Min]; OVa[2]=Vb[3-b12Min];
        OVb[0]=Va[0]; OVb[1]=Va[a12Min]; OVb[2]=Va[3-a12Min];
        return 1;
      };
   }
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  else // ( ncv == 0 )
   { 
     if ( VLT(Va[0], Va[1]) )
      { iMina=0; iMaxa=1; }
     else
      { iMina=1; iMaxa=0; }

     if ( VLT(Va[2], Va[iMina]) )
      { iMeda=iMina; iMina=2; }
     else if ( VLT(Va[2], Va[iMaxa]) )
      { iMeda=2; }
     else
      { iMeda=iMaxa; iMaxa=2; }

     if ( VLT(Vb[0], Vb[1]) )
      { iMinb=0; iMaxb=1; }
     else
      { iMinb=1; iMaxb=0; }

     if ( VLT(Vb[2], Vb[iMinb]) )
      { iMedb=iMinb; iMinb=2; }
     else if ( VLT(Vb[2], Vb[iMaxb]) )
      { iMedb=2; }
     else
      { iMedb=iMaxb; iMaxb=2; }

     if ( VLT(Va[iMina], Va[iMinb]) )
      { OVa[0] = Va[iMina]; OVa[1] = Va[iMeda]; OVa[2] = Va[iMaxa]; 
        OVb[0] = Vb[iMina]; OVb[1] = Vb[iMeda]; OVb[2] = Vb[iMaxa]; 
        return 0;
      }
     else
      { OVb[0] = Va[iMina]; OVb[1] = Va[iMeda]; OVb[2] = Va[iMaxa]; 
        OVa[0] = Vb[iMina]; OVa[1] = Vb[iMeda]; OVa[2] = Vb[iMaxa]; 
        return 1;
      };

   }; // if ncv ... else
  
}
