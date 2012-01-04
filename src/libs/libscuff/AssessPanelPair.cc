/*
 * AssessPanelPair.cc 
 * 
 * homer reid -- 11/2005 -- 10/2011
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
  double *OVb[3]={Vb[0], Vb[1], Vb[2]}; // 'original vertices, b

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

  double rMax=fmax(Pa->Radius, Pb->Radius);
  *rRel=VecDistance(Pa->Centroid, Pb->Centroid) / rMax;

  if ( *rRel > 2.0 ) // there can be no common vertices in this case 
   return 0;

  return AssessPanelPair(Va, Vb, rMax);

}

/***************************************************************/
/* alternate entry points in which the caller wants only a     */
/* subset of the full information returned by AssessPanelPair  */
/***************************************************************/
int AssessPanelPair(RWGObject *Oa, int npa, 
                    RWGObject *Ob, int npb, 
                    double *rRel)
{ 
  double *Va[3], *Vb[3];
  return AssessPanelPair(Oa, npa, Ob, npb, rRel, Va, Vb);
} 

int NumCommonVertices(RWGObject *Oa, int npa, 
                      RWGObject *Ob, int npb)
{ 
  double rRel;
  double *Va[3], *Vb[3];
  return AssessPanelPair(Oa, npa, Ob, npb, &rRel, Va, Vb);
} 
