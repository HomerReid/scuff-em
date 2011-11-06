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

  *rRel = VecDistance(Pa->Centroid, Pb->Centroid) / fmax(Pa->Radius, Pb->Radius);
  if ( *rRel > 2.0 ) // there can be no common vertices in this case 
   return 0;

  /***************************************************************/
  /* look for common vertices.                                   */
  /* NOTE: in an earlier incarnation of this code, i looked for  */
  /*       common vertices by simple integer comparisons         */
  /*       (comparing indices within a table of vertices), but   */
  /*       i specifically DON'T want to do that here for several */
  /*       reasons. ultimately it would be nice to avoid doing   */
  /*       9 separate comparisons here, although in practice it  */
  /*       won't matter much since most cases will be caught by  */
  /*       the rRel>2 check above.                               */
  /***************************************************************/
  int CVIa[3], CVIb[3];
  int ncv=0;
  int ia, ib;
  double ThresholdDistance2=CVTHRESHOLD*CVTHRESHOLD*rMax*rMax;
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
   { // vertices are already in acceptable order
   } 
  else if (ncv==1)
   { Va[0] = Oa->Vertices + 3*Pa->VI[  CVIa[0] ];
     Va[1] = Oa->Vertices + 3*Pa->VI[ (CVIa[0]+1)%3 ];
     Va[2] = Oa->Vertices + 3*Pa->VI[ (CVIa[0]+2)%3 ];
     Vb[0] = Ob->Vertices + 3*Pb->VI[  CVIb[0] ];
     Vb[1] = Ob->Vertices + 3*Pb->VI[ (CVIb[0]+1)%3 ];
     Vb[2] = Ob->Vertices + 3*Pb->VI[ (CVIb[0]+2)%3 ];
   }
  else if (ncv==2)
   { Va[0] = Oa->Vertices + 3*Pa->VI[  CVIa[0] ];
     Va[1] = Oa->Vertices + 3*Pa->VI[  CVIa[1] ];
     Va[2] = Oa->Vertices + 3*Pa->VI[  3-CVIa[0]-CVIa[1] ];
     Vb[0] = Ob->Vertices + 3*Pb->VI[  CVIb[0] ];
     Vb[1] = Ob->Vertices + 3*Pb->VI[  CVIb[1] ];
     Vb[2] = Ob->Vertices + 3*Pb->VI[  3-CVIb[0]-CVIb[1] ];
   }
  else if (ncv==3)
   { Va[0] = Oa->Vertices + 3*Pa->VI[  CVIa[0] ];
     Va[1] = Oa->Vertices + 3*Pa->VI[  CVIa[1] ];
     Va[2] = Oa->Vertices + 3*Pa->VI[  CVIa[2] ];
     Vb[0] = Ob->Vertices + 3*Pb->VI[  CVIb[0] ];
     Vb[1] = Ob->Vertices + 3*Pb->VI[  CVIb[1] ];
     Vb[2] = Ob->Vertices + 3*Pb->VI[  CVIb[2] ];
   };

  return ncv;
 
}
