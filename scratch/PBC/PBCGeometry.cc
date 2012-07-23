/*
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <libhrutil.h>
#include <libMDInterp.h>
#include <libscuff.h>
#include "PBCGeometry.h"

namespace scuff{

/***************************************************************/
/* do something about me please ********************************/
/***************************************************************/
int PBCGeometry::TriangleCubatureOrder=4;
double PBCGeometry::DeltaInterp=0.025;

/***************************************************************/
/* get the maximum and minimum cartesian coordinates obtained  */
/* by points on an RWGObject. (we do this by looking at the    */
/* panel vertices, which suffices because the panels are       */
/* flat and thus convex)                                       */
/***************************************************************/
void GetObjectExtents(RWGObject *O, double RMax[3], double RMin[3])
{ 
  memcpy(RMax, O->Vertices + 0, 3*sizeof(double));
  memcpy(RMin, O->Vertices + 0, 3*sizeof(double));
  for(int nv=1; nv<O->NumVertices; nv++)
   { 
     RMax[0] = fmax(RMax[0], O->Vertices[3*nv + 0]);
     RMin[0] = fmin(RMin[0], O->Vertices[3*nv + 0]);
     RMax[1] = fmax(RMax[1], O->Vertices[3*nv + 1]);
     RMin[1] = fmin(RMin[1], O->Vertices[3*nv + 1]);
     RMax[2] = fmax(RMax[2], O->Vertices[3*nv + 2]);
     RMin[2] = fmin(RMin[2], O->Vertices[3*nv + 2]);
   };
}

/***************************************************************/
/* PBCGeometry class constructor *******************************/
/***************************************************************/
PBCGeometry::PBCGeometry(RWGGeometry *pG, double **pLBV)
{
  /*--------------------------------------------------------------*/
  /*- initialize simple class fields -----------------------------*/
  /*--------------------------------------------------------------*/
  G=pG;
  LBV[0][0]=pLBV[0][0];
  LBV[0][1]=pLBV[0][1];
  LBV[1][0]=pLBV[1][0];
  LBV[1][1]=pLBV[1][1];
  CurrentOmega=-1.0;
  EpsTF = (cdouble *)mallocEC( (G->NumObjects+1)*sizeof(cdouble));
  MuTF  = (cdouble *)mallocEC( (G->NumObjects+1)*sizeof(cdouble));

  Log("Creating PBC geometry: unit cell geometry %s, lattice (%g,%g) x (%g,%g)",
       G->GeoFileName,LBV[0][0],LBV[0][1],LBV[1][0],LBV[1][1]); 

  /*--------------------------------------------------------------*/
  /*- add straddlers to all objects ... note that this modifies  -*/
  /*- G and its constituent objects. it would be better program  -*/
  /*- design to make a copy of G?                                -*/
  /*--------------------------------------------------------------*/
  NumStraddlers=(int *)mallocEC(2*G->NumObjects*sizeof(int));
  G->TotalBFs=G->TotalPanels=0;
  double *MyLBVs[2];
  MyLBVs[0]=LBV[0];
  MyLBVs[1]=LBV[1];
  for(int no=0; no<G->NumObjects; no++)
   {
     RWGObject *O = G->Objects[no];
     AddStraddlers(O, MyLBVs, NumStraddlers + 2*no);
 
     // FIXME
     if ( (NumStraddlers[2*no+0]==0) != (NumStraddlers[2*no+1]==0) )
      ErrExit("object %s: 1D straddling periodicity is not yet supported",O->Label); 

     G->TotalBFs    += O->NumBFs;
     G->TotalPanels += O->NumPanels;
     if ( no+1 < G->NumObjects )
      { G->BFIndexOffset[no+1]=G->BFIndexOffset[no] + O->NumBFs;
        G->PanelIndexOffset[no+1]=G->PanelIndexOffset[no] + O->NumPanels;
      };

     Log(" Detected (%i,%i) straddlers for object %s", NumStraddlers[2*no+0],NumStraddlers[2*no+1],O->Label);

   };

  /*--------------------------------------------------------------*/
  /*- allocate memory for the contributions of the innermost     -*/
  /*- lattice cells to the BEM matrix.                           -*/
  /*- note: P, M, Z stand for 'plus 1, minus 1, zero.'           -*/
  /*- Mab is the BEM interaction matrix between the unit-cell    -*/
  /*- geometry and a copy of itself translated through vector    -*/
  /*- a*LBV[0] + b*LBV[1].                                       -*/
  /*--------------------------------------------------------------*/
  MPP=new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX);
  MPM=new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX);
  MPZ=new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX);
  MZP=new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX);
  MZZ=new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX); // this one could be symmetric ...

  /*--------------------------------------------------------------*/
  /*- allocate interpolators for each object interior             */
  /*--------------------------------------------------------------*/
  GBarAB9_Interior=(Interp3D **)mallocEC(G->NumObjects * sizeof(Interp3D *));
  double RMaxTO[3], RMinTO[3]; // max/min coordinates values for this object
  //double RMax[3], RMin[3];   // max/min coord values for entire geometry (now class fields)
  double DeltaR[3];
  int NPoints[3];
  RWGObject *O;
  RMax[0] = RMax[1] = RMax[2] = -1.0e+9;
  RMin[0] = RMin[1] = RMin[2] = +1.0e+9;
  for(int no=0; no<G->NumObjects; no++)
   { 
     O=G->Objects[no];
     GetObjectExtents(O, RMaxTO, RMinTO);
     for(int i=0; i<3; i++)
      {  DeltaR[i]   = RMaxTO[i] - RMinTO[i];
         NPoints[i]  = 1 + (2.0*DeltaR[i] / PBCGeometry::DeltaInterp );
         if (NPoints[i] < 2)
          NPoints[i]=2;
         RMax[i]     = fmax(RMax[i], RMaxTO[i]);
         RMin[i]     = fmin(RMin[i], RMinTO[i]);
      };

     // FIXME to handle 1D periodicity
     if ( O->IsPEC )
      { 
        Log(" Object %s is PEC; no interpolation table needed",O->Label);
        GBarAB9_Interior[no]=0;
      }
     else if ( NumStraddlers[2*no+0]==0 && NumStraddlers[2*no+1]==0 )
      { 
        Log(" Object %s straddles no unit cell boundaries; no interpolation table needed",O->Label);
        GBarAB9_Interior[no]=0;
      }
     else
      { 
        Log("Creating %ix%ix%i i-table for object %s",NPoints[0],NPoints[1],NPoints[2],O->Label);
        GBarAB9_Interior[no]=new Interp3D( -DeltaR[0], DeltaR[0], NPoints[0],
                                           -DeltaR[1], DeltaR[1], NPoints[1],
                                                  0.0, DeltaR[2], 1 + NPoints[2]/2,
                                           2, 0, 0, 0);
      };

   };

  /*--------------------------------------------------------------*/
  /*- allocate interpolator for exterior medium ------------------*/
  /*--------------------------------------------------------------*/
  for(int i=0; i<3; i++)
   { DeltaR[i]   = RMax[i] - RMin[i];
     NPoints[i]  = 1 + (2.0*DeltaR[i] / PBCGeometry::DeltaInterp );
     if (NPoints[i] < 2)
      NPoints[i]=2;
   };
  GBarAB9_Exterior=new Interp3D( -DeltaR[0], DeltaR[0], NPoints[0],
                                 -DeltaR[1], DeltaR[1], NPoints[1],
                                        0.0, DeltaR[2], 1 + NPoints[2]/2,
                                 2, 0, 0, 0);

}

} // namespace scuff
