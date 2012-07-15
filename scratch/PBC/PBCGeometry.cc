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
int PBCGeometry::TriangleCubatureOrder=7;
double PBCGeometry::DeltaInterp=0.05;

/***************************************************************/
/* get the maximum and minimum cartesian coordinates obtained  */
/* by points on an RWGObject                                   */
/***************************************************************/
void GetXYZMaxMin(RWGObject *O, double XYZMax[3], double XYZMin[3])
{ 
  XYZMax[0] = XYZMax[1] = XYZMax[2] = -1.0e+9;
  XYZMin[0] = XYZMin[1] = XYZMin[2] = +1.0e+9;
  for(int nv=0; nv<O->NumVertices; nv++)
   { 
     XYZMax[0] = fmax(XYZMax[0], O->Vertices[3*nv + 0]);
     XYZMin[0] = fmin(XYZMin[0], O->Vertices[3*nv + 0]);
     XYZMax[1] = fmax(XYZMax[1], O->Vertices[3*nv + 1]);
     XYZMin[1] = fmin(XYZMin[1], O->Vertices[3*nv + 1]);
     XYZMax[2] = fmax(XYZMax[2], O->Vertices[3*nv + 2]);
     XYZMin[2] = fmin(XYZMin[2], O->Vertices[3*nv + 2]);
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
  /*- a*LBV[0] + b*LBV[1].                                      -*/
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
  double XYZMaxTO[3], XYZMinTO[3]; // 'x, y, z max/min, this object'
  //double XYZMax[3], XYZMin[3];  // these are now class data fields
  XYZMax[0] = XYZMax[1] = XYZMax[2] = -1.0e+9;
  XYZMin[0] = XYZMin[1] = XYZMin[2] = +1.0e+9;
  int no;
  int NXPoints, NYPoints, NZPoints;
  RWGObject *O;
  for(no=0; no<G->NumObjects; no++)
   { 
     O=G->Objects[no];
     GetXYZMaxMin(O, XYZMaxTO, XYZMinTO);

     XYZMax[0] = fmax(XYZMax[0], XYZMaxTO[0]);
     XYZMax[1] = fmax(XYZMax[1], XYZMaxTO[1]);
     XYZMax[2] = fmax(XYZMax[2], XYZMaxTO[2]);
     XYZMin[0] = fmin(XYZMin[0], XYZMinTO[0]);
     XYZMin[1] = fmin(XYZMin[1], XYZMinTO[1]);
     XYZMin[2] = fmin(XYZMin[2], XYZMinTO[2]);

     // FIXME to handle 1D periodicity
     if ( O->MP->IsPEC() )
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
      { NXPoints = (XYZMaxTO[0] - XYZMinTO[0]) / PBCGeometry::DeltaInterp;
        NYPoints = (XYZMaxTO[1] - XYZMinTO[1]) / PBCGeometry::DeltaInterp;
        NZPoints = (XYZMaxTO[2] - XYZMinTO[2]) / PBCGeometry::DeltaInterp;
        Log("Creating %ix%ix%i interpolation table for object %s",NXPoints,NYPoints,NZPoints,O->Label);
        GBarAB9_Interior[no]=new Interp3D( XYZMinTO[0], XYZMaxTO[0], NXPoints+1,
                                           XYZMinTO[1], XYZMaxTO[1], NYPoints+1,
                                           XYZMinTO[2], XYZMaxTO[2], NZPoints+1,
                                           2, 0, 0, 0);
      };

   };

  /*--------------------------------------------------------------*/
  /*- allocate interpolator for exterior medium ------------------*/
  /*--------------------------------------------------------------*/
  NXPoints = (XYZMax[0] - XYZMin[0]) / PBCGeometry::DeltaInterp; 
  NYPoints = (XYZMax[1] - XYZMin[1]) / PBCGeometry::DeltaInterp; 
  NZPoints = (XYZMax[2] - XYZMin[2]) / PBCGeometry::DeltaInterp; 
  if (NZPoints<2) 
   NZPoints=2;
  GBarAB9_Exterior=new Interp3D( XYZMin[0], XYZMax[0], NXPoints+1,
                                 XYZMin[1], XYZMax[1], NYPoints+1,
                                 XYZMin[2], XYZMax[2], NZPoints+1,
                                 2, 0, 0, 0);

}

} // namespace scuff
