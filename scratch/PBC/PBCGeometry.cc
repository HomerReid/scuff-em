/*
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <libscuff.h>
#include <PBC.h>

namespace scuff{

#define DELTAINTERP 0.05 // FIXME

/***************************************************************/
/***************************************************************/
/***************************************************************/
void GetXYZMaxMin(RWGObject *O, double XYZMax[3], XYZMin[3])
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

/***************************************************************/
/* PBCGeometry class constructor *******************************/
/***************************************************************/
PBCGeometry::PBCGeometry(RWGGeometry *pG, double pLBV[2][2])
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
  EpsTF = (cdouble *)mallocSE( (G->NumObjects+1)*sizeof(cdouble));
  MuTF  = (cdouble *)mallocSE( (G->NumObjects+1)*sizeof(cdouble));

  /*--------------------------------------------------------------*/
  /*- add straddlers to all objects ... maybe FIXME to prevent   -*/
  /*- G and its constituent objects from being modified?         -*/
  /*--------------------------------------------------------------*/
  NumStraddlers=(int *)mallocEC(2*G->NumObjects*sizeof(int));
  G->TotalBFs=G->TotalPanels=0;
  for(int no=0; no<G->NumObjects; no++)
   { RWGObject *O = G->Objects[no];
     AddStraddlers(O, LBV, NumStraddlers + 2*no);
     G->TotalBFs    += O->NumBFs;
     G->TotalPanels += O->NumPanels;
     if ( no+1 < G->NumObjects )
      { G->BFIndexOffset[no+1]=G->BFIndexOffset[no] + O->NumBFs;
        G->PanelIndexOffset[no+1]=G->PanelIndexOffset[no] + O->NumPanels;
      };
   };

  /*--------------------------------------------------------------*/
  /*- allocate memory for the contributions of the innermost     -*/
  /*- lattice cells to the BEM matrix                            -*/
  /*--------------------------------------------------------------*/
  MPP=new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX);
  MPM=new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX);
  MPZ=new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX);
  MZP=new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX);
  MZZ=new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX); // this one could be symmetric ...

  /*--------------------------------------------------------------*/
  /*- allocate interpolators for each object interior             */
  /*--------------------------------------------------------------*/
  GBarAB9_Interior=(Interp3D **)mallocSE(NumObject * sizeof(Interp3D *));
  double XYZMaxTO[3], XYZMinTO[3]; // 'x, y, z max/min, this object'
  double XYZMax[3], XYZMin[3];     // x,y,z max/min for geometry overall
  XYZMax[0] = XYZMax[1] = XYZMax[2] = -1.0e+9;
  XYZMin[0] = XYZMin[1] = XYZMin[2] = +1.0e+9;
  int no;
  int NXPoints, NYPoints, NZPoints;
  RWGObject *O;
  for(no=0; no<G->NumObjects; no++)
   { 
     O=G->Objects[no];
     GetMaxMinCoords(O, XYZMaxTO, XYZMinTO);
     XYZMax[0] = fmax(XYZMax[0], XYZMaxTO[0]);
     XYZMax[1] = fmax(XYZMax[1], XYZMaxTO[1]);
     XYZMax[2] = fmax(XYZMax[2], XYZMaxTO[2]);
     XYZMin[0] = fmin(XYZMin[0], XYZMinTO[0]);
     XYZMin[1] = fmin(XYZMin[1], XYZMinTO[1]);
     XYZMin[2] = fmin(XYZMin[2], XYZMinTO[2]);

     if ( O->MP->IsPEC() || (NumStraddlers[2*no+0]==0 && NumStraddlers[2*no+1]==0) )
      GBarAB9_Interior[no]=0;
     else
      { NXPoints = (XYZMaxTO[0] - XYZMinTO[0]) / DELTAINTERP; 
        NYPoints = (XYZMaxTO[1] - XYZMinTO[1]) / DELTAINTERP; 
        NZPoints = (XYZMaxTO[2] - XYZMinTO[2]) / DELTAINTERP; 
        GBarAB9_Interior[no]=new Interp3D( XYZMinTO[0], XYZMaxTO[0], NXPoints+1,
                                           XYZMinTO[1], XYZMaxTO[1], NXPoints+1,
                                           XYZMinTO[2], XYZMaxTO[2], NXPoints+1,
                                           1, 0, 0, 0);
      };

   };

  /*--------------------------------------------------------------*/
  /*- allocate interpolator for exterior medium ------------------*/
  /*--------------------------------------------------------------*/
  NXPoints = (XYZMax[0] - XYZMin[0]) / DELTAINTERP; 
  NYPoints = (XYZMax[1] - XYZMin[1]) / DELTAINTERP; 
  NZPoints = (XYZMax[2] - XYZMin[2]) / DELTAINTERP; 
  GBarAB9_Exterior=new Interp3D( XYZMin[0], XYZMax[0], NXPoints+1,
                                 XYZMin[1], XYZMax[1], NXPoints+1,
                                 XYZMin[2], XYZMax[2], NXPoints+1,
                                 1, 0, 0, 0);

}

/***************************************************************/
/***************************************************************/
/***************************************************************/


} // namespace scuff
