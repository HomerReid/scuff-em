/*
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <libscuff.h>
#include <PBC.h>

namespace scuff{

PBCAccelerator *CreatePBCAccelerator(RWGGeometry *G)
{
  /*--------------------------------------------------------------*/
  /*- figure out the maximum height differential between any      */
  /*- two points in the geometry so we know how large to make the */
  /*- Z-range of the interpolation table                          */
  /*--------------------------------------------------------------*/
  double ZMax=-1e9, ZMin=1e9;
  int no, nv;
  RWGObject *O;
  for(no=0; no<G->NumObjects; no++)
   { 
      O=G->Objects[no];
      for(nv=0; nv<O->NumVertices; nv++)
       { ZMax = fmax(ZMax, O->Vertices[3*nv + 2]);
         ZMin = fmin(ZMin, O->Vertices[3*nv + 2]);
       };
   };

  /*--------------------------------------------------------------*/
  /*- initialize the acceleration table --------------------------*/
  /*--------------------------------------------------------------*/

}

} // namespace scuff
