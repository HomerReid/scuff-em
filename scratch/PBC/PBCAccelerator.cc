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
  /*--------------------------------------------------------------*/ 
  /*--------------------------------------------------------------*/ 
  double ZMax=-1e9, ZMin=1e9;
  int no, nv;
  RWGObject *O;
  for(no=0; no<G->NumObjects; no++)
   { 
      O=G->Objects
      for(nv=0; nv<
    
   };

}



}
