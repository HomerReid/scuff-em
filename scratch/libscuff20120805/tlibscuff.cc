#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <libhrutil.h>
#include "libscuff.h"

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[]) 
{  
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  char *GeoFileName=0;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"geometry", PA_STRING,  1, 1, (void *)&GeoFileName,  0,  ".scuffgeo file"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (GeoFileName==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGGeometry *G = new RWGGeometry(GeoFileName);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  printf("\nRegions: \n");
  for(int nr=0; nr<G->NumRegions; nr++)
   printf("%i: label=%s, material=%s\n", nr, G->RegionLabels[nr], G->RegionMPs[nr]->Name);
  printf("\n\n");

  printf("Surfaces: \n");
  RWGSurface *S;
  for(int ns=0; ns<G->NumSurfaces; ns++)
   { S=G->Surfaces[ns];
     printf("Surface %i: label=%s (%s)\n", ns, S->Label, S->IsClosed ? "closed" : "open");
     printf(" %i panels /%i interior edges / %i exterior edges \n",
              S->NumPanels, S->NumEdges, S->NumExteriorEdges);
     printf(" regions %s (%i), %s (%i)\n", S->RegionLabels[0], S->RegionIndices[0], 
                                           S->RegionLabels[1], S->RegionIndices[1]);
     printf("\n");
   };
  printf("\n\n");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  double X[3];
  X[0]=0.1; 
  X[1]=0.1; 
  int nr;
  for( X[2]=-1.5; X[2]<=+1.5; X[2]+=1.0 )
   { nr=G->GetRegionIndex(X);
     printf("\n");
     printf("point %g %g %g: in region %i (%s) \n",X[0],X[1],X[2],nr,G->RegionLabels[nr]);
     printf("\n");
     for(nr=0; nr<G->NumRegions; nr++)
      printf(" %s region %i (%s)\n", G->PointInRegion(nr,X) ? "in" : "not in", nr, G->RegionLabels[nr]);
   };

}   
