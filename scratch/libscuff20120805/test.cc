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

  RWGGeometry *G = new RWGGeometry(GeoFileName);

}
