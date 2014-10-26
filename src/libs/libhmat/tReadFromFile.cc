#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "libhrutil.h"
#include "libhmat.h"

int main(int argc, char *argv[])
{ 
  HMatrix *M;
  int nr, nc, narg;
  char Options[1000];

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (argc<2)
   { fprintf(stderr,"usage: %s DataFile [options]\n",argv[0]);
     exit(1);
   };
  Options[0]=0;
  for(narg=2; narg<argc; narg++)
   { strcat(Options," ");
     strcat(Options,argv[narg]);
   };

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  M=new HMatrix(argv[1],LHM_TEXT,Options);
  if (M->ErrMsg)
   ErrExit(M->ErrMsg);

  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  printf("Your matrix is %i x %i: \n",M->NR, M->NC);
  for(nr=0; nr<M->NR; nr++)
   for(nc=0; nc<M->NC; nc++)
    printf("%+7.2e %c",M->GetEntryD(nr,nc),nc==(M->NC-1) ? '\n' : ' ');

  printf("\n\nThank you for your support.\n");

}
