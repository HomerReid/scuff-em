/*
 * tQR.cc -- test libhmat's QR routine
 *
 * usage: tQR --MatrixFile MyFile.txt
 *
 * where MyFile.txt is a text file containing your matrix 
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "libhrutil.h"
#include "libhmat.h"

int main(int argc, char *argv[])
{  
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  char *MatrixFile=0;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"MatrixFile", PA_STRING,  1, 1, (void *)&MatrixFile,  0,  "matrix text file"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);
  if (MatrixFile==0)
   OSUsage(argv[0],OSArray,"--MatrixFile option is mandatory");

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HMatrix *M = new HMatrix(MatrixFile);
  if (M->ErrMsg)
   ErrExit(M->ErrMsg);
  
  HMatrix *Q=0, *R=0;

  M->QR(&Q, &R);

  int nr, nrMax = M->NR > 5 ? 5 : M->NR;
  int nc, ncMax = M->NC > 5 ? 5 : M->NC;
  if (M->RealComplex == LHM_REAL)
   { 
     printf("Q: \n");
     for (nr=0; nr<nrMax; nr++)
      { for (nc=0; nc<ncMax; nc++)
         printf("%5.2f ",Q->GetEntryD(nr,nc));
        if (nc!=M->NC) printf(" ... ");
        printf("\n");
      };

     printf("R: \n");
     for (nr=0; nr<nrMax; nr++)
      { for (nc=0; nc<ncMax; nc++)
         printf("%5.2f ",R->GetEntryD(nr,nc));
        if (nc!=M->NC) printf(" ... ");
        printf("\n");
      };

   }
  else
   { 
     SetDefaultCD2SFormat("(%5.2f,%5.2f)");
     printf("Q: \n");
     for (nr=0; nr<nrMax; nr++)
      { for (nc=0; nc<ncMax; nc++)
         printf("%s ",CD2S(Q->GetEntry(nr,nc)));
        if (nc!=M->NC) printf(" ... ");
        printf("\n");
      };

     printf("R: \n");
     for (nr=0; nr<nrMax; nr++)
      { for (nc=0; nc<ncMax; nc++)
          printf("%s ",CD2S(R->GetEntry(nr,nc)));
        if (nc!=M->NC) printf(" ... ");
        printf("\n");
      };
   };

}
