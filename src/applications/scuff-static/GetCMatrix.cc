/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * GetCMatrix.cc
 *              
 *
 * homer reid    
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>

#include <libscuff.h>
#include <libSpherical.h>
#include <SSSolver.h>
#include <cmatheval.h>

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
#define MAXEPF   10    // max number of evaluation-point files
#define MAXCACHE 10    // max number of cache files for preload

#define MAXSTR   1000

/***************************************************************/
/* spherical-wave external field *******************************/
/***************************************************************/
typedef struct PESData 
 { 
   int l, m;

 } PESData;

// 1/sqrt(4*pi)
#define OORT4PI 0.28209479177387814347 
void PhiESpherical(double *x, void *UserData, double PhiE[4])
{ 
  PESData *PESD = (PESData *)UserData;

  int l = PESD->l;
  int m = PESD->m;
  double r, Theta, Phi;

  CoordinateC2S(x, &r, &Theta, &Phi);

  memset(PhiE, 0, 4.0*sizeof(double));
  if (l==0)
   { 
     PhiE[0] = -1.0 / OORT4PI;
   }
  else if (l==1 && m==1)
   { 
     PhiE[0] = -x[0] / OORT4PI;
     PhiE[1] = +1.0 / OORT4PI;
   }
  else if (l==1 && m==-1)
   { 
     PhiE[0] = -x[1] / OORT4PI;
     PhiE[2] = +1.0 / OORT4PI;
   }
  else if (l==1 && m==0)
   { 
     PhiE[0] = -x[2] / OORT4PI;
     PhiE[3] = +1.0 / OORT4PI;
   }
  else
   ErrExit("l>1 multipoles not supported");
 
} 


/***************************************************************/
/***************************************************************/
/***************************************************************/
HMatrix *GetCMatrix(SSSolver *SSS, HMatrix *M, HVector *Sigma, HMatrix *C)
{
  PESData PESDBuffer, *PESD = &PESDBuffer;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  int lMax = 1;
  int NAlpha = (lMax+1)*(lMax+1);

  if (C && ( C->NR!=NAlpha || C->NC!=NAlpha) )
   { Warn("incorrect matrix in GetCMatrix (reallocating)...");
     delete C;
     C=0;
   };
  if (!C)
   C = new HMatrix(NAlpha, NAlpha);
  C->Zero();

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  HVector *Moments = new HVector(NAlpha);
  for(int Alpha=0, l=0; l<=lMax; l++)
   for(int m=-l; m<=l; m++, Alpha++)
    { 
      /*--------------------------------------------------------------*/
      /*- solve the electrostatic problem with 'incident' field given */
      /*- by an (l,m) spherical wave                                  */
      /*--------------------------------------------------------------*/
      PESD->l = l;
      PESD->m = m;
      SSS->AssembleRHSVector(0, PhiESpherical, (void *)PESD, Sigma);
      M->LUSolve(Sigma);

      /*--------------------------------------------------------------*/
      /*- compute the spherical moments of the induced surface charge */
      /*- distribution and insert them as one full row of the C matrix*/
      /*--------------------------------------------------------------*/
      SSS->GetSphericalMoments(Sigma, 0, lMax, Moments);
      for(int AlphaP=0, lP=0; lP<=lMax; lP++)
       for(int mP=-lP; mP<=lP; mP++, AlphaP++)
        C->SetEntry(Alpha, AlphaP, Moments->GetEntry(AlphaP));
    };

  delete Moments;
  return C;

}

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{
  /***************************************************************/  
  /* process options *********************************************/
  /***************************************************************/
  InstallHRSignalHandler();
  char *GeoFile = 0;
  char *LambdaFile = 0;
  int lMax = 1;
  char *Cache=0;
  /* name               type    #args  max_instances  storage           count         description*/
  OptStruct OSArray[]=
   { 
     {"geometry",       PA_STRING,  1, 1,       (void *)&GeoFile,    0,             "geometry file"},
/**/
     {"LambdaFile",     PA_STRING,  1, 1,       (void *)&LambdaFile, 0,             "file of lambda values"},
/**/
     {"lMax",           PA_INT,     1, 1,       (void *)&lMax,       0,             "maximum L value"},
/**/
     {"Cache",          PA_STRING,  1, 1,       (void *)&Cache,      0,             "read/write cache"},
/**/
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  if (GeoFile==0)
   OSUsage(argv[0], OSArray, "--geometry option is mandatory");

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  SetLogFileName("GetCMatrix.log");
  Log("GetCMatrix running on %s",GetHostName());

  /*******************************************************************/
  /* figure out what the user wants to do about Lambda values ********/
  /*******************************************************************/
  HVector *LambdaVector;
  if (LambdaFile)
   { LambdaVector = new HVector(LambdaFile);
     if (LambdaVector->ErrMsg) 
      ErrExit(LambdaVector->ErrMsg);
   }
  else
   { LambdaVector = new HVector(1);
     LambdaVector->SetEntry(0,0.0);
   };

  /*******************************************************************/
  /* create the ScuffStaticGeometry **********************************/
  /*******************************************************************/
  SSSolver *SSS   = new SSSolver(GeoFile);
  RWGGeometry *G  = SSS->G;
  HMatrix *M      = SSS->AllocateBEMMatrix();
  HVector *Sigma  = SSS->AllocateRHSVector();

  /*******************************************************************/
  /* preload the scuff cache with any cache preload files the user   */
  /* may have specified                                              */
  /*******************************************************************/
  if (Cache)
   PreloadCache( Cache );

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  char *OutFileName=vstrdup("%s.CMatrix",GetFileBase(GeoFile));
  unlink(OutFileName);

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  int NAlpha = (lMax+1)*(lMax+1);
  HMatrix *C = new HMatrix(NAlpha, NAlpha);
  for(int n=0; n<LambdaVector->N; n++)
   { 
     /*******************************************************************/
     /* assemble and factorize the BEM matrix                           */
     /*******************************************************************/
     double Lambda = LambdaVector->GetEntryD(n);
     SSS->G->RegionMPs[1]->SetEps( cdouble(0.0,-Lambda) );
     SSS->AssembleBEMMatrix(M);
     M->LUFactorize();
     if (Cache)
      StoreCache( Cache );

     /*******************************************************************/
     /*******************************************************************/
     /*******************************************************************/
     GetCMatrix(SSS, M, Sigma, C);

     FILE *f=fopen(OutFileName,"a");
     fprintf(f,"%e ",Lambda);
     for(int nr=0; nr<NAlpha; nr++)
      for(int nc=0; nc<NAlpha; nc++)
       fprintf(f,"%e ",C->GetEntryD(nr,nc));
     fprintf(f,"\n");
     fclose(f);
   };

}
