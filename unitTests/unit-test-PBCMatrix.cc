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
 * unit-test-compactMatrix.cc -- SCUFF-EM unit test for the full BEM 
 *                               matrix of a geometry consisting of 
 *                               extended objects (in this case, an array of 
 *                               gold spheres atop a dielectric substrate)
 * 
 * homer reid        -- 11/2005 -- 10/2011
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <libhrutil.h>
#include "libscuff.h"
#include "libscuffInternals.h"

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{ 
  SetLogFileName("scuff-unit-tests.log");
  Log("SCUFF-EM PBC matrix unit test running on %s",GetHostName());
 
  /*--------------------------------------------------------------*/
  /*- use the --WriteFiles option to create the .hdf5 files that  */
  /*- are used as the comparison for subsequent runs.             */
  /*-                                                             */
  /*- when --WriteFiles is specified, the code automatically      */
  /*- writes out the scuff cache to a file named TwoSpheres.cache.*/
  /*-                                                             */
  /*- on subsequent runs, you can use the --UseCache option to    */
  /*- tell the code to use the cache; if this option is absent,   */
  /*- the code will not use the cache file even if it is present  */
  /*- in the working directory.                                   */
  /*--------------------------------------------------------------*/
  bool WriteFiles=0;
  bool UseCache=0;
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   { {"WriteFiles", PA_BOOL, 0, 1, (void *)&WriteFiles,  0,  "write .hdf5 files"},
     {"UseCache",   PA_BOOL, 0, 1, (void *)&UseCache,    0,  "use cache file"},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGGeometry *G = new RWGGeometry("SphereSubstrateArray_Coarse.scuffgeo");
  G->SetLogLevel(SCUFF_VERBOSELOGGING);
  HMatrix *M = G->AllocateBEMMatrix();
  double Elapsed;
<<<<<<< HEAD
=======
  double AvgMag;
>>>>>>> PBCV2P0
  double kBloch[2];

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (WriteFiles)
   { 
     void *pCC = HMatrix::OpenHDF5Context("SphereSubstrateArray_Coarse.hdf5");
   
     printf("Assembling BEM matrix at (Omega,kx,ky) = (0.05, 0.0, 0.0)...\n");
     kBloch[0]=kBloch[1]=0.0;
     Tic();
     G->AssembleBEMMatrix( cdouble(0.05, 0.0), kBloch, M );
     Elapsed=Toc();
     printf("...done in %.1f s; exporting to HDF5\n",Elapsed);
     M->ExportToHDF5(pCC,"M1");
     StoreCache("SphereSubstrateArray_Coarse.cache");
   
     printf("Assembling BEM matrix at (Omega,kx,ky) = (0.9,0.3,1.1)...\n");
     kBloch[0]=0.3;
     kBloch[1]=1.1;
     Tic();
     G->AssembleBEMMatrix( cdouble(0.9, 0.0), kBloch, M );
     Elapsed=Toc();
     printf("...done in %.1f s; exporting to HDF5\n",Elapsed);
     M->ExportToHDF5(pCC,"M2");

     HMatrix::CloseHDF5Context(pCC);

     StoreCache("SphereSubstrateArray_Coarse.cache");

   }
  else
   {
     if (UseCache)
      PreloadCache("SphereSubstrateArray_Coarse.cache");

     /***************************************************************/
     /* assemble and compare in first case   ***********************/
     /***************************************************************/
     HMatrix *MRef = new HMatrix("SphereSubstrateArray_Coarse.hdf5",LHM_HDF5,"M1");
     if (MRef->ErrMsg)
      ErrExit(MRef->ErrMsg);

     printf("Assembling BEM matrix at Omega = 0.05...\n");
     kBloch[0]=kBloch[1]=0.0;
     Tic();
     G->AssembleBEMMatrix( cdouble(0.05, 0.0), kBloch, M );
     Elapsed=Toc();
     printf("...done in %.0f s.\n",Elapsed);

<<<<<<< HEAD
=======
     AvgMag = 0.0;
     for(int nr=0; nr<M->NR; nr++)
      for(int nc=0; nc<M->NC; nc++)
       AvgMag += abs(M->GetEntry(nr,nc));
     AvgMag /= ((double)(M->NR * M->NC) );

>>>>>>> PBCV2P0
     printf("Comparing to reference...");
     int Mismatches=0;
     for(int nr=0; nr<M->NR; nr++)
      for(int nc=nr; nc<M->NC; nc++)
<<<<<<< HEAD
       if ( !EqualFloat(M->GetEntry(nr,nc), MRef->GetEntry(nr,nc)) )
=======
       if (    abs(M->GetEntry(nr,nc)) > (1.0e-6 * AvgMag)
            && RD(M->GetEntry(nr,nc), MRef->GetEntry(nr,nc)) > 1.0e-6 
          )
>>>>>>> PBCV2P0
        { 
          Mismatches++;
          Warn("MISMATCH(%i,%i): (%.8e,%.8e)--(%.8e,%.8e) (%.1e)\n",nr,nc,
                real(M->GetEntry(nr,nc)),    imag(M->GetEntry(nr,nc)),
                real(MRef->GetEntry(nr,nc)), imag(MRef->GetEntry(nr,nc)), 
                RD(M->GetEntry(nr,nc),MRef->GetEntry(nr,nc)));
          if (Mismatches>10)
           ErrExit("too many mismatches");
        };
      printf("...%i mismatches.\n",Mismatches);
      printf("\n");
      

     /***************************************************************/
     /* assemble and compare in second case *************************/
     /***************************************************************/
     MRef = new HMatrix("SphereSubstrateArray_Coarse.hdf5",LHM_HDF5,"M2");
     if (MRef->ErrMsg)
      ErrExit(MRef->ErrMsg);

     printf("Assembling BEM matrix at (Omega,kx,ky) = (0.9,0.3,1.1)...\n");
     kBloch[0]=0.3;
     kBloch[1]=1.1;
     Tic();
     G->AssembleBEMMatrix( cdouble(0.90, 0.00), kBloch, M );
     Elapsed=Toc();
     printf("...done in %.0f s.\n",Elapsed);

<<<<<<< HEAD
     printf("Comparing to reference...");
     for(int nr=0; nr<M->NR; nr++)
      for(int nc=nr; nc<M->NC; nc++)
       if ( !EqualFloat(M->GetEntry(nr,nc), MRef->GetEntry(nr,nc)) )
=======
     AvgMag = 0.0;
     for(int nr=0; nr<M->NR; nr++)
      for(int nc=0; nc<M->NC; nc++)
       AvgMag += abs(M->GetEntry(nr,nc));
     AvgMag /= ((double)(M->NR * M->NC) );

     printf("Comparing to reference...");
     for(int nr=0; nr<M->NR; nr++)
      for(int nc=nr; nc<M->NC; nc++)
       if (    abs(M->GetEntry(nr,nc)) > (1.0e-6 * AvgMag)
            && RD(M->GetEntry(nr,nc), MRef->GetEntry(nr,nc)) > 1.0e-6 
          )
>>>>>>> PBCV2P0
        { 
          Mismatches++;
          Warn("MISMATCH(%i,%i): (%.8e,%.8e)--(%.8e,%.8e) (%.1e)\n",nr,nc,
                real(M->GetEntry(nr,nc)),    imag(M->GetEntry(nr,nc)),
                real(MRef->GetEntry(nr,nc)), imag(MRef->GetEntry(nr,nc)), 
                RD(M->GetEntry(nr,nc),MRef->GetEntry(nr,nc)));
          if (Mismatches>10)
           ErrExit("too many mismatches");
        };
      printf("%i mismatches.\n",Mismatches);

   }; // if (WriteFiles) ... else 

  printf("All tests successfully passed.\n");
  return 0;

}
