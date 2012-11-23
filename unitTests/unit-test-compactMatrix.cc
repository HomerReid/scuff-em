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
 *                               compact objects (two spheres in    
 *                               this case)
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
  Log("SCUFF-EM compact matrix unit test running on %s",GetHostName);
 
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
   { {"WriteFiles", PA_BOOL, 0, 1, (void *)WriteFiles,  0,  "write .hdf5 files"},
     {"UseCache",   PA_BOOL, 0, 1, (void *)UseCache,    0,  ""},
     {0,0,0,0,0,0,0}
   };
  ProcessOptions(argc, argv, OSArray);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  RWGGeometry *G = new RWGGeometry("TwoSpheres.scuffgeo");
  HMatrix *M = G->AllocateBEMMatrix();

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (WriteFiles)
   { 
     void *pCC = HMatrix::CreateHDF5Context("TwoSpheres.hdf5");
   
     printf("Assembling BEM matrix at Omega = 0.05...\n");
     Tic();
     M->AssembleBEMMatrix( cdouble(0.05, 0.0) );
     Elapsed=Toc();
     printf("...done in %.1 s; exporting to HDF5\n",Elapsed);
     M->ExportToHDF5(pCC,"M0P05");
   
     printf("Assembling BEM matrix at Omega = 5+5I...\n");
     Tic();
     M->AssembleBEMMatrix( cdouble(5.00, 5.00) );
     Elapsed=Toc();
     printf("...done in %.1 s; exporting to HDF5\n",Elapsed);
     M->ExportToHDF5(pCC,"M5P5I");

     HMatrix::CloseHDF5Context(pCC);

     StoreCache("TwoSpheres.cache");

   }
  else
   {
     /***************************************************************/
     /***************************************************************/
     /***************************************************************/
     MRef = new HMatrix("TwoSpheres.hdf5",LHM_HDF5,"M0P05");
     if (MRef->ErrMsg)
      ErrExit(MRef->ErrMsg);

     if (UseCache)
      PreloadCache("TwoSpheres.cache");

     printf("Assembling BEM matrix at Omega = 0.05...\n");
     Tic();
     M->AssembleBEMMatrix( cdouble(0.05, 0.0) );
     Elapsed=Toc();
     printf("...done in %.1 s);

     printf("Comparing to reference...");
     int BadEntries=0;
     for(int nr=0; nr<M->NR; nr++)
      for(int nc=0; nc<M->NC; nc++)
       if ( !VecEqualFloat(

   };

}
