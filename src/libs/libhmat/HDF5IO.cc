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
 * HDF5IO.cc   -- libhmat class methods for exporting/importing HMatrices 
 *             -- and HVectors to/from binary HDF5 data files 
 *
 * homer reid  -- 12/2009
 *
 * --------------------------------------------------------------
 *
 *  1. exporting: 
 *
 * there are two basic ways to write HMatrices and HVectors to 
 * binary HDF5 files:
 *
 * (1) if you want each matrix and vector to be written to a 
 *     separate binary file, you can simply make one call to 
 *     ExportToHDF5() for each matrix/vector:
 * 
 *       M1->ExportToHDF5("MyFile1.hdf5","M1");
 *       M2->ExportToHDF5("MyFile2.hdf5","M2");
 *       V->ExportToHDF5("MyFile3.hdf5","V");
 *     
 *     this way, will get three separate binary HDF5 files, each 
 *     containing a single HDF5 dataset.
 *
 *
 * (2) if you want multiple matrices and/or vectors to be 
 *     written to the same binary file, you first create an 
 *     'HDF5 context' (referenced by an opaque pointer). 
 *     the routine that creates an HDF5 context is a 
 *     static class method inside the HMatrix class:
 * 
 *       void *pHC = HMatrix::OpenHDF5Context("MyFile.hdf5");
 * 
 *     as above, you make one call to ExportToHDF5() for each 
 *     matrix/vector you want to write to the file:
 *      
 *       M1->ExportToHDF5(pHC, "M1");
 *       M2->ExportToHDF5(pHC, "M2");
 *       V->ExportToHDF5(pHC, "V ");
 *     
 *     when you are done, go like this to close up the file:
 *     
 *       HMatrix::CloseHDF5Context(pHC);
 *     
 *     this way, you will have just a single binary HDF5 file 
 *     containing three HDF5 datasets.
 *
 * --------------------------------------------------------------
 *
 *  2. importing: 
 *
 *     HMatrix *M1 = new HMatrix("MyFile.hdf5", LHM_HDF5, "M1");
 *     HMatrix *M2 = new HMatrix("MyFile.hdf5", LHM_HDF5, "M2");
 *     HVector *V  = new HVector("MyFile.hdf5", LHM_HDF5, "V");
 *
 * --------------------------------------------------------------
 *
 *  3. error checking: 
 *
 *     in general, if anything fails, the ErrMsg field in the 
 *     structure will be non-null on return and will point to 
 *     a string that explains what happened. 
 *     in particular, if you try to create an HMatrix or HVector
 *     from an HDF5 that doesn't exist or doesn't contain the
 *     entity you are looking for, you will get back a class
 *     instance in which all fields are invalid except ErrMsg. 
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include <libhrutil.h>
#include "libhmat.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

// almost all of this file requires HDF5; dummy versions of 
// functions for use when compiling without HDF5 start down
// around line 500
#ifdef HAVE_HDF5

#include <hdf5.h>
#include <H5LTpublic.h>
#include <H5Epublic.h>


/***************************************************************/
/* this is a data structure containing everything that i need  */
/* to keep tratk of regarding an open HDF5 file.               */
/***************************************************************/
typedef struct HDF5Context
 { 
   hid_t file_id;

 } HDF5Context;

/***************************************************************/
/***************************************************************/
/***************************************************************/
void *HMatrix::OpenHDF5Context(const char *format, ...)
{ 
  va_list ap;
  char FileName[1000];
  va_start(ap,format);
  vsnprintfEC(FileName,1000,format,ap);
  va_end(ap);

  // turn off the HDF5 console error messages
  H5Eset_auto2( 0, 0, 0 );

  hid_t file_id = H5Fcreate(FileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  if (file_id>=0)
   { HDF5Context *HC=(HDF5Context *)mallocEC(sizeof(*HC));
     HC->file_id=file_id; 
     return HC;
   };
  return 0;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
void HMatrix::CloseHDF5Context(void *pHC)
{
  HDF5Context *HC=(HDF5Context *)pHC;
  H5Fclose(HC->file_id);
  free(HC);
}

/***************************************************************************/
/* export matrix to HDF5 file.                                             */
/*                                                                         */
/* how it works:                                                           */
/*                                                                         */
/*  1. for real-valued non-packed matrices (StorageType==LHM_NORMAL and    */
/*     RealComplex==LHM_REAL) the matrix is stored as an ordinary HDF5     */
/*     matrix and can be read in by other applications.  Note, however,    */
/*     that an HMatrix is stored in column-major (contiguous columns)      */
/*     format, whereas HDF5 uses the row-major(contiguous rows) format,    */
/*     and as a consequence the dimensions in the HDF5 file are reversed:  */
/*     an NR x NC matrix is stored as an NC x NR dataset.   [To read it    */
/*     as NR x NC in Matlab, you can either transpose the matrix after     */
/*     reading it or use hdf5read(..., 'V71Dimensions', true).]            */
/*                                                                         */
/*  2. for complex-valued non-packed matrices (StorageType==LHM_NORMAL and */
/*     RealComplex==LHM_COMPLEX) the matrix is stored as an ordinary HDF5  */ 
/*     matrix, but with 2x the rows dimension. The real and imaginary      */
/*     parts of the (M,N) entry in the original matrix are stored in the   */
/*     (N,2*M) and (N,2*M+1) entries of the matrix stored in the HDF5      */
/*     file [note the dimension reversal from above], and hence the        */
/*     original matrix can be recovered in MATLAB by reading the matrix    */
/*     stored in the HDF5 file and extracting the even/odd entries.        */
/*                                                                         */
/*   3. for packed matrices (StorageType==LHM_SYMMETRIC or LHM_HERMITIAN), */
/*      the content of the internal data buffer within the HMat is stored  */
/*      as the first row of a 1x(NE) (LHM_REAL) or 1x(2*NE) (LHM_COMPLEX)  */
/*      matrix, where NE is the number of elements in the reduced storage  */
/*      buffer.                                                            */
/*                                                                         */
/*   4. in all cases, StorageType and RealComplex are written to the HDF5  */
/*      file as integer attributes associated with the dataset.            */
/***************************************************************************/
void HMatrix::ExportToHDF5(void *pHC, const char *format, ...)
{ 
  HDF5Context *HC=(HDF5Context *)pHC;
  if (!HC) return;
  herr_t Status;
  hsize_t dims[2];

  va_list ap;
  char Name[1000];
  va_start(ap,format);
  vsnprintfEC(Name,1000,format,ap);
  va_end(ap);

  // turn off the HDF5 console error messages
  H5Eset_auto2( 0, 0, 0 );

  /*--------------------------------------------------------------*/
  /*- write the data ---------------------------------------------*/
  /*--------------------------------------------------------------*/
  Status=-1;
  if (StorageType==LHM_NORMAL && RealComplex==LHM_REAL )
   { 
     dims[0]=NC;
     dims[1]=NR;
     Status=H5LTmake_dataset_double(HC->file_id, Name, 2, dims, DM);
   }
  else if (StorageType==LHM_NORMAL  && RealComplex==LHM_COMPLEX )
   { dims[0]=NC;
     dims[1]=2*NR;
     Status=H5LTmake_dataset_double(HC->file_id, Name, 2, dims, (double *)ZM);
   }
  else if ( StorageType!=LHM_NORMAL && RealComplex==LHM_REAL )
   { 
     dims[0]=1;
     dims[1]=NumEntries();
     Status=H5LTmake_dataset_double(HC->file_id, Name, 2, dims, (double *)DM);
   }
  else if ( StorageType!=LHM_NORMAL && RealComplex==LHM_COMPLEX)
   { dims[0]=1;
     dims[1]=2*NumEntries();
     Status=H5LTmake_dataset_double(HC->file_id, Name, 2, dims, (double *)ZM);
   };

  if (Status!=0)
   ErrExit("%s:%i: internal error",__FILE__,__LINE__);
  
  /*--------------------------------------------------------------*/
  /*- write the attributes ---------------------------------------*/
  /*--------------------------------------------------------------*/
  H5LTset_attribute_int(HC->file_id, Name, "Storage_Type", &StorageType, 1);
  H5LTset_attribute_int(HC->file_id, Name, "RealComplex", &RealComplex, 1);
  
} 

/***************************************************************/
/* alternative version of ExportToHDF5 that creates an HDF5    */
/* file for this dataset alone                                 */
/***************************************************************/
void HMatrix::ExportToHDF5(const char *FileName, const char *format, ...)
{ 
  va_list ap;
  char Name[1000];
  va_start(ap,format);
  vsnprintfEC(Name,1000,format,ap);
  va_end(ap);

  void *pHC=OpenHDF5Context(FileName);
  if (!pHC) 
   { ErrMsg=vstrdup("could not open file %s",FileName);
     return;
   };
  ExportToHDF5(pHC,Name);
  CloseHDF5Context(pHC);
}

/***************************************************************/
/* HMatrix constructor that attempts to construct an HMatrix   */
/* from an HDF5 file.                                          */
/* On success, ErrMsg will be NULL on return.                  */
/* Otherwise ErrMsg says what went wrong.                      */
/* If Init==true, then we assume this routine is being called  */
/* from the HMatrix() constructor and that the matrix          */
/* dimensions, etc. are to be read from the data file.         */
/* Otherwise we assume this routine is being called on an      */
/* existing HMatrix(), in which case we consider it an error   */
/* if the matrix described by the data file doesn't have the   */
/* same dimensions, etc as the existing matrix.                */
/***************************************************************/
void HMatrix::ImportFromHDF5(const char *FileName, const char *Name, bool Init)
{ 
  hid_t file_id;
  herr_t status;

  // turn off the HDF5 console error messages
  H5Eset_auto2( 0, 0, 0 );

  /*--------------------------------------------------------------*/
  /*- try to open the file ---------------------------------------*/
  /*--------------------------------------------------------------*/
  file_id = H5Fopen(FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0 )
   { ErrMsg=vstrdup("could not open file %s",FileName); 
     return; 
   };

  /*--------------------------------------------------------------*/
  /*- look for a dataset with the given name and make sure it is -*/
  /*- a matrix (two dimensions)                                  -*/
  /*--------------------------------------------------------------*/
  status=H5LTfind_dataset(file_id, Name);
  if (status==0)
   { ErrMsg=vstrdup("file %s does not contain a dataset named %s",FileName,Name); 
     return; 
   };

  int Rank;
  status=H5LTget_dataset_ndims(file_id, Name, &Rank);
  if (status<0 || Rank!=2)
   { ErrMsg=vstrdup("file %s: dataset %s is not a matrix",FileName,Name); 
     return; 
   };

  /*--------------------------------------------------------------*/
  /*- read dataset dimensions ------------------------------------*/
  /*--------------------------------------------------------------*/
  H5T_class_t class_id;
  size_t type_size;
  hsize_t dims[2];
  status=H5LTget_dataset_info(file_id, Name, dims, &class_id, &type_size);
  if (status<0)
   { ErrMsg=vstrdup("file %s: dataset %s has incorrect type",FileName,Name); 
     return; 
   };
 
  /*--------------------------------------------------------------*/
  /*- try to read attributes from file. note that these will only-*/
  /*- be present if the HDF5 file was created by the libhmat     -*/
  /*- export routines; otherwise, we suppose the file was created-*/
  /*- by another application and assume it describes a real-valued*/
  /*- non-packed matrix.                                          */
  /*--------------------------------------------------------------*/
  int FileStorageType=LHM_NORMAL;
  int FileRealComplex=LHM_REAL;
  status=H5LTget_attribute_int(file_id, Name, "Storage_Type", &FileStorageType);
  status=H5LTget_attribute_int(file_id, Name, "RealComplex", &FileRealComplex);

  int FileNR=0, FileNC=0;
  if (FileStorageType==LHM_NORMAL && FileRealComplex==LHM_REAL)
   { FileNR=dims[1];
     FileNC=dims[0];
   }
  else if (FileStorageType==LHM_NORMAL && FileRealComplex==LHM_COMPLEX)
   { FileNR=dims[1]/2;
     FileNC=dims[0];
   }
  else if (FileStorageType!=LHM_NORMAL && FileRealComplex==LHM_REAL)
   FileNR=FileNC=int( floor( sqrt( 2.0*dims[1] ) ) );
  else if (FileStorageType!=LHM_NORMAL && FileRealComplex==LHM_COMPLEX)
   FileNR=FileNC=int( floor( sqrt( 2.0*(dims[1]/2) ) ) );

  if (Init)
   { StorageType=FileStorageType;
     RealComplex=FileRealComplex;
   }
  else 
   { 
     // make sure the matrix in the file matches the size, etc
     // of the existing matrix
     if (StorageType!=FileStorageType)
      ErrMsg=vstrdup("file %s, dataset %s: storage type mismatch",FileName,Name);
     else if (RealComplex!=FileRealComplex)
      ErrMsg=vstrdup("file %s, dataset %s: data type mismatch",FileName,Name);
     else if (NR!=FileNR || NC!=FileNC)
      ErrMsg=vstrdup("file %s, dataset %s: dimension mismatch (%ix%i) != (%ix%i)",FileName,Name,NR,NC,FileNR,FileNC);
     if (ErrMsg)
      return;
   };

  /*--------------------------------------------------------------*/
  /*- now switch off to handle the various cases for reading the  */
  /*- data in                                                     */
  /*--------------------------------------------------------------*/
  if (StorageType==LHM_NORMAL && RealComplex==LHM_REAL )
   { 
     if (Init) InitHMatrix(dims[1], dims[0], LHM_REAL);
     H5LTread_dataset_double(file_id, Name, DM);
   }
  else if (StorageType==LHM_NORMAL  && RealComplex==LHM_COMPLEX )
   { 
     if (Init) InitHMatrix(dims[1]/2, dims[0], LHM_COMPLEX);
     H5LTread_dataset_double(file_id, Name, (double *)ZM);
   }
  else if ( StorageType!=LHM_NORMAL && RealComplex==LHM_REAL )
   { 
     // recover NR=NC from the NumEntries=NR*(NR+1)/2 
     int rows = int( floor( sqrt( 2.0*dims[1] ) ) );
     if (Init) InitHMatrix(rows, rows, LHM_REAL, StorageType);
     H5LTread_dataset_double(file_id, Name, DM);
   }
  else if ( StorageType!=LHM_NORMAL && RealComplex==LHM_COMPLEX )
   {
     // recover NR=NC from the NumEntries=NR*(NR+1)/2 
     int rows = int( floor( sqrt( 2.0*(dims[1]/2) ) ) );
     InitHMatrix(rows, rows, LHM_COMPLEX, StorageType);
     H5LTread_dataset_double(file_id, Name, (double *)ZM);
   };

  /*--------------------------------------------------------------*/
  /*- close up the data file -------------------------------------*/
  /*--------------------------------------------------------------*/
  H5Fclose(file_id);
  
}  

/***************************************************************/
/* Export an HVector to an HDF5 file                           */
/***************************************************************/
void HVector::ExportToHDF5(void *pHC, const char *format, ...)
{ 
  va_list ap;
  char Name[1000];
  va_start(ap,format);
  vsnprintfEC(Name,1000,format,ap);
  va_end(ap);

  // turn off the HDF5 console error messages
  H5Eset_auto2( 0, 0, 0 );

  HDF5Context *HC=(HDF5Context *)pHC;
  if (HC==0) 
   { ErrMsg=vstrdup("%s:%i: ExportToHDF5 called with pHC=void",__FILE__,__LINE__);
     return;
   };
  hsize_t dims[1];

  /*--------------------------------------------------------------*/
  /*- write the data ---------------------------------------------*/
  /*--------------------------------------------------------------*/
  if (RealComplex==LHM_REAL )
   { 
     dims[0]=N;
     H5LTmake_dataset_double(HC->file_id, Name, 1, dims, DV);
   }
  else 
   { dims[0]=2*N;
     H5LTmake_dataset_double(HC->file_id, Name, 1, dims, (double *)ZV);
   };

  /*--------------------------------------------------------------*/
  /*- write the attribute  ---------------------------------------*/
  /*--------------------------------------------------------------*/
  H5LTset_attribute_int(HC->file_id, Name, "RealComplex", &RealComplex, 1);
  
} 

/***************************************************************/
/* alternative version of HVector::ExportToHDF5 that creates   */
/* an HDF5 file for this dataset alone                         */
/***************************************************************/
void HVector::ExportToHDF5(const char *FileName, const char *format, ...)
{ 
  va_list ap;
  char Name[1000];
  va_start(ap,format);
  vsnprintfEC(Name,1000,format,ap);
  va_end(ap);

  void *pHC=HMatrix::OpenHDF5Context(FileName);
  if (!pHC) 
   { ErrMsg=vstrdup("could not open file %s",FileName);
     return; 
   };
  ExportToHDF5(pHC,Name);
  HMatrix::CloseHDF5Context(pHC);
}

/***************************************************************/
/* HVector constructor that attempts to construct an HVector   */
/* from an HDF5 file.                                          */
/***************************************************************/
void HVector::ImportFromHDF5(const char *FileName, const char *Name)
{ 
  hid_t file_id;
  herr_t status;

  // turn off the HDF5 console error messages
  H5Eset_auto2( 0, 0, 0 );

  /*--------------------------------------------------------------*/
  /*- try to open the file ---------------------------------------*/
  /*--------------------------------------------------------------*/
  file_id = H5Fopen(FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file_id < 0 )
   { ErrMsg=vstrdup("could not open file %s",FileName); 
     return; 
   };

  /*--------------------------------------------------------------*/
  /*- look for a dataset with the given name and make sure it is -*/
  /*- a vector (one dimension)                                   -*/
  /*--------------------------------------------------------------*/
  status=H5LTfind_dataset(file_id, Name);
  if (status==0)
   { ErrMsg=vstrdup("file %s does not contain a dataset named %s",FileName,Name); 
     return; 
   };

  int Rank;
  status=H5LTget_dataset_ndims(file_id, Name, &Rank);
  if (status<0 || Rank!=1)
   { ErrMsg=vstrdup("file %s: dataset %s is not a vector",FileName,Name); 
     return; 
   };

  /*--------------------------------------------------------------*/
  /*- read dataset dimensions ------------------------------------*/
  /*--------------------------------------------------------------*/
  H5T_class_t class_id;
  size_t type_size;
  hsize_t dims[1];
  status=H5LTget_dataset_info(file_id, Name, dims, &class_id, &type_size);
  if (status<0)
   { ErrMsg=vstrdup("file %s: dataset %s has incorrect type",FileName,Name); 
     return; 
   };
 
  /*--------------------------------------------------------------*/
  /*- try to read attributes from file. note that these will only-*/
  /*- be present if the HDF5 file was created by the libhmat     -*/
  /*- export routines; otherwise, we suppose the file was created-*/
  /*- by another application and assume it describes a real-valued*/
  /*- vector.                                                     */
  /*--------------------------------------------------------------*/
  status=H5LTget_attribute_int(file_id, Name, "RealComplex", &RealComplex);
  if (status<0)
   RealComplex=LHM_REAL; 

  /*--------------------------------------------------------------*/
  /*- now switch off to handle the various cases for reading the  */
  /*- data in                                                     */
  /*--------------------------------------------------------------*/
  if (RealComplex==LHM_REAL)
   { 
     InitHVector(dims[0], LHM_REAL);
     H5LTread_dataset_double(file_id, Name, DV);
   }
  else if (RealComplex==LHM_COMPLEX)
   { 
     InitHVector(dims[0]/2, LHM_COMPLEX);
     H5LTread_dataset_double(file_id, Name, (double *)ZV);
   }

  /*--------------------------------------------------------------*/
  /*- close up the data file -------------------------------------*/
  /*--------------------------------------------------------------*/
  H5Fclose(file_id);
  
}  

#else // HAVE_HDF5

void WarnNoHDF5()
 { fprintf(stderr,"**\n");
   fprintf(stderr,"** warning: compiled without HDF5 support \n");
   fprintf(stderr,"**          (skipping HDF5 operations...) \n");
   fprintf(stderr,"**\n");
 }

void *HMatrix::OpenHDF5Context(const char *format, ...)
{ WarnNoHDF5();
  return 0;
}

void HMatrix::CloseHDF5Context(void *pHC)
{ WarnNoHDF5(); }

void HMatrix::ExportToHDF5(void *pHC, const char *format, ...)
{ WarnNoHDF5(); }

void HMatrix::ExportToHDF5(const char *FileName, const char *format, ...)
{ WarnNoHDF5(); }

void HMatrix::ImportFromHDF5(const char *FileName, const char *Name, bool Init)
{ WarnNoHDF5(); }

void HVector::ExportToHDF5(void *pHC, const char *format, ...)
{ WarnNoHDF5(); }

void HVector::ExportToHDF5(const char *FileName, const char *format, ...)
{ WarnNoHDF5(); }

void HVector::ImportFromHDF5(const char *FileName, const char *Name)
{ WarnNoHDF5(); }


#endif // HAVE_HDF5
