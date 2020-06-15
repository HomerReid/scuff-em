/* Copyright (C) 2005-2011 M. T. Homer Reid};
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
 * libhmat.h   -- header file for libhmat
 *
 * homer reid  -- 3/2007 -- 11/2009
 */

#ifndef LIBHMAT_H
#define LIBHMAT_H

#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <cmath>
#include <complex>
#include <vector>

/***************************************************************/
/***************************************************************/
/***************************************************************/
#ifndef cdouble
typedef std::complex<double> cdouble;
#endif 
#ifndef ivec
typedef std::vector <int> ivec;
#endif

/***************************************************************/
/***************************************************************/
/***************************************************************/
// values for the RealComplex argument to functions
#define LHM_REAL    0
#define LHM_COMPLEX 1

// values for the StorageType argument to functions
#define LHM_NORMAL    0    /* normal storage, no symmetry */
#define LHM_SYMMETRIC 1    /* non-hermitian symmetric, M_{ij} = M_{ji} */
#define LHM_HERMITIAN 2    /* hermitian symmetric, M_{ij} = M_{ji}^* */

// values for the FileType argument to functions
#define LHM_HDF5    0
#define LHM_TEXT    1
#define LHM_AUTO    -1 // infer from filename
int LHM_AUTO_FileType(const char *FileName);

// values for the How argument to Concat 
#define LHM_HORIZONTAL 0
#define LHM_VERTICAL 1

/***************************************************************/
/* HVector class definition ************************************/
/***************************************************************/
class HVector
 { 
  public:  

   /* constructor entry points to create vectors of known length */
   /* and data type (optionally with user-specified payload)     */
   HVector(int N, int RealComplex = LHM_REAL, void *data = NULL);
   HVector(int N, double *data);
   HVector(int N, cdouble *data);
   void InitHVector(int N, int RealComplex = LHM_REAL, void *data = NULL);

   /* constructors that attempt to read vector in from a data file */
   HVector(const char *FileName, int FileType = LHM_AUTO, const char *Options = "");
   void ReadFromFile(const char *FileName, int FileType = LHM_AUTO, const char *Options = "");

   /* HDF5 file IO */
   void ImportFromHDF5(const char *FileName, const char *Name);
   void ExportToHDF5(const char *FileName, const char *format, ...);
   void ExportToHDF5(void *pHC, const char *format, ...);

   /* MATLAB file IO */
   void ExportToMATLAB(void *pCC, const char *format, ... );

   /* ASCII text file IO */
   void ImportFromText(const char *FileName, const char *Options = "");

   void ExportToText(const char *FileName, const char *Options = "");

   /* copy constructor */
   HVector(HVector *V);

   /* destructor */
   ~HVector();

   /* get the value of an entry */
   cdouble GetEntry(int n);
   double GetEntryD(int n);

   /* copying another HVector of the same size */
   void Copy(HVector *V);
   HVector *Copy() { return new HVector(this); }

   /* set the value of an entry */
   void SetEntry(int n, double Entry);
   void SetEntry(int n, cdouble Entry);

   /* augment the value of an entry */
   void AddEntry(int n, double Entry);
   void AddEntry(int n, cdouble Entry);

   /* scale the entire vector by a scalar multiple */
   void Scale(double Alpha);
   void Scale(cdouble Alpha);

   cdouble Dot(HVector *B);   /* hermitian dot product    */
   cdouble DotU(HVector *B);  /* non-conjugated dot product */
   double  DotD(HVector *B);  /* real part of Dot()       */

   void Zero();

 // private:
   int N;
   int RealComplex;

   // pointers to the actual data storage. only one of these is 
   // used in a given instance so if i wanted to save 8 bytes i 
   // could put them into a union
   double *DV;
   cdouble *ZV;

   // flag to indicate whether we "own" the DV/ZV data & should free it
   bool ownsV;

   // if this field is nonzero on return from one of the 
   // constructors, it means the constructor failed and  
   // ErrMsg explains why 
   char *ErrMsg;
 };

/***************************************************************/
/* a couple of non-class methods that return newly allocated   */
/* HVectors                                                    */
/***************************************************************/
HVector *LinSpace(double Min, double Max, int Num);
HVector *LogSpace(double Min, double Max, int Num);
HVector *Concat(HVector *V1, HVector *V2);
HVector *GetOmegaList(char *OmegaFile,
                      cdouble *OmegaVals,  int nOmegaVals,
                      char *LambdaFile,
                      cdouble *LambdaVals, int nLambdaVals);

/***************************************************************/
/* HMatrix class definition ************************************/
/***************************************************************/
class SMatrix; // forward declaration used by a couple of HMatrix routines
class HMatrix
 { 
  public:  

   /* constructors that initialize the matrix with known sizes and */
   /* storage parameters                                           */
   HMatrix(int NRows, int NCols, int RealComplex = LHM_REAL,
	   int StorageType = LHM_NORMAL, void *data = NULL);
   HMatrix(int NRows, int NCols, int pRealComplex, void *data);
   HMatrix(int NRows, int NCols, cdouble *data);
   HMatrix(int NRows, int NCols, double *data);
   HMatrix(HMatrix *M, bool takedatandownership=false); // copy constructor
   void InitHMatrix(int NRows, int NCols, int RealComplex = LHM_REAL,
		    int StorageType = LHM_NORMAL, void *data = NULL);

   /* constructors that attempt to read matrix in from a data file */
   HMatrix(const char *FileName, int FileType = LHM_AUTO, const char *Options = "");
   void ReadFromFile(const char *FileName, int FileType = LHM_AUTO, const char *Options = "");

   /* constructor that converts an SMatrix to an HMatrix */
   HMatrix(SMatrix *S);

   /* HDF5 file IO */
   void ImportFromHDF5(const char *FileName, const char *MatrixName, bool Init=true);
   void ExportToHDF5(const char *FileName, const char *format, ...);
   void ExportToHDF5(void *pHC, const char *format, ...);

   /* MATLAB file IO */
   void ExportToMATLAB(void *pCC, const char *format, ...);

   /* ASCII text file IO */
   void ImportFromText(const char *FileName, const char *Options = "");
   void ExportToText(const char *FileName, const char *Options = "");

   /* destructor */
   ~HMatrix(); 
   /* get the value of an entry */
   cdouble GetEntry(size_t nr, size_t nc);
   double GetEntryD(size_t nr, size_t nc);

   /* routines for setting, augmenting, scaling matrix entries */
   void SetEntry(size_t nr, size_t nc, double Entry);
   void SetEntry(size_t nr, size_t nc, cdouble Entry);

   void AddEntry(size_t nr, size_t nc, double Entry);
   void AddEntry(size_t nr, size_t nc, cdouble Entry);

   void ScaleEntry(size_t nr, size_t nc, cdouble ScaleFactor);

   void *GetColumnPointer(size_t nc);

   /* matlab-style fetching or setting of swaths of a matrix */
   double *GetEntriesD(const char *RowString, int Col, double *Entries=0);
   cdouble *GetEntries(const char *RowString, int Col, cdouble *Entries=0);
   double *GetEntriesD(int Row, const char *ColString, double *Entries=0);
   cdouble *GetEntries(int Row, const char *ColString, cdouble *Entries=0);
   HMatrix *ExtractEntries(const char *RowColString, HMatrix *B=0);
   HMatrix *DoGetEntries(int RowStart, int RowStop, int RowInc, int RowLen,
                         int ColStart, int ColStop, int ColInc, int ColLen,
                         HMatrix *B, int RowOffset, int ColOffset);

   void SetEntriesD(const char *RowString, int Col, 
                    double *Entries, double Entry=0.0);
   void SetEntries(const char *RowString, int Col, 
                   cdouble *Entries, cdouble Entry=0.0);
   void SetEntriesD(int Row, const char *ColString, 
                    double *Entries, double Entry=0.0);
   void SetEntries(int Row, const char *ColString, 
                   cdouble *Entries, cdouble Entry=0.0);
   void SetEntriesD(const char *RowString, int Col, double Entry)
    { SetEntriesD(RowString, Col, 0, Entry); }
   void SetEntries(const char *RowString, int Col, cdouble Entry)
    { SetEntries(RowString, Col, 0, Entry); }
   void SetEntriesD(int Row, const char *ColString, double Entry)
    { SetEntriesD(Row, ColString, 0, Entry); }
   void SetEntries(int Row, const char *ColString, cdouble Entry)
    { SetEntries(Row, ColString, 0, Entry); }
   void DoSetEntries(int RowStart, int RowStop, int RowInc, int RowLen,
                     int ColStart, int ColStop, int ColInc, int ColLen,
                     HMatrix *B, cdouble Entry=0.0);

   /* routine for copying another HMatrix of the same size*/
   void Copy(HMatrix *M);
   HMatrix *Copy() { return new HMatrix(this); }

   /* get the trace */
   cdouble GetTrace();
   double GetTraceD();

   /* scale the entire matrix by a scalar multiple */
   void Scale(double Alpha);
   void Scale(cdouble Alpha);

   void Zero();
   void ZeroBlock(int RowOffset, int NumRows, int ColOffset, int NumCols);
   void Adjoint();   // conjugate transpose 
   void Transpose(); // non-conjugate transpose

   // routines for inserting a smaller HMatrix into a larger one
   void InsertBlock(HMatrix *B, int RowOffset, int ColOffset, double ScaleFactor=1.0);
   void InsertBlockAdjoint(HMatrix *B, int RowOffset, int ColOffset, double ScaleFactor=1.0);
   void InsertBlockTranspose(HMatrix *B, int RowOffset, int ColOffset, double ScaleFactor=1.0);
   // insert a subblock of B
   void InsertBlock(HMatrix *B, int RowOffset, int ColOffset,
                    int NRB, int NCB, int BRowOffset, int BColOffset, double ScaleFactor=1.0);

   // like InsertBlock, but addition rather than replacement
   void AddBlock(HMatrix *B, int RowOffset, int ColOffset, cdouble ScaleFactor=1.0);
   void AddBlockAdjoint(HMatrix *B, int RowOffset, int ColOffset, cdouble ScaleFactor=1.0);
   void AddBlock(SMatrix *B, int RowOffset, int ColOffset, cdouble ScaleFactor=1.0);
   void Add(HMatrix *B, cdouble ScaleFactor=1.0);
   void Add(SMatrix *B, cdouble ScaleFactor=1.0);

   // sort of the inverse of InsertBlock
   void ExtractBlock(int RowOffset, int ColOffset, HMatrix *B);

   /* sort the rows by the values of one or more columns */
   void Sort(std::vector<int> SortColumns,
             std::vector<char> SortType);
   void Sort(std::vector<int> SortColumns);
   void Sort(int SortColumn, const char *Options=0);

   // compute vector-matrix-vector product X'*this*Y
   cdouble BilinearProduct(HVector *X, HVector *Y=0);
   double BilinearProductD(HVector *X, HVector *Y=0);

   /*--------------------------------------------------------------*/
   /*- the following routines are wrappers around LAPACK or BLAS   */
   /*- calls, possibly with some pre- or post-processing           */
   /*--------------------------------------------------------------*/
   // matrix-matrix multiply, this * B -> C (xgemm)
   // options example: "--transA T" or "--transA T --transB C"
   //  (to use non-hermitian transpose of A and hermitian adjoint of B)
   void Multiply(HMatrix *B, HMatrix *C, const char *Options=0); 

   /* matrix-matrix multiply with only the diagonal elements */
   /* of the product matrix computed (using ddot or zdotu)   */
   void GetMatrixProductDiagonal(HMatrix *B, HVector *DAB);

   /* matrix-vector multiplication, this*X = Y (xgemv))      */
   /* Trans='C' or 'T' for hermitian / non-hermitian adjoint */
   void Apply(HVector *X, HVector *Y, char Trans=0);
   
   /* routines for LU-factorizing, solving, inverting */
   /* (xgetrf, xgetrs, xgetri) */
   int LUFactorize();
   int LUSolve(HVector *X);
   int LUSolve(HMatrix *X);
   int LUSolve(HMatrix *X, int nrhs);
   int LUSolve(HMatrix *X, char Trans);
   int LUSolve(HMatrix *X, char Trans, int nrhs);
   int LUInvert();

   /* routines for cholesky-factorizing, solving, inverting */
   /* (xpotrf, xpotrs, xpotri) */
   int CholFactorize();
   int CholSolve(HVector *X);
   int CholSolve(HMatrix *X);
   int CholSolve(HMatrix *X, int nrhs);

   /* routine for qr-factorizing (xgeqrf) */
   //int QR(HMatrix **R);
   int QR(HMatrix **Q, HMatrix **R);

   /* routine for eigenvalues and optionally vectors of */
   /* symmetric/hermitian matrices (dsyevr / zheevr)    */
   HVector *Eig(HVector *Lambda=0, HMatrix *U=0);

   /* routine for eigenvalues and optionally vectors of */
   /* non-symmetric matrices (dgeevx / zgeevx )         */
   HVector *NSEig(HVector *Lambda=0, HMatrix *U=0);

   /* singular-value decomposition */
   HVector *SVD(HVector *Sigma=0, HMatrix *U=0, HMatrix *VT=0);

   /* matrix norm*/
   double GetNorm(bool UseInfinityNorm=false);

   /* routine for estimating reciprocal of condition    */
   /* number (assumes LUFactorize() has been called)    */
   double GetRCond(double ANorm, bool UseInfinityNorm=false);

   size_t inline NumEntries() 
    { return StorageType==LHM_NORMAL ? ((size_t)NR)*NC : ((size_t)NR)*(NR+1)/2; }

 // private:
   int NR, NC;
   int RealComplex;
   int StorageType;
   int *ipiv;

   // pointers to the actual data storage. only one of these is 
   // used in a given instance so if i wanted to save 8 bytes i 
   // could put them into a union
   double *DM;
   cdouble *ZM;

   // internally-stored workspaces that are automatically
   // allocated as necessary for various LAPACK routines
   int lwork; // size currently allocated for work in doubles or cdoubles
   void *work;
   int liwork; // size currently allocated for iwork in ints
   int *iwork;

   // flag to indicate whether we "own" the DM/ZM data & should free it
   bool ownsM; 
   // if this field is nonzero on return from one of the 
   // constructors, it means the constructor failed and  
   // ErrMsg explains why 
   char *ErrMsg;

   // static class methods for opening and closing HDF5 contexts
   static void *OpenHDF5Context(const char *format, ... );
   static void CloseHDF5Context(void *pHC);

   // static class methods for opening and closing MATLAB contexts
   static void *OpenMATLABContext(const char *format, ... );
   static void CloseMATLABContext(void *pCC);

   // static class variable that, if true, yields a call to
   // ErrExit() whenever the HMatrix(DataFile) or HVector(DataFile)
   // constructor fails. 
   // This is set to true by default, but you can set it to
   // false if you like; in this case, failed constructors
   // will return an HMatrix/HVector with no data and with a 
   // non-NULL value of its ErrMsg field.
   static bool AbortOnIOError;

 };

// get the real part matrix
HMatrix *GetReal(HMatrix *A);

// make an unpacked copy of a symmetric/Hermitian matrix
HMatrix *CopyHMatrixUnpacked(HMatrix *Mpacked);

// contatenate A and B to create a new HMatrix
HMatrix *Concat(HMatrix *A, HMatrix *B, int How=LHM_HORIZONTAL);

// make a sorted copy of a matrix
HMatrix *SortHMatrix(HMatrix *M, std::vector<int> SortColumns, 
                                 std::vector<char> SortType);

// return a matrix with given properties, using an optional
// user-supplied candidate if it matches, or allocating a new one
// if not (in which case the candidate is deallocated if non-null)
HMatrix *CheckHMatrix(HMatrix *M, int NR, int NC, int RealComplex, const char *Context=0);
HMatrix *CheckHMatrix(HMatrix *M, HMatrix *M0, const char *Context=0);

/***************************************************************/
/* SMatrix class definition ************************************/
/***************************************************************/
class SMatrix
 { 
  public:  

    // constructor / destructor 
    // NR x NC matrix, nnz_row = approximate mean #nonzeros/row (0 is ok)
    SMatrix(int NR, int NC, int RealComplex = LHM_REAL);
    ~SMatrix();

    SMatrix(const char *FileName, char *MatrixName);

    void ExportToHDF5(void *pHC, const char *format, ...);
    void ExportToHDF5(const char *FileName, const char *format, ...);
    void ImportFromHDF5(const char *FileName, const char *Name);
    void ReadFromFile(const char *FileName, int FileType, const char *Options);

    // return the number of nonzero entries in a given row
    int GetNNZ(int nr);

    bool isNZ(int nr, int nc);

    void Zero();

    // return pointers to arrays of column indices and nonzero 
    // entries for a given row
    int GetRow(int nr, int **pColIndices, void **NZEntries);

    // like GetRow() but actually fills in user-allocated arrays
    int GetRowEntries(int nr, int *ColIndices, void *NZEntries);

    cdouble GetEntry(int nr, int nc);
    double GetEntryD(int nr, int nc) { return real(GetEntry(nr,nc)); }

    // MX = S * X
    cdouble Apply(HVector *X, HVector *MX);
    double ApplyD(HVector *X, HVector *MX);
    HVector *Apply(HVector *X); // return new MX
    void Apply(HMatrix *X, HMatrix *MX);
    HMatrix *Apply(HMatrix *X); // return new MX

    // compute vector-matrix-vector product X'*this*Y
    cdouble BilinearProduct(HVector *X, HVector *Y=0);
    double BilinearProductD(HVector *X, HVector *Y=0);

    // matrix assembly:
    //     BeginAssembly(estimate of # nonzero entries)
    //        ... SetEntry/AddEntry is fastest in order of rows!
    //     EndAssembly()
    void BeginAssembly(int est_nnz = 0 /* default: allocate on fly */);
    void EndAssembly();
    void SetEntry(int nr, int nc, cdouble Entry);
    void AddEntry(int nr, int nc, cdouble Entry,
		  bool compress = true);
    // (set compress = false to be faster (much faster for dense rows),
    //  but wastes matrix space if called multiple times for same nr & nc)

 // private:
    int NR, NC;       // number of rows, cols
    int nnz, nnz_alloc;      // total # entries, and total allocated
    int RealComplex;
    int cur_nr; // last-added row, during assembly, or NR+1 otherwise

    // data storage is CSR (compressed sparse row) format:
    //   for nr-th row, we get i = RowStart[nr], iend = RowStart[nr+1],
    //   and entries are M[i..iend-1] for col indices ColIndices[i..iend-1]
    int *RowStart; // if non-NULL, array of length NR+1
    double *DM; cdouble *ZM; // if non-NULL, array of length nnz_alloc
    int *ColIndices; // if non-NULL, array of length nnz_alloc

    char *ErrMsg;

 private:
    void Reallocate(int nnz_alloc); // internal allocation function
    int MakeEntry(int nr, int nc, bool force_new); // internal function to allocate entries
 };

#endif
