#  <span class=SC>libscuff</span> API Documentation: Utility routines
 
<span class=SC>libscuff</span> provides some simple classes 
called `HMatrix` and `HVector` for working with numerical matrices and vectors.

This page documents a representative *subset* of the 
matrix and vector functions provided; for the full API,  
consult the file `include/scuff-em/libhmat.h`
that comes with the <span class="SC">scuff-em</span>
distribution.

# 1. Creating from scratch
Create real-valued matrices and vectors of known dimensions:


```C++
  HVector *V = new HVector( 13);   // create a 13-component real-valued vector

  HMatrix *M = new HMatrix( 5, 6); // create a 5x6 real-valued matrix
```


Create complex-valued matrices and vectors:


```C++
  HVector *V = new HVector( 13, LHM_COMPLEX); 

  HMatrix *M = new HMatrix( 5, 6, LHM_COMPLEX);
```


For symmetric or Hermitian matrices, you have the option of using
*packed storage*; this roughly halves the RAM needed to store 
your matrix, but has the drawback of slowing down linear algebra 
operations.
 

```C++
  // real symmetric matrix, M_{ij} = M_{ji}
  HMatrix *M = new HMatrix( 5, 5, LHM_REAL, LHM_SYMMETRIC );

  // complex hermitian matrix, M_{ij} = M^*_{ji} 
  HMatrix *M = new HMatrix( 5, 5, LHM_COMPLEX, LHM_HERMITIAN );

  // complex symmetric matrix, M_{ij} = M_{ji} 
  HMatrix *M = new HMatrix( 5, 5, LHM_COMPLEX, LHM_SYMMETRIC );
```

# 2. Importing from text or binary files 

## Text files:

Create a new vector or matrix by importing a list of numbers in an ASCII text file:


```C++
  HVector *V = new HVector("MyVector.dat", LHM_TEXT );
  if (V->ErrMsg) 
   ErrExit(V->ErrMsg);

  HVector *M = new HMatrix("MyMatrix.dat", LHM_TEXT );
  if (M->ErrMsg) 
   ErrExit(M->ErrMsg);
```


Here `MyVector.dat` should be a file with 
a single number per line, while `MyMatrix.dat`
may have multiple numbers per line.
[Complex numbers](scuff-em/reference/scuffEMMisc.shtml#Complex)
are allowed. Blank lines and comments (lines beginning with a 
pound sign `#`) are skipped. 
The dimension of the resulting `HVector`, 
and the number of rows of the resulting `HMatrix`,
will be the number of non-blank non-comment lines in the file, 
while the number of columns of the `HMatrix` will
be the largest number of numbers read from any one line.


Note that, if the file import operation fails,
the constructor returns an `HVector` or
`HMatrix` whose `ErrMsg` field 
points to a nonempty error message. (In this case, all 
other class fields in the object should be assumed to be 
invalid.)


If the operation is successful, then `ErrMsg`
will be `NULL` on return.

### HDF5 files:

Create a new vector or matrix by importing from an HDF5 
binary data file:


```C++
   HMatrix *M1 = new HMatrix("MyFile.hdf5", LHM_HDF5, "M1");
   HMatrix *M2 = new HMatrix("MyFile.hdf5", LHM_HDF5, "M2");
   HVector *V  = new HVector("MyFile.hdf5", LHM_HDF5, "V");

   if (M1->ErrMsg) ErrExit(M1->ErrMsg);
   if (M2->ErrMsg) ErrExit(M2->ErrMsg);
   if (V->ErrMsg)  ErrExit(V->ErrMsg);
```


Note that the third parameter to the constructor here is 
the label of the dataset to within the HDF5 file.

<!---------------------------------------------------------------------->
<!---------------------------------------------------------------------->
<!---------------------------------------------------------------------->
# 3. Exporting to text or binary files

### Text files: 


Write the contents of a vector or matrix to an ASCII text file:

```C++
  M->ExportToText("MyMatrix.dat");
  V->ExportToText("MyVector.dat");
```


The argument to `ExportToText` supports `printf-`like
semantics for inserting numbers, etc. into the file name:

```C++
  int N = 3;
  M->ExportToText("Matrix_%i.dat", N);
```


### HDF5 files: 


There are two calling conventions for exporting matrices and vectors
to HDF5 files. 

If you want each matrix and vector to be exported to a separate 
HDF5 file, you can say simply 


```C++
  M->ExportToHDF5("MyMatrix.hdf5","M");
  V->ExportToHDF5("MyVector.hdf5","V");
```

Alternatively, you can write multiple matrices and vectors to 
a single HDF5 file.


```C++
  void *HC = HMatrix::OpenHDF5Context("MyFile.hdf5");
  M->ExportToHDF5(HC, "M");
  V->ExportToHDF5(HC, "V");
  HMatrix::CloseHDF5Context(HC);
```

# 4. Simple manipulations
 
** Note: All indices are zero-based.**
 
Get or set single entries:
 
```C++
   HMatrix *M = new HMatrix(M, N, LHM_COMPLEX);
   HVector *V = new HVector(N);

   M->SetEntry( 3, 4, 5.6 );
   M->SetEntry( 3, 4, cdouble(5.6+7.8) );

   V->SetEntry( 7, 8.009 );

   double D;
   cdouble Z;

   Z = M->GetEntry( 3, 4 );
   D = M->GetEntryD( 3, 4 ); // discard any imaginary part 

   Z = V->GetEntry(7);
   D = V->GetEntryD(8);
```

Get or set blocks of entries in an `HMatrix:`
```

  HMatrix *M = new HMatrix(10, 5, LHM_COMPLEX);

  cdouble Col3[10];
  M->GetEntriesD(":",3,Col3);

  cdouble M2_14[4];
  M->GetEntriesD(2,"14",M2_14);

```
 
Augment individual entries:

 
```C++
   M->AddEntry( 3, 4, 5.6 );  // the (3,4) entry gets increased by 5.6
   V->AddEntry( 0, cdouble(0.0, 2.3) ); // the 0th entry gets increased by 2.3i
```
 
Replace a matrix with its conjugate or non-conjugate tranpose 
(no distinction for real-valued matrices):

 
```C++
   M->Adjoint();    // conjugate transpose
   M->Transpose();  // non-conjugate transpose
```

# 5. Numerical linear algebra: <span class="SC">LAPACK/BLAS</span>
 
Multiply two matrices:

 
```C++
  int P, Q, R;
  A=new HMatrix(P, Q);
  B=new HMatrix(Q, R);
  C=new HMatrix(P, R); 
  ... 
  A->Multiply(B, C); // set C = A*B 
  A->Multiply(C, B, "--TransA T"); // set B = transpose(A) * C
```

 
Replace a matrix with its LU factorization:

 
```C++
   M->LUFactorize();
```

 
Solve a single linear system using an LU-factorized matrix:
 
```C++
   M=new HMatrix(N, N);
   // insert code to fill in M 
   M->LUFactorize();

   V=new HVector(N);
   // insert code to fill in V 
   M->LUSolve(V);   // replaces V with M^{-1} * V 
```

 
Of course, we only need to call `LUFactorize` once for
a given matrix, after which we can make any number of calls to 
`LUSolve()` with different vectors. 

 
We can also solve multiple simultaneous systems by passing
an `HMatrix` instead of an `HVector` to  
`LUSolve:`
 
```C++
   M=new HMatrix(N, N);
   // insert code to fill in M 
   M->LUFactorize();

   R=new HMatrix(N, M); 
   // insert code to fill in the M columns of R 
   M->LUSolve(R);   // replaces R with M^{-1} * R
```

{!Links.md!}
