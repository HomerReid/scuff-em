// -*- C++ -*-


%{
#define SWIG_FILE_WITH_INIT // to build as Python module
%}

%include "std_complex.i" // conversion for std::complex
%apply std::complex<double> { cdouble };

%include "std_vector.i"
namespace std {
   %template(IntVector) vector<int>;
   %template(DoubleVector) vector<double>;
}

%include "numpy.i" // numpy array conversions

// some additional backward compatibility declarations for supporting numpy < 1.7.0
%{
#if NPY_API_VERSION < 0x00000007
#define NPY_ARRAY_C_CONTIGUOUS NPY_C_CONTIGUOUS
#define NPY_ARRAY_ALIGNED  NPY_ALIGNED
#endif
%}

%init %{
  import_array();
%}
%fragment("NumPy_Fragments");
%fragment("NumPy_Macros");
%numpy_typemaps(double, NPY_DOUBLE, size_t)
%numpy_typemaps(cdouble, NPY_CDOUBLE, size_t)
%apply double IN_ARRAY1[ANY] { const double [3] };
%apply cdouble IN_ARRAY1[ANY] { const cdouble [3] };
%apply double INPLACE_ARRAY1[ANY] { double [6], double [3] };
%apply cdouble INPLACE_ARRAY1[ANY] { cdouble [6], cdouble [3] };
