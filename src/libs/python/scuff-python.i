// -*- C++ -*-


%{
#define SWIG_FILE_WITH_INIT // to build as Python module
%}

%include "std_complex.i" // conversion for std::complex
%include "numpy.i" // numpy array conversions
%init %{
  import_array();
%}
%fragment("NumPy_Fragments");
%fragment("NumPy_Macros");
%numpy_typemaps(double, NPY_DOUBLE, size_t)
%numpy_typemaps(cdouble, NPY_CDOUBLE, size_t)

//////////////////////////////////////////////////////////////////////////////
// Wrapper for EHFuncType

// wrap Python function EH = f(R)
%{
static void EHFunc_python(double *R, void *f, cdouble *EH) {
  npy_intp sz3 = 3, stride1 = sizeof(double);
  PyObject *Rpy = PyArray_New(&PyArray_Type, 1, &sz3, NPY_DOUBLE, &stride1,
                              R, // not NPY_WRITEABLE, effectively const
                              0, NPY_C_CONTIGUOUS | NPY_ALIGNED, NULL);
  PyObject *arglist = Py_BuildValue("O", Rpy);
  PyObject *result = PyEval_CallObject((PyObject *) f, arglist);
  PyArrayObject *EHpy;
  int is_new;
  Py_DECREF(arglist);
  Py_DECREF(Rpy);
  if (PyErr_Occurred()) {
    ErrExit("python exception handling not supported in libscuff");
  }
  else if ((EHpy = obj_to_array_contiguous_allow_conversion(result,
                                                            NPY_CDOUBLE,
                                                            &is_new))) {
    if (array_numdims(EHpy) != 1 || array_size(EHpy,0) != 0)
      ErrExit("6-component EH field must be returned to libscuff");
    memcpy(EH, array_data(EHpy), sizeof(cdouble) * 6);
    if (is_new) Py_DECREF((PyObject *) EHpy);
  }
  else {
    ErrExit("invalid EH field returned to libscuff");
  }
  Py_XDECREF(result);
}
%}

%typemap(in)(scuff::EHFuncType EHFunc, void *EHFuncUD) {
  $1 = EHFunc_python;
  Py_INCREF($input);
  $2 = (void *) $input;
}
%typemap(freearg)(scuff::EHFuncType EHFunc, void *EHFuncUD) {
  Py_XDECREF((PyObject *) $2);
}
%typecheck(SWIG_TYPECHECK_POINTER)(scuff::EHFuncType EHFunc, void *EHFuncUD) {
  $1 = PyCallable_Check($input);
}

//////////////////////////////////////////////////////////////////////////////
// HMatrix and HVector typemaps: these are mapped directly onto
// NumPy matrices, ideally with no data copying.
//
// The exception is the case of a packed symmetric/Hermitian HMatrix,
// for which there is no corresponding NumPy array type.  In this case
// we return an opaque HMatrix*, and correspondingly all of the routines
// that take HMatrix* arguments accept either a NumPy array (unpacked storage)
// or an opaque HMatrix*.
//
// We also return opaque HVector*/HMatrix* pointers if the user explicitly
// constructs an HVector/HMatrix object in Python.  This is both to honor
// the user's request and also to correctly work with the Python wrapper
// classes.  (SWIG doesn't provide an elegant way of using a different
// typemap in "new" constructors, but we can do so by checking $symname
// at runtime.)

%typemap(out)(HVector *) {
  if (!strcmp("$symname", "new_HVector")) {
    $result = SWIG_NewPointerObj((void*) $1, SWIGTYPE_p_HVector, 0);
  }
  else {
    npy_intp sz = npy_intp($1->N);
    $result = (PyObject*) PyArray_SimpleNewFromData
      (1, &sz, $1->RealComplex == LHM_COMPLEX ? NPY_CDOUBLE : NPY_DOUBLE,
       $1->RealComplex == LHM_COMPLEX ? (void*) $1->ZV : (void*) $1->DV);
    $1->ownsV = false; // prevent deallocation
    delete $1;
  }
}

%{
  static bool is_HVector_array(PyObject *o) {
    return is_array(o) && array_numdims(o) == 1
      && (array_is_contiguous(o) || array_is_fortran(o))
      && array_is_native(o)
      && PyArray_ISALIGNED(o)
      && (array_type(o) == NPY_DOUBLE || array_type(o) == NPY_CDOUBLE);
  }
%}

%typemap(in)(HVector *) (PyArrayObject* conv=NULL) {
  if (SWIG_IsOK(SWIG_ConvertPtr($input, (void**)&$1, SWIGTYPE_p_HVector, 0))) {
    // passed opaque HVector*
  }
  else if (is_HVector_array($input)) {
    // create HVector without making a copy of the data
    $1 = new HVector(array_size($input, 0), 
		     array_type($input) == NPY_DOUBLE ? LHM_REAL : LHM_COMPLEX,
		     array_data($input));
  }
  else {
    int is_new_object;
    conv = obj_to_array_contiguous_allow_conversion($input, NPY_DOUBLE, 
						     &is_new_object);
    if (!conv)
      conv = obj_to_array_contiguous_allow_conversion($input, NPY_CDOUBLE, 
						       &is_new_object);
    if (!conv)
      SWIG_exception_fail(SWIG_TypeError,
			  "in method '" "$symname" "', argument " "$argnum"
			  " of type '" "$type" "': expecting array");
    if (!is_new_object) Py_INCREF(conv); // prevent deallocation
    $1 = new HVector(array_size(conv, 0),
		     array_type(conv) == NPY_DOUBLE ? LHM_REAL : LHM_COMPLEX,
		     array_data(conv));
  }
}
%typemap(freearg)(HVector *) {
  if ($1 && !$1->ownsV) // only deallocate if wrapper around numpy array
    delete $1;
  Py_XDECREF(conv$argnum); // deallocate (if new object was created)
}
%typecheck(SWIG_TYPECHECK_POINTER)(HVector *) {
  $1 = SWIG_IsOK(SWIG_ConvertPtr($input, NULL, SWIGTYPE_p_HVector, 0))
    || is_HVector_array($input) || PySequence_Check($input);
}

%typemap(out)(HMatrix *) {
  if ($1->StorageType != LHM_NORMAL ||
      !strcmp("$symname", "new_HMatrix")) { // return opaque HMatrix*
    $result = SWIG_NewPointerObj((void*) $1, SWIGTYPE_p_HMatrix, 0);
  }
  else { // return numpy.array
    npy_intp dims[2];
    dims[0] = npy_intp($1->NR);
    dims[1] = npy_intp($1->NR);
    $result = (PyObject*) PyArray_New // TODO: how to return numpy.matrix?
      (&PyArray_Type, 2, dims,
       $1->RealComplex == LHM_COMPLEX ? NPY_CDOUBLE : NPY_DOUBLE, NULL,
       $1->RealComplex == LHM_COMPLEX ? (void*) $1->ZM : (void*) $1->DM, 0,
       NPY_FARRAY, NULL);
    $1->ownsM = false; // prevent deallocation
    delete $1;
  }
}

%{
  static bool is_HMatrix_array(PyObject *o) {
    return is_array(o) && array_numdims(o) == 2
      && (array_is_contiguous(o) || array_is_fortran(o))
      && array_is_native(o)
      && PyArray_ISALIGNED(o)
      && (array_type(o) == NPY_DOUBLE || array_type(o) == NPY_CDOUBLE);
  }
%}

%typemap(in)(HMatrix *) {
  if (SWIG_IsOK(SWIG_ConvertPtr($input, (void**)&$1, SWIGTYPE_p_HMatrix, 0))) {
    // passed opaque HMatrix*
  }
  else if (is_HMatrix_array($input)) {
    // create HMatrix without making a copy of the data
    // ... note transposition of dimensions for C-ordered array
    //     (FIXME: transpose the data in this case?)
    $1 = new HMatrix(array_size($input, !array_is_fortran($input)),
		     array_size($input, array_is_fortran($input)),
		     array_type($input) == NPY_DOUBLE ? LHM_REAL : LHM_COMPLEX,
		     LHM_NORMAL,
		     array_data($input));
  }
  else {
    // TODO: allow conversions (with a copy) from other types?
    SWIG_exception_fail(SWIG_TypeError,
			"in method '" "$symname" "', argument " "$argnum"
			" of type '" "$type" "': expecting numpy array");
  }
}
%typemap(freearg)(HMatrix *) {
  if ($1 && !$1->ownsM) // only deallocate if wrapper around numpy array
    delete $1; 
}
%typecheck(SWIG_TYPECHECK_POINTER)(HMatrix *) {
  $1 = SWIG_IsOK(SWIG_ConvertPtr($input, NULL, SWIGTYPE_p_HMatrix, 0))
    || is_HMatrix_array($input);
}
