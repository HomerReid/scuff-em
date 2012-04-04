/*
   FieldGrid.cc -- classes and functions so that we can compute
                   arbitrary functions of the fields on an arbitrary 2d grid 
		   for convenient visualization and other analyses

   SGJ, 3/2012
*/

#include "libscuff.h"
#include "cmatheval.h"

#include <string.h>

static const double MU0 = 4e-7 * 3.14159265358979323846; // magnetic constant
static const double C0 = 299792458.0;                  // vacuum speed of light
static const double EPS0 = 1.0 / (C0*C0*MU0);        // electric constant

namespace scuff {

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

// variables for use in field expressions; note that that if you change this
// list you mst change some of the hard-coded indices (10,20,24) below.
static const char *PFF_names[24] = {
  "x", "y", "z", 
  "dAx", "dAy", "dAz", "dA", "nx", "ny", "nz",
  "Ex", "Ey", "Ez", "Hx", "Hy", "Hz",
  "Eps","Mu","eps", "mu",

  "eps0","mu0","c","Z0" // constants
};
  
ParsedFieldFunc::ParsedFieldFunc(const char *ExprString) {
  if (!ExprString) ErrExit("invalid NULL string passed to ParsedFieldFunc");
  if (!(Expr = cevaluator_create(const_cast<char*>(ExprString))))
    ErrExit("error parsing FieldFunc expression: %s", ExprString);

  // check to make sure that only PFF_names are used as variables
  char **names;
  int count;
  cevaluator_get_variables(Expr, &names, &count);
  for (int i = 0; i < count; ++i) {
    int j;
    for (j = 0; j < 24 && strcmp(names[i], PFF_names[j]); ++j)
      ;
    if (j == 24)
      ErrExit("unrecognized variable %s in FieldFunc expression %s",
	      names[i], ExprString);
  }

  // set constants
  cevaluator_set_var(Expr, "eps0", EPS0);
  cevaluator_set_var(Expr, "mu0", MU0);
  cevaluator_set_var(Expr, "c", C0);
  cevaluator_set_var(Expr, "Z0", ZVAC);

  // we have to tell cmatheval which vars may be complex for is_real to work
  // (variables 10..19)
  cdouble cvals[10];
  for (int i = 0; i < 10; ++i) cvals[i] = cdouble(0,1);

  isReal = cevaluator_is_real(Expr, 10, const_cast<char**>(&PFF_names[10]), 
			      cvals);

  // Make first 20 variables thread-safe (using fixed index into vals array)
  for (int i = 0; i < 20; ++i)
    cevaluator_set_var_index(Expr, const_cast<char*>(PFF_names[i]), i);
}

ParsedFieldFunc::~ParsedFieldFunc() {
  cevaluator_destroy(Expr);
}

char *ParsedFieldFunc::String() const {
  return strdup(cevaluator_get_string(Expr));
}

cdouble ParsedFieldFunc::Eval(const double X[3],
			      const double dA[3],
			      const cdouble EH[3],
			      cdouble Eps, cdouble Mu) const {
  double dAi;
  cdouble vals[20]; // values corresponding to PFF_names variables
  vals[0] = X[0];
  vals[1] = X[1];
  vals[2] = X[2];
  vals[3] = dA[0];
  vals[4] = dA[1];
  vals[5] = dA[2];
  vals[6] = (dAi = VecNorm(dA)); dAi = 1.0 / dAi;
  vals[7] = (dA[0] * dAi);
  vals[8] = (dA[1] * dAi);
  vals[9] = (dA[2] * dAi);
  vals[10] = EH[0];
  vals[11] = EH[1];
  vals[12] = EH[2];
  vals[13] = EH[3];
  vals[14] = EH[4];
  vals[15] = EH[5];
  vals[16] = Eps;
  vals[17] = Mu;
  vals[18] = Eps*EPS0;
  vals[19] = Mu*MU0;
  // note: thread-safe since these 20 variables are "indexed" vars
  //       and hence the symbol table in Expr is not modified
  return cevaluator_evaluate(Expr, 20, const_cast<char**>(PFF_names), vals);
}

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

PlaneGrid::PlaneGrid(int n1, int n2, const double c0[3],
		     const double s1[3], const double s2[3]) {
  N1 = n1; N2 = n2;
  VecPlusEquals(VecPlusEquals(VecCopy(c0, X0), -0.5, s1), -0.5, s2);
  VecScale(VecCopy(s1, S1), 1.0 / N1);
  VecScale(VecCopy(s2, S2), 1.0 / N2);
  VecPlusEquals(VecPlusEquals(X0, 0.5, S1), 0.5, S2);
  VecCross(S1, S2, dA0);
}

PlaneGrid::PlaneGrid(int n1, int n2, const double c0[3],
		     double s1, double s2, CartesianDirection normal) {
  N1 = n1; N2 = n2;
  VecCopy(c0, X0);
  VecZero(S1); VecZero(S2);
  X0[(normal + 1) % 3] -= 0.5 * s1;
  X0[(normal + 2) % 3] -= 0.5 * s2;
  S1[(normal + 1) % 3] = s1 / N1;
  S2[(normal + 2) % 3] = s2 / N2;
  VecPlusEquals(VecPlusEquals(X0, 0.5, S1), 0.5, S2);
  VecCross(S1, S2, dA0);
}

void PlaneGrid::GetPoint(int n1, int n2, double X[3], double dA[3]) const {
  X[0] = X0[0] + n1 * S1[0] + n2 * S2[0];
  X[1] = X0[1] + n1 * S1[1] + n2 * S2[1];
  X[2] = X0[2] + n1 * S1[2] + n2 * S2[2];
  dA[0] = dA0[0]; dA[1] = dA0[1]; dA[2] = dA0[2];
}

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/
/* Routine for evaluating arbitrary functions of the fields on a 2d
   surface grid.  Returns a NULL-terminated (malloc'ed) array of
   HMatrix pointers.  If a non-NULL incident field function is
   supplied, then the total of scattered plus incident fields is
   used.  If KN == NULL, then the scattered fields are set to zero. */

HMatrix **RWGGeometry::GetFieldsGrids(SurfaceGrid &grid,
				      int nfuncs, FieldFunc **funcs,
				      cdouble Omega, HVector *KN,
				      EHFuncType2 EHFunc, void *EHFuncUD,
				      int nThread) {
  HMatrix **Ms = (HMatrix **) mallocEC(sizeof(HMatrix*) * (nfuncs + 1));

  for (int i = 0; i < nfuncs; ++i)
    Ms[i] = grid.AllocateGrid(funcs[i]->IsReal());

  // NULL-terminate the returned array; this makes it easier for SWIG
  // to automatically convert the return value;
  Ms[nfuncs] = NULL;
  if (!nfuncs) return Ms;

  if (nThread <= 0) nThread = GetNumThreads();

  // number of threads to use in GetFields
#ifdef USE_OPENMP // grid loop is parallelized
  int nThreadFields = nThread / (grid.N1 * grid.N2);
  if (nThreadFields < 1) nThreadFields = 1;
#else // grid loop not parallelized
  int nThreadFields = nThread;
#endif

  // make sure kdPanels is initialized outside loop (this is not thread-safe)
  for (int no = 0; no < NumObjects; ++no)
    Objects[no]->InitkdPanels(false, LogLevel);

#ifdef USE_OPENMP
#pragma omp parallel for firstprivate(grid,nFuncs,funcs,Omega,KN,EHFunc,EHFuncUD,nThreadFields,Ms,this), schedule(static), num_threads(nThread), collapse(2)
#endif
  for (int n1 = 0; n1 < grid.N1; ++n1)
    for (int n2 = 0; n2 < grid.N2; ++n2) {
      if (LogLevel >= SCUFF_VERBOSELOGGING) {
	int count = n1*grid.N2 + n2;
	for (int PerCent=0; PerCent<9; PerCent++)
	  if (count == (PerCent*grid.N1*grid.N2)/10)
	    MutexLog("%i0 %% (%i/%i)...",PerCent,count,grid.N1*grid.N2);
      }

      double X[3], dA[3];
      grid.GetPoint(n1, n2, X, dA);

      int no = GetObjectIndex(X);
      cdouble EH[6];
      if (KN)
	GetFields(X, no, Omega, KN, EH, nThreadFields);
      else
	memset(EH, 0, sizeof(cdouble) * 6);

      if (EHFunc) {
	cdouble EHi[6];
	EHFunc(X, EHFuncUD, EHi, -2, no);
	for (int j = 0; j < 6; ++j) EH[j] += EHi[j];
      }

      // in theory this could be done outside the grid loop and
      // cached somewhere per-object
      cdouble Eps, Mu;
      if (no >= 0)
	Objects[no]->MP->GetEpsMu(Omega, &Eps, &Mu);
      else
	ExteriorMP->GetEpsMu(Omega, &Eps, &Mu);

      // Call the field functions and store in the output matrices
      for (int i = 0; i < nfuncs; ++i)
	Ms[i]->SetEntry(n1,n2, funcs[i]->Eval(X, dA, EH, Eps, Mu));
    }

  return Ms;
}

/***************************************************************************/
// Convenience wrappers of GetFieldsGrids with simpler arguments.

typedef struct {
  EHFuncType EHFunc; void *EHFuncUD;
} EHFunc_wrap_data;

// EHFuncType2 wrapper for EHFunc, assuming exterior incident fields
static void EHFunc_wrap(const double R[3], void *UserData, cdouble EH[6],
			  int exterior_index, int interior_index) {
  (void) exterior_index; // unused
  if (interior_index == -1) {
    EHFunc_wrap_data *d = (EHFunc_wrap_data *) UserData;
    d->EHFunc(R, d->EHFuncUD, EH);
  }
  else
    memset(EH, 0, sizeof(cdouble) * 6);
}

HMatrix **RWGGeometry::GetFieldsGrids(SurfaceGrid &grid, 
				      int nfuncs, FieldFunc **funcs,
				      cdouble Omega, HVector *KN,
				      EHFuncType EHFunc, void *EHFuncUD,
                                      int nThread) {
  EHFunc_wrap_data d;
  d.EHFunc = EHFunc; d.EHFuncUD = EHFuncUD;
  return GetFieldsGrids(grid, nfuncs, funcs, Omega, KN,
                        EHFunc ? EHFunc_wrap : NULL, (void*) &d,
                        nThread);  
}

// EHFuncType2 wrapper for IncField
static void IncField_wrap(const double R[3], void *UserData, cdouble EH[6],
			  int exterior_index, int interior_index) {
  (void) exterior_index; // unused
  IncField *inc = (IncField *) UserData;
  memset(EH, 0, sizeof(cdouble) * 6);
  for (; inc; inc = inc->Next) 
    if (inc->ObjectIndex == interior_index) {
      cdouble EHi[6];
      inc->GetFields(R, EHi);
      for (int j = 0; j < 6; ++j) EH[j] += EHi[j];
    }
}

HMatrix **RWGGeometry::GetFieldsGrids(SurfaceGrid &grid, 
				      int nfuncs, FieldFunc **funcs,
				    cdouble Omega, HVector *KN, IncField *inc,
                                      int nThread) {
  UpdateIncFields(inc, Omega); // update ObjectIndex/Eps/Mu/Omega of inc
  return GetFieldsGrids(grid, nfuncs, funcs, Omega, KN,
			inc ? IncField_wrap : NULL, (void*) inc,
			nThread);
}

HMatrix **RWGGeometry::GetFieldsGrids(SurfaceGrid &grid, const char *exprs_,
			  cdouble Omega, HVector *KN, IncField *inc,
				      int nThread) {
  FieldFunc **f;
  int nf = 1;

  // count the number of comma-separated expressions
  for (const char *c = strchr(exprs_,','); c; c = strchr(c+1,','))
    ++nf;
  
  f = (FieldFunc**) mallocEC(sizeof(FieldFunc*) * nf);

  // extract each comma-separated expression and parse it
  char *exprs = strdup(exprs_);
  char *expr = exprs;
  for (int i = 0; i < nf; ++i) {
    char *c = strchr(expr,',');
    if (c) *c = 0;
    f[i] = new ParsedFieldFunc(expr);
    expr = c+1;
  }
  free(exprs);

  if (LogLevel >= SCUFF_VERBOSELOGGING)
    Log(" GetFieldsGrids %dx%d of: %s ...", grid.N1,grid.N2, exprs_);

  HMatrix **Ms = GetFieldsGrids(grid, nf, f, Omega, KN, inc, nThread);

  for (int i; i < nf; ++i) delete f[i];
  free(f);
  
  return Ms;
}

HMatrix *RWGGeometry::GetFieldsGrid(SurfaceGrid &grid, FieldFunc &func,
				    cdouble Omega, HVector *KN, IncField *inc,
				    int nThread) {
  FieldFunc *funcs = &func;
  HMatrix *M, **Ms = GetFieldsGrids(grid, 1, &funcs, Omega, KN, inc, nThread);
  M = *Ms;
  free(Ms);
  return M;
}

HMatrix *RWGGeometry::GetFieldsGrid(SurfaceGrid &grid, const char *expr,
				    cdouble Omega, HVector *KN, IncField *inc,
				    int nThread) {
  ParsedFieldFunc f(expr);
  if (LogLevel >= SCUFF_VERBOSELOGGING)
    Log(" GetFieldsGrid %dx%d of: %s ...", grid.N1,grid.N2, expr);
  return GetFieldsGrid(grid, f, Omega, KN, inc, nThread);
}

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

} // namespace scuff
