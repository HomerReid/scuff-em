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
				      IncField *IF) {

  if (IF)
   UpdateIncFields(IF, Omega);

  HMatrix **Ms = (HMatrix **) mallocEC(sizeof(HMatrix*) * (nfuncs + 1));

  for (int i = 0; i < nfuncs; ++i)
    Ms[i] = grid.AllocateGrid(funcs[i]->IsReal());

  // NULL-terminate the returned array; this makes it easier for SWIG
  // to automatically convert the return value;
  Ms[nfuncs] = NULL;
  if (!nfuncs) return Ms;

  int nThread = GetNumThreads();

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
#pragma omp parallel for firstprivate(grid,nFuncs,funcs,Omega,KN,IF,nThreadFields,Ms,this), schedule(static), num_threads(nThread), collapse(2)
#endif
  for (int n1 = 0; n1 < grid.N1; ++n1)
    for (int n2 = 0; n2 < grid.N2; ++n2) {
      int count = n1*grid.N2 + n2;
      if (LogLevel >= SCUFF_VERBOSELOGGING)
	for (int PerCent=0; PerCent<9; PerCent++)
	  if (count == (PerCent*grid.N1*grid.N2)/10)
	    MutexLog("%i0 %% (%i/%i)...",PerCent,count,grid.N1*grid.N2);

      double X[3], dA[3];
      grid.GetPoint(n1, n2, X, dA);

      int no = GetObjectIndex(X);
      cdouble EH[6];
      if (KN)
	GetFields(0, KN, Omega, X, EH);
      else
	memset(EH, 0, sizeof(cdouble) * 6);

     if (IF)
      { cdouble EHi[6];
        IncField *IFNode;
        for(IFNode=IF; IFNode; IFNode=IFNode->Next)
         if ( IFNode->ObjectIndex == no )
          { IFNode->GetFields(X, EHi);
            SixVecPlusEquals(EH, 1.0, EHi);
          };
      };

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

HMatrix **RWGGeometry::GetFieldsGrids(SurfaceGrid &grid, const char *exprs_,
                                      cdouble Omega, HVector *KN, IncField *inc) {
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

  HMatrix **Ms = GetFieldsGrids(grid, nf, f, Omega, KN, inc);

  for (int i=0; i < nf; ++i) delete f[i];
  free(f);
  
  return Ms;
}

HMatrix *RWGGeometry::GetFieldsGrid(SurfaceGrid &grid, FieldFunc &func,
				    cdouble Omega, HVector *KN, IncField *inc) {
  FieldFunc *funcs = &func;
  HMatrix *M, **Ms = GetFieldsGrids(grid, 1, &funcs, Omega, KN, inc);
  M = *Ms;
  free(Ms);
  return M;
}

HMatrix *RWGGeometry::GetFieldsGrid(SurfaceGrid &grid, const char *expr,
				    cdouble Omega, HVector *KN, IncField *inc) {
  ParsedFieldFunc f(expr);
  if (LogLevel >= SCUFF_VERBOSELOGGING)
    Log(" GetFieldsGrid %dx%d of: %s ...", grid.N1,grid.N2, expr);
  return GetFieldsGrid(grid, f, Omega, KN, inc);
}

/***************************************************************************/
/***************************************************************************/
/***************************************************************************/

} // namespace scuff
