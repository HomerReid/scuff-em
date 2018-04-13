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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#else
#	define F77_FUNC(name,NAME) name ## _
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifndef __CLAPACK_H
#define __CLAPACK_H

#include "lapack_names.h"

typedef int logical;
typedef int L_fp;
typedef int ftnlen;

typedef std::complex<double> cdouble;
typedef std::complex<float>  cfloat;

/* following prototype added by homer 6/2005 */
double zlange_(char *Norm, int *M, int *N, cdouble *A, int *LDA, double *Work);
/* following prototype added by homer 3/2015 */
double dlange_(char *Norm, int *M, int *N, double *A, int *LDA, double *Work);

/* following prototype added by homer 11/2008 */
double dlamch_(const char *which);

/* following prototypes added by homer 3/2012 */
double ddot_(int *N, double *X, int *IncX, double *Y, int *IncY);
cdouble zdotc_(int *N, cdouble *X, int *IncX, cdouble *Y, int *IncY);
cdouble zdotu_(int *N, cdouble *X, int *IncX, cdouble *Y, int *IncY);
 
/* Subroutine */ int cbdsqr_(const char *uplo, int *n, int *ncvt, int *
	nru, int *ncc, float *d__, float *e, cfloat *vt, int *ldvt, 
	cfloat *u, int *ldu, cfloat *c__, int *ldc, float *rwork, 
	int *info);
 
/* Subroutine */ int cgbbrd_(const char *vect, int *m, int *n, int *ncc,
	 int *kl, int *ku, cfloat *ab, int *ldab, float *d__, 
	float *e, cfloat *q, int *ldq, cfloat *pt, int *ldpt, 
	cfloat *c__, int *ldc, cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int cgbcon_(const char *norm, int *n, int *kl, int *ku,
	 cfloat *ab, int *ldab, int *ipiv, float *anorm, float *rcond, 
	cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int cgbequ_(int *m, int *n, int *kl, int *ku,
	 cfloat *ab, int *ldab, float *r__, float *c__, float *rowcnd, float 
	*colcnd, float *amax, int *info);
 
/* Subroutine */ int cgbrfs_(const char *trans, int *n, int *kl, int *
	ku, int *nrhs, cfloat *ab, int *ldab, cfloat *afb, int *
	ldafb, int *ipiv, cfloat *b, int *ldb, cfloat *x, int *
	ldx, float *ferr, float *berr, cfloat *work, float *rwork, int *
	info);
 
/* Subroutine */ int cgbsv_(int *n, int *kl, int *ku, int *
	nrhs, cfloat *ab, int *ldab, int *ipiv, cfloat *b, int *
	ldb, int *info);
 
/* Subroutine */ int cgbsvx_(const char *fact, const char *trans, int *n, int *kl,
	 int *ku, int *nrhs, cfloat *ab, int *ldab, cfloat *afb,
	 int *ldafb, int *ipiv, const char *equed, float *r__, float *c__, 
	cfloat *b, int *ldb, cfloat *x, int *ldx, float *rcond, float 
	*ferr, float *berr, cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int cgbtf2_(int *m, int *n, int *kl, int *ku,
	 cfloat *ab, int *ldab, int *ipiv, int *info);
 
/* Subroutine */ int cgbtrf_(int *m, int *n, int *kl, int *ku,
	 cfloat *ab, int *ldab, int *ipiv, int *info);
 
/* Subroutine */ int cgbtrs_(const char *trans, int *n, int *kl, int *
	ku, int *nrhs, cfloat *ab, int *ldab, int *ipiv, cfloat 
	*b, int *ldb, int *info);
 
/* Subroutine */ int cgebak_(const char *job, const char *side, int *n, int *ilo, 
	int *ihi, float *scale, int *m, cfloat *v, int *ldv, 
	int *info);
 
/* Subroutine */ int cgebal_(const char *job, int *n, cfloat *a, int *lda, 
	int *ilo, int *ihi, float *scale, int *info);
 
/* Subroutine */ int cgebd2_(int *m, int *n, cfloat *a, int *lda,
	 float *d__, float *e, cfloat *tauq, cfloat *taup, cfloat *work, 
	int *info);
 
/* Subroutine */ int cgebrd_(int *m, int *n, cfloat *a, int *lda,
	 float *d__, float *e, cfloat *tauq, cfloat *taup, cfloat *work, 
	int *lwork, int *info);
 
/* Subroutine */ int cgecon_(const char *norm, int *n, cfloat *a, int *lda,
	 float *anorm, float *rcond, cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int cgeequ_(int *m, int *n, cfloat *a, int *lda,
	 float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, 
	int *info);
 
/* Subroutine */ int cgees_(const char *jobvs, const char *sort, L_fp select, int *n, 
	cfloat *a, int *lda, int *sdim, cfloat *w, cfloat *vs, 
	int *ldvs, cfloat *work, int *lwork, float *rwork, logical *
	bwork, int *info);
 
/* Subroutine */ int cgeesx_(const char *jobvs, const char *sort, L_fp select, const char *
	sense, int *n, cfloat *a, int *lda, int *sdim, cfloat *
	w, cfloat *vs, int *ldvs, float *rconde, float *rcondv, cfloat *
	work, int *lwork, float *rwork, logical *bwork, int *info);
 
/* Subroutine */ int cgeev_(const char *jobvl, const char *jobvr, int *n, cfloat *a, 
	int *lda, cfloat *w, cfloat *vl, int *ldvl, cfloat *vr, 
	int *ldvr, cfloat *work, int *lwork, float *rwork, int *
	info);
 
/* Subroutine */ int cgeevx_(const char *balanc, const char *jobvl, const char *jobvr, const char *
	sense, int *n, cfloat *a, int *lda, cfloat *w, cfloat *vl, 
	int *ldvl, cfloat *vr, int *ldvr, int *ilo, int *ihi,
	 float *scale, float *abnrm, float *rconde, float *rcondv, cfloat *work, 
	int *lwork, float *rwork, int *info);
 
/* Subroutine */ int cgegs_(const char *jobvsl, const char *jobvsr, int *n, cfloat *
	a, int *lda, cfloat *b, int *ldb, cfloat *alpha, cfloat *
	beta, cfloat *vsl, int *ldvsl, cfloat *vsr, int *ldvsr, 
	cfloat *work, int *lwork, float *rwork, int *info);
 
/* Subroutine */ int cgegv_(const char *jobvl, const char *jobvr, int *n, cfloat *a, 
	int *lda, cfloat *b, int *ldb, cfloat *alpha, cfloat *beta,
	 cfloat *vl, int *ldvl, cfloat *vr, int *ldvr, cfloat *
	work, int *lwork, float *rwork, int *info);
 
/* Subroutine */ int cgehd2_(int *n, int *ilo, int *ihi, cfloat *
	a, int *lda, cfloat *tau, cfloat *work, int *info);
 
/* Subroutine */ int cgehrd_(int *n, int *ilo, int *ihi, cfloat *
	a, int *lda, cfloat *tau, cfloat *work, int *lwork, int 
	*info);
 
/* Subroutine */ int cgelq2_(int *m, int *n, cfloat *a, int *lda,
	 cfloat *tau, cfloat *work, int *info);
 
/* Subroutine */ int cgelqf_(int *m, int *n, cfloat *a, int *lda,
	 cfloat *tau, cfloat *work, int *lwork, int *info);
 
/* Subroutine */ int cgels_(const char *trans, int *m, int *n, int *
	nrhs, cfloat *a, int *lda, cfloat *b, int *ldb, cfloat *
	work, int *lwork, int *info);
 
/* Subroutine */ int cgelsx_(int *m, int *n, int *nrhs, cfloat *
	a, int *lda, cfloat *b, int *ldb, int *jpvt, float *rcond,
	 int *rank, cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int cgelsy_(int *m, int *n, int *nrhs, cfloat *
	a, int *lda, cfloat *b, int *ldb, int *jpvt, float *rcond,
	 int *rank, cfloat *work, int *lwork, float *rwork, int *
	info);
 
/* Subroutine */ int cgeql2_(int *m, int *n, cfloat *a, int *lda,
	 cfloat *tau, cfloat *work, int *info);
 
/* Subroutine */ int cgeqlf_(int *m, int *n, cfloat *a, int *lda,
	 cfloat *tau, cfloat *work, int *lwork, int *info);
 
/* Subroutine */ int cgeqp3_(int *m, int *n, cfloat *a, int *lda,
	 int *jpvt, cfloat *tau, cfloat *work, int *lwork, float *
	rwork, int *info);
 
/* Subroutine */ int cgeqpf_(int *m, int *n, cfloat *a, int *lda,
	 int *jpvt, cfloat *tau, cfloat *work, float *rwork, int *
	info);
 
/* Subroutine */ int cgeqr2_(int *m, int *n, cfloat *a, int *lda,
	 cfloat *tau, cfloat *work, int *info);
 
/* Subroutine */ int cgeqrf_(int *m, int *n, cfloat *a, int *lda,
	 cfloat *tau, cfloat *work, int *lwork, int *info);
 
/* Subroutine */ int cgerfs_(const char *trans, int *n, int *nrhs, cfloat *
	a, int *lda, cfloat *af, int *ldaf, int *ipiv, cfloat *
	b, int *ldb, cfloat *x, int *ldx, float *ferr, float *berr, 
	cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int cgerq2_(int *m, int *n, cfloat *a, int *lda,
	 cfloat *tau, cfloat *work, int *info);
 
/* Subroutine */ int cgerqf_(int *m, int *n, cfloat *a, int *lda,
	 cfloat *tau, cfloat *work, int *lwork, int *info);
 
/* Subroutine */ int cgesc2_(int *n, cfloat *a, int *lda, cfloat *
	rhs, int *ipiv, int *jpiv, float *scale);
 
/* Subroutine */ int cgesv_(int *n, int *nrhs, cfloat *a, int *
	lda, int *ipiv, cfloat *b, int *ldb, int *info);
 
/* Subroutine */ int cgesvx_(const char *fact, const char *trans, int *n, int *
	nrhs, cfloat *a, int *lda, cfloat *af, int *ldaf, int *
	ipiv, const char *equed, float *r__, float *c__, cfloat *b, int *ldb, 
	cfloat *x, int *ldx, float *rcond, float *ferr, float *berr, 
	cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int cgetc2_(int *n, cfloat *a, int *lda, int *
	ipiv, int *jpiv, int *info);
 
/* Subroutine */ int cgetf2_(int *m, int *n, cfloat *a, int *lda,
	 int *ipiv, int *info);
 
/* Subroutine */ int cgetrf_(int *m, int *n, cfloat *a, int *lda,
	 int *ipiv, int *info);
 
/* Subroutine */ int cgetri_(int *n, cfloat *a, int *lda, int *
	ipiv, cfloat *work, int *lwork, int *info);
 
/* Subroutine */ int cgetrs_(const char *trans, int *n, int *nrhs, cfloat *
	a, int *lda, int *ipiv, cfloat *b, int *ldb, int *
	info);
 
/* Subroutine */ int cggbak_(const char *job, const char *side, int *n, int *ilo, 
	int *ihi, float *lscale, float *rscale, int *m, cfloat *v, 
	int *ldv, int *info);
 
/* Subroutine */ int cggbal_(const char *job, int *n, cfloat *a, int *lda, 
	cfloat *b, int *ldb, int *ilo, int *ihi, float *lscale, 
	float *rscale, float *work, int *info);
 
/* Subroutine */ int cgges_(const char *jobvsl, const char *jobvsr, const char *sort, L_fp 
	selctg, int *n, cfloat *a, int *lda, cfloat *b, int *
	ldb, int *sdim, cfloat *alpha, cfloat *beta, cfloat *vsl, 
	int *ldvsl, cfloat *vsr, int *ldvsr, cfloat *work, int *
	lwork, float *rwork, logical *bwork, int *info);
 
/* Subroutine */ int cggesx_(const char *jobvsl, const char *jobvsr, const char *sort, L_fp 
	selctg, const char *sense, int *n, cfloat *a, int *lda, cfloat *b,
	 int *ldb, int *sdim, cfloat *alpha, cfloat *beta, cfloat *
	vsl, int *ldvsl, cfloat *vsr, int *ldvsr, float *rconde, float 
	*rcondv, cfloat *work, int *lwork, float *rwork, int *iwork, 
	int *liwork, logical *bwork, int *info);
 
/* Subroutine */ int cggev_(const char *jobvl, const char *jobvr, int *n, cfloat *a, 
	int *lda, cfloat *b, int *ldb, cfloat *alpha, cfloat *beta,
	 cfloat *vl, int *ldvl, cfloat *vr, int *ldvr, cfloat *
	work, int *lwork, float *rwork, int *info);
 
/* Subroutine */ int cggevx_(const char *balanc, const char *jobvl, const char *jobvr, const char *
	sense, int *n, cfloat *a, int *lda, cfloat *b, int *ldb,
	 cfloat *alpha, cfloat *beta, cfloat *vl, int *ldvl, cfloat *
	vr, int *ldvr, int *ilo, int *ihi, float *lscale, float *
	rscale, float *abnrm, float *bbnrm, float *rconde, float *rcondv, cfloat 
	*work, int *lwork, float *rwork, int *iwork, logical *bwork, 
	int *info);
 
/* Subroutine */ int cggglm_(int *n, int *m, int *p, cfloat *a, 
	int *lda, cfloat *b, int *ldb, cfloat *d__, cfloat *x, 
	cfloat *y, cfloat *work, int *lwork, int *info);
 
/* Subroutine */ int cgghrd_(const char *compq, const char *compz, int *n, int *
	ilo, int *ihi, cfloat *a, int *lda, cfloat *b, int *ldb,
	 cfloat *q, int *ldq, cfloat *z__, int *ldz, int *info);
 
/* Subroutine */ int cgglse_(int *m, int *n, int *p, cfloat *a, 
	int *lda, cfloat *b, int *ldb, cfloat *c__, cfloat *d__, 
	cfloat *x, cfloat *work, int *lwork, int *info);
 
/* Subroutine */ int cggqrf_(int *n, int *m, int *p, cfloat *a, 
	int *lda, cfloat *taua, cfloat *b, int *ldb, cfloat *taub, 
	cfloat *work, int *lwork, int *info);
 
/* Subroutine */ int cggrqf_(int *m, int *p, int *n, cfloat *a, 
	int *lda, cfloat *taua, cfloat *b, int *ldb, cfloat *taub, 
	cfloat *work, int *lwork, int *info);
 
/* Subroutine */ int cggsvd_(const char *jobu, const char *jobv, const char *jobq, int *m, 
	int *n, int *p, int *k, int *l, cfloat *a, int *
	lda, cfloat *b, int *ldb, float *alpha, float *beta, cfloat *u, 
	int *ldu, cfloat *v, int *ldv, cfloat *q, int *ldq, 
	cfloat *work, float *rwork, int *iwork, int *info);
 
/* Subroutine */ int cggsvp_(const char *jobu, const char *jobv, const char *jobq, int *m, 
	int *p, int *n, cfloat *a, int *lda, cfloat *b, int 
	*ldb, float *tola, float *tolb, int *k, int *l, cfloat *u, 
	int *ldu, cfloat *v, int *ldv, cfloat *q, int *ldq, 
	int *iwork, float *rwork, cfloat *tau, cfloat *work, int *
	info);
 
/* Subroutine */ int cgtcon_(const char *norm, int *n, cfloat *dl, cfloat *
	d__, cfloat *du, cfloat *du2, int *ipiv, float *anorm, float *
	rcond, cfloat *work, int *info);
 
/* Subroutine */ int cgtrfs_(const char *trans, int *n, int *nrhs, cfloat *
	dl, cfloat *d__, cfloat *du, cfloat *dlf, cfloat *df, cfloat *
	duf, cfloat *du2, int *ipiv, cfloat *b, int *ldb, cfloat *
	x, int *ldx, float *ferr, float *berr, cfloat *work, float *rwork, 
	int *info);
 
/* Subroutine */ int cgtsv_(int *n, int *nrhs, cfloat *dl, cfloat *
	d__, cfloat *du, cfloat *b, int *ldb, int *info);
 
/* Subroutine */ int cgtsvx_(const char *fact, const char *trans, int *n, int *
	nrhs, cfloat *dl, cfloat *d__, cfloat *du, cfloat *dlf, cfloat *
	df, cfloat *duf, cfloat *du2, int *ipiv, cfloat *b, int *
	ldb, cfloat *x, int *ldx, float *rcond, float *ferr, float *berr, 
	cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int cgttrf_(int *n, cfloat *dl, cfloat *d__, cfloat *
	du, cfloat *du2, int *ipiv, int *info);
 
/* Subroutine */ int cgttrs_(const char *trans, int *n, int *nrhs, cfloat *
	dl, cfloat *d__, cfloat *du, cfloat *du2, int *ipiv, cfloat *
	b, int *ldb, int *info);
 
/* Subroutine */ int cgtts2_(int *itrans, int *n, int *nrhs, 
	cfloat *dl, cfloat *d__, cfloat *du, cfloat *du2, int *ipiv, 
	cfloat *b, int *ldb);
 
/* Subroutine */ int chbev_(const char *jobz, const char *uplo, int *n, int *kd, 
	cfloat *ab, int *ldab, float *w, cfloat *z__, int *ldz, 
	cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int chbevd_(const char *jobz, const char *uplo, int *n, int *kd, 
	cfloat *ab, int *ldab, float *w, cfloat *z__, int *ldz, 
	cfloat *work, int *lwork, float *rwork, int *lrwork, int *
	iwork, int *liwork, int *info);
 
/* Subroutine */ int chbevx_(const char *jobz, const char *range, const char *uplo, int *n, 
	int *kd, cfloat *ab, int *ldab, cfloat *q, int *ldq, 
	float *vl, float *vu, int *il, int *iu, float *abstol, int *
	m, float *w, cfloat *z__, int *ldz, cfloat *work, float *rwork, 
	int *iwork, int *ifail, int *info);
 
/* Subroutine */ int chbgst_(const char *vect, const char *uplo, int *n, int *ka, 
	int *kb, cfloat *ab, int *ldab, cfloat *bb, int *ldbb, 
	cfloat *x, int *ldx, cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int chbgv_(const char *jobz, const char *uplo, int *n, int *ka, 
	int *kb, cfloat *ab, int *ldab, cfloat *bb, int *ldbb, 
	float *w, cfloat *z__, int *ldz, cfloat *work, float *rwork, 
	int *info);
 
/* Subroutine */ int chbgvx_(const char *jobz, const char *range, const char *uplo, int *n, 
	int *ka, int *kb, cfloat *ab, int *ldab, cfloat *bb, 
	int *ldbb, cfloat *q, int *ldq, float *vl, float *vu, int *
	il, int *iu, float *abstol, int *m, float *w, cfloat *z__, 
	int *ldz, cfloat *work, float *rwork, int *iwork, int *
	ifail, int *info);
 
/* Subroutine */ int chbtrd_(const char *vect, const char *uplo, int *n, int *kd, 
	cfloat *ab, int *ldab, float *d__, float *e, cfloat *q, int *
	ldq, cfloat *work, int *info);
 
/* Subroutine */ int checon_(const char *uplo, int *n, cfloat *a, int *lda,
	 int *ipiv, float *anorm, float *rcond, cfloat *work, int *
	info);
 
/* Subroutine */ int cheev_(const char *jobz, const char *uplo, int *n, cfloat *a, 
	int *lda, float *w, cfloat *work, int *lwork, float *rwork, 
	int *info);
 
/* Subroutine */ int cheevd_(const char *jobz, const char *uplo, int *n, cfloat *a, 
	int *lda, float *w, cfloat *work, int *lwork, float *rwork, 
	int *lrwork, int *iwork, int *liwork, int *info);
 
/* Subroutine */ int cheevr_(const char *jobz, const char *range, const char *uplo, int *n, 
	cfloat *a, int *lda, float *vl, float *vu, int *il, int *
	iu, float *abstol, int *m, float *w, cfloat *z__, int *ldz, 
	int *isuppz, cfloat *work, int *lwork, float *rwork, int *
	lrwork, int *iwork, int *liwork, int *info);
 
/* Subroutine */ int cheevx_(const char *jobz, const char *range, const char *uplo, int *n, 
	cfloat *a, int *lda, float *vl, float *vu, int *il, int *
	iu, float *abstol, int *m, float *w, cfloat *z__, int *ldz, 
	cfloat *work, int *lwork, float *rwork, int *iwork, int *
	ifail, int *info);
 
/* Subroutine */ int chegs2_(int *itype, const char *uplo, int *n, cfloat *
	a, int *lda, cfloat *b, int *ldb, int *info);
 
/* Subroutine */ int chegst_(int *itype, const char *uplo, int *n, cfloat *
	a, int *lda, cfloat *b, int *ldb, int *info);
 
/* Subroutine */ int chegv_(int *itype, const char *jobz, const char *uplo, int *
	n, cfloat *a, int *lda, cfloat *b, int *ldb, float *w, 
	cfloat *work, int *lwork, float *rwork, int *info);
 
/* Subroutine */ int chegvd_(int *itype, const char *jobz, const char *uplo, int *
	n, cfloat *a, int *lda, cfloat *b, int *ldb, float *w, 
	cfloat *work, int *lwork, float *rwork, int *lrwork, int *
	iwork, int *liwork, int *info);
 
/* Subroutine */ int chegvx_(int *itype, const char *jobz, const char *range, const char *
	uplo, int *n, cfloat *a, int *lda, cfloat *b, int *ldb, 
	float *vl, float *vu, int *il, int *iu, float *abstol, int *
	m, float *w, cfloat *z__, int *ldz, cfloat *work, int *lwork,
	 float *rwork, int *iwork, int *ifail, int *info);
 
/* Subroutine */ int cherfs_(const char *uplo, int *n, int *nrhs, cfloat *
	a, int *lda, cfloat *af, int *ldaf, int *ipiv, cfloat *
	b, int *ldb, cfloat *x, int *ldx, float *ferr, float *berr, 
	cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int chesv_(const char *uplo, int *n, int *nrhs, cfloat *a,
	 int *lda, int *ipiv, cfloat *b, int *ldb, cfloat *work,
	 int *lwork, int *info);
 
/* Subroutine */ int chesvx_(const char *fact, const char *uplo, int *n, int *
	nrhs, cfloat *a, int *lda, cfloat *af, int *ldaf, int *
	ipiv, cfloat *b, int *ldb, cfloat *x, int *ldx, float *rcond,
	 float *ferr, float *berr, cfloat *work, int *lwork, float *rwork, 
	int *info);
 
/* Subroutine */ int chetf2_(const char *uplo, int *n, cfloat *a, int *lda,
	 int *ipiv, int *info);
 
/* Subroutine */ int chetrd_(const char *uplo, int *n, cfloat *a, int *lda,
	 float *d__, float *e, cfloat *tau, cfloat *work, int *lwork, 
	int *info);
 
/* Subroutine */ int chetrf_(const char *uplo, int *n, cfloat *a, int *lda,
	 int *ipiv, cfloat *work, int *lwork, int *info);
 
/* Subroutine */ int chetri_(const char *uplo, int *n, cfloat *a, int *lda,
	 int *ipiv, cfloat *work, int *info);
 
/* Subroutine */ int chetrs_(const char *uplo, int *n, int *nrhs, cfloat *
	a, int *lda, int *ipiv, cfloat *b, int *ldb, int *
	info);
 
/* Subroutine */ int chgeqz_(const char *job, const char *compq, const char *compz, int *n, 
	int *ilo, int *ihi, cfloat *a, int *lda, cfloat *b, 
	int *ldb, cfloat *alpha, cfloat *beta, cfloat *q, int *ldq,
	 cfloat *z__, int *ldz, cfloat *work, int *lwork, float *
	rwork, int *info);
 
/* Subroutine */ int chpcon_(const char *uplo, int *n, cfloat *ap, int *
	ipiv, float *anorm, float *rcond, cfloat *work, int *info);
 
/* Subroutine */ int chpev_(const char *jobz, const char *uplo, int *n, cfloat *ap, 
	float *w, cfloat *z__, int *ldz, cfloat *work, float *rwork, 
	int *info);
 
/* Subroutine */ int chpevd_(const char *jobz, const char *uplo, int *n, cfloat *ap, 
	float *w, cfloat *z__, int *ldz, cfloat *work, int *lwork, 
	float *rwork, int *lrwork, int *iwork, int *liwork, 
	int *info);
 
/* Subroutine */ int chpevx_(const char *jobz, const char *range, const char *uplo, int *n, 
	cfloat *ap, float *vl, float *vu, int *il, int *iu, float *
	abstol, int *m, float *w, cfloat *z__, int *ldz, cfloat *
	work, float *rwork, int *iwork, int *ifail, int *info);
 
/* Subroutine */ int chpgst_(int *itype, const char *uplo, int *n, cfloat *
	ap, cfloat *bp, int *info);
 
/* Subroutine */ int chpgv_(int *itype, const char *jobz, const char *uplo, int *
	n, cfloat *ap, cfloat *bp, float *w, cfloat *z__, int *ldz, 
	cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int chpgvd_(int *itype, const char *jobz, const char *uplo, int *
	n, cfloat *ap, cfloat *bp, float *w, cfloat *z__, int *ldz, 
	cfloat *work, int *lwork, float *rwork, int *lrwork, int *
	iwork, int *liwork, int *info);
 
/* Subroutine */ int chpgvx_(int *itype, const char *jobz, const char *range, const char *
	uplo, int *n, cfloat *ap, cfloat *bp, float *vl, float *vu, 
	int *il, int *iu, float *abstol, int *m, float *w, cfloat *
	z__, int *ldz, cfloat *work, float *rwork, int *iwork, 
	int *ifail, int *info);
 
/* Subroutine */ int chprfs_(const char *uplo, int *n, int *nrhs, cfloat *
	ap, cfloat *afp, int *ipiv, cfloat *b, int *ldb, cfloat *x,
	 int *ldx, float *ferr, float *berr, cfloat *work, float *rwork, 
	int *info);
 
/* Subroutine */ int chpsv_(const char *uplo, int *n, int *nrhs, cfloat *
	ap, int *ipiv, cfloat *b, int *ldb, int *info);
 
/* Subroutine */ int chpsvx_(const char *fact, const char *uplo, int *n, int *
	nrhs, cfloat *ap, cfloat *afp, int *ipiv, cfloat *b, int *
	ldb, cfloat *x, int *ldx, float *rcond, float *ferr, float *berr, 
	cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int chptrd_(const char *uplo, int *n, cfloat *ap, float *d__, 
	float *e, cfloat *tau, int *info);
 
/* Subroutine */ int chptrf_(const char *uplo, int *n, cfloat *ap, int *
	ipiv, int *info);
 
/* Subroutine */ int chptri_(const char *uplo, int *n, cfloat *ap, int *
	ipiv, cfloat *work, int *info);
 
/* Subroutine */ int chptrs_(const char *uplo, int *n, int *nrhs, cfloat *
	ap, int *ipiv, cfloat *b, int *ldb, int *info);
 
/* Subroutine */ int chsein_(const char *side, const char *eigsrc, const char *initv, logical *
	select, int *n, cfloat *h__, int *ldh, cfloat *w, cfloat *
	vl, int *ldvl, cfloat *vr, int *ldvr, int *mm, int *
	m, cfloat *work, float *rwork, int *ifaill, int *ifailr, 
	int *info);
 
/* Subroutine */ int chseqr_(const char *job, const char *compz, int *n, int *ilo,
	 int *ihi, cfloat *h__, int *ldh, cfloat *w, cfloat *z__, 
	int *ldz, cfloat *work, int *lwork, int *info);
 
/* Subroutine */ int clabrd_(int *m, int *n, int *nb, cfloat *a, 
	int *lda, float *d__, float *e, cfloat *tauq, cfloat *taup, 
	cfloat *x, int *ldx, cfloat *y, int *ldy);
 
/* Subroutine */ int clacgv_(int *n, cfloat *x, int *incx);
 
/* Subroutine */ int clacon_(int *n, cfloat *v, cfloat *x, float *est, 
	int *kase);
 
/* Subroutine */ int clacp2_(const char *uplo, int *m, int *n, float *a, 
	int *lda, cfloat *b, int *ldb);
 
/* Subroutine */ int clacpy_(const char *uplo, int *m, int *n, cfloat *a, 
	int *lda, cfloat *b, int *ldb);
 
/* Subroutine */ int clacrm_(int *m, int *n, cfloat *a, int *lda,
	 float *b, int *ldb, cfloat *c__, int *ldc, float *rwork);
 
/* Subroutine */ int clacrt_(int *n, cfloat *cx, int *incx, cfloat *
	cy, int *incy, cfloat *c__, cfloat *s);
 
/* Subroutine */ int claed0_(int *qsiz, int *n, float *d__, float *e, 
	cfloat *q, int *ldq, cfloat *qstore, int *ldqs, float *rwork,
	 int *iwork, int *info);
 
/* Subroutine */ int claed7_(int *n, int *cutpnt, int *qsiz, 
	int *tlvls, int *curlvl, int *curpbm, float *d__, cfloat *
	q, int *ldq, float *rho, int *indxq, float *qstore, int *
	qptr, int *prmptr, int *perm, int *givptr, int *
	givcol, float *givnum, cfloat *work, float *rwork, int *iwork, 
	int *info);
 
/* Subroutine */ int claed8_(int *k, int *n, int *qsiz, cfloat *
	q, int *ldq, float *d__, float *rho, int *cutpnt, float *z__, 
	float *dlamda, cfloat *q2, int *ldq2, float *w, int *indxp, 
	int *indx, int *indxq, int *perm, int *givptr, 
	int *givcol, float *givnum, int *info);
 
/* Subroutine */ int claein_(logical *rightv, logical *noinit, int *n, 
	cfloat *h__, int *ldh, cfloat *w, cfloat *v, cfloat *b, 
	int *ldb, float *rwork, float *eps3, float *smlnum, int *info);
 
/* Subroutine */ int claesy_(cfloat *a, cfloat *b, cfloat *c__, cfloat *
	rt1, cfloat *rt2, cfloat *evscal, cfloat *cs1, cfloat *sn1);
 
/* Subroutine */ int claev2_(cfloat *a, cfloat *b, cfloat *c__, float *rt1, 
	float *rt2, float *cs1, cfloat *sn1);
 
/* Subroutine */ int clags2_(logical *upper, float *a1, cfloat *a2, float *a3, 
	float *b1, cfloat *b2, float *b3, float *csu, cfloat *snu, float *csv, 
	cfloat *snv, float *csq, cfloat *snq);
 
/* Subroutine */ int clagtm_(const char *trans, int *n, int *nrhs, float *
	alpha, cfloat *dl, cfloat *d__, cfloat *du, cfloat *x, int *
	ldx, float *beta, cfloat *b, int *ldb);
 
/* Subroutine */ int clahef_(const char *uplo, int *n, int *nb, int *kb,
	 cfloat *a, int *lda, int *ipiv, cfloat *w, int *ldw, 
	int *info);
 
/* Subroutine */ int clahqr_(logical *wantt, logical *wantz, int *n, 
	int *ilo, int *ihi, cfloat *h__, int *ldh, cfloat *w, 
	int *iloz, int *ihiz, cfloat *z__, int *ldz, int *
	info);
 
/* Subroutine */ int clahrd_(int *n, int *k, int *nb, cfloat *a, 
	int *lda, cfloat *tau, cfloat *t, int *ldt, cfloat *y, 
	int *ldy);
 
/* Subroutine */ int claic1_(int *job, int *j, cfloat *x, float *sest,
	 cfloat *w, cfloat *gamma, float *sestpr, cfloat *s, cfloat *c__);
 
/* Subroutine */ int clals0_(int *icompq, int *nl, int *nr, 
	int *sqre, int *nrhs, cfloat *b, int *ldb, cfloat *bx, 
	int *ldbx, int *perm, int *givptr, int *givcol, 
	int *ldgcol, float *givnum, int *ldgnum, float *poles, float *
	difl, float *difr, float *z__, int *k, float *c__, float *s, float *
	rwork, int *info);
 
/* Subroutine */ int clalsa_(int *icompq, int *smlsiz, int *n, 
	int *nrhs, cfloat *b, int *ldb, cfloat *bx, int *ldbx, 
	float *u, int *ldu, float *vt, int *k, float *difl, float *difr, 
	float *z__, float *poles, int *givptr, int *givcol, int *
	ldgcol, int *perm, float *givnum, float *c__, float *s, float *rwork, 
	int *iwork, int *info);
 
/* Subroutine */ int clapll_(int *n, cfloat *x, int *incx, cfloat *
	y, int *incy, float *ssmin);
 
/* Subroutine */ int clapmt_(logical *forwrd, int *m, int *n, cfloat 
	*x, int *ldx, int *k);
 
/* Subroutine */ int claqgb_(int *m, int *n, int *kl, int *ku,
	 cfloat *ab, int *ldab, float *r__, float *c__, float *rowcnd, float 
	*colcnd, float *amax, const char *equed);
 
/* Subroutine */ int claqge_(int *m, int *n, cfloat *a, int *lda,
	 float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, const char *
	equed);
 
/* Subroutine */ int claqhb_(const char *uplo, int *n, int *kd, cfloat *ab,
	 int *ldab, float *s, float *scond, float *amax, const char *equed);
 
/* Subroutine */ int claqhe_(const char *uplo, int *n, cfloat *a, int *lda,
	 float *s, float *scond, float *amax, const char *equed);
 
/* Subroutine */ int claqhp_(const char *uplo, int *n, cfloat *ap, float *s, 
	float *scond, float *amax, const char *equed);
 
/* Subroutine */ int claqp2_(int *m, int *n, int *offset, cfloat 
	*a, int *lda, int *jpvt, cfloat *tau, float *vn1, float *vn2, 
	cfloat *work);
 
/* Subroutine */ int claqps_(int *m, int *n, int *offset, int 
	*nb, int *kb, cfloat *a, int *lda, int *jpvt, cfloat *
	tau, float *vn1, float *vn2, cfloat *auxv, cfloat *f, int *ldf);
 
/* Subroutine */ int claqsb_(const char *uplo, int *n, int *kd, cfloat *ab,
	 int *ldab, float *s, float *scond, float *amax, const char *equed);
 
/* Subroutine */ int claqsp_(const char *uplo, int *n, cfloat *ap, float *s, 
	float *scond, float *amax, const char *equed);
 
/* Subroutine */ int claqsy_(const char *uplo, int *n, cfloat *a, int *lda,
	 float *s, float *scond, float *amax, const char *equed);
 
/* Subroutine */ int clar1v_(int *n, int *b1, int *bn, float *
	sigma, float *d__, float *l, float *ld, float *lld, float *gersch, cfloat 
	*z__, float *ztz, float *mingma, int *r__, int *isuppz, float *
	work);
 
/* Subroutine */ int clar2v_(int *n, cfloat *x, cfloat *y, cfloat *z__,
	 int *incx, float *c__, cfloat *s, int *incc);
 
/* Subroutine */ int clarcm_(int *m, int *n, float *a, int *lda, 
	cfloat *b, int *ldb, cfloat *c__, int *ldc, float *rwork);
 
/* Subroutine */ int clarf_(const char *side, int *m, int *n, cfloat *v, 
	int *incv, cfloat *tau, cfloat *c__, int *ldc, cfloat *
	work);
 
/* Subroutine */ int clarfb_(const char *side, const char *trans, const char *direct, const char *
	storev, int *m, int *n, int *k, cfloat *v, int *ldv, 
	cfloat *t, int *ldt, cfloat *c__, int *ldc, cfloat *work, 
	int *ldwork);
 
/* Subroutine */ int clarfg_(int *n, cfloat *alpha, cfloat *x, int *
	incx, cfloat *tau);
 
/* Subroutine */ int clarft_(const char *direct, const char *storev, int *n, int *
	k, cfloat *v, int *ldv, cfloat *tau, cfloat *t, int *ldt);
 
/* Subroutine */ int clarfx_(const char *side, int *m, int *n, cfloat *v, 
	cfloat *tau, cfloat *c__, int *ldc, cfloat *work);
 
/* Subroutine */ int clargv_(int *n, cfloat *x, int *incx, cfloat *
	y, int *incy, float *c__, int *incc);
 
/* Subroutine */ int clarnv_(int *idist, int *iseed, int *n, 
	cfloat *x);
 
/* Subroutine */ int clarrv_(int *n, float *d__, float *l, int *isplit, 
	int *m, float *w, int *iblock, float *gersch, float *tol, 
	cfloat *z__, int *ldz, int *isuppz, float *work, int *
	iwork, int *info);
 
/* Subroutine */ int clartg_(cfloat *f, cfloat *g, float *cs, cfloat *sn, 
	cfloat *r__);
 
/* Subroutine */ int clartv_(int *n, cfloat *x, int *incx, cfloat *
	y, int *incy, float *c__, cfloat *s, int *incc);
 
/* Subroutine */ int clarz_(const char *side, int *m, int *n, int *l, 
	cfloat *v, int *incv, cfloat *tau, cfloat *c__, int *ldc, 
	cfloat *work);
 
/* Subroutine */ int clarzb_(const char *side, const char *trans, const char *direct, const char *
	storev, int *m, int *n, int *k, int *l, cfloat *v, 
	int *ldv, cfloat *t, int *ldt, cfloat *c__, int *ldc, 
	cfloat *work, int *ldwork);
 
/* Subroutine */ int clarzt_(const char *direct, const char *storev, int *n, int *
	k, cfloat *v, int *ldv, cfloat *tau, cfloat *t, int *ldt);
 
/* Subroutine */ int clascl_(const char *type__, int *kl, int *ku, float *
	cfrom, float *cto, int *m, int *n, cfloat *a, int *lda, 
	int *info);
 
/* Subroutine */ int claset_(const char *uplo, int *m, int *n, cfloat *
	alpha, cfloat *beta, cfloat *a, int *lda);
 
/* Subroutine */ int clasr_(const char *side, const char *pivot, const char *direct, int *m,
	 int *n, float *c__, float *s, cfloat *a, int *lda);
 
/* Subroutine */ int classq_(int *n, cfloat *x, int *incx, float *
	scale, float *sumsq);
 
/* Subroutine */ int claswp_(int *n, cfloat *a, int *lda, int *
	k1, int *k2, int *ipiv, int *incx);
 
/* Subroutine */ int clasyf_(const char *uplo, int *n, int *nb, int *kb,
	 cfloat *a, int *lda, int *ipiv, cfloat *w, int *ldw, 
	int *info);
 
/* Subroutine */ int clatbs_(const char *uplo, const char *trans, const char *diag, const char *
	normin, int *n, int *kd, cfloat *ab, int *ldab, cfloat *
	x, float *scale, float *cnorm, int *info);
 
/* Subroutine */ int clatdf_(int *ijob, int *n, cfloat *z__, int 
	*ldz, cfloat *rhs, float *rdsum, float *rdscal, int *ipiv, int 
	*jpiv);
 
/* Subroutine */ int clatps_(const char *uplo, const char *trans, const char *diag, const char *
	normin, int *n, cfloat *ap, cfloat *x, float *scale, float *cnorm,
	 int *info);
 
/* Subroutine */ int clatrd_(const char *uplo, int *n, int *nb, cfloat *a, 
	int *lda, float *e, cfloat *tau, cfloat *w, int *ldw);
 
/* Subroutine */ int clatrs_(const char *uplo, const char *trans, const char *diag, const char *
	normin, int *n, cfloat *a, int *lda, cfloat *x, float *scale,
	 float *cnorm, int *info);
 
/* Subroutine */ int clatrz_(int *m, int *n, int *l, cfloat *a, 
	int *lda, cfloat *tau, cfloat *work);
 
/* Subroutine */ int clatzm_(const char *side, int *m, int *n, cfloat *v, 
	int *incv, cfloat *tau, cfloat *c1, cfloat *c2, int *ldc, 
	cfloat *work);
 
/* Subroutine */ int clauu2_(const char *uplo, int *n, cfloat *a, int *lda,
	 int *info);
 
/* Subroutine */ int clauum_(const char *uplo, int *n, cfloat *a, int *lda,
	 int *info);
 
/* Subroutine */ int cpbcon_(const char *uplo, int *n, int *kd, cfloat *ab,
	 int *ldab, float *anorm, float *rcond, cfloat *work, float *rwork, 
	int *info);
 
/* Subroutine */ int cpbequ_(const char *uplo, int *n, int *kd, cfloat *ab,
	 int *ldab, float *s, float *scond, float *amax, int *info);
 
/* Subroutine */ int cpbrfs_(const char *uplo, int *n, int *kd, int *
	nrhs, cfloat *ab, int *ldab, cfloat *afb, int *ldafb, 
	cfloat *b, int *ldb, cfloat *x, int *ldx, float *ferr, float *
	berr, cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int cpbstf_(const char *uplo, int *n, int *kd, cfloat *ab,
	 int *ldab, int *info);
 
/* Subroutine */ int cpbsv_(const char *uplo, int *n, int *kd, int *
	nrhs, cfloat *ab, int *ldab, cfloat *b, int *ldb, int *
	info);
 
/* Subroutine */ int cpbsvx_(const char *fact, const char *uplo, int *n, int *kd, 
	int *nrhs, cfloat *ab, int *ldab, cfloat *afb, int *
	ldafb, const char *equed, float *s, cfloat *b, int *ldb, cfloat *x, 
	int *ldx, float *rcond, float *ferr, float *berr, cfloat *work, 
	float *rwork, int *info);
 
/* Subroutine */ int cpbtf2_(const char *uplo, int *n, int *kd, cfloat *ab,
	 int *ldab, int *info);
 
/* Subroutine */ int cpbtrf_(const char *uplo, int *n, int *kd, cfloat *ab,
	 int *ldab, int *info);
 
/* Subroutine */ int cpbtrs_(const char *uplo, int *n, int *kd, int *
	nrhs, cfloat *ab, int *ldab, cfloat *b, int *ldb, int *
	info);
 
/* Subroutine */ int cpocon_(const char *uplo, int *n, cfloat *a, int *lda,
	 float *anorm, float *rcond, cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int cpoequ_(int *n, cfloat *a, int *lda, float *s, 
	float *scond, float *amax, int *info);
 
/* Subroutine */ int cporfs_(const char *uplo, int *n, int *nrhs, cfloat *
	a, int *lda, cfloat *af, int *ldaf, cfloat *b, int *ldb,
	 cfloat *x, int *ldx, float *ferr, float *berr, cfloat *work, 
	float *rwork, int *info);
 
/* Subroutine */ int cposv_(const char *uplo, int *n, int *nrhs, cfloat *a,
	 int *lda, cfloat *b, int *ldb, int *info);
 
/* Subroutine */ int cposvx_(const char *fact, const char *uplo, int *n, int *
	nrhs, cfloat *a, int *lda, cfloat *af, int *ldaf, const char *
	equed, float *s, cfloat *b, int *ldb, cfloat *x, int *ldx, 
	float *rcond, float *ferr, float *berr, cfloat *work, float *rwork, 
	int *info);
 
/* Subroutine */ int cpotf2_(const char *uplo, int *n, cfloat *a, int *lda,
	 int *info);
 
/* Subroutine */ int cpotrf_(const char *uplo, int *n, cfloat *a, int *lda,
	 int *info);
 
/* Subroutine */ int cpotri_(const char *uplo, int *n, cfloat *a, int *lda,
	 int *info);
 
/* Subroutine */ int cpotrs_(const char *uplo, int *n, int *nrhs, cfloat *
	a, int *lda, cfloat *b, int *ldb, int *info);
 
/* Subroutine */ int cppcon_(const char *uplo, int *n, cfloat *ap, float *anorm,
	 float *rcond, cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int cppequ_(const char *uplo, int *n, cfloat *ap, float *s, 
	float *scond, float *amax, int *info);
 
/* Subroutine */ int cpprfs_(const char *uplo, int *n, int *nrhs, cfloat *
	ap, cfloat *afp, cfloat *b, int *ldb, cfloat *x, int *ldx, 
	float *ferr, float *berr, cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int cppsv_(const char *uplo, int *n, int *nrhs, cfloat *
	ap, cfloat *b, int *ldb, int *info);
 
/* Subroutine */ int cppsvx_(const char *fact, const char *uplo, int *n, int *
	nrhs, cfloat *ap, cfloat *afp, const char *equed, float *s, cfloat *b, 
	int *ldb, cfloat *x, int *ldx, float *rcond, float *ferr, float 
	*berr, cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int cpptrf_(const char *uplo, int *n, cfloat *ap, int *
	info);
 
/* Subroutine */ int cpptri_(const char *uplo, int *n, cfloat *ap, int *
	info);
 
/* Subroutine */ int cpptrs_(const char *uplo, int *n, int *nrhs, cfloat *
	ap, cfloat *b, int *ldb, int *info);
 
/* Subroutine */ int cptcon_(int *n, float *d__, cfloat *e, float *anorm, 
	float *rcond, float *rwork, int *info);
 
/* Subroutine */ int cptrfs_(const char *uplo, int *n, int *nrhs, float *d__,
	 cfloat *e, float *df, cfloat *ef, cfloat *b, int *ldb, cfloat 
	*x, int *ldx, float *ferr, float *berr, cfloat *work, float *rwork, 
	int *info);
 
/* Subroutine */ int cptsv_(int *n, int *nrhs, float *d__, cfloat *e, 
	cfloat *b, int *ldb, int *info);
 
/* Subroutine */ int cptsvx_(const char *fact, int *n, int *nrhs, float *d__,
	 cfloat *e, float *df, cfloat *ef, cfloat *b, int *ldb, cfloat 
	*x, int *ldx, float *rcond, float *ferr, float *berr, cfloat *work, 
	float *rwork, int *info);
 
/* Subroutine */ int cpttrf_(int *n, float *d__, cfloat *e, int *info);
 
/* Subroutine */ int cpttrs_(const char *uplo, int *n, int *nrhs, float *d__,
	 cfloat *e, cfloat *b, int *ldb, int *info);
 
/* Subroutine */ int cptts2_(int *iuplo, int *n, int *nrhs, float *
	d__, cfloat *e, cfloat *b, int *ldb);
 
/* Subroutine */ int crot_(int *n, cfloat *cx, int *incx, cfloat *
	cy, int *incy, float *c__, cfloat *s);
 
/* Subroutine */ int cspcon_(const char *uplo, int *n, cfloat *ap, int *
	ipiv, float *anorm, float *rcond, cfloat *work, int *info);
 
/* Subroutine */ int cspmv_(const char *uplo, int *n, cfloat *alpha, cfloat *
	ap, cfloat *x, int *incx, cfloat *beta, cfloat *y, int *
	incy);
 
/* Subroutine */ int cspr_(const char *uplo, int *n, cfloat *alpha, cfloat *x,
	 int *incx, cfloat *ap);
 
/* Subroutine */ int csprfs_(const char *uplo, int *n, int *nrhs, cfloat *
	ap, cfloat *afp, int *ipiv, cfloat *b, int *ldb, cfloat *x,
	 int *ldx, float *ferr, float *berr, cfloat *work, float *rwork, 
	int *info);
 
/* Subroutine */ int cspsv_(const char *uplo, int *n, int *nrhs, cfloat *
	ap, int *ipiv, cfloat *b, int *ldb, int *info);
 
/* Subroutine */ int cspsvx_(const char *fact, const char *uplo, int *n, int *
	nrhs, cfloat *ap, cfloat *afp, int *ipiv, cfloat *b, int *
	ldb, cfloat *x, int *ldx, float *rcond, float *ferr, float *berr, 
	cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int csptrf_(const char *uplo, int *n, cfloat *ap, int *
	ipiv, int *info);
 
/* Subroutine */ int csptri_(const char *uplo, int *n, cfloat *ap, int *
	ipiv, cfloat *work, int *info);
 
/* Subroutine */ int csptrs_(const char *uplo, int *n, int *nrhs, cfloat *
	ap, int *ipiv, cfloat *b, int *ldb, int *info);
 
/* Subroutine */ int csrot_(int *n, cfloat *cx, int *incx, cfloat *
	cy, int *incy, float *c__, float *s);
 
/* Subroutine */ int csrscl_(int *n, float *sa, cfloat *sx, int *incx);
 
/* Subroutine */ int cstedc_(const char *compz, int *n, float *d__, float *e, 
	cfloat *z__, int *ldz, cfloat *work, int *lwork, float *
	rwork, int *lrwork, int *iwork, int *liwork, int *
	info);
 
/* Subroutine */ int cstein_(int *n, float *d__, float *e, int *m, float 
	*w, int *iblock, int *isplit, cfloat *z__, int *ldz, 
	float *work, int *iwork, int *ifail, int *info);
 
/* Subroutine */ int csteqr_(const char *compz, int *n, float *d__, float *e, 
	cfloat *z__, int *ldz, float *work, int *info);
 
/* Subroutine */ int csycon_(const char *uplo, int *n, cfloat *a, int *lda,
	 int *ipiv, float *anorm, float *rcond, cfloat *work, int *
	info);
 
/* Subroutine */ int csymv_(const char *uplo, int *n, cfloat *alpha, cfloat *
	a, int *lda, cfloat *x, int *incx, cfloat *beta, cfloat *y,
	 int *incy);
 
/* Subroutine */ int csyr_(const char *uplo, int *n, cfloat *alpha, cfloat *x,
	 int *incx, cfloat *a, int *lda);
 
/* Subroutine */ int csyrfs_(const char *uplo, int *n, int *nrhs, cfloat *
	a, int *lda, cfloat *af, int *ldaf, int *ipiv, cfloat *
	b, int *ldb, cfloat *x, int *ldx, float *ferr, float *berr, 
	cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int csysv_(const char *uplo, int *n, int *nrhs, cfloat *a,
	 int *lda, int *ipiv, cfloat *b, int *ldb, cfloat *work,
	 int *lwork, int *info);
 
/* Subroutine */ int csysvx_(const char *fact, const char *uplo, int *n, int *
	nrhs, cfloat *a, int *lda, cfloat *af, int *ldaf, int *
	ipiv, cfloat *b, int *ldb, cfloat *x, int *ldx, float *rcond,
	 float *ferr, float *berr, cfloat *work, int *lwork, float *rwork, 
	int *info);
 
/* Subroutine */ int csytf2_(const char *uplo, int *n, cfloat *a, int *lda,
	 int *ipiv, int *info);
 
/* Subroutine */ int csytrf_(const char *uplo, int *n, cfloat *a, int *lda,
	 int *ipiv, cfloat *work, int *lwork, int *info);
 
/* Subroutine */ int csytri_(const char *uplo, int *n, cfloat *a, int *lda,
	 int *ipiv, cfloat *work, int *info);
 
/* Subroutine */ int csytrs_(const char *uplo, int *n, int *nrhs, cfloat *
	a, int *lda, int *ipiv, cfloat *b, int *ldb, int *
	info);
 
/* Subroutine */ int ctbcon_(const char *norm, const char *uplo, const char *diag, int *n, 
	int *kd, cfloat *ab, int *ldab, float *rcond, cfloat *work, 
	float *rwork, int *info);
 
/* Subroutine */ int ctbrfs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *kd, int *nrhs, cfloat *ab, int *ldab, cfloat *b, 
	int *ldb, cfloat *x, int *ldx, float *ferr, float *berr, 
	cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int ctbtrs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *kd, int *nrhs, cfloat *ab, int *ldab, cfloat *b, 
	int *ldb, int *info);
 
/* Subroutine */ int ctgevc_(const char *side, const char *howmny, logical *select, 
	int *n, cfloat *a, int *lda, cfloat *b, int *ldb, 
	cfloat *vl, int *ldvl, cfloat *vr, int *ldvr, int *mm, 
	int *m, cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int ctgex2_(logical *wantq, logical *wantz, int *n, 
	cfloat *a, int *lda, cfloat *b, int *ldb, cfloat *q, 
	int *ldq, cfloat *z__, int *ldz, int *j1, int *info);
 
/* Subroutine */ int ctgexc_(logical *wantq, logical *wantz, int *n, 
	cfloat *a, int *lda, cfloat *b, int *ldb, cfloat *q, 
	int *ldq, cfloat *z__, int *ldz, int *ifst, int *
	ilst, int *info);
 
/* Subroutine */ int ctgsen_(int *ijob, logical *wantq, logical *wantz, 
	logical *select, int *n, cfloat *a, int *lda, cfloat *b, 
	int *ldb, cfloat *alpha, cfloat *beta, cfloat *q, int *ldq,
	 cfloat *z__, int *ldz, int *m, float *pl, float *pr, float *
	dif, cfloat *work, int *lwork, int *iwork, int *liwork, 
	int *info);
 
/* Subroutine */ int ctgsja_(const char *jobu, const char *jobv, const char *jobq, int *m, 
	int *p, int *n, int *k, int *l, cfloat *a, int *
	lda, cfloat *b, int *ldb, float *tola, float *tolb, float *alpha, 
	float *beta, cfloat *u, int *ldu, cfloat *v, int *ldv, 
	cfloat *q, int *ldq, cfloat *work, int *ncycle, int *
	info);
 
/* Subroutine */ int ctgsna_(const char *job, const char *howmny, logical *select, 
	int *n, cfloat *a, int *lda, cfloat *b, int *ldb, 
	cfloat *vl, int *ldvl, cfloat *vr, int *ldvr, float *s, float 
	*dif, int *mm, int *m, cfloat *work, int *lwork, int 
	*iwork, int *info);
 
/* Subroutine */ int ctgsy2_(const char *trans, int *ijob, int *m, int *
	n, cfloat *a, int *lda, cfloat *b, int *ldb, cfloat *c__, 
	int *ldc, cfloat *d__, int *ldd, cfloat *e, int *lde, 
	cfloat *f, int *ldf, float *scale, float *rdsum, float *rdscal, 
	int *info);
 
/* Subroutine */ int ctgsyl_(const char *trans, int *ijob, int *m, int *
	n, cfloat *a, int *lda, cfloat *b, int *ldb, cfloat *c__, 
	int *ldc, cfloat *d__, int *ldd, cfloat *e, int *lde, 
	cfloat *f, int *ldf, float *scale, float *dif, cfloat *work, 
	int *lwork, int *iwork, int *info);
 
/* Subroutine */ int ctpcon_(const char *norm, const char *uplo, const char *diag, int *n, 
	cfloat *ap, float *rcond, cfloat *work, float *rwork, int *info);
 
/* Subroutine */ int ctprfs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *nrhs, cfloat *ap, cfloat *b, int *ldb, cfloat *x, 
	int *ldx, float *ferr, float *berr, cfloat *work, float *rwork, 
	int *info);
 
/* Subroutine */ int ctptri_(const char *uplo, const char *diag, int *n, cfloat *ap, 
	int *info);
 
/* Subroutine */ int ctptrs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *nrhs, cfloat *ap, cfloat *b, int *ldb, int *info);
 
/* Subroutine */ int ctrcon_(const char *norm, const char *uplo, const char *diag, int *n, 
	cfloat *a, int *lda, float *rcond, cfloat *work, float *rwork, 
	int *info);
 
/* Subroutine */ int ctrevc_(const char *side, const char *howmny, logical *select, 
	int *n, cfloat *t, int *ldt, cfloat *vl, int *ldvl, 
	cfloat *vr, int *ldvr, int *mm, int *m, cfloat *work, 
	float *rwork, int *info);
 
/* Subroutine */ int ctrexc_(const char *compq, int *n, cfloat *t, int *
	ldt, cfloat *q, int *ldq, int *ifst, int *ilst, int *
	info);
 
/* Subroutine */ int ctrrfs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *nrhs, cfloat *a, int *lda, cfloat *b, int *ldb, 
	cfloat *x, int *ldx, float *ferr, float *berr, cfloat *work, float 
	*rwork, int *info);
 
/* Subroutine */ int ctrsen_(const char *job, const char *compq, logical *select, int 
	*n, cfloat *t, int *ldt, cfloat *q, int *ldq, cfloat *w, 
	int *m, float *s, float *sep, cfloat *work, int *lwork, 
	int *info);
 
/* Subroutine */ int ctrsna_(const char *job, const char *howmny, logical *select, 
	int *n, cfloat *t, int *ldt, cfloat *vl, int *ldvl, 
	cfloat *vr, int *ldvr, float *s, float *sep, int *mm, int *
	m, cfloat *work, int *ldwork, float *rwork, int *info);
 
/* Subroutine */ int ctrsyl_(const char *trana, const char *tranb, int *isgn, int 
	*m, int *n, cfloat *a, int *lda, cfloat *b, int *ldb, 
	cfloat *c__, int *ldc, float *scale, int *info);
 
/* Subroutine */ int ctrti2_(const char *uplo, const char *diag, int *n, cfloat *a, 
	int *lda, int *info);
 
/* Subroutine */ int ctrtri_(const char *uplo, const char *diag, int *n, cfloat *a, 
	int *lda, int *info);
 
/* Subroutine */ int ctrtrs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *nrhs, cfloat *a, int *lda, cfloat *b, int *ldb, 
	int *info);
 
/* Subroutine */ int ctzrqf_(int *m, int *n, cfloat *a, int *lda,
	 cfloat *tau, int *info);
 
/* Subroutine */ int ctzrzf_(int *m, int *n, cfloat *a, int *lda,
	 cfloat *tau, cfloat *work, int *lwork, int *info);
 
/* Subroutine */ int cung2l_(int *m, int *n, int *k, cfloat *a, 
	int *lda, cfloat *tau, cfloat *work, int *info);
 
/* Subroutine */ int cung2r_(int *m, int *n, int *k, cfloat *a, 
	int *lda, cfloat *tau, cfloat *work, int *info);
 
/* Subroutine */ int cungbr_(const char *vect, int *m, int *n, int *k, 
	cfloat *a, int *lda, cfloat *tau, cfloat *work, int *lwork,
	 int *info);
 
/* Subroutine */ int cunghr_(int *n, int *ilo, int *ihi, cfloat *
	a, int *lda, cfloat *tau, cfloat *work, int *lwork, int 
	*info);
 
/* Subroutine */ int cungl2_(int *m, int *n, int *k, cfloat *a, 
	int *lda, cfloat *tau, cfloat *work, int *info);
 
/* Subroutine */ int cunglq_(int *m, int *n, int *k, cfloat *a, 
	int *lda, cfloat *tau, cfloat *work, int *lwork, int *
	info);
 
/* Subroutine */ int cungql_(int *m, int *n, int *k, cfloat *a, 
	int *lda, cfloat *tau, cfloat *work, int *lwork, int *
	info);
 
/* Subroutine */ int cungqr_(int *m, int *n, int *k, cfloat *a, 
	int *lda, cfloat *tau, cfloat *work, int *lwork, int *
	info);
 
/* Subroutine */ int cungr2_(int *m, int *n, int *k, cfloat *a, 
	int *lda, cfloat *tau, cfloat *work, int *info);
 
/* Subroutine */ int cungrq_(int *m, int *n, int *k, cfloat *a, 
	int *lda, cfloat *tau, cfloat *work, int *lwork, int *
	info);
 
/* Subroutine */ int cungtr_(const char *uplo, int *n, cfloat *a, int *lda,
	 cfloat *tau, cfloat *work, int *lwork, int *info);
 
/* Subroutine */ int cunm2l_(const char *side, const char *trans, int *m, int *n, 
	int *k, cfloat *a, int *lda, cfloat *tau, cfloat *c__, 
	int *ldc, cfloat *work, int *info);
 
/* Subroutine */ int cunm2r_(const char *side, const char *trans, int *m, int *n, 
	int *k, cfloat *a, int *lda, cfloat *tau, cfloat *c__, 
	int *ldc, cfloat *work, int *info);
 
/* Subroutine */ int cunmbr_(const char *vect, const char *side, const char *trans, int *m, 
	int *n, int *k, cfloat *a, int *lda, cfloat *tau, 
	cfloat *c__, int *ldc, cfloat *work, int *lwork, int *
	info);
 
/* Subroutine */ int cunmhr_(const char *side, const char *trans, int *m, int *n, 
	int *ilo, int *ihi, cfloat *a, int *lda, cfloat *tau, 
	cfloat *c__, int *ldc, cfloat *work, int *lwork, int *
	info);
 
/* Subroutine */ int cunml2_(const char *side, const char *trans, int *m, int *n, 
	int *k, cfloat *a, int *lda, cfloat *tau, cfloat *c__, 
	int *ldc, cfloat *work, int *info);
 
/* Subroutine */ int cunmlq_(const char *side, const char *trans, int *m, int *n, 
	int *k, cfloat *a, int *lda, cfloat *tau, cfloat *c__, 
	int *ldc, cfloat *work, int *lwork, int *info);
 
/* Subroutine */ int cunmql_(const char *side, const char *trans, int *m, int *n, 
	int *k, cfloat *a, int *lda, cfloat *tau, cfloat *c__, 
	int *ldc, cfloat *work, int *lwork, int *info);
 
/* Subroutine */ int cunmqr_(const char *side, const char *trans, int *m, int *n, 
	int *k, cfloat *a, int *lda, cfloat *tau, cfloat *c__, 
	int *ldc, cfloat *work, int *lwork, int *info);
 
/* Subroutine */ int cunmr2_(const char *side, const char *trans, int *m, int *n, 
	int *k, cfloat *a, int *lda, cfloat *tau, cfloat *c__, 
	int *ldc, cfloat *work, int *info);
 
/* Subroutine */ int cunmr3_(const char *side, const char *trans, int *m, int *n, 
	int *k, int *l, cfloat *a, int *lda, cfloat *tau, 
	cfloat *c__, int *ldc, cfloat *work, int *info);
 
/* Subroutine */ int cunmrq_(const char *side, const char *trans, int *m, int *n, 
	int *k, cfloat *a, int *lda, cfloat *tau, cfloat *c__, 
	int *ldc, cfloat *work, int *lwork, int *info);
 
/* Subroutine */ int cunmrz_(const char *side, const char *trans, int *m, int *n, 
	int *k, int *l, cfloat *a, int *lda, cfloat *tau, 
	cfloat *c__, int *ldc, cfloat *work, int *lwork, int *
	info);
 
/* Subroutine */ int cunmtr_(const char *side, const char *uplo, const char *trans, int *m, 
	int *n, cfloat *a, int *lda, cfloat *tau, cfloat *c__, 
	int *ldc, cfloat *work, int *lwork, int *info);
 
/* Subroutine */ int cupgtr_(const char *uplo, int *n, cfloat *ap, cfloat *
	tau, cfloat *q, int *ldq, cfloat *work, int *info);
 
/* Subroutine */ int cupmtr_(const char *side, const char *uplo, const char *trans, int *m, 
	int *n, cfloat *ap, cfloat *tau, cfloat *c__, int *ldc, 
	cfloat *work, int *info);
 
/* Subroutine */ int dbdsdc_(const char *uplo, const char *compq, int *n, double *
	d__, double *e, double *u, int *ldu, double *vt, 
	int *ldvt, double *q, int *iq, double *work, int *
	iwork, int *info);
 
/* Subroutine */ int dbdsqr_(const char *uplo, int *n, int *ncvt, int *
	nru, int *ncc, double *d__, double *e, double *vt, 
	int *ldvt, double *u, int *ldu, double *c__, int *
	ldc, double *work, int *info);
 
/* Subroutine */ int ddisna_(const char *job, int *m, int *n, double *
	d__, double *sep, int *info);
 
/* Subroutine */ int dgbbrd_(const char *vect, int *m, int *n, int *ncc,
	 int *kl, int *ku, double *ab, int *ldab, double *
	d__, double *e, double *q, int *ldq, double *pt, 
	int *ldpt, double *c__, int *ldc, double *work, 
	int *info);
 
/* Subroutine */ int dgbcon_(const char *norm, int *n, int *kl, int *ku,
	 double *ab, int *ldab, int *ipiv, double *anorm, 
	double *rcond, double *work, int *iwork, int *info);
 
/* Subroutine */ int dgbequ_(int *m, int *n, int *kl, int *ku,
	 double *ab, int *ldab, double *r__, double *c__, 
	double *rowcnd, double *colcnd, double *amax, int *
	info);
 
/* Subroutine */ int dgbrfs_(const char *trans, int *n, int *kl, int *
	ku, int *nrhs, double *ab, int *ldab, double *afb, 
	int *ldafb, int *ipiv, double *b, int *ldb, 
	double *x, int *ldx, double *ferr, double *berr, 
	double *work, int *iwork, int *info);
 
/* Subroutine */ int dgbsv_(int *n, int *kl, int *ku, int *
	nrhs, double *ab, int *ldab, int *ipiv, double *b, 
	int *ldb, int *info);
 
/* Subroutine */ int dgbsvx_(const char *fact, const char *trans, int *n, int *kl,
	 int *ku, int *nrhs, double *ab, int *ldab, 
	double *afb, int *ldafb, int *ipiv, const char *equed, 
	double *r__, double *c__, double *b, int *ldb, 
	double *x, int *ldx, double *rcond, double *ferr, 
	double *berr, double *work, int *iwork, int *info);
 
/* Subroutine */ int dgbtf2_(int *m, int *n, int *kl, int *ku,
	 double *ab, int *ldab, int *ipiv, int *info);
 
/* Subroutine */ int dgbtrf_(int *m, int *n, int *kl, int *ku,
	 double *ab, int *ldab, int *ipiv, int *info);
 
/* Subroutine */ int dgbtrs_(const char *trans, int *n, int *kl, int *
	ku, int *nrhs, double *ab, int *ldab, int *ipiv, 
	double *b, int *ldb, int *info);
 
/* Subroutine */ int dgebak_(const char *job, const char *side, int *n, int *ilo, 
	int *ihi, double *scale, int *m, double *v, int *
	ldv, int *info);
 
/* Subroutine */ int dgebal_(const char *job, int *n, double *a, int *
	lda, int *ilo, int *ihi, double *scale, int *info);
 
/* Subroutine */ int dgebd2_(int *m, int *n, double *a, int *
	lda, double *d__, double *e, double *tauq, double *
	taup, double *work, int *info);
 
/* Subroutine */ int dgebrd_(int *m, int *n, double *a, int *
	lda, double *d__, double *e, double *tauq, double *
	taup, double *work, int *lwork, int *info);
 
/* Subroutine */ int dgecon_(const char *norm, int *n, double *a, int *
	lda, double *anorm, double *rcond, double *work, int *
	iwork, int *info);
 
/* Subroutine */ int dgeequ_(int *m, int *n, double *a, int *
	lda, double *r__, double *c__, double *rowcnd, double 
	*colcnd, double *amax, int *info);
 
/* Subroutine */ int dgees_(const char *jobvs, const char *sort, L_fp select, int *n, 
	double *a, int *lda, int *sdim, double *wr, 
	double *wi, double *vs, int *ldvs, double *work, 
	int *lwork, logical *bwork, int *info);
 
/* Subroutine */ int dgeesx_(const char *jobvs, const char *sort, L_fp select, const char *
	sense, int *n, double *a, int *lda, int *sdim, 
	double *wr, double *wi, double *vs, int *ldvs, 
	double *rconde, double *rcondv, double *work, int *
	lwork, int *iwork, int *liwork, logical *bwork, int *info);
 
/* Subroutine */ int dgeev_(const char *jobvl, const char *jobvr, int *n, double *
	a, int *lda, double *wr, double *wi, double *vl, 
	int *ldvl, double *vr, int *ldvr, double *work, 
	int *lwork, int *info);
 
/* Subroutine */ int dgeevx_(const char *balanc, const char *jobvl, const char *jobvr, const char *
	sense, int *n, double *a, int *lda, double *wr, 
	double *wi, double *vl, int *ldvl, double *vr, 
	int *ldvr, int *ilo, int *ihi, double *scale, 
	double *abnrm, double *rconde, double *rcondv, double 
	*work, int *lwork, int *iwork, int *info);
 
/* Subroutine */ int dgegs_(const char *jobvsl, const char *jobvsr, int *n, 
	double *a, int *lda, double *b, int *ldb, double *
	alphar, double *alphai, double *beta, double *vsl, 
	int *ldvsl, double *vsr, int *ldvsr, double *work, 
	int *lwork, int *info);
 
/* Subroutine */ int dgegv_(const char *jobvl, const char *jobvr, int *n, double *
	a, int *lda, double *b, int *ldb, double *alphar, 
	double *alphai, double *beta, double *vl, int *ldvl, 
	double *vr, int *ldvr, double *work, int *lwork, 
	int *info);
 
/* Subroutine */ int dgehd2_(int *n, int *ilo, int *ihi, 
	double *a, int *lda, double *tau, double *work, 
	int *info);
 
/* Subroutine */ int dgehrd_(int *n, int *ilo, int *ihi, 
	double *a, int *lda, double *tau, double *work, 
	int *lwork, int *info);
 
/* Subroutine */ int dgelq2_(int *m, int *n, double *a, int *
	lda, double *tau, double *work, int *info);
 
/* Subroutine */ int dgelqf_(int *m, int *n, double *a, int *
	lda, double *tau, double *work, int *lwork, int *info);
 
/* Subroutine */ int dgels_(const char *trans, int *m, int *n, int *
	nrhs, double *a, int *lda, double *b, int *ldb, 
	double *work, int *lwork, int *info);
 
/* Subroutine */ int dgelsd_(int *m, int *n, int *nrhs, 
	double *a, int *lda, double *b, int *ldb, double *
	s, double *rcond, int *rank, double *work, int *lwork,
	 int *iwork, int *info);
 
/* Subroutine */ int dgelss_(int *m, int *n, int *nrhs, 
	double *a, int *lda, double *b, int *ldb, double *
	s, double *rcond, int *rank, double *work, int *lwork,
	 int *info);
 
/* Subroutine */ int dgelsx_(int *m, int *n, int *nrhs, 
	double *a, int *lda, double *b, int *ldb, int *
	jpvt, double *rcond, int *rank, double *work, int *
	info);
 
/* Subroutine */ int dgelsy_(int *m, int *n, int *nrhs, 
	double *a, int *lda, double *b, int *ldb, int *
	jpvt, double *rcond, int *rank, double *work, int *
	lwork, int *info);
 
/* Subroutine */ int dgeql2_(int *m, int *n, double *a, int *
	lda, double *tau, double *work, int *info);
 
/* Subroutine */ int dgeqlf_(int *m, int *n, double *a, int *
	lda, double *tau, double *work, int *lwork, int *info);
 
/* Subroutine */ int dgeqp3_(int *m, int *n, double *a, int *
	lda, int *jpvt, double *tau, double *work, int *lwork,
	 int *info);
 
/* Subroutine */ int dgeqpf_(int *m, int *n, double *a, int *
	lda, int *jpvt, double *tau, double *work, int *info);
 
/* Subroutine */ int dgeqr2_(int *m, int *n, double *a, int *
	lda, double *tau, double *work, int *info);
 
/* Subroutine */ int dgeqrf_(int *m, int *n, double *a, int *
	lda, double *tau, double *work, int *lwork, int *info);
 
/* Subroutine */ int dgerfs_(const char *trans, int *n, int *nrhs, 
	double *a, int *lda, double *af, int *ldaf, int *
	ipiv, double *b, int *ldb, double *x, int *ldx, 
	double *ferr, double *berr, double *work, int *iwork, 
	int *info);
 
/* Subroutine */ int dgerq2_(int *m, int *n, double *a, int *
	lda, double *tau, double *work, int *info);
 
/* Subroutine */ int dgerqf_(int *m, int *n, double *a, int *
	lda, double *tau, double *work, int *lwork, int *info);
 
/* Subroutine */ int dgesc2_(int *n, double *a, int *lda, 
	double *rhs, int *ipiv, int *jpiv, double *scale);
 
/* Subroutine */ int dgesdd_(const char *jobz, int *m, int *n, double *
	a, int *lda, double *s, double *u, int *ldu, 
	double *vt, int *ldvt, double *work, int *lwork, 
	int *iwork, int *info);
 
/* Subroutine */ int dgesv_(int *n, int *nrhs, double *a, int 
	*lda, int *ipiv, double *b, int *ldb, int *info);
 
/* Subroutine */ int dgesvd_(const char *jobu, const char *jobvt, int *m, int *n, 
	double *a, int *lda, double *s, double *u, int *
	ldu, double *vt, int *ldvt, double *work, int *lwork, 
	int *info);
 
/* Subroutine */ int dgesvx_(const char *fact, const char *trans, int *n, int *
	nrhs, double *a, int *lda, double *af, int *ldaf, 
	int *ipiv, const char *equed, double *r__, double *c__, 
	double *b, int *ldb, double *x, int *ldx, double *
	rcond, double *ferr, double *berr, double *work, int *
	iwork, int *info);
 
/* Subroutine */ int dgetc2_(int *n, double *a, int *lda, int 
	*ipiv, int *jpiv, int *info);
 
/* Subroutine */ int dgetf2_(int *m, int *n, double *a, int *
	lda, int *ipiv, int *info);
 
/* Subroutine */ int dgetrf_(int *m, int *n, double *a, int *
	lda, int *ipiv, int *info);
 
/* Subroutine */ int dgetri_(int *n, double *a, int *lda, int 
	*ipiv, double *work, int *lwork, int *info);
 
/* Subroutine */ int dgetrs_(const char *trans, int *n, int *nrhs, 
	double *a, int *lda, int *ipiv, double *b, int *
	ldb, int *info);
 
/* Subroutine */ int dggbak_(const char *job, const char *side, int *n, int *ilo, 
	int *ihi, double *lscale, double *rscale, int *m, 
	double *v, int *ldv, int *info);
 
/* Subroutine */ int dggbal_(const char *job, int *n, double *a, int *
	lda, double *b, int *ldb, int *ilo, int *ihi, 
	double *lscale, double *rscale, double *work, int *
	info);
 
/* Subroutine */ int dgges_(const char *jobvsl, const char *jobvsr, const char *sort, L_fp 
	delctg, int *n, double *a, int *lda, double *b, 
	int *ldb, int *sdim, double *alphar, double *alphai, 
	double *beta, double *vsl, int *ldvsl, double *vsr, 
	int *ldvsr, double *work, int *lwork, logical *bwork, 
	int *info);
 
/* Subroutine */ int dggesx_(const char *jobvsl, const char *jobvsr, const char *sort, L_fp 
	delctg, const char *sense, int *n, double *a, int *lda, 
	double *b, int *ldb, int *sdim, double *alphar, 
	double *alphai, double *beta, double *vsl, int *ldvsl,
	 double *vsr, int *ldvsr, double *rconde, double *
	rcondv, double *work, int *lwork, int *iwork, int *
	liwork, logical *bwork, int *info);
 
/* Subroutine */ int dggev_(const char *jobvl, const char *jobvr, int *n, double *
	a, int *lda, double *b, int *ldb, double *alphar, 
	double *alphai, double *beta, double *vl, int *ldvl, 
	double *vr, int *ldvr, double *work, int *lwork, 
	int *info);
 
/* Subroutine */ int dggevx_(const char *balanc, const char *jobvl, const char *jobvr, const char *
	sense, int *n, double *a, int *lda, double *b, 
	int *ldb, double *alphar, double *alphai, double *
	beta, double *vl, int *ldvl, double *vr, int *ldvr, 
	int *ilo, int *ihi, double *lscale, double *rscale, 
	double *abnrm, double *bbnrm, double *rconde, double *
	rcondv, double *work, int *lwork, int *iwork, logical *
	bwork, int *info);
 
/* Subroutine */ int dggglm_(int *n, int *m, int *p, double *
	a, int *lda, double *b, int *ldb, double *d__, 
	double *x, double *y, double *work, int *lwork, 
	int *info);
 
/* Subroutine */ int dgghrd_(const char *compq, const char *compz, int *n, int *
	ilo, int *ihi, double *a, int *lda, double *b, 
	int *ldb, double *q, int *ldq, double *z__, int *
	ldz, int *info);
 
/* Subroutine */ int dgglse_(int *m, int *n, int *p, double *
	a, int *lda, double *b, int *ldb, double *c__, 
	double *d__, double *x, double *work, int *lwork, 
	int *info);
 
/* Subroutine */ int dggqrf_(int *n, int *m, int *p, double *
	a, int *lda, double *taua, double *b, int *ldb, 
	double *taub, double *work, int *lwork, int *info);
 
/* Subroutine */ int dggrqf_(int *m, int *p, int *n, double *
	a, int *lda, double *taua, double *b, int *ldb, 
	double *taub, double *work, int *lwork, int *info);
 
/* Subroutine */ int dggsvd_(const char *jobu, const char *jobv, const char *jobq, int *m, 
	int *n, int *p, int *k, int *l, double *a, 
	int *lda, double *b, int *ldb, double *alpha, 
	double *beta, double *u, int *ldu, double *v, int 
	*ldv, double *q, int *ldq, double *work, int *iwork, 
	int *info);
 
/* Subroutine */ int dggsvp_(const char *jobu, const char *jobv, const char *jobq, int *m, 
	int *p, int *n, double *a, int *lda, double *b, 
	int *ldb, double *tola, double *tolb, int *k, int 
	*l, double *u, int *ldu, double *v, int *ldv, 
	double *q, int *ldq, int *iwork, double *tau, 
	double *work, int *info);
 
/* Subroutine */ int dgtcon_(const char *norm, int *n, double *dl, 
	double *d__, double *du, double *du2, int *ipiv, 
	double *anorm, double *rcond, double *work, int *
	iwork, int *info);
 
/* Subroutine */ int dgtrfs_(const char *trans, int *n, int *nrhs, 
	double *dl, double *d__, double *du, double *dlf, 
	double *df, double *duf, double *du2, int *ipiv, 
	double *b, int *ldb, double *x, int *ldx, double *
	ferr, double *berr, double *work, int *iwork, int *
	info);
 
/* Subroutine */ int dgtsv_(int *n, int *nrhs, double *dl, 
	double *d__, double *du, double *b, int *ldb, int 
	*info);
 
/* Subroutine */ int dgtsvx_(const char *fact, const char *trans, int *n, int *
	nrhs, double *dl, double *d__, double *du, double *
	dlf, double *df, double *duf, double *du2, int *ipiv, 
	double *b, int *ldb, double *x, int *ldx, double *
	rcond, double *ferr, double *berr, double *work, int *
	iwork, int *info);
 
/* Subroutine */ int dgttrf_(int *n, double *dl, double *d__, 
	double *du, double *du2, int *ipiv, int *info);
 
/* Subroutine */ int dgttrs_(const char *trans, int *n, int *nrhs, 
	double *dl, double *d__, double *du, double *du2, 
	int *ipiv, double *b, int *ldb, int *info);
 
/* Subroutine */ int dgtts2_(int *itrans, int *n, int *nrhs, 
	double *dl, double *d__, double *du, double *du2, 
	int *ipiv, double *b, int *ldb);
 
/* Subroutine */ int dhgeqz_(const char *job, const char *compq, const char *compz, int *n, 
	int *ilo, int *ihi, double *a, int *lda, double *
	b, int *ldb, double *alphar, double *alphai, double *
	beta, double *q, int *ldq, double *z__, int *ldz, 
	double *work, int *lwork, int *info);
 
/* Subroutine */ int dhsein_(const char *side, const char *eigsrc, const char *initv, logical *
	select, int *n, double *h__, int *ldh, double *wr, 
	double *wi, double *vl, int *ldvl, double *vr, 
	int *ldvr, int *mm, int *m, double *work, int *
	ifaill, int *ifailr, int *info);
 
/* Subroutine */ int dhseqr_(const char *job, const char *compz, int *n, int *ilo,
	 int *ihi, double *h__, int *ldh, double *wr, 
	double *wi, double *z__, int *ldz, double *work, 
	int *lwork, int *info);
 
/* Subroutine */ int dlabad_(double *small, double *large);
 
/* Subroutine */ int dlabrd_(int *m, int *n, int *nb, double *
	a, int *lda, double *d__, double *e, double *tauq, 
	double *taup, double *x, int *ldx, double *y, int 
	*ldy);
 
/* Subroutine */ int dlacon_(int *n, double *v, double *x, 
	int *isgn, double *est, int *kase);
 
/* Subroutine */ int dlacpy_(const char *uplo, int *m, int *n, double *
	a, int *lda, double *b, int *ldb);
 
/* Subroutine */ int dladiv_(double *a, double *b, double *c__, 
	double *d__, double *p, double *q);
 
/* Subroutine */ int dlae2_(double *a, double *b, double *c__, 
	double *rt1, double *rt2);
 
/* Subroutine */ int dlaebz_(int *ijob, int *nitmax, int *n, 
	int *mmax, int *minp, int *nbmin, double *abstol, 
	double *reltol, double *pivmin, double *d__, double *
	e, double *e2, int *nval, double *ab, double *c__, 
	int *mout, int *nab, double *work, int *iwork, 
	int *info);
 
/* Subroutine */ int dlaed0_(int *icompq, int *qsiz, int *n, 
	double *d__, double *e, double *q, int *ldq, 
	double *qstore, int *ldqs, double *work, int *iwork, 
	int *info);
 
/* Subroutine */ int dlaed1_(int *n, double *d__, double *q, 
	int *ldq, int *indxq, double *rho, int *cutpnt, 
	double *work, int *iwork, int *info);
 
/* Subroutine */ int dlaed2_(int *k, int *n, int *n1, double *
	d__, double *q, int *ldq, int *indxq, double *rho, 
	double *z__, double *dlamda, double *w, double *q2, 
	int *indx, int *indxc, int *indxp, int *coltyp, 
	int *info);
 
/* Subroutine */ int dlaed3_(int *k, int *n, int *n1, double *
	d__, double *q, int *ldq, double *rho, double *dlamda,
	 double *q2, int *indx, int *ctot, double *w, 
	double *s, int *info);
 
/* Subroutine */ int dlaed4_(int *n, int *i__, double *d__, 
	double *z__, double *delta, double *rho, double *dlam,
	 int *info);
 
/* Subroutine */ int dlaed5_(int *i__, double *d__, double *z__, 
	double *delta, double *rho, double *dlam);
 
/* Subroutine */ int dlaed6_(int *kniter, logical *orgati, double *
	rho, double *d__, double *z__, double *finit, double *
	tau, int *info);
 
/* Subroutine */ int dlaed7_(int *icompq, int *n, int *qsiz, 
	int *tlvls, int *curlvl, int *curpbm, double *d__, 
	double *q, int *ldq, int *indxq, double *rho, int 
	*cutpnt, double *qstore, int *qptr, int *prmptr, int *
	perm, int *givptr, int *givcol, double *givnum, 
	double *work, int *iwork, int *info);
 
/* Subroutine */ int dlaed8_(int *icompq, int *k, int *n, int 
	*qsiz, double *d__, double *q, int *ldq, int *indxq, 
	double *rho, int *cutpnt, double *z__, double *dlamda,
	 double *q2, int *ldq2, double *w, int *perm, int 
	*givptr, int *givcol, double *givnum, int *indxp, int 
	*indx, int *info);
 
/* Subroutine */ int dlaed9_(int *k, int *kstart, int *kstop, 
	int *n, double *d__, double *q, int *ldq, double *
	rho, double *dlamda, double *w, double *s, int *lds, 
	int *info);
 
/* Subroutine */ int dlaeda_(int *n, int *tlvls, int *curlvl, 
	int *curpbm, int *prmptr, int *perm, int *givptr, 
	int *givcol, double *givnum, double *q, int *qptr, 
	double *z__, double *ztemp, int *info);
 
/* Subroutine */ int dlaein_(logical *rightv, logical *noinit, int *n, 
	double *h__, int *ldh, double *wr, double *wi, 
	double *vr, double *vi, double *b, int *ldb, 
	double *work, double *eps3, double *smlnum, double *
	bignum, int *info);
 
/* Subroutine */ int dlaev2_(double *a, double *b, double *c__, 
	double *rt1, double *rt2, double *cs1, double *sn1);
 
/* Subroutine */ int dlaexc_(logical *wantq, int *n, double *t, 
	int *ldt, double *q, int *ldq, int *j1, int *n1, 
	int *n2, double *work, int *info);
 
/* Subroutine */ int dlag2_(double *a, int *lda, double *b, 
	int *ldb, double *safmin, double *scale1, double *
	scale2, double *wr1, double *wr2, double *wi);
 
/* Subroutine */ int dlags2_(logical *upper, double *a1, double *a2, 
	double *a3, double *b1, double *b2, double *b3, 
	double *csu, double *snu, double *csv, double *snv, 
	double *csq, double *snq);
 
/* Subroutine */ int dlagtf_(int *n, double *a, double *lambda, 
	double *b, double *c__, double *tol, double *d__, 
	int *in, int *info);
 
/* Subroutine */ int dlagtm_(const char *trans, int *n, int *nrhs, 
	double *alpha, double *dl, double *d__, double *du, 
	double *x, int *ldx, double *beta, double *b, int 
	*ldb);
 
/* Subroutine */ int dlagts_(int *job, int *n, double *a, 
	double *b, double *c__, double *d__, int *in, 
	double *y, double *tol, int *info);
 
/* Subroutine */ int dlagv2_(double *a, int *lda, double *b, 
	int *ldb, double *alphar, double *alphai, double *
	beta, double *csl, double *snl, double *csr, double *
	snr);
 
/* Subroutine */ int dlahqr_(logical *wantt, logical *wantz, int *n, 
	int *ilo, int *ihi, double *h__, int *ldh, double 
	*wr, double *wi, int *iloz, int *ihiz, double *z__, 
	int *ldz, int *info);
 
/* Subroutine */ int dlahrd_(int *n, int *k, int *nb, double *
	a, int *lda, double *tau, double *t, int *ldt, 
	double *y, int *ldy);
 
/* Subroutine */ int dlaic1_(int *job, int *j, double *x, 
	double *sest, double *w, double *gamma, double *
	sestpr, double *s, double *c__);
 
/* Subroutine */ int dlaln2_(logical *ltrans, int *na, int *nw, 
	double *smin, double *ca, double *a, int *lda, 
	double *d1, double *d2, double *b, int *ldb, 
	double *wr, double *wi, double *x, int *ldx, 
	double *scale, double *xnorm, int *info);
 
/* Subroutine */ int dlals0_(int *icompq, int *nl, int *nr, 
	int *sqre, int *nrhs, double *b, int *ldb, double 
	*bx, int *ldbx, int *perm, int *givptr, int *givcol, 
	int *ldgcol, double *givnum, int *ldgnum, double *
	poles, double *difl, double *difr, double *z__, int *
	k, double *c__, double *s, double *work, int *info);
 
/* Subroutine */ int dlalsa_(int *icompq, int *smlsiz, int *n, 
	int *nrhs, double *b, int *ldb, double *bx, int *
	ldbx, double *u, int *ldu, double *vt, int *k, 
	double *difl, double *difr, double *z__, double *
	poles, int *givptr, int *givcol, int *ldgcol, int *
	perm, double *givnum, double *c__, double *s, double *
	work, int *iwork, int *info);
 
/* Subroutine */ int dlalsd_(const char *uplo, int *smlsiz, int *n, int 
	*nrhs, double *d__, double *e, double *b, int *ldb, 
	double *rcond, int *rank, double *work, int *iwork, 
	int *info);
 
/* Subroutine */ int dlamc1_(int *beta, int *t, logical *rnd, logical 
	*ieee1);
 
/* Subroutine */ int dlamc2_(int *beta, int *t, logical *rnd, 
	double *eps, int *emin, double *rmin, int *emax, 
	double *rmax);
 
/* Subroutine */ int dlamc4_(int *emin, double *start, int *base);
 
/* Subroutine */ int dlamc5_(int *beta, int *p, int *emin, 
	logical *ieee, int *emax, double *rmax);
 
/* Subroutine */ int dlamrg_(int *n1, int *n2, double *a, int 
	*dtrd1, int *dtrd2, int *index);
 
/* Subroutine */ int dlanv2_(double *a, double *b, double *c__, 
	double *d__, double *rt1r, double *rt1i, double *rt2r,
	 double *rt2i, double *cs, double *sn);
 
/* Subroutine */ int dlapll_(int *n, double *x, int *incx, 
	double *y, int *incy, double *ssmin);
 
/* Subroutine */ int dlapmt_(logical *forwrd, int *m, int *n, 
	double *x, int *ldx, int *k);
 
/* Subroutine */ int dlaqgb_(int *m, int *n, int *kl, int *ku,
	 double *ab, int *ldab, double *r__, double *c__, 
	double *rowcnd, double *colcnd, double *amax, const char *equed);
 
/* Subroutine */ int dlaqge_(int *m, int *n, double *a, int *
	lda, double *r__, double *c__, double *rowcnd, double 
	*colcnd, double *amax, const char *equed);
 
/* Subroutine */ int dlaqp2_(int *m, int *n, int *offset, 
	double *a, int *lda, int *jpvt, double *tau, 
	double *vn1, double *vn2, double *work);
 
/* Subroutine */ int dlaqps_(int *m, int *n, int *offset, int 
	*nb, int *kb, double *a, int *lda, int *jpvt, 
	double *tau, double *vn1, double *vn2, double *auxv, 
	double *f, int *ldf);
 
/* Subroutine */ int dlaqsb_(const char *uplo, int *n, int *kd, double *
	ab, int *ldab, double *s, double *scond, double *amax,
	 const char *equed);
 
/* Subroutine */ int dlaqsp_(const char *uplo, int *n, double *ap, 
	double *s, double *scond, double *amax, const char *equed);
 
/* Subroutine */ int dlaqsy_(const char *uplo, int *n, double *a, int *
	lda, double *s, double *scond, double *amax, const char *equed);
 
/* Subroutine */ int dlaqtr_(logical *ltran, logical *lfloat, int *n, 
	double *t, int *ldt, double *b, double *w, double 
	*scale, double *x, double *work, int *info);
 
/* Subroutine */ int dlar1v_(int *n, int *b1, int *bn, double 
	*sigma, double *d__, double *l, double *ld, double *
	lld, double *gersch, double *z__, double *ztz, double 
	*mingma, int *r__, int *isuppz, double *work);
 
/* Subroutine */ int dlar2v_(int *n, double *x, double *y, 
	double *z__, int *incx, double *c__, double *s, 
	int *incc);
 
/* Subroutine */ int dlarf_(const char *side, int *m, int *n, double *v,
	 int *incv, double *tau, double *c__, int *ldc, 
	double *work);
 
/* Subroutine */ int dlarfb_(const char *side, const char *trans, const char *direct, const char *
	storev, int *m, int *n, int *k, double *v, int *
	ldv, double *t, int *ldt, double *c__, int *ldc, 
	double *work, int *ldwork);
 
/* Subroutine */ int dlarfg_(int *n, double *alpha, double *x, 
	int *incx, double *tau);
 
/* Subroutine */ int dlarft_(const char *direct, const char *storev, int *n, int *
	k, double *v, int *ldv, double *tau, double *t, 
	int *ldt);
 
/* Subroutine */ int dlarfx_(const char *side, int *m, int *n, double *
	v, double *tau, double *c__, int *ldc, double *work);
 
/* Subroutine */ int dlargv_(int *n, double *x, int *incx, 
	double *y, int *incy, double *c__, int *incc);
 
/* Subroutine */ int dlarnv_(int *idist, int *iseed, int *n, 
	double *x);
 
/* Subroutine */ int dlarrb_(int *n, double *d__, double *l, 
	double *ld, double *lld, int *ifirst, int *ilast, 
	double *sigma, double *reltol, double *w, double *
	wgap, double *werr, double *work, int *iwork, int *
	info);
 
/* Subroutine */ int dlarre_(int *n, double *d__, double *e, 
	double *tol, int *nsplit, int *isplit, int *m, 
	double *w, double *woff, double *gersch, double *work,
	 int *info);
 
/* Subroutine */ int dlarrf_(int *n, double *d__, double *l, 
	double *ld, double *lld, int *ifirst, int *ilast, 
	double *w, double *dplus, double *lplus, double *work,
	 int *iwork, int *info);
 
/* Subroutine */ int dlarrv_(int *n, double *d__, double *l, 
	int *isplit, int *m, double *w, int *iblock, 
	double *gersch, double *tol, double *z__, int *ldz, 
	int *isuppz, double *work, int *iwork, int *info);
 
/* Subroutine */ int dlartg_(double *f, double *g, double *cs, 
	double *sn, double *r__);
 
/* Subroutine */ int dlartv_(int *n, double *x, int *incx, 
	double *y, int *incy, double *c__, double *s, int 
	*incc);
 
/* Subroutine */ int dlaruv_(int *iseed, int *n, double *x);
 
/* Subroutine */ int dlarz_(const char *side, int *m, int *n, int *l, 
	double *v, int *incv, double *tau, double *c__, 
	int *ldc, double *work);
 
/* Subroutine */ int dlarzb_(const char *side, const char *trans, const char *direct, const char *
	storev, int *m, int *n, int *k, int *l, double *v,
	 int *ldv, double *t, int *ldt, double *c__, int *
	ldc, double *work, int *ldwork);
 
/* Subroutine */ int dlarzt_(const char *direct, const char *storev, int *n, int *
	k, double *v, int *ldv, double *tau, double *t, 
	int *ldt);
 
/* Subroutine */ int dlas2_(double *f, double *g, double *h__, 
	double *ssmin, double *ssmax);
 
/* Subroutine */ int dlascl_(const char *type__, int *kl, int *ku, 
	double *cfrom, double *cto, int *m, int *n, 
	double *a, int *lda, int *info);
 
/* Subroutine */ int dlasd0_(int *n, int *sqre, double *d__, 
	double *e, double *u, int *ldu, double *vt, int *
	ldvt, int *smlsiz, int *iwork, double *work, int *
	info);
 
/* Subroutine */ int dlasd1_(int *nl, int *nr, int *sqre, 
	double *d__, double *alpha, double *beta, double *u, 
	int *ldu, double *vt, int *ldvt, int *idxq, int *
	iwork, double *work, int *info);
 
/* Subroutine */ int dlasd2_(int *nl, int *nr, int *sqre, int 
	*k, double *d__, double *z__, double *alpha, double *
	beta, double *u, int *ldu, double *vt, int *ldvt, 
	double *dsigma, double *u2, int *ldu2, double *vt2, 
	int *ldvt2, int *idxp, int *idx, int *idxc, int *
	idxq, int *coltyp, int *info);
 
/* Subroutine */ int dlasd3_(int *nl, int *nr, int *sqre, int 
	*k, double *d__, double *q, int *ldq, double *dsigma, 
	double *u, int *ldu, double *u2, int *ldu2, 
	double *vt, int *ldvt, double *vt2, int *ldvt2, 
	int *idxc, int *ctot, double *z__, int *info);
 
/* Subroutine */ int dlasd4_(int *n, int *i__, double *d__, 
	double *z__, double *delta, double *rho, double *
	sigma, double *work, int *info);
 
/* Subroutine */ int dlasd5_(int *i__, double *d__, double *z__, 
	double *delta, double *rho, double *dsigma, double *
	work);
 
/* Subroutine */ int dlasd6_(int *icompq, int *nl, int *nr, 
	int *sqre, double *d__, double *vf, double *vl, 
	double *alpha, double *beta, int *idxq, int *perm, 
	int *givptr, int *givcol, int *ldgcol, double *givnum,
	 int *ldgnum, double *poles, double *difl, double *
	difr, double *z__, int *k, double *c__, double *s, 
	double *work, int *iwork, int *info);
 
/* Subroutine */ int dlasd7_(int *icompq, int *nl, int *nr, 
	int *sqre, int *k, double *d__, double *z__, 
	double *zw, double *vf, double *vfw, double *vl, 
	double *vlw, double *alpha, double *beta, double *
	dsigma, int *idx, int *idxp, int *idxq, int *perm, 
	int *givptr, int *givcol, int *ldgcol, double *givnum,
	 int *ldgnum, double *c__, double *s, int *info);
 
/* Subroutine */ int dlasd8_(int *icompq, int *k, double *d__, 
	double *z__, double *vf, double *vl, double *difl, 
	double *difr, int *lddifr, double *dsigma, double *
	work, int *info);
 
/* Subroutine */ int dlasd9_(int *icompq, int *ldu, int *k, 
	double *d__, double *z__, double *vf, double *vl, 
	double *difl, double *difr, double *dsigma, double *
	work, int *info);
 
/* Subroutine */ int dlasda_(int *icompq, int *smlsiz, int *n, 
	int *sqre, double *d__, double *e, double *u, int 
	*ldu, double *vt, int *k, double *difl, double *difr, 
	double *z__, double *poles, int *givptr, int *givcol, 
	int *ldgcol, int *perm, double *givnum, double *c__, 
	double *s, double *work, int *iwork, int *info);
 
/* Subroutine */ int dlasdq_(const char *uplo, int *sqre, int *n, int *
	ncvt, int *nru, int *ncc, double *d__, double *e, 
	double *vt, int *ldvt, double *u, int *ldu, 
	double *c__, int *ldc, double *work, int *info);
 
/* Subroutine */ int dlasdt_(int *n, int *lvl, int *nd, int *
	inode, int *ndiml, int *ndimr, int *msub);
 
/* Subroutine */ int dlaset_(const char *uplo, int *m, int *n, double *
	alpha, double *beta, double *a, int *lda);
 
/* Subroutine */ int dlasq1_(int *n, double *d__, double *e, 
	double *work, int *info);
 
/* Subroutine */ int dlasq2_(int *n, double *z__, int *info);
 
/* Subroutine */ int dlasq3_(int *i0, int *n0, double *z__, 
	int *pp, double *dmin__, double *sigma, double *desig,
	 double *qmax, int *nfail, int *iter, int *ndiv, 
	logical *ieee);
 
/* Subroutine */ int dlasq4_(int *i0, int *n0, double *z__, 
	int *pp, int *n0in, double *dmin__, double *dmin1, 
	double *dmin2, double *dn, double *dn1, double *dn2, 
	double *tau, int *ttype);
 
/* Subroutine */ int dlasq5_(int *i0, int *n0, double *z__, 
	int *pp, double *tau, double *dmin__, double *dmin1, 
	double *dmin2, double *dn, double *dnm1, double *dnm2,
	 logical *ieee);
 
/* Subroutine */ int dlasq6_(int *i0, int *n0, double *z__, 
	int *pp, double *dmin__, double *dmin1, double *dmin2,
	 double *dn, double *dnm1, double *dnm2);
 
/* Subroutine */ int dlasr_(const char *side, const char *pivot, const char *direct, int *m,
	 int *n, double *c__, double *s, double *a, int *
	lda);
 
/* Subroutine */ int dlasrt_(const char *id, int *n, double *d__, int *
	info);
 
/* Subroutine */ int dlassq_(int *n, double *x, int *incx, 
	double *scale, double *sumsq);
 
/* Subroutine */ int dlasv2_(double *f, double *g, double *h__, 
	double *ssmin, double *ssmax, double *snr, double *
	csr, double *snl, double *csl);
 
/* Subroutine */ int dlaswp_(int *n, double *a, int *lda, int 
	*k1, int *k2, int *ipiv, int *incx);
 
/* Subroutine */ int dlasy2_(logical *ltranl, logical *ltranr, int *isgn, 
	int *n1, int *n2, double *tl, int *ldtl, double *
	tr, int *ldtr, double *b, int *ldb, double *scale, 
	double *x, int *ldx, double *xnorm, int *info);
 
/* Subroutine */ int dlasyf_(const char *uplo, int *n, int *nb, int *kb,
	 double *a, int *lda, int *ipiv, double *w, int *
	ldw, int *info);
 
/* Subroutine */ int dlatbs_(const char *uplo, const char *trans, const char *diag, const char *
	normin, int *n, int *kd, double *ab, int *ldab, 
	double *x, double *scale, double *cnorm, int *info);
 
/* Subroutine */ int dlatdf_(int *ijob, int *n, double *z__, 
	int *ldz, double *rhs, double *rdsum, double *rdscal, 
	int *ipiv, int *jpiv);
 
/* Subroutine */ int dlatps_(const char *uplo, const char *trans, const char *diag, const char *
	normin, int *n, double *ap, double *x, double *scale, 
	double *cnorm, int *info);
 
/* Subroutine */ int dlatrd_(const char *uplo, int *n, int *nb, double *
	a, int *lda, double *e, double *tau, double *w, 
	int *ldw);
 
/* Subroutine */ int dlatrs_(const char *uplo, const char *trans, const char *diag, const char *
	normin, int *n, double *a, int *lda, double *x, 
	double *scale, double *cnorm, int *info);
 
/* Subroutine */ int dlatrz_(int *m, int *n, int *l, double *
	a, int *lda, double *tau, double *work);
 
/* Subroutine */ int dlatzm_(const char *side, int *m, int *n, double *
	v, int *incv, double *tau, double *c1, double *c2, 
	int *ldc, double *work);
 
/* Subroutine */ int dlauu2_(const char *uplo, int *n, double *a, int *
	lda, int *info);
 
/* Subroutine */ int dlauum_(const char *uplo, int *n, double *a, int *
	lda, int *info);
 
/* Subroutine */ int dopgtr_(const char *uplo, int *n, double *ap, 
	double *tau, double *q, int *ldq, double *work, 
	int *info);
 
/* Subroutine */ int dopmtr_(const char *side, const char *uplo, const char *trans, int *m, 
	int *n, double *ap, double *tau, double *c__, int 
	*ldc, double *work, int *info);
 
/* Subroutine */ int dorg2l_(int *m, int *n, int *k, double *
	a, int *lda, double *tau, double *work, int *info);
 
/* Subroutine */ int dorg2r_(int *m, int *n, int *k, double *
	a, int *lda, double *tau, double *work, int *info);
 
/* Subroutine */ int dorgbr_(const char *vect, int *m, int *n, int *k, 
	double *a, int *lda, double *tau, double *work, 
	int *lwork, int *info);
 
/* Subroutine */ int dorghr_(int *n, int *ilo, int *ihi, 
	double *a, int *lda, double *tau, double *work, 
	int *lwork, int *info);
 
/* Subroutine */ int dorgl2_(int *m, int *n, int *k, double *
	a, int *lda, double *tau, double *work, int *info);
 
/* Subroutine */ int dorglq_(int *m, int *n, int *k, double *
	a, int *lda, double *tau, double *work, int *lwork, 
	int *info);
 
/* Subroutine */ int dorgql_(int *m, int *n, int *k, double *
	a, int *lda, double *tau, double *work, int *lwork, 
	int *info);
 
/* Subroutine */ int dorgqr_(int *m, int *n, int *k, double *
	a, int *lda, double *tau, double *work, int *lwork, 
	int *info);
 
/* Subroutine */ int dorgr2_(int *m, int *n, int *k, double *
	a, int *lda, double *tau, double *work, int *info);
 
/* Subroutine */ int dorgrq_(int *m, int *n, int *k, double *
	a, int *lda, double *tau, double *work, int *lwork, 
	int *info);
 
/* Subroutine */ int dorgtr_(const char *uplo, int *n, double *a, int *
	lda, double *tau, double *work, int *lwork, int *info);
 
/* Subroutine */ int dorm2l_(const char *side, const char *trans, int *m, int *n, 
	int *k, double *a, int *lda, double *tau, double *
	c__, int *ldc, double *work, int *info);
 
/* Subroutine */ int dorm2r_(const char *side, const char *trans, int *m, int *n, 
	int *k, double *a, int *lda, double *tau, double *
	c__, int *ldc, double *work, int *info);
 
/* Subroutine */ int dormbr_(const char *vect, const char *side, const char *trans, int *m, 
	int *n, int *k, double *a, int *lda, double *tau, 
	double *c__, int *ldc, double *work, int *lwork, 
	int *info);
 
/* Subroutine */ int dormhr_(const char *side, const char *trans, int *m, int *n, 
	int *ilo, int *ihi, double *a, int *lda, double *
	tau, double *c__, int *ldc, double *work, int *lwork, 
	int *info);
 
/* Subroutine */ int dorml2_(const char *side, const char *trans, int *m, int *n, 
	int *k, double *a, int *lda, double *tau, double *
	c__, int *ldc, double *work, int *info);
 
/* Subroutine */ int dormlq_(const char *side, const char *trans, int *m, int *n, 
	int *k, double *a, int *lda, double *tau, double *
	c__, int *ldc, double *work, int *lwork, int *info);
 
/* Subroutine */ int dormql_(const char *side, const char *trans, int *m, int *n, 
	int *k, double *a, int *lda, double *tau, double *
	c__, int *ldc, double *work, int *lwork, int *info);
 
/* Subroutine */ int dormqr_(const char *side, const char *trans, int *m, int *n, 
	int *k, double *a, int *lda, double *tau, double *
	c__, int *ldc, double *work, int *lwork, int *info);
 
/* Subroutine */ int dormr2_(const char *side, const char *trans, int *m, int *n, 
	int *k, double *a, int *lda, double *tau, double *
	c__, int *ldc, double *work, int *info);
 
/* Subroutine */ int dormr3_(const char *side, const char *trans, int *m, int *n, 
	int *k, int *l, double *a, int *lda, double *tau, 
	double *c__, int *ldc, double *work, int *info);
 
/* Subroutine */ int dormrq_(const char *side, const char *trans, int *m, int *n, 
	int *k, double *a, int *lda, double *tau, double *
	c__, int *ldc, double *work, int *lwork, int *info);
 
/* Subroutine */ int dormrz_(const char *side, const char *trans, int *m, int *n, 
	int *k, int *l, double *a, int *lda, double *tau, 
	double *c__, int *ldc, double *work, int *lwork, 
	int *info);
 
/* Subroutine */ int dormtr_(const char *side, const char *uplo, const char *trans, int *m, 
	int *n, double *a, int *lda, double *tau, double *
	c__, int *ldc, double *work, int *lwork, int *info);
 
/* Subroutine */ int dpbcon_(const char *uplo, int *n, int *kd, double *
	ab, int *ldab, double *anorm, double *rcond, double *
	work, int *iwork, int *info);
 
/* Subroutine */ int dpbequ_(const char *uplo, int *n, int *kd, double *
	ab, int *ldab, double *s, double *scond, double *amax,
	 int *info);
 
/* Subroutine */ int dpbrfs_(const char *uplo, int *n, int *kd, int *
	nrhs, double *ab, int *ldab, double *afb, int *ldafb, 
	double *b, int *ldb, double *x, int *ldx, double *
	ferr, double *berr, double *work, int *iwork, int *
	info);
 
/* Subroutine */ int dpbstf_(const char *uplo, int *n, int *kd, double *
	ab, int *ldab, int *info);
 
/* Subroutine */ int dpbsv_(const char *uplo, int *n, int *kd, int *
	nrhs, double *ab, int *ldab, double *b, int *ldb, 
	int *info);
 
/* Subroutine */ int dpbsvx_(const char *fact, const char *uplo, int *n, int *kd, 
	int *nrhs, double *ab, int *ldab, double *afb, 
	int *ldafb, const char *equed, double *s, double *b, int *
	ldb, double *x, int *ldx, double *rcond, double *ferr,
	 double *berr, double *work, int *iwork, int *info);
 
/* Subroutine */ int dpbtf2_(const char *uplo, int *n, int *kd, double *
	ab, int *ldab, int *info);
 
/* Subroutine */ int dpbtrf_(const char *uplo, int *n, int *kd, double *
	ab, int *ldab, int *info);
 
/* Subroutine */ int dpbtrs_(const char *uplo, int *n, int *kd, int *
	nrhs, double *ab, int *ldab, double *b, int *ldb, 
	int *info);
 
/* Subroutine */ int dpocon_(const char *uplo, int *n, double *a, int *
	lda, double *anorm, double *rcond, double *work, int *
	iwork, int *info);
 
/* Subroutine */ int dpoequ_(int *n, double *a, int *lda, 
	double *s, double *scond, double *amax, int *info);
 
/* Subroutine */ int dporfs_(const char *uplo, int *n, int *nrhs, 
	double *a, int *lda, double *af, int *ldaf, 
	double *b, int *ldb, double *x, int *ldx, double *
	ferr, double *berr, double *work, int *iwork, int *
	info);
 
/* Subroutine */ int dposv_(const char *uplo, int *n, int *nrhs, double 
	*a, int *lda, double *b, int *ldb, int *info);
 
/* Subroutine */ int dposvx_(const char *fact, const char *uplo, int *n, int *
	nrhs, double *a, int *lda, double *af, int *ldaf, 
	const char *equed, double *s, double *b, int *ldb, double *
	x, int *ldx, double *rcond, double *ferr, double *
	berr, double *work, int *iwork, int *info);
 
/* Subroutine */ int dpotf2_(const char *uplo, int *n, double *a, int *
	lda, int *info);
 
/* Subroutine */ int dpotrf_(const char *uplo, int *n, double *a, int *
	lda, int *info);
 
/* Subroutine */ int dpotri_(const char *uplo, int *n, double *a, int *
	lda, int *info);
 
/* Subroutine */ int dpotrs_(const char *uplo, int *n, int *nrhs, 
	double *a, int *lda, double *b, int *ldb, int *
	info);
 
/* Subroutine */ int dppcon_(const char *uplo, int *n, double *ap, 
	double *anorm, double *rcond, double *work, int *
	iwork, int *info);
 
/* Subroutine */ int dppequ_(const char *uplo, int *n, double *ap, 
	double *s, double *scond, double *amax, int *info);
 
/* Subroutine */ int dpprfs_(const char *uplo, int *n, int *nrhs, 
	double *ap, double *afp, double *b, int *ldb, 
	double *x, int *ldx, double *ferr, double *berr, 
	double *work, int *iwork, int *info);
 
/* Subroutine */ int dppsv_(const char *uplo, int *n, int *nrhs, double 
	*ap, double *b, int *ldb, int *info);
 
/* Subroutine */ int dppsvx_(const char *fact, const char *uplo, int *n, int *
	nrhs, double *ap, double *afp, const char *equed, double *s, 
	double *b, int *ldb, double *x, int *ldx, double *
	rcond, double *ferr, double *berr, double *work, int *
	iwork, int *info);
 
/* Subroutine */ int dpptrf_(const char *uplo, int *n, double *ap, int *
	info);
 
/* Subroutine */ int dpptri_(const char *uplo, int *n, double *ap, int *
	info);
 
/* Subroutine */ int dpptrs_(const char *uplo, int *n, int *nrhs, 
	double *ap, double *b, int *ldb, int *info);
 
/* Subroutine */ int dptcon_(int *n, double *d__, double *e, 
	double *anorm, double *rcond, double *work, int *info);
 
/* Subroutine */ int dpteqr_(const char *compz, int *n, double *d__, 
	double *e, double *z__, int *ldz, double *work, 
	int *info);
 
/* Subroutine */ int dptrfs_(int *n, int *nrhs, double *d__, 
	double *e, double *df, double *ef, double *b, int 
	*ldb, double *x, int *ldx, double *ferr, double *berr,
	 double *work, int *info);
 
/* Subroutine */ int dptsv_(int *n, int *nrhs, double *d__, 
	double *e, double *b, int *ldb, int *info);
 
/* Subroutine */ int dptsvx_(const char *fact, int *n, int *nrhs, 
	double *d__, double *e, double *df, double *ef, 
	double *b, int *ldb, double *x, int *ldx, double *
	rcond, double *ferr, double *berr, double *work, int *
	info);
 
/* Subroutine */ int dpttrf_(int *n, double *d__, double *e, 
	int *info);
 
/* Subroutine */ int dpttrs_(int *n, int *nrhs, double *d__, 
	double *e, double *b, int *ldb, int *info);
 
/* Subroutine */ int dptts2_(int *n, int *nrhs, double *d__, 
	double *e, double *b, int *ldb);
 
/* Subroutine */ int drscl_(int *n, double *sa, double *sx, 
	int *incx);
 
/* Subroutine */ int dsbev_(const char *jobz, const char *uplo, int *n, int *kd, 
	double *ab, int *ldab, double *w, double *z__, 
	int *ldz, double *work, int *info);
 
/* Subroutine */ int dsbevd_(const char *jobz, const char *uplo, int *n, int *kd, 
	double *ab, int *ldab, double *w, double *z__, 
	int *ldz, double *work, int *lwork, int *iwork, 
	int *liwork, int *info);
 
/* Subroutine */ int dsbevx_(const char *jobz, const char *range, const char *uplo, int *n, 
	int *kd, double *ab, int *ldab, double *q, int *
	ldq, double *vl, double *vu, int *il, int *iu, 
	double *abstol, int *m, double *w, double *z__, 
	int *ldz, double *work, int *iwork, int *ifail, 
	int *info);
 
/* Subroutine */ int dsbgst_(const char *vect, const char *uplo, int *n, int *ka, 
	int *kb, double *ab, int *ldab, double *bb, int *
	ldbb, double *x, int *ldx, double *work, int *info);
 
/* Subroutine */ int dsbgv_(const char *jobz, const char *uplo, int *n, int *ka, 
	int *kb, double *ab, int *ldab, double *bb, int *
	ldbb, double *w, double *z__, int *ldz, double *work, 
	int *info);
 
/* Subroutine */ int dsbgvd_(const char *jobz, const char *uplo, int *n, int *ka, 
	int *kb, double *ab, int *ldab, double *bb, int *
	ldbb, double *w, double *z__, int *ldz, double *work, 
	int *lwork, int *iwork, int *liwork, int *info);
 
/* Subroutine */ int dsbgvx_(const char *jobz, const char *range, const char *uplo, int *n, 
	int *ka, int *kb, double *ab, int *ldab, double *
	bb, int *ldbb, double *q, int *ldq, double *vl, 
	double *vu, int *il, int *iu, double *abstol, int 
	*m, double *w, double *z__, int *ldz, double *work, 
	int *iwork, int *ifail, int *info);
 
/* Subroutine */ int dsbtrd_(const char *vect, const char *uplo, int *n, int *kd, 
	double *ab, int *ldab, double *d__, double *e, 
	double *q, int *ldq, double *work, int *info);
 
/* Subroutine */ int dspcon_(const char *uplo, int *n, double *ap, int *
	ipiv, double *anorm, double *rcond, double *work, int 
	*iwork, int *info);
 
/* Subroutine */ int dspev_(const char *jobz, const char *uplo, int *n, double *
	ap, double *w, double *z__, int *ldz, double *work, 
	int *info);
 
/* Subroutine */ int dspevd_(const char *jobz, const char *uplo, int *n, double *
	ap, double *w, double *z__, int *ldz, double *work, 
	int *lwork, int *iwork, int *liwork, int *info);
 
/* Subroutine */ int dspevx_(const char *jobz, const char *range, const char *uplo, int *n, 
	double *ap, double *vl, double *vu, int *il, int *
	iu, double *abstol, int *m, double *w, double *z__, 
	int *ldz, double *work, int *iwork, int *ifail, 
	int *info);
 
/* Subroutine */ int dspgst_(int *itype, const char *uplo, int *n, 
	double *ap, double *bp, int *info);
 
/* Subroutine */ int dspgv_(int *itype, const char *jobz, const char *uplo, int *
	n, double *ap, double *bp, double *w, double *z__, 
	int *ldz, double *work, int *info);
 
/* Subroutine */ int dspgvd_(int *itype, const char *jobz, const char *uplo, int *
	n, double *ap, double *bp, double *w, double *z__, 
	int *ldz, double *work, int *lwork, int *iwork, 
	int *liwork, int *info);
 
/* Subroutine */ int dspgvx_(int *itype, const char *jobz, const char *range, const char *
	uplo, int *n, double *ap, double *bp, double *vl, 
	double *vu, int *il, int *iu, double *abstol, int 
	*m, double *w, double *z__, int *ldz, double *work, 
	int *iwork, int *ifail, int *info);
 
/* Subroutine */ int dsprfs_(const char *uplo, int *n, int *nrhs, 
	double *ap, double *afp, int *ipiv, double *b, 
	int *ldb, double *x, int *ldx, double *ferr, 
	double *berr, double *work, int *iwork, int *info);
 
/* Subroutine */ int dspsv_(const char *uplo, int *n, int *nrhs, double 
	*ap, int *ipiv, double *b, int *ldb, int *info);
 
/* Subroutine */ int dspsvx_(const char *fact, const char *uplo, int *n, int *
	nrhs, double *ap, double *afp, int *ipiv, double *b, 
	int *ldb, double *x, int *ldx, double *rcond, 
	double *ferr, double *berr, double *work, int *iwork, 
	int *info);
 
/* Subroutine */ int dsptrd_(const char *uplo, int *n, double *ap, 
	double *d__, double *e, double *tau, int *info);
 
/* Subroutine */ int dsptrf_(const char *uplo, int *n, double *ap, int *
	ipiv, int *info);
 
/* Subroutine */ int dsptri_(const char *uplo, int *n, double *ap, int *
	ipiv, double *work, int *info);
 
/* Subroutine */ int dsptrs_(const char *uplo, int *n, int *nrhs, 
	double *ap, int *ipiv, double *b, int *ldb, int *
	info);
 
/* Subroutine */ int dstebz_(const char *range, const char *order, int *n, double 
	*vl, double *vu, int *il, int *iu, double *abstol, 
	double *d__, double *e, int *m, int *nsplit, 
	double *w, int *iblock, int *isplit, double *work, 
	int *iwork, int *info);
 
/* Subroutine */ int dstedc_(const char *compz, int *n, double *d__, 
	double *e, double *z__, int *ldz, double *work, 
	int *lwork, int *iwork, int *liwork, int *info);
 
/* Subroutine */ int dstegr_(const char *jobz, const char *range, int *n, double *
	d__, double *e, double *vl, double *vu, int *il, 
	int *iu, double *abstol, int *m, double *w, 
	double *z__, int *ldz, int *isuppz, double *work, 
	int *lwork, int *iwork, int *liwork, int *info);
 
/* Subroutine */ int dstein_(int *n, double *d__, double *e, 
	int *m, double *w, int *iblock, int *isplit, 
	double *z__, int *ldz, double *work, int *iwork, 
	int *ifail, int *info);
 
/* Subroutine */ int dsteqr_(const char *compz, int *n, double *d__, 
	double *e, double *z__, int *ldz, double *work, 
	int *info);
 
/* Subroutine */ int dsterf_(int *n, double *d__, double *e, 
	int *info);
 
/* Subroutine */ int dstev_(const char *jobz, int *n, double *d__, 
	double *e, double *z__, int *ldz, double *work, 
	int *info);
 
/* Subroutine */ int dstevd_(const char *jobz, int *n, double *d__, 
	double *e, double *z__, int *ldz, double *work, 
	int *lwork, int *iwork, int *liwork, int *info);
 
/* Subroutine */ int dstevr_(const char *jobz, const char *range, int *n, double *
	d__, double *e, double *vl, double *vu, int *il, 
	int *iu, double *abstol, int *m, double *w, 
	double *z__, int *ldz, int *isuppz, double *work, 
	int *lwork, int *iwork, int *liwork, int *info);
 
/* Subroutine */ int dstevx_(const char *jobz, const char *range, int *n, double *
	d__, double *e, double *vl, double *vu, int *il, 
	int *iu, double *abstol, int *m, double *w, 
	double *z__, int *ldz, double *work, int *iwork, 
	int *ifail, int *info);
 
/* Subroutine */ int dsycon_(const char *uplo, int *n, double *a, int *
	lda, int *ipiv, double *anorm, double *rcond, double *
	work, int *iwork, int *info);
 
/* Subroutine */ int dsyev_(const char *jobz, const char *uplo, int *n, double *a,
	 int *lda, double *w, double *work, int *lwork, 
	int *info);
 
/* Subroutine */ int dsyevd_(const char *jobz, const char *uplo, int *n, double *
	a, int *lda, double *w, double *work, int *lwork, 
	int *iwork, int *liwork, int *info);
 
/* Subroutine */ int dsyevr_(const char *jobz, const char *range, const char *uplo, int *n, 
	double *a, int *lda, double *vl, double *vu, int *
	il, int *iu, double *abstol, int *m, double *w, 
	double *z__, int *ldz, int *isuppz, double *work, 
	int *lwork, int *iwork, int *liwork, int *info);
 
/* Subroutine */ int dsyevx_(const char *jobz, const char *range, const char *uplo, int *n, 
	double *a, int *lda, double *vl, double *vu, int *
	il, int *iu, double *abstol, int *m, double *w, 
	double *z__, int *ldz, double *work, int *lwork, 
	int *iwork, int *ifail, int *info);
 
/* Subroutine */ int dsygs2_(int *itype, const char *uplo, int *n, 
	double *a, int *lda, double *b, int *ldb, int *
	info);
 
/* Subroutine */ int dsygst_(int *itype, const char *uplo, int *n, 
	double *a, int *lda, double *b, int *ldb, int *
	info);
 
/* Subroutine */ int dsygv_(int *itype, const char *jobz, const char *uplo, int *
	n, double *a, int *lda, double *b, int *ldb, 
	double *w, double *work, int *lwork, int *info);
 
/* Subroutine */ int dsygvd_(int *itype, const char *jobz, const char *uplo, int *
	n, double *a, int *lda, double *b, int *ldb, 
	double *w, double *work, int *lwork, int *iwork, 
	int *liwork, int *info);
 
/* Subroutine */ int dsygvx_(int *itype, const char *jobz, const char *range, const char *
	uplo, int *n, double *a, int *lda, double *b, int 
	*ldb, double *vl, double *vu, int *il, int *iu, 
	double *abstol, int *m, double *w, double *z__, 
	int *ldz, double *work, int *lwork, int *iwork, 
	int *ifail, int *info);
 
/* Subroutine */ int dsyrfs_(const char *uplo, int *n, int *nrhs, 
	double *a, int *lda, double *af, int *ldaf, int *
	ipiv, double *b, int *ldb, double *x, int *ldx, 
	double *ferr, double *berr, double *work, int *iwork, 
	int *info);
 
/* Subroutine */ int dsysv_(const char *uplo, int *n, int *nrhs, double 
	*a, int *lda, int *ipiv, double *b, int *ldb, 
	double *work, int *lwork, int *info);
 
/* Subroutine */ int dsysvx_(const char *fact, const char *uplo, int *n, int *
	nrhs, double *a, int *lda, double *af, int *ldaf, 
	int *ipiv, double *b, int *ldb, double *x, int *
	ldx, double *rcond, double *ferr, double *berr, 
	double *work, int *lwork, int *iwork, int *info);
 
/* Subroutine */ int dsytd2_(const char *uplo, int *n, double *a, int *
	lda, double *d__, double *e, double *tau, int *info);
 
/* Subroutine */ int dsytf2_(const char *uplo, int *n, double *a, int *
	lda, int *ipiv, int *info);
 
/* Subroutine */ int dsytrd_(const char *uplo, int *n, double *a, int *
	lda, double *d__, double *e, double *tau, double *
	work, int *lwork, int *info);
 
/* Subroutine */ int dsytrf_(const char *uplo, int *n, double *a, int *
	lda, int *ipiv, double *work, int *lwork, int *info);
 
/* Subroutine */ int dsytri_(const char *uplo, int *n, double *a, int *
	lda, int *ipiv, double *work, int *info);
 
/* Subroutine */ int dsytrs_(const char *uplo, int *n, int *nrhs, 
	double *a, int *lda, int *ipiv, double *b, int *
	ldb, int *info);
 
/* Subroutine */ int dtbcon_(const char *norm, const char *uplo, const char *diag, int *n, 
	int *kd, double *ab, int *ldab, double *rcond, 
	double *work, int *iwork, int *info);
 
/* Subroutine */ int dtbrfs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *kd, int *nrhs, double *ab, int *ldab, double 
	*b, int *ldb, double *x, int *ldx, double *ferr, 
	double *berr, double *work, int *iwork, int *info);
 
/* Subroutine */ int dtbtrs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *kd, int *nrhs, double *ab, int *ldab, double 
	*b, int *ldb, int *info);
 
/* Subroutine */ int dtgevc_(const char *side, const char *howmny, logical *select, 
	int *n, double *a, int *lda, double *b, int *ldb, 
	double *vl, int *ldvl, double *vr, int *ldvr, int 
	*mm, int *m, double *work, int *info);
 
/* Subroutine */ int dtgex2_(logical *wantq, logical *wantz, int *n, 
	double *a, int *lda, double *b, int *ldb, double *
	q, int *ldq, double *z__, int *ldz, int *j1, int *
	n1, int *n2, double *work, int *lwork, int *info);
 
/* Subroutine */ int dtgexc_(logical *wantq, logical *wantz, int *n, 
	double *a, int *lda, double *b, int *ldb, double *
	q, int *ldq, double *z__, int *ldz, int *ifst, 
	int *ilst, double *work, int *lwork, int *info);
 
/* Subroutine */ int dtgsen_(int *ijob, logical *wantq, logical *wantz, 
	logical *select, int *n, double *a, int *lda, double *
	b, int *ldb, double *alphar, double *alphai, double *
	beta, double *q, int *ldq, double *z__, int *ldz, 
	int *m, double *pl, double *pr, double *dif, 
	double *work, int *lwork, int *iwork, int *liwork, 
	int *info);
 
/* Subroutine */ int dtgsja_(const char *jobu, const char *jobv, const char *jobq, int *m, 
	int *p, int *n, int *k, int *l, double *a, 
	int *lda, double *b, int *ldb, double *tola, 
	double *tolb, double *alpha, double *beta, double *u, 
	int *ldu, double *v, int *ldv, double *q, int *
	ldq, double *work, int *ncycle, int *info);
 
/* Subroutine */ int dtgsna_(const char *job, const char *howmny, logical *select, 
	int *n, double *a, int *lda, double *b, int *ldb, 
	double *vl, int *ldvl, double *vr, int *ldvr, 
	double *s, double *dif, int *mm, int *m, double *
	work, int *lwork, int *iwork, int *info);
 
/* Subroutine */ int dtgsy2_(const char *trans, int *ijob, int *m, int *
	n, double *a, int *lda, double *b, int *ldb, 
	double *c__, int *ldc, double *d__, int *ldd, 
	double *e, int *lde, double *f, int *ldf, double *
	scale, double *rdsum, double *rdscal, int *iwork, int 
	*pq, int *info);
 
/* Subroutine */ int dtgsyl_(const char *trans, int *ijob, int *m, int *
	n, double *a, int *lda, double *b, int *ldb, 
	double *c__, int *ldc, double *d__, int *ldd, 
	double *e, int *lde, double *f, int *ldf, double *
	scale, double *dif, double *work, int *lwork, int *
	iwork, int *info);
 
/* Subroutine */ int dtpcon_(const char *norm, const char *uplo, const char *diag, int *n, 
	double *ap, double *rcond, double *work, int *iwork, 
	int *info);
 
/* Subroutine */ int dtprfs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *nrhs, double *ap, double *b, int *ldb, 
	double *x, int *ldx, double *ferr, double *berr, 
	double *work, int *iwork, int *info);
 
/* Subroutine */ int dtptri_(const char *uplo, const char *diag, int *n, double *
	ap, int *info);
 
/* Subroutine */ int dtptrs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *nrhs, double *ap, double *b, int *ldb, int *
	info);
 
/* Subroutine */ int dtrcon_(const char *norm, const char *uplo, const char *diag, int *n, 
	double *a, int *lda, double *rcond, double *work, 
	int *iwork, int *info);
 
/* Subroutine */ int dtrevc_(const char *side, const char *howmny, logical *select, 
	int *n, double *t, int *ldt, double *vl, int *
	ldvl, double *vr, int *ldvr, int *mm, int *m, 
	double *work, int *info);
 
/* Subroutine */ int dtrexc_(const char *compq, int *n, double *t, int *
	ldt, double *q, int *ldq, int *ifst, int *ilst, 
	double *work, int *info);
 
/* Subroutine */ int dtrrfs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *nrhs, double *a, int *lda, double *b, int *
	ldb, double *x, int *ldx, double *ferr, double *berr, 
	double *work, int *iwork, int *info);
 
/* Subroutine */ int dtrsen_(const char *job, const char *compq, logical *select, int 
	*n, double *t, int *ldt, double *q, int *ldq, 
	double *wr, double *wi, int *m, double *s, double 
	*sep, double *work, int *lwork, int *iwork, int *
	liwork, int *info);
 
/* Subroutine */ int dtrsna_(const char *job, const char *howmny, logical *select, 
	int *n, double *t, int *ldt, double *vl, int *
	ldvl, double *vr, int *ldvr, double *s, double *sep, 
	int *mm, int *m, double *work, int *ldwork, int *
	iwork, int *info);
 
/* Subroutine */ int dtrsyl_(const char *trana, const char *tranb, int *isgn, int 
	*m, int *n, double *a, int *lda, double *b, int *
	ldb, double *c__, int *ldc, double *scale, int *info);
 
/* Subroutine */ int dtrti2_(const char *uplo, const char *diag, int *n, double *
	a, int *lda, int *info);
 
/* Subroutine */ int dtrtri_(const char *uplo, const char *diag, int *n, double *
	a, int *lda, int *info);
 
/* Subroutine */ int dtrtrs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *nrhs, double *a, int *lda, double *b, int *
	ldb, int *info);
 
/* Subroutine */ int dtzrqf_(int *m, int *n, double *a, int *
	lda, double *tau, int *info);
 
/* Subroutine */ int dtzrzf_(int *m, int *n, double *a, int *
	lda, double *tau, double *work, int *lwork, int *info);
 
int icmax1_(int *n, cfloat *cx, int *incx);
 
int ieeeck_(int *ispec, float *zero, float *one);
 
int ilaenv_(int *ispec, const char *name__, const char *opts, int *n1, 
	int *n2, int *n3, int *n4, ftnlen name_len, ftnlen 
	opts_len);
 
int izmax1_(int *n, cdouble *cx, int *incx);
 
/* Subroutine */ int sbdsdc_(const char *uplo, const char *compq, int *n, float *d__, 
	float *e, float *u, int *ldu, float *vt, int *ldvt, float *q, 
	int *iq, float *work, int *iwork, int *info);
 
/* Subroutine */ int sbdsqr_(const char *uplo, int *n, int *ncvt, int *
	nru, int *ncc, float *d__, float *e, float *vt, int *ldvt, float *
	u, int *ldu, float *c__, int *ldc, float *work, int *info);
 
/* Subroutine */ int sdisna_(const char *job, int *m, int *n, float *d__, 
	float *sep, int *info);
 
/* Subroutine */ int sgbbrd_(const char *vect, int *m, int *n, int *ncc,
	 int *kl, int *ku, float *ab, int *ldab, float *d__, float *
	e, float *q, int *ldq, float *pt, int *ldpt, float *c__, int 
	*ldc, float *work, int *info);
 
/* Subroutine */ int sgbcon_(const char *norm, int *n, int *kl, int *ku,
	 float *ab, int *ldab, int *ipiv, float *anorm, float *rcond, 
	float *work, int *iwork, int *info);
 
/* Subroutine */ int sgbequ_(int *m, int *n, int *kl, int *ku,
	 float *ab, int *ldab, float *r__, float *c__, float *rowcnd, float *
	colcnd, float *amax, int *info);
 
/* Subroutine */ int sgbrfs_(const char *trans, int *n, int *kl, int *
	ku, int *nrhs, float *ab, int *ldab, float *afb, int *ldafb,
	 int *ipiv, float *b, int *ldb, float *x, int *ldx, float *
	ferr, float *berr, float *work, int *iwork, int *info);
 
/* Subroutine */ int sgbsv_(int *n, int *kl, int *ku, int *
	nrhs, float *ab, int *ldab, int *ipiv, float *b, int *ldb, 
	int *info);
 
/* Subroutine */ int sgbsvx_(const char *fact, const char *trans, int *n, int *kl,
	 int *ku, int *nrhs, float *ab, int *ldab, float *afb, 
	int *ldafb, int *ipiv, const char *equed, float *r__, float *c__, 
	float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr,
	 float *berr, float *work, int *iwork, int *info);
 
/* Subroutine */ int sgbtf2_(int *m, int *n, int *kl, int *ku,
	 float *ab, int *ldab, int *ipiv, int *info);
 
/* Subroutine */ int sgbtrf_(int *m, int *n, int *kl, int *ku,
	 float *ab, int *ldab, int *ipiv, int *info);
 
/* Subroutine */ int sgbtrs_(const char *trans, int *n, int *kl, int *
	ku, int *nrhs, float *ab, int *ldab, int *ipiv, float *b, 
	int *ldb, int *info);
 
/* Subroutine */ int sgebak_(const char *job, const char *side, int *n, int *ilo, 
	int *ihi, float *scale, int *m, float *v, int *ldv, int 
	*info);
 
/* Subroutine */ int sgebal_(const char *job, int *n, float *a, int *lda, 
	int *ilo, int *ihi, float *scale, int *info);
 
/* Subroutine */ int sgebd2_(int *m, int *n, float *a, int *lda, 
	float *d__, float *e, float *tauq, float *taup, float *work, int *info);
 
/* Subroutine */ int sgebrd_(int *m, int *n, float *a, int *lda, 
	float *d__, float *e, float *tauq, float *taup, float *work, int *
	lwork, int *info);
 
/* Subroutine */ int sgecon_(const char *norm, int *n, float *a, int *lda, 
	float *anorm, float *rcond, float *work, int *iwork, int *info);
 
/* Subroutine */ int sgeequ_(int *m, int *n, float *a, int *lda, 
	float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, int 
	*info);
 
/* Subroutine */ int sgees_(const char *jobvs, const char *sort, L_fp select, int *n, 
	float *a, int *lda, int *sdim, float *wr, float *wi, float *vs, 
	int *ldvs, float *work, int *lwork, logical *bwork, int *
	info);
 
/* Subroutine */ int sgeesx_(const char *jobvs, const char *sort, L_fp select, const char *
	sense, int *n, float *a, int *lda, int *sdim, float *wr, 
	float *wi, float *vs, int *ldvs, float *rconde, float *rcondv, float *
	work, int *lwork, int *iwork, int *liwork, logical *bwork,
	 int *info);
 
/* Subroutine */ int sgeev_(const char *jobvl, const char *jobvr, int *n, float *a, 
	int *lda, float *wr, float *wi, float *vl, int *ldvl, float *vr, 
	int *ldvr, float *work, int *lwork, int *info);
 
/* Subroutine */ int sgeevx_(const char *balanc, const char *jobvl, const char *jobvr, const char *
	sense, int *n, float *a, int *lda, float *wr, float *wi, float *
	vl, int *ldvl, float *vr, int *ldvr, int *ilo, int *
	ihi, float *scale, float *abnrm, float *rconde, float *rcondv, float *work,
	 int *lwork, int *iwork, int *info);
 
/* Subroutine */ int sgegs_(const char *jobvsl, const char *jobvsr, int *n, float *a, 
	int *lda, float *b, int *ldb, float *alphar, float *alphai, float 
	*beta, float *vsl, int *ldvsl, float *vsr, int *ldvsr, float *
	work, int *lwork, int *info);
 
/* Subroutine */ int sgegv_(const char *jobvl, const char *jobvr, int *n, float *a, 
	int *lda, float *b, int *ldb, float *alphar, float *alphai, float 
	*beta, float *vl, int *ldvl, float *vr, int *ldvr, float *work, 
	int *lwork, int *info);
 
/* Subroutine */ int sgehd2_(int *n, int *ilo, int *ihi, float *a, 
	int *lda, float *tau, float *work, int *info);
 
/* Subroutine */ int sgehrd_(int *n, int *ilo, int *ihi, float *a, 
	int *lda, float *tau, float *work, int *lwork, int *info);
 
/* Subroutine */ int sgelq2_(int *m, int *n, float *a, int *lda, 
	float *tau, float *work, int *info);
 
/* Subroutine */ int sgelqf_(int *m, int *n, float *a, int *lda, 
	float *tau, float *work, int *lwork, int *info);
 
/* Subroutine */ int sgels_(const char *trans, int *m, int *n, int *
	nrhs, float *a, int *lda, float *b, int *ldb, float *work, 
	int *lwork, int *info);
 
/* Subroutine */ int sgelsd_(int *m, int *n, int *nrhs, float *a, 
	int *lda, float *b, int *ldb, float *s, float *rcond, int *
	rank, float *work, int *lwork, int *iwork, int *info);
 
/* Subroutine */ int sgelss_(int *m, int *n, int *nrhs, float *a, 
	int *lda, float *b, int *ldb, float *s, float *rcond, int *
	rank, float *work, int *lwork, int *info);
 
/* Subroutine */ int sgelsx_(int *m, int *n, int *nrhs, float *a, 
	int *lda, float *b, int *ldb, int *jpvt, float *rcond, 
	int *rank, float *work, int *info);
 
/* Subroutine */ int sgelsy_(int *m, int *n, int *nrhs, float *a, 
	int *lda, float *b, int *ldb, int *jpvt, float *rcond, 
	int *rank, float *work, int *lwork, int *info);
 
/* Subroutine */ int sgeql2_(int *m, int *n, float *a, int *lda, 
	float *tau, float *work, int *info);
 
/* Subroutine */ int sgeqlf_(int *m, int *n, float *a, int *lda, 
	float *tau, float *work, int *lwork, int *info);
 
/* Subroutine */ int sgeqp3_(int *m, int *n, float *a, int *lda, 
	int *jpvt, float *tau, float *work, int *lwork, int *info);
 
/* Subroutine */ int sgeqpf_(int *m, int *n, float *a, int *lda, 
	int *jpvt, float *tau, float *work, int *info);
 
/* Subroutine */ int sgeqr2_(int *m, int *n, float *a, int *lda, 
	float *tau, float *work, int *info);
 
/* Subroutine */ int sgeqrf_(int *m, int *n, float *a, int *lda, 
	float *tau, float *work, int *lwork, int *info);
 
/* Subroutine */ int sgerfs_(const char *trans, int *n, int *nrhs, float *a, 
	int *lda, float *af, int *ldaf, int *ipiv, float *b, 
	int *ldb, float *x, int *ldx, float *ferr, float *berr, float *
	work, int *iwork, int *info);
 
/* Subroutine */ int sgerq2_(int *m, int *n, float *a, int *lda, 
	float *tau, float *work, int *info);
 
/* Subroutine */ int sgerqf_(int *m, int *n, float *a, int *lda, 
	float *tau, float *work, int *lwork, int *info);
 
/* Subroutine */ int sgesc2_(int *n, float *a, int *lda, float *rhs, 
	int *ipiv, int *jpiv, float *scale);
 
/* Subroutine */ int sgesdd_(const char *jobz, int *m, int *n, float *a, 
	int *lda, float *s, float *u, int *ldu, float *vt, int *ldvt,
	 float *work, int *lwork, int *iwork, int *info);
 
/* Subroutine */ int sgesv_(int *n, int *nrhs, float *a, int *lda, 
	int *ipiv, float *b, int *ldb, int *info);
 
/* Subroutine */ int sgesvd_(const char *jobu, const char *jobvt, int *m, int *n, 
	float *a, int *lda, float *s, float *u, int *ldu, float *vt, 
	int *ldvt, float *work, int *lwork, int *info);
 
/* Subroutine */ int sgesvx_(const char *fact, const char *trans, int *n, int *
	nrhs, float *a, int *lda, float *af, int *ldaf, int *ipiv, 
	const char *equed, float *r__, float *c__, float *b, int *ldb, float *x, 
	int *ldx, float *rcond, float *ferr, float *berr, float *work, 
	int *iwork, int *info);
 
/* Subroutine */ int sgetc2_(int *n, float *a, int *lda, int *ipiv,
	 int *jpiv, int *info);
 
/* Subroutine */ int sgetf2_(int *m, int *n, float *a, int *lda, 
	int *ipiv, int *info);
 
/* Subroutine */ int sgetrf_(int *m, int *n, float *a, int *lda, 
	int *ipiv, int *info);
 
/* Subroutine */ int sgetri_(int *n, float *a, int *lda, int *ipiv,
	 float *work, int *lwork, int *info);
 
/* Subroutine */ int sgetrs_(const char *trans, int *n, int *nrhs, float *a, 
	int *lda, int *ipiv, float *b, int *ldb, int *info);
 
/* Subroutine */ int sggbak_(const char *job, const char *side, int *n, int *ilo, 
	int *ihi, float *lscale, float *rscale, int *m, float *v, 
	int *ldv, int *info);
 
/* Subroutine */ int sggbal_(const char *job, int *n, float *a, int *lda, 
	float *b, int *ldb, int *ilo, int *ihi, float *lscale, float 
	*rscale, float *work, int *info);
 
/* Subroutine */ int sgges_(const char *jobvsl, const char *jobvsr, const char *sort, L_fp 
	selctg, int *n, float *a, int *lda, float *b, int *ldb, 
	int *sdim, float *alphar, float *alphai, float *beta, float *vsl, 
	int *ldvsl, float *vsr, int *ldvsr, float *work, int *lwork,
	 logical *bwork, int *info);
 
/* Subroutine */ int sggesx_(const char *jobvsl, const char *jobvsr, const char *sort, L_fp 
	selctg, const char *sense, int *n, float *a, int *lda, float *b, 
	int *ldb, int *sdim, float *alphar, float *alphai, float *beta, 
	float *vsl, int *ldvsl, float *vsr, int *ldvsr, float *rconde, 
	float *rcondv, float *work, int *lwork, int *iwork, int *
	liwork, logical *bwork, int *info);
 
/* Subroutine */ int sggev_(const char *jobvl, const char *jobvr, int *n, float *a, 
	int *lda, float *b, int *ldb, float *alphar, float *alphai, float 
	*beta, float *vl, int *ldvl, float *vr, int *ldvr, float *work, 
	int *lwork, int *info);
 
/* Subroutine */ int sggevx_(const char *balanc, const char *jobvl, const char *jobvr, const char *
	sense, int *n, float *a, int *lda, float *b, int *ldb, float 
	*alphar, float *alphai, float *beta, float *vl, int *ldvl, float *vr, 
	int *ldvr, int *ilo, int *ihi, float *lscale, float *rscale,
	 float *abnrm, float *bbnrm, float *rconde, float *rcondv, float *work, 
	int *lwork, int *iwork, logical *bwork, int *info);
 
/* Subroutine */ int sggglm_(int *n, int *m, int *p, float *a, 
	int *lda, float *b, int *ldb, float *d__, float *x, float *y, 
	float *work, int *lwork, int *info);
 
/* Subroutine */ int sgghrd_(const char *compq, const char *compz, int *n, int *
	ilo, int *ihi, float *a, int *lda, float *b, int *ldb, float 
	*q, int *ldq, float *z__, int *ldz, int *info);
 
/* Subroutine */ int sgglse_(int *m, int *n, int *p, float *a, 
	int *lda, float *b, int *ldb, float *c__, float *d__, float *x, 
	float *work, int *lwork, int *info);
 
/* Subroutine */ int sggqrf_(int *n, int *m, int *p, float *a, 
	int *lda, float *taua, float *b, int *ldb, float *taub, float *
	work, int *lwork, int *info);
 
/* Subroutine */ int sggrqf_(int *m, int *p, int *n, float *a, 
	int *lda, float *taua, float *b, int *ldb, float *taub, float *
	work, int *lwork, int *info);
 
/* Subroutine */ int sggsvd_(const char *jobu, const char *jobv, const char *jobq, int *m, 
	int *n, int *p, int *k, int *l, float *a, int *lda,
	 float *b, int *ldb, float *alpha, float *beta, float *u, int *
	ldu, float *v, int *ldv, float *q, int *ldq, float *work, 
	int *iwork, int *info);
 
/* Subroutine */ int sggsvp_(const char *jobu, const char *jobv, const char *jobq, int *m, 
	int *p, int *n, float *a, int *lda, float *b, int *ldb, 
	float *tola, float *tolb, int *k, int *l, float *u, int *ldu,
	 float *v, int *ldv, float *q, int *ldq, int *iwork, float *
	tau, float *work, int *info);
 
/* Subroutine */ int sgtcon_(const char *norm, int *n, float *dl, float *d__, 
	float *du, float *du2, int *ipiv, float *anorm, float *rcond, float *
	work, int *iwork, int *info);
 
/* Subroutine */ int sgtrfs_(const char *trans, int *n, int *nrhs, float *dl,
	 float *d__, float *du, float *dlf, float *df, float *duf, float *du2, 
	int *ipiv, float *b, int *ldb, float *x, int *ldx, float *
	ferr, float *berr, float *work, int *iwork, int *info);
 
/* Subroutine */ int sgtsv_(int *n, int *nrhs, float *dl, float *d__, 
	float *du, float *b, int *ldb, int *info);
 
/* Subroutine */ int sgtsvx_(const char *fact, const char *trans, int *n, int *
	nrhs, float *dl, float *d__, float *du, float *dlf, float *df, float *duf, 
	float *du2, int *ipiv, float *b, int *ldb, float *x, int *
	ldx, float *rcond, float *ferr, float *berr, float *work, int *iwork, 
	int *info);
 
/* Subroutine */ int sgttrf_(int *n, float *dl, float *d__, float *du, float *
	du2, int *ipiv, int *info);
 
/* Subroutine */ int sgttrs_(const char *trans, int *n, int *nrhs, float *dl,
	 float *d__, float *du, float *du2, int *ipiv, float *b, int *ldb,
	 int *info);
 
/* Subroutine */ int sgtts2_(int *itrans, int *n, int *nrhs, float 
	*dl, float *d__, float *du, float *du2, int *ipiv, float *b, int *
	ldb);
 
/* Subroutine */ int shgeqz_(const char *job, const char *compq, const char *compz, int *n, 
	int *ilo, int *ihi, float *a, int *lda, float *b, int *
	ldb, float *alphar, float *alphai, float *beta, float *q, int *ldq, 
	float *z__, int *ldz, float *work, int *lwork, int *info);
 
/* Subroutine */ int shsein_(const char *side, const char *eigsrc, const char *initv, logical *
	select, int *n, float *h__, int *ldh, float *wr, float *wi, float 
	*vl, int *ldvl, float *vr, int *ldvr, int *mm, int *m, 
	float *work, int *ifaill, int *ifailr, int *info);
 
/* Subroutine */ int shseqr_(const char *job, const char *compz, int *n, int *ilo,
	 int *ihi, float *h__, int *ldh, float *wr, float *wi, float *z__,
	 int *ldz, float *work, int *lwork, int *info);
 
/* Subroutine */ int slabad_(float *small, float *large);
 
/* Subroutine */ int slabrd_(int *m, int *n, int *nb, float *a, 
	int *lda, float *d__, float *e, float *tauq, float *taup, float *x, 
	int *ldx, float *y, int *ldy);
 
/* Subroutine */ int slacon_(int *n, float *v, float *x, int *isgn, 
	float *est, int *kase);
 
/* Subroutine */ int slacpy_(const char *uplo, int *m, int *n, float *a, 
	int *lda, float *b, int *ldb);
 
/* Subroutine */ int sladiv_(float *a, float *b, float *c__, float *d__, float *p, 
	float *q);
 
/* Subroutine */ int slae2_(float *a, float *b, float *c__, float *rt1, float *rt2);
 
/* Subroutine */ int slaebz_(int *ijob, int *nitmax, int *n, 
	int *mmax, int *minp, int *nbmin, float *abstol, float *
	reltol, float *pivmin, float *d__, float *e, float *e2, int *nval, 
	float *ab, float *c__, int *mout, int *nab, float *work, int 
	*iwork, int *info);
 
/* Subroutine */ int slaed0_(int *icompq, int *qsiz, int *n, float 
	*d__, float *e, float *q, int *ldq, float *qstore, int *ldqs, 
	float *work, int *iwork, int *info);
 
/* Subroutine */ int slaed1_(int *n, float *d__, float *q, int *ldq, 
	int *indxq, float *rho, int *cutpnt, float *work, int *
	iwork, int *info);
 
/* Subroutine */ int slaed2_(int *k, int *n, int *n1, float *d__, 
	float *q, int *ldq, int *indxq, float *rho, float *z__, float *
	dlamda, float *w, float *q2, int *indx, int *indxc, int *
	indxp, int *coltyp, int *info);
 
/* Subroutine */ int slaed3_(int *k, int *n, int *n1, float *d__, 
	float *q, int *ldq, float *rho, float *dlamda, float *q2, int *
	indx, int *ctot, float *w, float *s, int *info);
 
/* Subroutine */ int slaed4_(int *n, int *i__, float *d__, float *z__, 
	float *delta, float *rho, float *dlam, int *info);
 
/* Subroutine */ int slaed5_(int *i__, float *d__, float *z__, float *delta, 
	float *rho, float *dlam);
 
/* Subroutine */ int slaed6_(int *kniter, logical *orgati, float *rho, 
	float *d__, float *z__, float *finit, float *tau, int *info);
 
/* Subroutine */ int slaed7_(int *icompq, int *n, int *qsiz, 
	int *tlvls, int *curlvl, int *curpbm, float *d__, float *q, 
	int *ldq, int *indxq, float *rho, int *cutpnt, float *
	qstore, int *qptr, int *prmptr, int *perm, int *
	givptr, int *givcol, float *givnum, float *work, int *iwork, 
	int *info);
 
/* Subroutine */ int slaed8_(int *icompq, int *k, int *n, int 
	*qsiz, float *d__, float *q, int *ldq, int *indxq, float *rho, 
	int *cutpnt, float *z__, float *dlamda, float *q2, int *ldq2, 
	float *w, int *perm, int *givptr, int *givcol, float *
	givnum, int *indxp, int *indx, int *info);
 
/* Subroutine */ int slaed9_(int *k, int *kstart, int *kstop, 
	int *n, float *d__, float *q, int *ldq, float *rho, float *dlamda,
	 float *w, float *s, int *lds, int *info);
 
/* Subroutine */ int slaeda_(int *n, int *tlvls, int *curlvl, 
	int *curpbm, int *prmptr, int *perm, int *givptr, 
	int *givcol, float *givnum, float *q, int *qptr, float *z__, 
	float *ztemp, int *info);
 
/* Subroutine */ int slaein_(logical *rightv, logical *noinit, int *n, 
	float *h__, int *ldh, float *wr, float *wi, float *vr, float *vi, float 
	*b, int *ldb, float *work, float *eps3, float *smlnum, float *bignum, 
	int *info);
 
/* Subroutine */ int slaev2_(float *a, float *b, float *c__, float *rt1, float *
	rt2, float *cs1, float *sn1);
 
/* Subroutine */ int slaexc_(logical *wantq, int *n, float *t, int *
	ldt, float *q, int *ldq, int *j1, int *n1, int *n2, 
	float *work, int *info);
 
/* Subroutine */ int slag2_(float *a, int *lda, float *b, int *ldb, 
	float *safmin, float *scale1, float *scale2, float *wr1, float *wr2, float *
	wi);
 
/* Subroutine */ int slags2_(logical *upper, float *a1, float *a2, float *a3, 
	float *b1, float *b2, float *b3, float *csu, float *snu, float *csv, float *
	snv, float *csq, float *snq);
 
/* Subroutine */ int slagtf_(int *n, float *a, float *lambda, float *b, float 
	*c__, float *tol, float *d__, int *in, int *info);
 
/* Subroutine */ int slagtm_(const char *trans, int *n, int *nrhs, float *
	alpha, float *dl, float *d__, float *du, float *x, int *ldx, float *
	beta, float *b, int *ldb);
 
/* Subroutine */ int slagts_(int *job, int *n, float *a, float *b, float 
	*c__, float *d__, int *in, float *y, float *tol, int *info);
 
/* Subroutine */ int slagv2_(float *a, int *lda, float *b, int *ldb, 
	float *alphar, float *alphai, float *beta, float *csl, float *snl, float *
	csr, float *snr);
 
/* Subroutine */ int slahqr_(logical *wantt, logical *wantz, int *n, 
	int *ilo, int *ihi, float *h__, int *ldh, float *wr, float *
	wi, int *iloz, int *ihiz, float *z__, int *ldz, int *
	info);
 
/* Subroutine */ int slahrd_(int *n, int *k, int *nb, float *a, 
	int *lda, float *tau, float *t, int *ldt, float *y, int *ldy);
 
/* Subroutine */ int slaic1_(int *job, int *j, float *x, float *sest, 
	float *w, float *gamma, float *sestpr, float *s, float *c__);
 
/* Subroutine */ int slaln2_(logical *ltrans, int *na, int *nw, float *
	smin, float *ca, float *a, int *lda, float *d1, float *d2, float *b, 
	int *ldb, float *wr, float *wi, float *x, int *ldx, float *scale, 
	float *xnorm, int *info);
 
/* Subroutine */ int slals0_(int *icompq, int *nl, int *nr, 
	int *sqre, int *nrhs, float *b, int *ldb, float *bx, 
	int *ldbx, int *perm, int *givptr, int *givcol, 
	int *ldgcol, float *givnum, int *ldgnum, float *poles, float *
	difl, float *difr, float *z__, int *k, float *c__, float *s, float *
	work, int *info);
 
/* Subroutine */ int slalsa_(int *icompq, int *smlsiz, int *n, 
	int *nrhs, float *b, int *ldb, float *bx, int *ldbx, float *
	u, int *ldu, float *vt, int *k, float *difl, float *difr, float *
	z__, float *poles, int *givptr, int *givcol, int *ldgcol, 
	int *perm, float *givnum, float *c__, float *s, float *work, int *
	iwork, int *info);
 
/* Subroutine */ int slalsd_(const char *uplo, int *smlsiz, int *n, int 
	*nrhs, float *d__, float *e, float *b, int *ldb, float *rcond, 
	int *rank, float *work, int *iwork, int *info);
 
/* Subroutine */ int slamc1_(int *beta, int *t, logical *rnd, logical 
	*ieee1);
 
/* Subroutine */ int slamc2_(int *beta, int *t, logical *rnd, float *
	eps, int *emin, float *rmin, int *emax, float *rmax);
 
/* Subroutine */ int slamc4_(int *emin, float *start, int *base);
 
/* Subroutine */ int slamc5_(int *beta, int *p, int *emin, 
	logical *ieee, int *emax, float *rmax);
 
/* Subroutine */ int slamrg_(int *n1, int *n2, float *a, int *
	strd1, int *strd2, int *index);
 
/* Subroutine */ int slanv2_(float *a, float *b, float *c__, float *d__, float *
	rt1r, float *rt1i, float *rt2r, float *rt2i, float *cs, float *sn);
 
/* Subroutine */ int slapll_(int *n, float *x, int *incx, float *y, 
	int *incy, float *ssmin);
 
/* Subroutine */ int slapmt_(logical *forwrd, int *m, int *n, float *x,
	 int *ldx, int *k);
 
/* Subroutine */ int slaqgb_(int *m, int *n, int *kl, int *ku,
	 float *ab, int *ldab, float *r__, float *c__, float *rowcnd, float *
	colcnd, float *amax, const char *equed);
 
/* Subroutine */ int slaqge_(int *m, int *n, float *a, int *lda, 
	float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, const char *
	equed);
 
/* Subroutine */ int slaqp2_(int *m, int *n, int *offset, float *a,
	 int *lda, int *jpvt, float *tau, float *vn1, float *vn2, float *
	work);
 
/* Subroutine */ int slaqps_(int *m, int *n, int *offset, int 
	*nb, int *kb, float *a, int *lda, int *jpvt, float *tau, 
	float *vn1, float *vn2, float *auxv, float *f, int *ldf);
 
/* Subroutine */ int slaqsb_(const char *uplo, int *n, int *kd, float *ab, 
	int *ldab, float *s, float *scond, float *amax, const char *equed);
 
/* Subroutine */ int slaqsp_(const char *uplo, int *n, float *ap, float *s, float *
	scond, float *amax, const char *equed);
 
/* Subroutine */ int slaqsy_(const char *uplo, int *n, float *a, int *lda, 
	float *s, float *scond, float *amax, const char *equed);
 
/* Subroutine */ int slaqtr_(logical *ltran, logical *lfloat, int *n, float 
	*t, int *ldt, float *b, float *w, float *scale, float *x, float *work, 
	int *info);
 
/* Subroutine */ int slar1v_(int *n, int *b1, int *bn, float *
	sigma, float *d__, float *l, float *ld, float *lld, float *gersch, float *
	z__, float *ztz, float *mingma, int *r__, int *isuppz, float *
	work);
 
/* Subroutine */ int slar2v_(int *n, float *x, float *y, float *z__, int 
	*incx, float *c__, float *s, int *incc);
 
/* Subroutine */ int slarf_(const char *side, int *m, int *n, float *v, 
	int *incv, float *tau, float *c__, int *ldc, float *work);
 
/* Subroutine */ int slarfb_(const char *side, const char *trans, const char *direct, const char *
	storev, int *m, int *n, int *k, float *v, int *ldv, 
	float *t, int *ldt, float *c__, int *ldc, float *work, int *
	ldwork);
 
/* Subroutine */ int slarfg_(int *n, float *alpha, float *x, int *incx, 
	float *tau);
 
/* Subroutine */ int slarft_(const char *direct, const char *storev, int *n, int *
	k, float *v, int *ldv, float *tau, float *t, int *ldt);
 
/* Subroutine */ int slarfx_(const char *side, int *m, int *n, float *v, 
	float *tau, float *c__, int *ldc, float *work);
 
/* Subroutine */ int slargv_(int *n, float *x, int *incx, float *y, 
	int *incy, float *c__, int *incc);
 
/* Subroutine */ int slarnv_(int *idist, int *iseed, int *n, float 
	*x);
 
/* Subroutine */ int slarrb_(int *n, float *d__, float *l, float *ld, float *
	lld, int *ifirst, int *ilast, float *sigma, float *reltol, float 
	*w, float *wgap, float *werr, float *work, int *iwork, int *info);
 
/* Subroutine */ int slarre_(int *n, float *d__, float *e, float *tol, 
	int *nsplit, int *isplit, int *m, float *w, float *woff, 
	float *gersch, float *work, int *info);
 
/* Subroutine */ int slarrf_(int *n, float *d__, float *l, float *ld, float *
	lld, int *ifirst, int *ilast, float *w, float *dplus, float *
	lplus, float *work, int *iwork, int *info);
 
/* Subroutine */ int slarrv_(int *n, float *d__, float *l, int *isplit, 
	int *m, float *w, int *iblock, float *gersch, float *tol, float *
	z__, int *ldz, int *isuppz, float *work, int *iwork, 
	int *info);
 
/* Subroutine */ int slartg_(float *f, float *g, float *cs, float *sn, float *r__);
 
/* Subroutine */ int slartv_(int *n, float *x, int *incx, float *y, 
	int *incy, float *c__, float *s, int *incc);
 
/* Subroutine */ int slaruv_(int *iseed, int *n, float *x);
 
/* Subroutine */ int slarz_(const char *side, int *m, int *n, int *l, 
	float *v, int *incv, float *tau, float *c__, int *ldc, float *
	work);
 
/* Subroutine */ int slarzb_(const char *side, const char *trans, const char *direct, const char *
	storev, int *m, int *n, int *k, int *l, float *v, 
	int *ldv, float *t, int *ldt, float *c__, int *ldc, float *
	work, int *ldwork);
 
/* Subroutine */ int slarzt_(const char *direct, const char *storev, int *n, int *
	k, float *v, int *ldv, float *tau, float *t, int *ldt);
 
/* Subroutine */ int slas2_(float *f, float *g, float *h__, float *ssmin, float *
	ssmax);
 
/* Subroutine */ int slascl_(const char *type__, int *kl, int *ku, float *
	cfrom, float *cto, int *m, int *n, float *a, int *lda, 
	int *info);
 
/* Subroutine */ int slasd0_(int *n, int *sqre, float *d__, float *e, 
	float *u, int *ldu, float *vt, int *ldvt, int *smlsiz, 
	int *iwork, float *work, int *info);
 
/* Subroutine */ int slasd1_(int *nl, int *nr, int *sqre, float *
	d__, float *alpha, float *beta, float *u, int *ldu, float *vt, 
	int *ldvt, int *idxq, int *iwork, float *work, int *
	info);
 
/* Subroutine */ int slasd2_(int *nl, int *nr, int *sqre, int 
	*k, float *d__, float *z__, float *alpha, float *beta, float *u, int *
	ldu, float *vt, int *ldvt, float *dsigma, float *u2, int *ldu2, 
	float *vt2, int *ldvt2, int *idxp, int *idx, int *idxc,
	 int *idxq, int *coltyp, int *info);
 
/* Subroutine */ int slasd3_(int *nl, int *nr, int *sqre, int 
	*k, float *d__, float *q, int *ldq, float *dsigma, float *u, int *
	ldu, float *u2, int *ldu2, float *vt, int *ldvt, float *vt2, 
	int *ldvt2, int *idxc, int *ctot, float *z__, int *
	info);
 
/* Subroutine */ int slasd4_(int *n, int *i__, float *d__, float *z__, 
	float *delta, float *rho, float *sigma, float *work, int *info);
 
/* Subroutine */ int slasd5_(int *i__, float *d__, float *z__, float *delta, 
	float *rho, float *dsigma, float *work);
 
/* Subroutine */ int slasd6_(int *icompq, int *nl, int *nr, 
	int *sqre, float *d__, float *vf, float *vl, float *alpha, float *beta,
	 int *idxq, int *perm, int *givptr, int *givcol, 
	int *ldgcol, float *givnum, int *ldgnum, float *poles, float *
	difl, float *difr, float *z__, int *k, float *c__, float *s, float *
	work, int *iwork, int *info);
 
/* Subroutine */ int slasd7_(int *icompq, int *nl, int *nr, 
	int *sqre, int *k, float *d__, float *z__, float *zw, float *vf, 
	float *vfw, float *vl, float *vlw, float *alpha, float *beta, float *dsigma,
	 int *idx, int *idxp, int *idxq, int *perm, int *
	givptr, int *givcol, int *ldgcol, float *givnum, int *
	ldgnum, float *c__, float *s, int *info);
 
/* Subroutine */ int slasd8_(int *icompq, int *k, float *d__, float *
	z__, float *vf, float *vl, float *difl, float *difr, int *lddifr, 
	float *dsigma, float *work, int *info);
 
/* Subroutine */ int slasd9_(int *icompq, int *ldu, int *k, float *
	d__, float *z__, float *vf, float *vl, float *difl, float *difr, float *
	dsigma, float *work, int *info);
 
/* Subroutine */ int slasda_(int *icompq, int *smlsiz, int *n, 
	int *sqre, float *d__, float *e, float *u, int *ldu, float *vt, 
	int *k, float *difl, float *difr, float *z__, float *poles, int *
	givptr, int *givcol, int *ldgcol, int *perm, float *givnum,
	 float *c__, float *s, float *work, int *iwork, int *info);
 
/* Subroutine */ int slasdq_(const char *uplo, int *sqre, int *n, int *
	ncvt, int *nru, int *ncc, float *d__, float *e, float *vt, 
	int *ldvt, float *u, int *ldu, float *c__, int *ldc, float *
	work, int *info);
 
/* Subroutine */ int slasdt_(int *n, int *lvl, int *nd, int *
	inode, int *ndiml, int *ndimr, int *msub);
 
/* Subroutine */ int slaset_(const char *uplo, int *m, int *n, float *alpha, 
	float *beta, float *a, int *lda);
 
/* Subroutine */ int slasq1_(int *n, float *d__, float *e, float *work, 
	int *info);
 
/* Subroutine */ int slasq2_(int *n, float *z__, int *info);
 
/* Subroutine */ int slasq3_(int *i0, int *n0, float *z__, int *pp,
	 float *dmin__, float *sigma, float *desig, float *qmax, int *nfail, 
	int *iter, int *ndiv, logical *ieee);
 
/* Subroutine */ int slasq4_(int *i0, int *n0, float *z__, int *pp,
	 int *n0in, float *dmin__, float *dmin1, float *dmin2, float *dn, 
	float *dn1, float *dn2, float *tau, int *ttype);
 
/* Subroutine */ int slasq5_(int *i0, int *n0, float *z__, int *pp,
	 float *tau, float *dmin__, float *dmin1, float *dmin2, float *dn, float *
	dnm1, float *dnm2, logical *ieee);
 
/* Subroutine */ int slasq6_(int *i0, int *n0, float *z__, int *pp,
	 float *dmin__, float *dmin1, float *dmin2, float *dn, float *dnm1, float *
	dnm2);
 
/* Subroutine */ int slasr_(const char *side, const char *pivot, const char *direct, int *m,
	 int *n, float *c__, float *s, float *a, int *lda);
 
/* Subroutine */ int slasrt_(const char *id, int *n, float *d__, int *info);
 
/* Subroutine */ int slassq_(int *n, float *x, int *incx, float *scale, 
	float *sumsq);
 
/* Subroutine */ int slasv2_(float *f, float *g, float *h__, float *ssmin, float *
	ssmax, float *snr, float *csr, float *snl, float *csl);
 
/* Subroutine */ int slaswp_(int *n, float *a, int *lda, int *k1, 
	int *k2, int *ipiv, int *incx);
 
/* Subroutine */ int slasy2_(logical *ltranl, logical *ltranr, int *isgn, 
	int *n1, int *n2, float *tl, int *ldtl, float *tr, int *
	ldtr, float *b, int *ldb, float *scale, float *x, int *ldx, float 
	*xnorm, int *info);
 
/* Subroutine */ int slasyf_(const char *uplo, int *n, int *nb, int *kb,
	 float *a, int *lda, int *ipiv, float *w, int *ldw, int 
	*info);
 
/* Subroutine */ int slatbs_(const char *uplo, const char *trans, const char *diag, const char *
	normin, int *n, int *kd, float *ab, int *ldab, float *x, 
	float *scale, float *cnorm, int *info);
 
/* Subroutine */ int slatdf_(int *ijob, int *n, float *z__, int *
	ldz, float *rhs, float *rdsum, float *rdscal, int *ipiv, int *
	jpiv);
 
/* Subroutine */ int slatps_(const char *uplo, const char *trans, const char *diag, const char *
	normin, int *n, float *ap, float *x, float *scale, float *cnorm, 
	int *info);
 
/* Subroutine */ int slatrd_(const char *uplo, int *n, int *nb, float *a, 
	int *lda, float *e, float *tau, float *w, int *ldw);
 
/* Subroutine */ int slatrs_(const char *uplo, const char *trans, const char *diag, const char *
	normin, int *n, float *a, int *lda, float *x, float *scale, float 
	*cnorm, int *info);
 
/* Subroutine */ int slatrz_(int *m, int *n, int *l, float *a, 
	int *lda, float *tau, float *work);
 
/* Subroutine */ int slatzm_(const char *side, int *m, int *n, float *v, 
	int *incv, float *tau, float *c1, float *c2, int *ldc, float *
	work);
 
/* Subroutine */ int slauu2_(const char *uplo, int *n, float *a, int *lda, 
	int *info);
 
/* Subroutine */ int slauum_(const char *uplo, int *n, float *a, int *lda, 
	int *info);
 
/* Subroutine */ int sopgtr_(const char *uplo, int *n, float *ap, float *tau, 
	float *q, int *ldq, float *work, int *info);
 
/* Subroutine */ int sopmtr_(const char *side, const char *uplo, const char *trans, int *m, 
	int *n, float *ap, float *tau, float *c__, int *ldc, float *work, 
	int *info);
 
/* Subroutine */ int sorg2l_(int *m, int *n, int *k, float *a, 
	int *lda, float *tau, float *work, int *info);
 
/* Subroutine */ int sorg2r_(int *m, int *n, int *k, float *a, 
	int *lda, float *tau, float *work, int *info);
 
/* Subroutine */ int sorgbr_(const char *vect, int *m, int *n, int *k, 
	float *a, int *lda, float *tau, float *work, int *lwork, int 
	*info);
 
/* Subroutine */ int sorghr_(int *n, int *ilo, int *ihi, float *a, 
	int *lda, float *tau, float *work, int *lwork, int *info);
 
/* Subroutine */ int sorgl2_(int *m, int *n, int *k, float *a, 
	int *lda, float *tau, float *work, int *info);
 
/* Subroutine */ int sorglq_(int *m, int *n, int *k, float *a, 
	int *lda, float *tau, float *work, int *lwork, int *info);
 
/* Subroutine */ int sorgql_(int *m, int *n, int *k, float *a, 
	int *lda, float *tau, float *work, int *lwork, int *info);
 
/* Subroutine */ int sorgqr_(int *m, int *n, int *k, float *a, 
	int *lda, float *tau, float *work, int *lwork, int *info);
 
/* Subroutine */ int sorgr2_(int *m, int *n, int *k, float *a, 
	int *lda, float *tau, float *work, int *info);
 
/* Subroutine */ int sorgrq_(int *m, int *n, int *k, float *a, 
	int *lda, float *tau, float *work, int *lwork, int *info);
 
/* Subroutine */ int sorgtr_(const char *uplo, int *n, float *a, int *lda, 
	float *tau, float *work, int *lwork, int *info);
 
/* Subroutine */ int sorm2l_(const char *side, const char *trans, int *m, int *n, 
	int *k, float *a, int *lda, float *tau, float *c__, int *ldc,
	 float *work, int *info);
 
/* Subroutine */ int sorm2r_(const char *side, const char *trans, int *m, int *n, 
	int *k, float *a, int *lda, float *tau, float *c__, int *ldc,
	 float *work, int *info);
 
/* Subroutine */ int sormbr_(const char *vect, const char *side, const char *trans, int *m, 
	int *n, int *k, float *a, int *lda, float *tau, float *c__, 
	int *ldc, float *work, int *lwork, int *info);
 
/* Subroutine */ int sormhr_(const char *side, const char *trans, int *m, int *n, 
	int *ilo, int *ihi, float *a, int *lda, float *tau, float *
	c__, int *ldc, float *work, int *lwork, int *info);
 
/* Subroutine */ int sorml2_(const char *side, const char *trans, int *m, int *n, 
	int *k, float *a, int *lda, float *tau, float *c__, int *ldc,
	 float *work, int *info);
 
/* Subroutine */ int sormlq_(const char *side, const char *trans, int *m, int *n, 
	int *k, float *a, int *lda, float *tau, float *c__, int *ldc,
	 float *work, int *lwork, int *info);
 
/* Subroutine */ int sormql_(const char *side, const char *trans, int *m, int *n, 
	int *k, float *a, int *lda, float *tau, float *c__, int *ldc,
	 float *work, int *lwork, int *info);
 
/* Subroutine */ int sormqr_(const char *side, const char *trans, int *m, int *n, 
	int *k, float *a, int *lda, float *tau, float *c__, int *ldc,
	 float *work, int *lwork, int *info);
 
/* Subroutine */ int sormr2_(const char *side, const char *trans, int *m, int *n, 
	int *k, float *a, int *lda, float *tau, float *c__, int *ldc,
	 float *work, int *info);
 
/* Subroutine */ int sormr3_(const char *side, const char *trans, int *m, int *n, 
	int *k, int *l, float *a, int *lda, float *tau, float *c__, 
	int *ldc, float *work, int *info);
 
/* Subroutine */ int sormrq_(const char *side, const char *trans, int *m, int *n, 
	int *k, float *a, int *lda, float *tau, float *c__, int *ldc,
	 float *work, int *lwork, int *info);
 
/* Subroutine */ int sormrz_(const char *side, const char *trans, int *m, int *n, 
	int *k, int *l, float *a, int *lda, float *tau, float *c__, 
	int *ldc, float *work, int *lwork, int *info);
 
/* Subroutine */ int sormtr_(const char *side, const char *uplo, const char *trans, int *m, 
	int *n, float *a, int *lda, float *tau, float *c__, int *ldc,
	 float *work, int *lwork, int *info);
 
/* Subroutine */ int spbcon_(const char *uplo, int *n, int *kd, float *ab, 
	int *ldab, float *anorm, float *rcond, float *work, int *iwork, 
	int *info);
 
/* Subroutine */ int spbequ_(const char *uplo, int *n, int *kd, float *ab, 
	int *ldab, float *s, float *scond, float *amax, int *info);
 
/* Subroutine */ int spbrfs_(const char *uplo, int *n, int *kd, int *
	nrhs, float *ab, int *ldab, float *afb, int *ldafb, float *b, 
	int *ldb, float *x, int *ldx, float *ferr, float *berr, float *
	work, int *iwork, int *info);
 
/* Subroutine */ int spbstf_(const char *uplo, int *n, int *kd, float *ab, 
	int *ldab, int *info);
 
/* Subroutine */ int spbsv_(const char *uplo, int *n, int *kd, int *
	nrhs, float *ab, int *ldab, float *b, int *ldb, int *info);
 
/* Subroutine */ int spbsvx_(const char *fact, const char *uplo, int *n, int *kd, 
	int *nrhs, float *ab, int *ldab, float *afb, int *ldafb, 
	const char *equed, float *s, float *b, int *ldb, float *x, int *ldx, 
	float *rcond, float *ferr, float *berr, float *work, int *iwork, 
	int *info);
 
/* Subroutine */ int spbtf2_(const char *uplo, int *n, int *kd, float *ab, 
	int *ldab, int *info);
 
/* Subroutine */ int spbtrf_(const char *uplo, int *n, int *kd, float *ab, 
	int *ldab, int *info);
 
/* Subroutine */ int spbtrs_(const char *uplo, int *n, int *kd, int *
	nrhs, float *ab, int *ldab, float *b, int *ldb, int *info);
 
/* Subroutine */ int spocon_(const char *uplo, int *n, float *a, int *lda, 
	float *anorm, float *rcond, float *work, int *iwork, int *info);
 
/* Subroutine */ int spoequ_(int *n, float *a, int *lda, float *s, float 
	*scond, float *amax, int *info);
 
/* Subroutine */ int sporfs_(const char *uplo, int *n, int *nrhs, float *a, 
	int *lda, float *af, int *ldaf, float *b, int *ldb, float *x,
	 int *ldx, float *ferr, float *berr, float *work, int *iwork, 
	int *info);
 
/* Subroutine */ int sposv_(const char *uplo, int *n, int *nrhs, float *a, 
	int *lda, float *b, int *ldb, int *info);
 
/* Subroutine */ int sposvx_(const char *fact, const char *uplo, int *n, int *
	nrhs, float *a, int *lda, float *af, int *ldaf, const char *equed, 
	float *s, float *b, int *ldb, float *x, int *ldx, float *rcond, 
	float *ferr, float *berr, float *work, int *iwork, int *info);
 
/* Subroutine */ int spotf2_(const char *uplo, int *n, float *a, int *lda, 
	int *info);
 
/* Subroutine */ int spotrf_(const char *uplo, int *n, float *a, int *lda, 
	int *info);
 
/* Subroutine */ int spotri_(const char *uplo, int *n, float *a, int *lda, 
	int *info);
 
/* Subroutine */ int spotrs_(const char *uplo, int *n, int *nrhs, float *a, 
	int *lda, float *b, int *ldb, int *info);
 
/* Subroutine */ int sppcon_(const char *uplo, int *n, float *ap, float *anorm, 
	float *rcond, float *work, int *iwork, int *info);
 
/* Subroutine */ int sppequ_(const char *uplo, int *n, float *ap, float *s, float *
	scond, float *amax, int *info);
 
/* Subroutine */ int spprfs_(const char *uplo, int *n, int *nrhs, float *ap, 
	float *afp, float *b, int *ldb, float *x, int *ldx, float *ferr, 
	float *berr, float *work, int *iwork, int *info);
 
/* Subroutine */ int sppsv_(const char *uplo, int *n, int *nrhs, float *ap, 
	float *b, int *ldb, int *info);
 
/* Subroutine */ int sppsvx_(const char *fact, const char *uplo, int *n, int *
	nrhs, float *ap, float *afp, const char *equed, float *s, float *b, int *
	ldb, float *x, int *ldx, float *rcond, float *ferr, float *berr, float 
	*work, int *iwork, int *info);
 
/* Subroutine */ int spptrf_(const char *uplo, int *n, float *ap, int *info);
 
/* Subroutine */ int spptri_(const char *uplo, int *n, float *ap, int *info);
 
/* Subroutine */ int spptrs_(const char *uplo, int *n, int *nrhs, float *ap, 
	float *b, int *ldb, int *info);
 
/* Subroutine */ int sptcon_(int *n, float *d__, float *e, float *anorm, 
	float *rcond, float *work, int *info);
 
/* Subroutine */ int spteqr_(const char *compz, int *n, float *d__, float *e, 
	float *z__, int *ldz, float *work, int *info);
 
/* Subroutine */ int sptrfs_(int *n, int *nrhs, float *d__, float *e, 
	float *df, float *ef, float *b, int *ldb, float *x, int *ldx, 
	float *ferr, float *berr, float *work, int *info);
 
/* Subroutine */ int sptsv_(int *n, int *nrhs, float *d__, float *e, 
	float *b, int *ldb, int *info);
 
/* Subroutine */ int sptsvx_(const char *fact, int *n, int *nrhs, float *d__,
	 float *e, float *df, float *ef, float *b, int *ldb, float *x, int 
	*ldx, float *rcond, float *ferr, float *berr, float *work, int *info);
 
/* Subroutine */ int spttrf_(int *n, float *d__, float *e, int *info);
 
/* Subroutine */ int spttrs_(int *n, int *nrhs, float *d__, float *e, 
	float *b, int *ldb, int *info);
 
/* Subroutine */ int sptts2_(int *n, int *nrhs, float *d__, float *e, 
	float *b, int *ldb);
 
/* Subroutine */ int srscl_(int *n, float *sa, float *sx, int *incx);
 
/* Subroutine */ int ssbev_(const char *jobz, const char *uplo, int *n, int *kd, 
	float *ab, int *ldab, float *w, float *z__, int *ldz, float *work,
	 int *info);
 
/* Subroutine */ int ssbevd_(const char *jobz, const char *uplo, int *n, int *kd, 
	float *ab, int *ldab, float *w, float *z__, int *ldz, float *work,
	 int *lwork, int *iwork, int *liwork, int *info);
 
/* Subroutine */ int ssbevx_(const char *jobz, const char *range, const char *uplo, int *n, 
	int *kd, float *ab, int *ldab, float *q, int *ldq, float *vl,
	 float *vu, int *il, int *iu, float *abstol, int *m, float *
	w, float *z__, int *ldz, float *work, int *iwork, int *
	ifail, int *info);
 
/* Subroutine */ int ssbgst_(const char *vect, const char *uplo, int *n, int *ka, 
	int *kb, float *ab, int *ldab, float *bb, int *ldbb, float *
	x, int *ldx, float *work, int *info);
 
/* Subroutine */ int ssbgv_(const char *jobz, const char *uplo, int *n, int *ka, 
	int *kb, float *ab, int *ldab, float *bb, int *ldbb, float *
	w, float *z__, int *ldz, float *work, int *info);
 
/* Subroutine */ int ssbgvd_(const char *jobz, const char *uplo, int *n, int *ka, 
	int *kb, float *ab, int *ldab, float *bb, int *ldbb, float *
	w, float *z__, int *ldz, float *work, int *lwork, int *
	iwork, int *liwork, int *info);
 
/* Subroutine */ int ssbgvx_(const char *jobz, const char *range, const char *uplo, int *n, 
	int *ka, int *kb, float *ab, int *ldab, float *bb, int *
	ldbb, float *q, int *ldq, float *vl, float *vu, int *il, int 
	*iu, float *abstol, int *m, float *w, float *z__, int *ldz, float 
	*work, int *iwork, int *ifail, int *info);
 
/* Subroutine */ int ssbtrd_(const char *vect, const char *uplo, int *n, int *kd, 
	float *ab, int *ldab, float *d__, float *e, float *q, int *ldq, 
	float *work, int *info);
 
/* Subroutine */ int sspcon_(const char *uplo, int *n, float *ap, int *ipiv, 
	float *anorm, float *rcond, float *work, int *iwork, int *info);
 
/* Subroutine */ int sspev_(const char *jobz, const char *uplo, int *n, float *ap, 
	float *w, float *z__, int *ldz, float *work, int *info);
 
/* Subroutine */ int sspevd_(const char *jobz, const char *uplo, int *n, float *ap, 
	float *w, float *z__, int *ldz, float *work, int *lwork, int 
	*iwork, int *liwork, int *info);
 
/* Subroutine */ int sspevx_(const char *jobz, const char *range, const char *uplo, int *n, 
	float *ap, float *vl, float *vu, int *il, int *iu, float *abstol, 
	int *m, float *w, float *z__, int *ldz, float *work, int *
	iwork, int *ifail, int *info);
 
/* Subroutine */ int sspgst_(int *itype, const char *uplo, int *n, float *ap,
	 float *bp, int *info);
 
/* Subroutine */ int sspgv_(int *itype, const char *jobz, const char *uplo, int *
	n, float *ap, float *bp, float *w, float *z__, int *ldz, float *work, 
	int *info);
 
/* Subroutine */ int sspgvd_(int *itype, const char *jobz, const char *uplo, int *
	n, float *ap, float *bp, float *w, float *z__, int *ldz, float *work, 
	int *lwork, int *iwork, int *liwork, int *info);
 
/* Subroutine */ int sspgvx_(int *itype, const char *jobz, const char *range, const char *
	uplo, int *n, float *ap, float *bp, float *vl, float *vu, int *il,
	 int *iu, float *abstol, int *m, float *w, float *z__, int *
	ldz, float *work, int *iwork, int *ifail, int *info);
 
/* Subroutine */ int ssprfs_(const char *uplo, int *n, int *nrhs, float *ap, 
	float *afp, int *ipiv, float *b, int *ldb, float *x, int *
	ldx, float *ferr, float *berr, float *work, int *iwork, int *
	info);
 
/* Subroutine */ int sspsv_(const char *uplo, int *n, int *nrhs, float *ap, 
	int *ipiv, float *b, int *ldb, int *info);
 
/* Subroutine */ int sspsvx_(const char *fact, const char *uplo, int *n, int *
	nrhs, float *ap, float *afp, int *ipiv, float *b, int *ldb, float 
	*x, int *ldx, float *rcond, float *ferr, float *berr, float *work, 
	int *iwork, int *info);
 
/* Subroutine */ int ssptrd_(const char *uplo, int *n, float *ap, float *d__, 
	float *e, float *tau, int *info);
 
/* Subroutine */ int ssptrf_(const char *uplo, int *n, float *ap, int *ipiv, 
	int *info);
 
/* Subroutine */ int ssptri_(const char *uplo, int *n, float *ap, int *ipiv, 
	float *work, int *info);
 
/* Subroutine */ int ssptrs_(const char *uplo, int *n, int *nrhs, float *ap, 
	int *ipiv, float *b, int *ldb, int *info);
 
/* Subroutine */ int sstebz_(const char *range, const char *order, int *n, float *vl, 
	float *vu, int *il, int *iu, float *abstol, float *d__, float *e, 
	int *m, int *nsplit, float *w, int *iblock, int *
	isplit, float *work, int *iwork, int *info);
 
/* Subroutine */ int sstedc_(const char *compz, int *n, float *d__, float *e, 
	float *z__, int *ldz, float *work, int *lwork, int *iwork, 
	int *liwork, int *info);
 
/* Subroutine */ int sstegr_(const char *jobz, const char *range, int *n, float *d__, 
	float *e, float *vl, float *vu, int *il, int *iu, float *abstol, 
	int *m, float *w, float *z__, int *ldz, int *isuppz, float *
	work, int *lwork, int *iwork, int *liwork, int *info);
 
/* Subroutine */ int sstein_(int *n, float *d__, float *e, int *m, float 
	*w, int *iblock, int *isplit, float *z__, int *ldz, float *
	work, int *iwork, int *ifail, int *info);
 
/* Subroutine */ int ssteqr_(const char *compz, int *n, float *d__, float *e, 
	float *z__, int *ldz, float *work, int *info);
 
/* Subroutine */ int ssterf_(int *n, float *d__, float *e, int *info);
 
/* Subroutine */ int sstev_(const char *jobz, int *n, float *d__, float *e, float *
	z__, int *ldz, float *work, int *info);
 
/* Subroutine */ int sstevd_(const char *jobz, int *n, float *d__, float *e, float 
	*z__, int *ldz, float *work, int *lwork, int *iwork, 
	int *liwork, int *info);
 
/* Subroutine */ int sstevr_(const char *jobz, const char *range, int *n, float *d__, 
	float *e, float *vl, float *vu, int *il, int *iu, float *abstol, 
	int *m, float *w, float *z__, int *ldz, int *isuppz, float *
	work, int *lwork, int *iwork, int *liwork, int *info);
 
/* Subroutine */ int sstevx_(const char *jobz, const char *range, int *n, float *d__, 
	float *e, float *vl, float *vu, int *il, int *iu, float *abstol, 
	int *m, float *w, float *z__, int *ldz, float *work, int *
	iwork, int *ifail, int *info);
 
/* Subroutine */ int ssycon_(const char *uplo, int *n, float *a, int *lda, 
	int *ipiv, float *anorm, float *rcond, float *work, int *iwork, 
	int *info);
 
/* Subroutine */ int ssyev_(const char *jobz, const char *uplo, int *n, float *a, 
	int *lda, float *w, float *work, int *lwork, int *info);
 
/* Subroutine */ int ssyevd_(const char *jobz, const char *uplo, int *n, float *a, 
	int *lda, float *w, float *work, int *lwork, int *iwork, 
	int *liwork, int *info);
 
/* Subroutine */ int ssyevr_(const char *jobz, const char *range, const char *uplo, int *n, 
	float *a, int *lda, float *vl, float *vu, int *il, int *iu, 
	float *abstol, int *m, float *w, float *z__, int *ldz, int *
	isuppz, float *work, int *lwork, int *iwork, int *liwork, 
	int *info);
 
/* Subroutine */ int ssyevx_(const char *jobz, const char *range, const char *uplo, int *n, 
	float *a, int *lda, float *vl, float *vu, int *il, int *iu, 
	float *abstol, int *m, float *w, float *z__, int *ldz, float *
	work, int *lwork, int *iwork, int *ifail, int *info);
 
/* Subroutine */ int ssygs2_(int *itype, const char *uplo, int *n, float *a, 
	int *lda, float *b, int *ldb, int *info);
 
/* Subroutine */ int ssygst_(int *itype, const char *uplo, int *n, float *a, 
	int *lda, float *b, int *ldb, int *info);
 
/* Subroutine */ int ssygv_(int *itype, const char *jobz, const char *uplo, int *
	n, float *a, int *lda, float *b, int *ldb, float *w, float *work, 
	int *lwork, int *info);
 
/* Subroutine */ int ssygvd_(int *itype, const char *jobz, const char *uplo, int *
	n, float *a, int *lda, float *b, int *ldb, float *w, float *work, 
	int *lwork, int *iwork, int *liwork, int *info);
 
/* Subroutine */ int ssygvx_(int *itype, const char *jobz, const char *range, const char *
	uplo, int *n, float *a, int *lda, float *b, int *ldb, float *
	vl, float *vu, int *il, int *iu, float *abstol, int *m, 
	float *w, float *z__, int *ldz, float *work, int *lwork, int 
	*iwork, int *ifail, int *info);
 
/* Subroutine */ int ssyrfs_(const char *uplo, int *n, int *nrhs, float *a, 
	int *lda, float *af, int *ldaf, int *ipiv, float *b, 
	int *ldb, float *x, int *ldx, float *ferr, float *berr, float *
	work, int *iwork, int *info);
 
/* Subroutine */ int ssysv_(const char *uplo, int *n, int *nrhs, float *a, 
	int *lda, int *ipiv, float *b, int *ldb, float *work, 
	int *lwork, int *info);
 
/* Subroutine */ int ssysvx_(const char *fact, const char *uplo, int *n, int *
	nrhs, float *a, int *lda, float *af, int *ldaf, int *ipiv, 
	float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr,
	 float *berr, float *work, int *lwork, int *iwork, int *
	info);
 
/* Subroutine */ int ssytd2_(const char *uplo, int *n, float *a, int *lda, 
	float *d__, float *e, float *tau, int *info);
 
/* Subroutine */ int ssytf2_(const char *uplo, int *n, float *a, int *lda, 
	int *ipiv, int *info);
 
/* Subroutine */ int ssytrd_(const char *uplo, int *n, float *a, int *lda, 
	float *d__, float *e, float *tau, float *work, int *lwork, int *
	info);
 
/* Subroutine */ int ssytrf_(const char *uplo, int *n, float *a, int *lda, 
	int *ipiv, float *work, int *lwork, int *info);
 
/* Subroutine */ int ssytri_(const char *uplo, int *n, float *a, int *lda, 
	int *ipiv, float *work, int *info);
 
/* Subroutine */ int ssytrs_(const char *uplo, int *n, int *nrhs, float *a, 
	int *lda, int *ipiv, float *b, int *ldb, int *info);
 
/* Subroutine */ int stbcon_(const char *norm, const char *uplo, const char *diag, int *n, 
	int *kd, float *ab, int *ldab, float *rcond, float *work, 
	int *iwork, int *info);
 
/* Subroutine */ int stbrfs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *kd, int *nrhs, float *ab, int *ldab, float *b, int 
	*ldb, float *x, int *ldx, float *ferr, float *berr, float *work, 
	int *iwork, int *info);
 
/* Subroutine */ int stbtrs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *kd, int *nrhs, float *ab, int *ldab, float *b, int 
	*ldb, int *info);
 
/* Subroutine */ int stgevc_(const char *side, const char *howmny, logical *select, 
	int *n, float *a, int *lda, float *b, int *ldb, float *vl, 
	int *ldvl, float *vr, int *ldvr, int *mm, int *m, float 
	*work, int *info);
 
/* Subroutine */ int stgex2_(logical *wantq, logical *wantz, int *n, float 
	*a, int *lda, float *b, int *ldb, float *q, int *ldq, float *
	z__, int *ldz, int *j1, int *n1, int *n2, float *work, 
	int *lwork, int *info);
 
/* Subroutine */ int stgexc_(logical *wantq, logical *wantz, int *n, float 
	*a, int *lda, float *b, int *ldb, float *q, int *ldq, float *
	z__, int *ldz, int *ifst, int *ilst, float *work, int *
	lwork, int *info);
 
/* Subroutine */ int stgsen_(int *ijob, logical *wantq, logical *wantz, 
	logical *select, int *n, float *a, int *lda, float *b, int *
	ldb, float *alphar, float *alphai, float *beta, float *q, int *ldq, 
	float *z__, int *ldz, int *m, float *pl, float *pr, float *dif, 
	float *work, int *lwork, int *iwork, int *liwork, int *
	info);
 
/* Subroutine */ int stgsja_(const char *jobu, const char *jobv, const char *jobq, int *m, 
	int *p, int *n, int *k, int *l, float *a, int *lda,
	 float *b, int *ldb, float *tola, float *tolb, float *alpha, float *
	beta, float *u, int *ldu, float *v, int *ldv, float *q, int *
	ldq, float *work, int *ncycle, int *info);
 
/* Subroutine */ int stgsna_(const char *job, const char *howmny, logical *select, 
	int *n, float *a, int *lda, float *b, int *ldb, float *vl, 
	int *ldvl, float *vr, int *ldvr, float *s, float *dif, int *
	mm, int *m, float *work, int *lwork, int *iwork, int *
	info);
 
/* Subroutine */ int stgsy2_(const char *trans, int *ijob, int *m, int *
	n, float *a, int *lda, float *b, int *ldb, float *c__, int *
	ldc, float *d__, int *ldd, float *e, int *lde, float *f, int 
	*ldf, float *scale, float *rdsum, float *rdscal, int *iwork, int 
	*pq, int *info);
 
/* Subroutine */ int stgsyl_(const char *trans, int *ijob, int *m, int *
	n, float *a, int *lda, float *b, int *ldb, float *c__, int *
	ldc, float *d__, int *ldd, float *e, int *lde, float *f, int 
	*ldf, float *scale, float *dif, float *work, int *lwork, int *
	iwork, int *info);
 
/* Subroutine */ int stpcon_(const char *norm, const char *uplo, const char *diag, int *n, 
	float *ap, float *rcond, float *work, int *iwork, int *info);
 
/* Subroutine */ int stprfs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *nrhs, float *ap, float *b, int *ldb, float *x, int *ldx,
	 float *ferr, float *berr, float *work, int *iwork, int *info);
 
/* Subroutine */ int stptri_(const char *uplo, const char *diag, int *n, float *ap, 
	int *info);
 
/* Subroutine */ int stptrs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *nrhs, float *ap, float *b, int *ldb, int *info);
 
/* Subroutine */ int strcon_(const char *norm, const char *uplo, const char *diag, int *n, 
	float *a, int *lda, float *rcond, float *work, int *iwork, 
	int *info);
 
/* Subroutine */ int strevc_(const char *side, const char *howmny, logical *select, 
	int *n, float *t, int *ldt, float *vl, int *ldvl, float *vr, 
	int *ldvr, int *mm, int *m, float *work, int *info);
 
/* Subroutine */ int strexc_(const char *compq, int *n, float *t, int *ldt, 
	float *q, int *ldq, int *ifst, int *ilst, float *work, 
	int *info);
 
/* Subroutine */ int strrfs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *nrhs, float *a, int *lda, float *b, int *ldb, float *x, 
	int *ldx, float *ferr, float *berr, float *work, int *iwork, 
	int *info);
 
/* Subroutine */ int strsen_(const char *job, const char *compq, logical *select, int 
	*n, float *t, int *ldt, float *q, int *ldq, float *wr, float *wi, 
	int *m, float *s, float *sep, float *work, int *lwork, int *
	iwork, int *liwork, int *info);
 
/* Subroutine */ int strsna_(const char *job, const char *howmny, logical *select, 
	int *n, float *t, int *ldt, float *vl, int *ldvl, float *vr, 
	int *ldvr, float *s, float *sep, int *mm, int *m, float *
	work, int *ldwork, int *iwork, int *info);
 
/* Subroutine */ int strsyl_(const char *trana, const char *tranb, int *isgn, int 
	*m, int *n, float *a, int *lda, float *b, int *ldb, float *
	c__, int *ldc, float *scale, int *info);
 
/* Subroutine */ int strti2_(const char *uplo, const char *diag, int *n, float *a, 
	int *lda, int *info);
 
/* Subroutine */ int strtri_(const char *uplo, const char *diag, int *n, float *a, 
	int *lda, int *info);
 
/* Subroutine */ int strtrs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *nrhs, float *a, int *lda, float *b, int *ldb, int *
	info);
 
/* Subroutine */ int stzrqf_(int *m, int *n, float *a, int *lda, 
	float *tau, int *info);
 
/* Subroutine */ int stzrzf_(int *m, int *n, float *a, int *lda, 
	float *tau, float *work, int *lwork, int *info);
 
/* Subroutine */ int xerbla_(char *srname, int *info);
 
/* Subroutine */ int zbdsqr_(const char *uplo, int *n, int *ncvt, int *
	nru, int *ncc, double *d__, double *e, cdouble *vt, 
	int *ldvt, cdouble *u, int *ldu, cdouble *c__, 
	int *ldc, double *rwork, int *info);
 
/* Subroutine */ int zdrot_(int *n, cdouble *cx, int *incx, 
	cdouble *cy, int *incy, double *c__, double *s);
 
/* Subroutine */ int zdrscl_(int *n, double *sa, cdouble *sx, 
	int *incx);
 
/* Subroutine */ int zgbbrd_(const char *vect, int *m, int *n, int *ncc,
	 int *kl, int *ku, cdouble *ab, int *ldab, 
	double *d__, double *e, cdouble *q, int *ldq, 
	cdouble *pt, int *ldpt, cdouble *c__, int *ldc, 
	cdouble *work, double *rwork, int *info);
 
/* Subroutine */ int zgbcon_(const char *norm, int *n, int *kl, int *ku,
	 cdouble *ab, int *ldab, int *ipiv, double *anorm, 
	double *rcond, cdouble *work, double *rwork, int *
	info);
 
/* Subroutine */ int zgbequ_(int *m, int *n, int *kl, int *ku,
	 cdouble *ab, int *ldab, double *r__, double *c__, 
	double *rowcnd, double *colcnd, double *amax, int *
	info);
 
/* Subroutine */ int zgbrfs_(const char *trans, int *n, int *kl, int *
	ku, int *nrhs, cdouble *ab, int *ldab, cdouble *
	afb, int *ldafb, int *ipiv, cdouble *b, int *ldb, 
	cdouble *x, int *ldx, double *ferr, double *berr, 
	cdouble *work, double *rwork, int *info);
 
/* Subroutine */ int zgbsv_(int *n, int *kl, int *ku, int *
	nrhs, cdouble *ab, int *ldab, int *ipiv, cdouble *
	b, int *ldb, int *info);
 
/* Subroutine */ int zgbsvx_(const char *fact, const char *trans, int *n, int *kl,
	 int *ku, int *nrhs, cdouble *ab, int *ldab, 
	cdouble *afb, int *ldafb, int *ipiv, const char *equed, 
	double *r__, double *c__, cdouble *b, int *ldb, 
	cdouble *x, int *ldx, double *rcond, double *ferr, 
	double *berr, cdouble *work, double *rwork, int *
	info);
 
/* Subroutine */ int zgbtf2_(int *m, int *n, int *kl, int *ku,
	 cdouble *ab, int *ldab, int *ipiv, int *info);
 
/* Subroutine */ int zgbtrf_(int *m, int *n, int *kl, int *ku,
	 cdouble *ab, int *ldab, int *ipiv, int *info);
 
/* Subroutine */ int zgbtrs_(const char *trans, int *n, int *kl, int *
	ku, int *nrhs, cdouble *ab, int *ldab, int *ipiv, 
	cdouble *b, int *ldb, int *info);
 
/* Subroutine */ int zgebak_(const char *job, const char *side, int *n, int *ilo, 
	int *ihi, double *scale, int *m, cdouble *v, 
	int *ldv, int *info);
 
/* Subroutine */ int zgebal_(const char *job, int *n, cdouble *a, int 
	*lda, int *ilo, int *ihi, double *scale, int *info);
 
/* Subroutine */ int zgebd2_(int *m, int *n, cdouble *a, 
	int *lda, double *d__, double *e, cdouble *tauq, 
	cdouble *taup, cdouble *work, int *info);
 
/* Subroutine */ int zgebrd_(int *m, int *n, cdouble *a, 
	int *lda, double *d__, double *e, cdouble *tauq, 
	cdouble *taup, cdouble *work, int *lwork, int *
	info);
 
/* Subroutine */ int zgecon_(const char *norm, int *n, cdouble *a, 
	int *lda, double *anorm, double *rcond, cdouble *
	work, double *rwork, int *info);
 
/* Subroutine */ int zgeequ_(int *m, int *n, cdouble *a, 
	int *lda, double *r__, double *c__, double *rowcnd, 
	double *colcnd, double *amax, int *info);
 
/* Subroutine */ int zgees_(const char *jobvs, const char *sort, L_fp select, int *n, 
	cdouble *a, int *lda, int *sdim, cdouble *w, 
	cdouble *vs, int *ldvs, cdouble *work, int *lwork,
	 double *rwork, logical *bwork, int *info);
 
/* Subroutine */ int zgeesx_(const char *jobvs, const char *sort, L_fp select, const char *
	sense, int *n, cdouble *a, int *lda, int *sdim, 
	cdouble *w, cdouble *vs, int *ldvs, double *
	rconde, double *rcondv, cdouble *work, int *lwork, 
	double *rwork, logical *bwork, int *info);
 
/* Subroutine */ int zgeev_(const char *jobvl, const char *jobvr, int *n, 
	cdouble *a, int *lda, cdouble *w, cdouble *vl, 
	int *ldvl, cdouble *vr, int *ldvr, cdouble *work, 
	int *lwork, double *rwork, int *info);
 
/* Subroutine */ int zgeevx_(const char *balanc, const char *jobvl, const char *jobvr, const char *
	sense, int *n, cdouble *a, int *lda, cdouble *w, 
	cdouble *vl, int *ldvl, cdouble *vr, int *ldvr, 
	int *ilo, int *ihi, double *scale, double *abnrm, 
	double *rconde, double *rcondv, cdouble *work, int *
	lwork, double *rwork, int *info);
 
/* Subroutine */ int zgegs_(const char *jobvsl, const char *jobvsr, int *n, 
	cdouble *a, int *lda, cdouble *b, int *ldb, 
	cdouble *alpha, cdouble *beta, cdouble *vsl, 
	int *ldvsl, cdouble *vsr, int *ldvsr, cdouble *
	work, int *lwork, double *rwork, int *info);
 
/* Subroutine */ int zgegv_(const char *jobvl, const char *jobvr, int *n, 
	cdouble *a, int *lda, cdouble *b, int *ldb, 
	cdouble *alpha, cdouble *beta, cdouble *vl, int 
	*ldvl, cdouble *vr, int *ldvr, cdouble *work, int 
	*lwork, double *rwork, int *info);
 
/* Subroutine */ int zgehd2_(int *n, int *ilo, int *ihi, 
	cdouble *a, int *lda, cdouble *tau, cdouble *
	work, int *info);
 
/* Subroutine */ int zgehrd_(int *n, int *ilo, int *ihi, 
	cdouble *a, int *lda, cdouble *tau, cdouble *
	work, int *lwork, int *info);
 
/* Subroutine */ int zgelq2_(int *m, int *n, cdouble *a, 
	int *lda, cdouble *tau, cdouble *work, int *info);
 
/* Subroutine */ int zgelqf_(int *m, int *n, cdouble *a, 
	int *lda, cdouble *tau, cdouble *work, int *lwork,
	 int *info);
 
/* Subroutine */ int zgels_(const char *trans, int *m, int *n, int *
	nrhs, cdouble *a, int *lda, cdouble *b, int *ldb, 
	cdouble *work, int *lwork, int *info);
 
/* Subroutine */ int zgelsx_(int *m, int *n, int *nrhs, 
	cdouble *a, int *lda, cdouble *b, int *ldb, 
	int *jpvt, double *rcond, int *rank, cdouble *work, 
	double *rwork, int *info);
 
/* Subroutine */ int zgelsy_(int *m, int *n, int *nrhs, 
	cdouble *a, int *lda, cdouble *b, int *ldb, 
	int *jpvt, double *rcond, int *rank, cdouble *work, 
	int *lwork, double *rwork, int *info);
 
/* Subroutine */ int zgeql2_(int *m, int *n, cdouble *a, 
	int *lda, cdouble *tau, cdouble *work, int *info);
 
/* Subroutine */ int zgeqlf_(int *m, int *n, cdouble *a, 
	int *lda, cdouble *tau, cdouble *work, int *lwork,
	 int *info);
 
/* Subroutine */ int zgeqp3_(int *m, int *n, cdouble *a, 
	int *lda, int *jpvt, cdouble *tau, cdouble *work, 
	int *lwork, double *rwork, int *info);
 
/* Subroutine */ int zgeqpf_(int *m, int *n, cdouble *a, 
	int *lda, int *jpvt, cdouble *tau, cdouble *work, 
	double *rwork, int *info);
 
/* Subroutine */ int zgeqr2_(int *m, int *n, cdouble *a, 
	int *lda, cdouble *tau, cdouble *work, int *info);
 
/* Subroutine */ int zgeqrf_(int *m, int *n, cdouble *a, 
	int *lda, cdouble *tau, cdouble *work, int *lwork,
	 int *info);
 
/* Subroutine */ int zgerfs_(const char *trans, int *n, int *nrhs, 
	cdouble *a, int *lda, cdouble *af, int *ldaf, 
	int *ipiv, cdouble *b, int *ldb, cdouble *x, 
	int *ldx, double *ferr, double *berr, cdouble *work,
	 double *rwork, int *info);
 
/* Subroutine */ int zgerq2_(int *m, int *n, cdouble *a, 
	int *lda, cdouble *tau, cdouble *work, int *info);
 
/* Subroutine */ int zgerqf_(int *m, int *n, cdouble *a, 
	int *lda, cdouble *tau, cdouble *work, int *lwork,
	 int *info);
 
/* Subroutine */ int zgesc2_(int *n, cdouble *a, int *lda, 
	cdouble *rhs, int *ipiv, int *jpiv, double *scale);
 
/* Subroutine */ int zgesv_(int *n, int *nrhs, cdouble *a, 
	int *lda, int *ipiv, cdouble *b, int *ldb, int *
	info);
 
/* Subroutine */ int zgesvd_(const char *jobu, const char *jobvt, 
        int *m, int *n, cdouble *a, int *lda, double *s, 
        cdouble *u, int * ldu, cdouble *vt, int *ldvt, 
        cdouble *work, int *lwork, double *rwork, int *info);
 
/* Subroutine */ int zgesvx_(const char *fact, const char *trans, int *n, int *
	nrhs, cdouble *a, int *lda, cdouble *af, int *
	ldaf, int *ipiv, const char *equed, double *r__, double *c__, 
	cdouble *b, int *ldb, cdouble *x, int *ldx, 
	double *rcond, double *ferr, double *berr, cdouble *
	work, double *rwork, int *info);
 
/* Subroutine */ int zgetc2_(int *n, cdouble *a, int *lda, 
	int *ipiv, int *jpiv, int *info);
 
/* Subroutine */ int zgetf2_(int *m, int *n, cdouble *a, 
	int *lda, int *ipiv, int *info);
 
/* Subroutine */ int zgetrf_(int *m, int *n, cdouble *a, 
	int *lda, int *ipiv, int *info);
 
/* Subroutine */ int zgetri_(int *n, cdouble *a, int *lda, 
	int *ipiv, cdouble *work, int *lwork, int *info);
 
/* Subroutine */ int zgetrs_(const char *trans, int *n, int *nrhs, 
	cdouble *a, int *lda, int *ipiv, cdouble *b, 
	int *ldb, int *info);
 
/* Subroutine */ int zggbak_(const char *job, const char *side, int *n, int *ilo, 
	int *ihi, double *lscale, double *rscale, int *m, 
	cdouble *v, int *ldv, int *info);
 
/* Subroutine */ int zggbal_(const char *job, int *n, cdouble *a, int 
	*lda, cdouble *b, int *ldb, int *ilo, int *ihi, 
	double *lscale, double *rscale, double *work, int *
	info);
 
/* Subroutine */ int zgges_(const char *jobvsl, const char *jobvsr, const char *sort, L_fp 
	delctg, int *n, cdouble *a, int *lda, cdouble *b, 
	int *ldb, int *sdim, cdouble *alpha, cdouble *
	beta, cdouble *vsl, int *ldvsl, cdouble *vsr, int 
	*ldvsr, cdouble *work, int *lwork, double *rwork, 
	logical *bwork, int *info);
 
/* Subroutine */ int zggesx_(const char *jobvsl, const char *jobvsr, const char *sort, L_fp 
	delctg, const char *sense, int *n, cdouble *a, int *lda, 
	cdouble *b, int *ldb, int *sdim, cdouble *alpha, 
	cdouble *beta, cdouble *vsl, int *ldvsl, 
	cdouble *vsr, int *ldvsr, double *rconde, double *
	rcondv, cdouble *work, int *lwork, double *rwork, 
	int *iwork, int *liwork, logical *bwork, int *info);
 
/* Subroutine */ int zggev_(const char *jobvl, const char *jobvr, int *n, 
	cdouble *a, int *lda, cdouble *b, int *ldb, 
	cdouble *alpha, cdouble *beta, cdouble *vl, int 
	*ldvl, cdouble *vr, int *ldvr, cdouble *work, int 
	*lwork, double *rwork, int *info);
 
/* Subroutine */ int zggevx_(const char *balanc, const char *jobvl, const char *jobvr, const char *
	sense, int *n, cdouble *a, int *lda, cdouble *b, 
	int *ldb, cdouble *alpha, cdouble *beta, 
	cdouble *vl, int *ldvl, cdouble *vr, int *ldvr, 
	int *ilo, int *ihi, double *lscale, double *rscale, 
	double *abnrm, double *bbnrm, double *rconde, double *
	rcondv, cdouble *work, int *lwork, double *rwork, 
	int *iwork, logical *bwork, int *info);
 
/* Subroutine */ int zggglm_(int *n, int *m, int *p, 
	cdouble *a, int *lda, cdouble *b, int *ldb, 
	cdouble *d__, cdouble *x, cdouble *y, cdouble 
	*work, int *lwork, int *info);
 
/* Subroutine */ int zgghrd_(const char *compq, const char *compz, int *n, int *
	ilo, int *ihi, cdouble *a, int *lda, cdouble *b, 
	int *ldb, cdouble *q, int *ldq, cdouble *z__, 
	int *ldz, int *info);
 
/* Subroutine */ int zgglse_(int *m, int *n, int *p, 
	cdouble *a, int *lda, cdouble *b, int *ldb, 
	cdouble *c__, cdouble *d__, cdouble *x, 
	cdouble *work, int *lwork, int *info);
 
/* Subroutine */ int zggqrf_(int *n, int *m, int *p, 
	cdouble *a, int *lda, cdouble *taua, cdouble *b,
	 int *ldb, cdouble *taub, cdouble *work, int *
	lwork, int *info);
 
/* Subroutine */ int zggrqf_(int *m, int *p, int *n, 
	cdouble *a, int *lda, cdouble *taua, cdouble *b,
	 int *ldb, cdouble *taub, cdouble *work, int *
	lwork, int *info);
 
/* Subroutine */ int zggsvd_(const char *jobu, const char *jobv, const char *jobq, int *m, 
	int *n, int *p, int *k, int *l, cdouble *a, 
	int *lda, cdouble *b, int *ldb, double *alpha, 
	double *beta, cdouble *u, int *ldu, cdouble *v, 
	int *ldv, cdouble *q, int *ldq, cdouble *work, 
	double *rwork, int *iwork, int *info);
 
/* Subroutine */ int zggsvp_(const char *jobu, const char *jobv, const char *jobq, int *m, 
	int *p, int *n, cdouble *a, int *lda, cdouble 
	*b, int *ldb, double *tola, double *tolb, int *k, 
	int *l, cdouble *u, int *ldu, cdouble *v, int 
	*ldv, cdouble *q, int *ldq, int *iwork, double *
	rwork, cdouble *tau, cdouble *work, int *info);
 
/* Subroutine */ int zgtcon_(const char *norm, int *n, cdouble *dl, 
	cdouble *d__, cdouble *du, cdouble *du2, int *
	ipiv, double *anorm, double *rcond, cdouble *work, 
	int *info);
 
/* Subroutine */ int zgtrfs_(const char *trans, int *n, int *nrhs, 
	cdouble *dl, cdouble *d__, cdouble *du, 
	cdouble *dlf, cdouble *df, cdouble *duf, 
	cdouble *du2, int *ipiv, cdouble *b, int *ldb, 
	cdouble *x, int *ldx, double *ferr, double *berr, 
	cdouble *work, double *rwork, int *info);
 
/* Subroutine */ int zgtsv_(int *n, int *nrhs, cdouble *dl, 
	cdouble *d__, cdouble *du, cdouble *b, int *ldb,
	 int *info);
 
/* Subroutine */ int zgtsvx_(const char *fact, const char *trans, int *n, int *
	nrhs, cdouble *dl, cdouble *d__, cdouble *du, 
	cdouble *dlf, cdouble *df, cdouble *duf, 
	cdouble *du2, int *ipiv, cdouble *b, int *ldb, 
	cdouble *x, int *ldx, double *rcond, double *ferr, 
	double *berr, cdouble *work, double *rwork, int *
	info);
 
/* Subroutine */ int zgttrf_(int *n, cdouble *dl, cdouble *
	d__, cdouble *du, cdouble *du2, int *ipiv, int *
	info);
 
/* Subroutine */ int zgttrs_(const char *trans, int *n, int *nrhs, 
	cdouble *dl, cdouble *d__, cdouble *du, 
	cdouble *du2, int *ipiv, cdouble *b, int *ldb, 
	int *info);
 
/* Subroutine */ int zgtts2_(int *itrans, int *n, int *nrhs, 
	cdouble *dl, cdouble *d__, cdouble *du, 
	cdouble *du2, int *ipiv, cdouble *b, int *ldb);
 
/* Subroutine */ int zhbev_(const char *jobz, const char *uplo, int *n, int *kd, 
	cdouble *ab, int *ldab, double *w, cdouble *z__, 
	int *ldz, cdouble *work, double *rwork, int *info);
 
/* Subroutine */ int zhbevd_(const char *jobz, const char *uplo, int *n, int *kd, 
	cdouble *ab, int *ldab, double *w, cdouble *z__, 
	int *ldz, cdouble *work, int *lwork, double *rwork, 
	int *lrwork, int *iwork, int *liwork, int *info);
 
/* Subroutine */ int zhbevx_(const char *jobz, const char *range, const char *uplo, int *n, 
	int *kd, cdouble *ab, int *ldab, cdouble *q, 
	int *ldq, double *vl, double *vu, int *il, int *
	iu, double *abstol, int *m, double *w, cdouble *z__,
	 int *ldz, cdouble *work, double *rwork, int *iwork,
	 int *ifail, int *info);
 
/* Subroutine */ int zhbgst_(const char *vect, const char *uplo, int *n, int *ka, 
	int *kb, cdouble *ab, int *ldab, cdouble *bb, 
	int *ldbb, cdouble *x, int *ldx, cdouble *work, 
	double *rwork, int *info);
 
/* Subroutine */ int zhbgv_(const char *jobz, const char *uplo, int *n, int *ka, 
	int *kb, cdouble *ab, int *ldab, cdouble *bb, 
	int *ldbb, double *w, cdouble *z__, int *ldz, 
	cdouble *work, double *rwork, int *info);
 
/* Subroutine */ int zhbgvx_(const char *jobz, const char *range, const char *uplo, int *n, 
	int *ka, int *kb, cdouble *ab, int *ldab, 
	cdouble *bb, int *ldbb, cdouble *q, int *ldq, 
	double *vl, double *vu, int *il, int *iu, double *
	abstol, int *m, double *w, cdouble *z__, int *ldz, 
	cdouble *work, double *rwork, int *iwork, int *
	ifail, int *info);
 
/* Subroutine */ int zhbtrd_(const char *vect, const char *uplo, int *n, int *kd, 
	cdouble *ab, int *ldab, double *d__, double *e, 
	cdouble *q, int *ldq, cdouble *work, int *info);
 
/* Subroutine */ int zhecon_(const char *uplo, int *n, cdouble *a, 
	int *lda, int *ipiv, double *anorm, double *rcond, 
	cdouble *work, int *info);
 
/* Subroutine */ int zheev_(const char *jobz, const char *uplo, int *n, cdouble 
	*a, int *lda, double *w, cdouble *work, int *lwork, 
	double *rwork, int *info);
 
/* Subroutine */ int zheevd_(const char *jobz, const char *uplo, int *n, 
	cdouble *a, int *lda, double *w, cdouble *work, 
	int *lwork, double *rwork, int *lrwork, int *iwork, 
	int *liwork, int *info);
 
/* Subroutine */ int zheevr_(const char *jobz, const char *range, const char *uplo, int *n, 
	cdouble *a, int *lda, double *vl, double *vu, 
	int *il, int *iu, double *abstol, int *m, double *
	w, cdouble *z__, int *ldz, int *isuppz, cdouble *
	work, int *lwork, double *rwork, int *lrwork, int *
	iwork, int *liwork, int *info);
 
/* Subroutine */ int zheevx_(const char *jobz, const char *range, const char *uplo, int *n, 
	cdouble *a, int *lda, double *vl, double *vu, 
	int *il, int *iu, double *abstol, int *m, double *
	w, cdouble *z__, int *ldz, cdouble *work, int *
	lwork, double *rwork, int *iwork, int *ifail, int *
	info);
 
/* Subroutine */ int zhegs2_(int *itype, const char *uplo, int *n, 
	cdouble *a, int *lda, cdouble *b, int *ldb, 
	int *info);
 
/* Subroutine */ int zhegst_(int *itype, const char *uplo, int *n, 
	cdouble *a, int *lda, cdouble *b, int *ldb, 
	int *info);
 
/* Subroutine */ int zhegv_(int *itype, const char *jobz, const char *uplo, int *
	n, cdouble *a, int *lda, cdouble *b, int *ldb, 
	double *w, cdouble *work, int *lwork, double *rwork,
	 int *info);
 
/* Subroutine */ int zhegvd_(int *itype, const char *jobz, const char *uplo, int *
	n, cdouble *a, int *lda, cdouble *b, int *ldb, 
	double *w, cdouble *work, int *lwork, double *rwork,
	 int *lrwork, int *iwork, int *liwork, int *info);
 
/* Subroutine */ int zhegvx_(int *itype, const char *jobz, const char *range, const char *
	uplo, int *n, cdouble *a, int *lda, cdouble *b, 
	int *ldb, double *vl, double *vu, int *il, int *
	iu, double *abstol, int *m, double *w, cdouble *z__,
	 int *ldz, cdouble *work, int *lwork, double *rwork,
	 int *iwork, int *ifail, int *info);
 
/* Subroutine */ int zherfs_(const char *uplo, int *n, int *nrhs, 
	cdouble *a, int *lda, cdouble *af, int *ldaf, 
	int *ipiv, cdouble *b, int *ldb, cdouble *x, 
	int *ldx, double *ferr, double *berr, cdouble *work,
	 double *rwork, int *info);
 
/* Subroutine */ int zhesv_(const char *uplo, int *n, int *nrhs, 
	cdouble *a, int *lda, int *ipiv, cdouble *b, 
	int *ldb, cdouble *work, int *lwork, int *info);
 
/* Subroutine */ int zhesvx_(const char *fact, const char *uplo, int *n, int *
	nrhs, cdouble *a, int *lda, cdouble *af, int *
	ldaf, int *ipiv, cdouble *b, int *ldb, cdouble *x,
	 int *ldx, double *rcond, double *ferr, double *berr, 
	cdouble *work, int *lwork, double *rwork, int *info);
 
/* Subroutine */ int zhetf2_(const char *uplo, int *n, cdouble *a, 
	int *lda, int *ipiv, int *info);
 
/* Subroutine */ int zhetrd_(const char *uplo, int *n, cdouble *a, 
	int *lda, double *d__, double *e, cdouble *tau, 
	cdouble *work, int *lwork, int *info);
 
/* Subroutine */ int zhetrf_(const char *uplo, int *n, cdouble *a, 
	int *lda, int *ipiv, cdouble *work, int *lwork, 
	int *info);
 
/* Subroutine */ int zhetri_(const char *uplo, int *n, cdouble *a, 
	int *lda, int *ipiv, cdouble *work, int *info);
 
/* Subroutine */ int zhetrs_(const char *uplo, int *n, int *nrhs, 
	cdouble *a, int *lda, int *ipiv, cdouble *b, 
	int *ldb, int *info);
 
/* Subroutine */ int zhgeqz_(const char *job, const char *compq, const char *compz, int *n, 
	int *ilo, int *ihi, cdouble *a, int *lda, 
	cdouble *b, int *ldb, cdouble *alpha, cdouble *
	beta, cdouble *q, int *ldq, cdouble *z__, int *
	ldz, cdouble *work, int *lwork, double *rwork, int *
	info);
 
/* Subroutine */ int zhpcon_(const char *uplo, int *n, cdouble *ap, 
	int *ipiv, double *anorm, double *rcond, cdouble *
	work, int *info);
 
/* Subroutine */ int zhpev_(const char *jobz, const char *uplo, int *n, cdouble 
	*ap, double *w, cdouble *z__, int *ldz, cdouble *
	work, double *rwork, int *info);
 
/* Subroutine */ int zhpevd_(const char *jobz, const char *uplo, int *n, 
	cdouble *ap, double *w, cdouble *z__, int *ldz, 
	cdouble *work, int *lwork, double *rwork, int *
	lrwork, int *iwork, int *liwork, int *info);
 
/* Subroutine */ int zhpevx_(const char *jobz, const char *range, const char *uplo, int *n, 
	cdouble *ap, double *vl, double *vu, int *il, 
	int *iu, double *abstol, int *m, double *w, 
	cdouble *z__, int *ldz, cdouble *work, double *
	rwork, int *iwork, int *ifail, int *info);
 
/* Subroutine */ int zhpgst_(int *itype, const char *uplo, int *n, 
	cdouble *ap, cdouble *bp, int *info);
 
/* Subroutine */ int zhpgv_(int *itype, const char *jobz, const char *uplo, int *
	n, cdouble *ap, cdouble *bp, double *w, cdouble 
	*z__, int *ldz, cdouble *work, double *rwork, int *
	info);
 
/* Subroutine */ int zhpgvd_(int *itype, const char *jobz, const char *uplo, int *
	n, cdouble *ap, cdouble *bp, double *w, cdouble 
	*z__, int *ldz, cdouble *work, int *lwork, double *
	rwork, int *lrwork, int *iwork, int *liwork, int *
	info);
 
/* Subroutine */ int zhpgvx_(int *itype, const char *jobz, const char *range, const char *
	uplo, int *n, cdouble *ap, cdouble *bp, double *
	vl, double *vu, int *il, int *iu, double *abstol, 
	int *m, double *w, cdouble *z__, int *ldz, 
	cdouble *work, double *rwork, int *iwork, int *
	ifail, int *info);
 
/* Subroutine */ int zhprfs_(const char *uplo, int *n, int *nrhs, 
	cdouble *ap, cdouble *afp, int *ipiv, cdouble *
	b, int *ldb, cdouble *x, int *ldx, double *ferr, 
	double *berr, cdouble *work, double *rwork, int *
	info);
 
/* Subroutine */ int zhpsv_(const char *uplo, int *n, int *nrhs, 
	cdouble *ap, int *ipiv, cdouble *b, int *ldb, 
	int *info);
 
/* Subroutine */ int zhpsvx_(const char *fact, const char *uplo, int *n, int *
	nrhs, cdouble *ap, cdouble *afp, int *ipiv, 
	cdouble *b, int *ldb, cdouble *x, int *ldx, 
	double *rcond, double *ferr, double *berr, cdouble *
	work, double *rwork, int *info);
 
/* Subroutine */ int zhptrd_(const char *uplo, int *n, cdouble *ap, 
	double *d__, double *e, cdouble *tau, int *info);
 
/* Subroutine */ int zhptrf_(const char *uplo, int *n, cdouble *ap, 
	int *ipiv, int *info);
 
/* Subroutine */ int zhptri_(const char *uplo, int *n, cdouble *ap, 
	int *ipiv, cdouble *work, int *info);
 
/* Subroutine */ int zhptrs_(const char *uplo, int *n, int *nrhs, 
	cdouble *ap, int *ipiv, cdouble *b, int *ldb, 
	int *info);
 
/* Subroutine */ int zhsein_(const char *side, const char *eigsrc, const char *initv, logical *
	select, int *n, cdouble *h__, int *ldh, cdouble *
	w, cdouble *vl, int *ldvl, cdouble *vr, int *ldvr,
	 int *mm, int *m, cdouble *work, double *rwork, 
	int *ifaill, int *ifailr, int *info);
 
/* Subroutine */ int zhseqr_(const char *job, const char *compz, int *n, int *ilo,
	 int *ihi, cdouble *h__, int *ldh, cdouble *w, 
	cdouble *z__, int *ldz, cdouble *work, int *lwork,
	 int *info);
 
/* Subroutine */ int zlabrd_(int *m, int *n, int *nb, 
	cdouble *a, int *lda, double *d__, double *e, 
	cdouble *tauq, cdouble *taup, cdouble *x, int *
	ldx, cdouble *y, int *ldy);
 
/* Subroutine */ int zlacgv_(int *n, cdouble *x, int *incx);
 
/* Subroutine */ int zlacon_(int *n, cdouble *v, cdouble *x, 
	double *est, int *kase);
 
/* Subroutine */ int zlacp2_(const char *uplo, int *m, int *n, double *
	a, int *lda, cdouble *b, int *ldb);
 
/* Subroutine */ int zlacpy_(const char *uplo, int *m, int *n, 
	cdouble *a, int *lda, cdouble *b, int *ldb);
 
/* Subroutine */ int zlacrm_(int *m, int *n, cdouble *a, 
	int *lda, double *b, int *ldb, cdouble *c__, 
	int *ldc, double *rwork);
 
/* Subroutine */ int zlacrt_(int *n, cdouble *cx, int *incx, 
	cdouble *cy, int *incy, cdouble *c__, cdouble *
	s);
 
/* Subroutine */ int zlaed0_(int *qsiz, int *n, double *d__, 
	double *e, cdouble *q, int *ldq, cdouble *qstore, 
	int *ldqs, double *rwork, int *iwork, int *info);
 
/* Subroutine */ int zlaed7_(int *n, int *cutpnt, int *qsiz, 
	int *tlvls, int *curlvl, int *curpbm, double *d__, 
	cdouble *q, int *ldq, double *rho, int *indxq, 
	double *qstore, int *qptr, int *prmptr, int *perm, 
	int *givptr, int *givcol, double *givnum, cdouble *
	work, double *rwork, int *iwork, int *info);
 
/* Subroutine */ int zlaed8_(int *k, int *n, int *qsiz, 
	cdouble *q, int *ldq, double *d__, double *rho, 
	int *cutpnt, double *z__, double *dlamda, cdouble *
	q2, int *ldq2, double *w, int *indxp, int *indx, 
	int *indxq, int *perm, int *givptr, int *givcol, 
	double *givnum, int *info);
 
/* Subroutine */ int zlaein_(logical *rightv, logical *noinit, int *n, 
	cdouble *h__, int *ldh, cdouble *w, cdouble *v, 
	cdouble *b, int *ldb, double *rwork, double *eps3, 
	double *smlnum, int *info);
 
/* Subroutine */ int zlaesy_(cdouble *a, cdouble *b, 
	cdouble *c__, cdouble *rt1, cdouble *rt2, 
	cdouble *evscal, cdouble *cs1, cdouble *sn1);
 
/* Subroutine */ int zlaev2_(cdouble *a, cdouble *b, 
	cdouble *c__, double *rt1, double *rt2, double *cs1,
	 cdouble *sn1);
 
/* Subroutine */ int zlags2_(logical *upper, double *a1, cdouble *
	a2, double *a3, double *b1, cdouble *b2, double *b3,
	 double *csu, cdouble *snu, double *csv, cdouble *
	snv, double *csq, cdouble *snq);
 
/* Subroutine */ int zlagtm_(const char *trans, int *n, int *nrhs, 
	double *alpha, cdouble *dl, cdouble *d__, 
	cdouble *du, cdouble *x, int *ldx, double *beta, 
	cdouble *b, int *ldb);
 
/* Subroutine */ int zlahef_(const char *uplo, int *n, int *nb, int *kb,
	 cdouble *a, int *lda, int *ipiv, cdouble *w, 
	int *ldw, int *info);
 
/* Subroutine */ int zlahqr_(logical *wantt, logical *wantz, int *n, 
	int *ilo, int *ihi, cdouble *h__, int *ldh, 
	cdouble *w, int *iloz, int *ihiz, cdouble *z__, 
	int *ldz, int *info);
 
/* Subroutine */ int zlahrd_(int *n, int *k, int *nb, 
	cdouble *a, int *lda, cdouble *tau, cdouble *t, 
	int *ldt, cdouble *y, int *ldy);
 
/* Subroutine */ int zlaic1_(int *job, int *j, cdouble *x, 
	double *sest, cdouble *w, cdouble *gamma, double *
	sestpr, cdouble *s, cdouble *c__);
 
/* Subroutine */ int zlals0_(int *icompq, int *nl, int *nr, 
	int *sqre, int *nrhs, cdouble *b, int *ldb, 
	cdouble *bx, int *ldbx, int *perm, int *givptr, 
	int *givcol, int *ldgcol, double *givnum, int *ldgnum,
	 double *poles, double *difl, double *difr, double *
	z__, int *k, double *c__, double *s, double *rwork, 
	int *info);
 
/* Subroutine */ int zlalsa_(int *icompq, int *smlsiz, int *n, 
	int *nrhs, cdouble *b, int *ldb, cdouble *bx, 
	int *ldbx, double *u, int *ldu, double *vt, int *
	k, double *difl, double *difr, double *z__, double *
	poles, int *givptr, int *givcol, int *ldgcol, int *
	perm, double *givnum, double *c__, double *s, double *
	rwork, int *iwork, int *info);
 
/* Subroutine */ int zlapll_(int *n, cdouble *x, int *incx, 
	cdouble *y, int *incy, double *ssmin);
 
/* Subroutine */ int zlapmt_(logical *forwrd, int *m, int *n, 
	cdouble *x, int *ldx, int *k);
 
/* Subroutine */ int zlaqgb_(int *m, int *n, int *kl, int *ku,
	 cdouble *ab, int *ldab, double *r__, double *c__, 
	double *rowcnd, double *colcnd, double *amax, const char *equed);
 
/* Subroutine */ int zlaqge_(int *m, int *n, cdouble *a, 
	int *lda, double *r__, double *c__, double *rowcnd, 
	double *colcnd, double *amax, const char *equed);
 
/* Subroutine */ int zlaqhb_(const char *uplo, int *n, int *kd, 
	cdouble *ab, int *ldab, double *s, double *scond, 
	double *amax, const char *equed);
 
/* Subroutine */ int zlaqhe_(const char *uplo, int *n, cdouble *a, 
	int *lda, double *s, double *scond, double *amax, 
	const char *equed);
 
/* Subroutine */ int zlaqhp_(const char *uplo, int *n, cdouble *ap, 
	double *s, double *scond, double *amax, const char *equed);
 
/* Subroutine */ int zlaqp2_(int *m, int *n, int *offset, 
	cdouble *a, int *lda, int *jpvt, cdouble *tau, 
	double *vn1, double *vn2, cdouble *work);
 
/* Subroutine */ int zlaqps_(int *m, int *n, int *offset, int 
	*nb, int *kb, cdouble *a, int *lda, int *jpvt, 
	cdouble *tau, double *vn1, double *vn2, cdouble *
	auxv, cdouble *f, int *ldf);
 
/* Subroutine */ int zlaqsb_(const char *uplo, int *n, int *kd, 
	cdouble *ab, int *ldab, double *s, double *scond, 
	double *amax, const char *equed);
 
/* Subroutine */ int zlaqsp_(const char *uplo, int *n, cdouble *ap, 
	double *s, double *scond, double *amax, const char *equed);
 
/* Subroutine */ int zlaqsy_(const char *uplo, int *n, cdouble *a, 
	int *lda, double *s, double *scond, double *amax, 
	const char *equed);
 
/* Subroutine */ int zlar1v_(int *n, int *b1, int *bn, double 
	*sigma, double *d__, double *l, double *ld, double *
	lld, double *gersch, cdouble *z__, double *ztz, 
	double *mingma, int *r__, int *isuppz, double *work);
 
/* Subroutine */ int zlar2v_(int *n, cdouble *x, cdouble *y, 
	cdouble *z__, int *incx, double *c__, cdouble *s, 
	int *incc);
 
/* Subroutine */ int zlarcm_(int *m, int *n, double *a, int *
	lda, cdouble *b, int *ldb, cdouble *c__, int *ldc,
	 double *rwork);
 
/* Subroutine */ int zlarf_(const char *side, int *m, int *n, cdouble 
	*v, int *incv, cdouble *tau, cdouble *c__, int *
	ldc, cdouble *work);
 
/* Subroutine */ int zlarfb_(const char *side, const char *trans, const char *direct, const char *
	storev, int *m, int *n, int *k, cdouble *v, int 
	*ldv, cdouble *t, int *ldt, cdouble *c__, int *
	ldc, cdouble *work, int *ldwork);
 
/* Subroutine */ int zlarfg_(int *n, cdouble *alpha, cdouble *
	x, int *incx, cdouble *tau);
 
/* Subroutine */ int zlarft_(const char *direct, const char *storev, int *n, int *
	k, cdouble *v, int *ldv, cdouble *tau, cdouble *
	t, int *ldt);
 
/* Subroutine */ int zlarfx_(const char *side, int *m, int *n, 
	cdouble *v, cdouble *tau, cdouble *c__, int *
	ldc, cdouble *work);
 
/* Subroutine */ int zlargv_(int *n, cdouble *x, int *incx, 
	cdouble *y, int *incy, double *c__, int *incc);
 
/* Subroutine */ int zlarnv_(int *idist, int *iseed, int *n, 
	cdouble *x);
 
/* Subroutine */ int zlarrv_(int *n, double *d__, double *l, 
	int *isplit, int *m, double *w, int *iblock, 
	double *gersch, double *tol, cdouble *z__, int *ldz,
	 int *isuppz, double *work, int *iwork, int *info);
 
/* Subroutine */ int zlartg_(cdouble *f, cdouble *g, double *
	cs, cdouble *sn, cdouble *r__);
 
/* Subroutine */ int zlartv_(int *n, cdouble *x, int *incx, 
	cdouble *y, int *incy, double *c__, cdouble *s, 
	int *incc);
 
/* Subroutine */ int zlarz_(const char *side, int *m, int *n, int *l, 
	cdouble *v, int *incv, cdouble *tau, cdouble *
	c__, int *ldc, cdouble *work);
 
/* Subroutine */ int zlarzb_(const char *side, const char *trans, const char *direct, const char *
	storev, int *m, int *n, int *k, int *l, cdouble 
	*v, int *ldv, cdouble *t, int *ldt, cdouble *c__, 
	int *ldc, cdouble *work, int *ldwork);
 
/* Subroutine */ int zlarzt_(const char *direct, const char *storev, int *n, int *
	k, cdouble *v, int *ldv, cdouble *tau, cdouble *
	t, int *ldt);
 
/* Subroutine */ int zlascl_(const char *type__, int *kl, int *ku, 
	double *cfrom, double *cto, int *m, int *n, 
	cdouble *a, int *lda, int *info);
 
/* Subroutine */ int zlaset_(const char *uplo, int *m, int *n, 
	cdouble *alpha, cdouble *beta, cdouble *a, int *
	lda);
 
/* Subroutine */ int zlasr_(const char *side, const char *pivot, const char *direct, int *m,
	 int *n, double *c__, double *s, cdouble *a, 
	int *lda);
 
/* Subroutine */ int zlassq_(int *n, cdouble *x, int *incx, 
	double *scale, double *sumsq);
 
/* Subroutine */ int zlaswp_(int *n, cdouble *a, int *lda, 
	int *k1, int *k2, int *ipiv, int *incx);
 
/* Subroutine */ int zlasyf_(const char *uplo, int *n, int *nb, int *kb,
	 cdouble *a, int *lda, int *ipiv, cdouble *w, 
	int *ldw, int *info);
 
/* Subroutine */ int zlatbs_(const char *uplo, const char *trans, const char *diag, const char *
	normin, int *n, int *kd, cdouble *ab, int *ldab, 
	cdouble *x, double *scale, double *cnorm, int *info);
 
/* Subroutine */ int zlatdf_(int *ijob, int *n, cdouble *z__, 
	int *ldz, cdouble *rhs, double *rdsum, double *
	rdscal, int *ipiv, int *jpiv);
 
/* Subroutine */ int zlatps_(const char *uplo, const char *trans, const char *diag, const char *
	normin, int *n, cdouble *ap, cdouble *x, double *
	scale, double *cnorm, int *info);
 
/* Subroutine */ int zlatrd_(const char *uplo, int *n, int *nb, 
	cdouble *a, int *lda, double *e, cdouble *tau, 
	cdouble *w, int *ldw);
 
/* Subroutine */ int zlatrs_(const char *uplo, const char *trans, const char *diag, const char *
	normin, int *n, cdouble *a, int *lda, cdouble *x, 
	double *scale, double *cnorm, int *info);
 
/* Subroutine */ int zlatrz_(int *m, int *n, int *l, 
	cdouble *a, int *lda, cdouble *tau, cdouble *
	work);
 
/* Subroutine */ int zlatzm_(const char *side, int *m, int *n, 
	cdouble *v, int *incv, cdouble *tau, cdouble *
	c1, cdouble *c2, int *ldc, cdouble *work);
 
/* Subroutine */ int zlauu2_(const char *uplo, int *n, cdouble *a, 
	int *lda, int *info);
 
/* Subroutine */ int zlauum_(const char *uplo, int *n, cdouble *a, 
	int *lda, int *info);
 
/* Subroutine */ int zpbcon_(const char *uplo, int *n, int *kd, 
	cdouble *ab, int *ldab, double *anorm, double *
	rcond, cdouble *work, double *rwork, int *info);
 
/* Subroutine */ int zpbequ_(const char *uplo, int *n, int *kd, 
	cdouble *ab, int *ldab, double *s, double *scond, 
	double *amax, int *info);
 
/* Subroutine */ int zpbrfs_(const char *uplo, int *n, int *kd, int *
	nrhs, cdouble *ab, int *ldab, cdouble *afb, int *
	ldafb, cdouble *b, int *ldb, cdouble *x, int *ldx,
	 double *ferr, double *berr, cdouble *work, double *
	rwork, int *info);
 
/* Subroutine */ int zpbstf_(const char *uplo, int *n, int *kd, 
	cdouble *ab, int *ldab, int *info);
 
/* Subroutine */ int zpbsv_(const char *uplo, int *n, int *kd, int *
	nrhs, cdouble *ab, int *ldab, cdouble *b, int *
	ldb, int *info);
 
/* Subroutine */ int zpbsvx_(const char *fact, const char *uplo, int *n, int *kd, 
	int *nrhs, cdouble *ab, int *ldab, cdouble *afb, 
	int *ldafb, const char *equed, double *s, cdouble *b, int 
	*ldb, cdouble *x, int *ldx, double *rcond, double *
	ferr, double *berr, cdouble *work, double *rwork, 
	int *info);
 
/* Subroutine */ int zpbtf2_(const char *uplo, int *n, int *kd, 
	cdouble *ab, int *ldab, int *info);
 
/* Subroutine */ int zpbtrf_(const char *uplo, int *n, int *kd, 
	cdouble *ab, int *ldab, int *info);
 
/* Subroutine */ int zpbtrs_(const char *uplo, int *n, int *kd, int *
	nrhs, cdouble *ab, int *ldab, cdouble *b, int *
	ldb, int *info);
 
/* Subroutine */ int zpocon_(const char *uplo, int *n, cdouble *a, 
	int *lda, double *anorm, double *rcond, cdouble *
	work, double *rwork, int *info);
 
/* Subroutine */ int zpoequ_(int *n, cdouble *a, int *lda, 
	double *s, double *scond, double *amax, int *info);
 
/* Subroutine */ int zporfs_(const char *uplo, int *n, int *nrhs, 
	cdouble *a, int *lda, cdouble *af, int *ldaf, 
	cdouble *b, int *ldb, cdouble *x, int *ldx, 
	double *ferr, double *berr, cdouble *work, double *
	rwork, int *info);
 
/* Subroutine */ int zposv_(const char *uplo, int *n, int *nrhs, 
	cdouble *a, int *lda, cdouble *b, int *ldb, 
	int *info);
 
/* Subroutine */ int zposvx_(const char *fact, const char *uplo, int *n, int *
	nrhs, cdouble *a, int *lda, cdouble *af, int *
	ldaf, const char *equed, double *s, cdouble *b, int *ldb, 
	cdouble *x, int *ldx, double *rcond, double *ferr, 
	double *berr, cdouble *work, double *rwork, int *
	info);
 
/* Subroutine */ int zpotf2_(const char *uplo, int *n, cdouble *a, 
	int *lda, int *info);
 
/* Subroutine */ int zpotrf_(const char *uplo, int *n, cdouble *a, 
	int *lda, int *info);
 
/* Subroutine */ int zpotri_(const char *uplo, int *n, cdouble *a, 
	int *lda, int *info);
 
/* Subroutine */ int zpotrs_(const char *uplo, int *n, int *nrhs, 
	cdouble *a, int *lda, cdouble *b, int *ldb, 
	int *info);
 
/* Subroutine */ int zppcon_(const char *uplo, int *n, cdouble *ap, 
	double *anorm, double *rcond, cdouble *work, double 
	*rwork, int *info);
 
/* Subroutine */ int zppequ_(const char *uplo, int *n, cdouble *ap, 
	double *s, double *scond, double *amax, int *info);
 
/* Subroutine */ int zpprfs_(const char *uplo, int *n, int *nrhs, 
	cdouble *ap, cdouble *afp, cdouble *b, int *ldb,
	 cdouble *x, int *ldx, double *ferr, double *berr, 
	cdouble *work, double *rwork, int *info);
 
/* Subroutine */ int zppsv_(const char *uplo, int *n, int *nrhs, 
	cdouble *ap, cdouble *b, int *ldb, int *info);
 
/* Subroutine */ int zppsvx_(const char *fact, const char *uplo, int *n, int *
	nrhs, cdouble *ap, cdouble *afp, const char *equed, double *
	s, cdouble *b, int *ldb, cdouble *x, int *ldx, 
	double *rcond, double *ferr, double *berr, cdouble *
	work, double *rwork, int *info);
 
/* Subroutine */ int zpptrf_(const char *uplo, int *n, cdouble *ap, 
	int *info);
 
/* Subroutine */ int zpptri_(const char *uplo, int *n, cdouble *ap, 
	int *info);
 
/* Subroutine */ int zpptrs_(const char *uplo, int *n, int *nrhs, 
	cdouble *ap, cdouble *b, int *ldb, int *info);
 
/* Subroutine */ int zptcon_(int *n, double *d__, cdouble *e, 
	double *anorm, double *rcond, double *rwork, int *
	info);
 
/* Subroutine */ int zptrfs_(const char *uplo, int *n, int *nrhs, 
	double *d__, cdouble *e, double *df, cdouble *ef, 
	cdouble *b, int *ldb, cdouble *x, int *ldx, 
	double *ferr, double *berr, cdouble *work, double *
	rwork, int *info);
 
/* Subroutine */ int zptsv_(int *n, int *nrhs, double *d__, 
	cdouble *e, cdouble *b, int *ldb, int *info);
 
/* Subroutine */ int zptsvx_(const char *fact, int *n, int *nrhs, 
	double *d__, cdouble *e, double *df, cdouble *ef, 
	cdouble *b, int *ldb, cdouble *x, int *ldx, 
	double *rcond, double *ferr, double *berr, cdouble *
	work, double *rwork, int *info);
 
/* Subroutine */ int zpttrf_(int *n, double *d__, cdouble *e, 
	int *info);
 
/* Subroutine */ int zpttrs_(const char *uplo, int *n, int *nrhs, 
	double *d__, cdouble *e, cdouble *b, int *ldb, 
	int *info);
 
/* Subroutine */ int zptts2_(int *iuplo, int *n, int *nrhs, 
	double *d__, cdouble *e, cdouble *b, int *ldb);
 
/* Subroutine */ int zrot_(int *n, cdouble *cx, int *incx, 
	cdouble *cy, int *incy, double *c__, cdouble *s);
 
/* Subroutine */ int zspcon_(const char *uplo, int *n, cdouble *ap, 
	int *ipiv, double *anorm, double *rcond, cdouble *
	work, int *info);
 
/* Subroutine */ int zspmv_(const char *uplo, int *n, cdouble *alpha, 
	cdouble *ap, cdouble *x, int *incx, cdouble *
	beta, cdouble *y, int *incy);
 
/* Subroutine */ int zspr_(const char *uplo, int *n, cdouble *alpha, 
	cdouble *x, int *incx, cdouble *ap);
 
/* Subroutine */ int zsprfs_(const char *uplo, int *n, int *nrhs, 
	cdouble *ap, cdouble *afp, int *ipiv, cdouble *
	b, int *ldb, cdouble *x, int *ldx, double *ferr, 
	double *berr, cdouble *work, double *rwork, int *
	info);
 
/* Subroutine */ int zspsv_(const char *uplo, int *n, int *nrhs, 
	cdouble *ap, int *ipiv, cdouble *b, int *ldb, 
	int *info);
 
/* Subroutine */ int zspsvx_(const char *fact, const char *uplo, int *n, int *
	nrhs, cdouble *ap, cdouble *afp, int *ipiv, 
	cdouble *b, int *ldb, cdouble *x, int *ldx, 
	double *rcond, double *ferr, double *berr, cdouble *
	work, double *rwork, int *info);
 
/* Subroutine */ int zsptrf_(const char *uplo, int *n, cdouble *ap, 
	int *ipiv, int *info);
 
/* Subroutine */ int zsptri_(const char *uplo, int *n, cdouble *ap, 
	int *ipiv, cdouble *work, int *info);
 
/* Subroutine */ int zsptrs_(const char *uplo, int *n, int *nrhs, 
	cdouble *ap, int *ipiv, cdouble *b, int *ldb, 
	int *info);
 
/* Subroutine */ int zstedc_(const char *compz, int *n, double *d__, 
	double *e, cdouble *z__, int *ldz, cdouble *work, 
	int *lwork, double *rwork, int *lrwork, int *iwork, 
	int *liwork, int *info);
 
/* Subroutine */ int zstein_(int *n, double *d__, double *e, 
	int *m, double *w, int *iblock, int *isplit, 
	cdouble *z__, int *ldz, double *work, int *iwork, 
	int *ifail, int *info);
 
/* Subroutine */ int zsteqr_(const char *compz, int *n, double *d__, 
	double *e, cdouble *z__, int *ldz, double *work, 
	int *info);
 
/* Subroutine */ int zsycon_(const char *uplo, int *n, cdouble *a, 
	int *lda, int *ipiv, double *anorm, double *rcond, 
	cdouble *work, int *info);
 
/* Subroutine */ int zsymv_(const char *uplo, int *n, cdouble *alpha, 
	cdouble *a, int *lda, cdouble *x, int *incx, 
	cdouble *beta, cdouble *y, int *incy);
 
/* Subroutine */ int zsyr_(const char *uplo, int *n, cdouble *alpha, 
	cdouble *x, int *incx, cdouble *a, int *lda);
 
/* Subroutine */ int zsyrfs_(const char *uplo, int *n, int *nrhs, 
	cdouble *a, int *lda, cdouble *af, int *ldaf, 
	int *ipiv, cdouble *b, int *ldb, cdouble *x, 
	int *ldx, double *ferr, double *berr, cdouble *work,
	 double *rwork, int *info);
 
/* Subroutine */ int zsysv_(const char *uplo, int *n, int *nrhs, 
	cdouble *a, int *lda, int *ipiv, cdouble *b, 
	int *ldb, cdouble *work, int *lwork, int *info);
 
/* Subroutine */ int zsysvx_(const char *fact, const char *uplo, int *n, int *
	nrhs, cdouble *a, int *lda, cdouble *af, int *
	ldaf, int *ipiv, cdouble *b, int *ldb, cdouble *x,
	 int *ldx, double *rcond, double *ferr, double *berr, 
	cdouble *work, int *lwork, double *rwork, int *info);
 
/* Subroutine */ int zsytf2_(const char *uplo, int *n, cdouble *a, 
	int *lda, int *ipiv, int *info);
 
/* Subroutine */ int zsytrf_(const char *uplo, int *n, cdouble *a, 
	int *lda, int *ipiv, cdouble *work, int *lwork, 
	int *info);
 
/* Subroutine */ int zsytri_(const char *uplo, int *n, cdouble *a, 
	int *lda, int *ipiv, cdouble *work, int *info);
 
/* Subroutine */ int zsytrs_(const char *uplo, int *n, int *nrhs, 
	cdouble *a, int *lda, int *ipiv, cdouble *b, 
	int *ldb, int *info);
 
/* Subroutine */ int ztbcon_(const char *norm, const char *uplo, const char *diag, int *n, 
	int *kd, cdouble *ab, int *ldab, double *rcond, 
	cdouble *work, double *rwork, int *info);
 
/* Subroutine */ int ztbrfs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *kd, int *nrhs, cdouble *ab, int *ldab, 
	cdouble *b, int *ldb, cdouble *x, int *ldx, 
	double *ferr, double *berr, cdouble *work, double *
	rwork, int *info);
 
/* Subroutine */ int ztbtrs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *kd, int *nrhs, cdouble *ab, int *ldab, 
	cdouble *b, int *ldb, int *info);
 
/* Subroutine */ int ztgevc_(const char *side, const char *howmny, logical *select, 
	int *n, cdouble *a, int *lda, cdouble *b, int 
	*ldb, cdouble *vl, int *ldvl, cdouble *vr, int *
	ldvr, int *mm, int *m, cdouble *work, double *rwork,
	 int *info);
 
/* Subroutine */ int ztgex2_(logical *wantq, logical *wantz, int *n, 
	cdouble *a, int *lda, cdouble *b, int *ldb, 
	cdouble *q, int *ldq, cdouble *z__, int *ldz, 
	int *j1, int *info);
 
/* Subroutine */ int ztgexc_(logical *wantq, logical *wantz, int *n, 
	cdouble *a, int *lda, cdouble *b, int *ldb, 
	cdouble *q, int *ldq, cdouble *z__, int *ldz, 
	int *ifst, int *ilst, int *info);
 
/* Subroutine */ int ztgsen_(int *ijob, logical *wantq, logical *wantz, 
	logical *select, int *n, cdouble *a, int *lda, 
	cdouble *b, int *ldb, cdouble *alpha, cdouble *
	beta, cdouble *q, int *ldq, cdouble *z__, int *
	ldz, int *m, double *pl, double *pr, double *dif, 
	cdouble *work, int *lwork, int *iwork, int *liwork, 
	int *info);
 
/* Subroutine */ int ztgsja_(const char *jobu, const char *jobv, const char *jobq, int *m, 
	int *p, int *n, int *k, int *l, cdouble *a, 
	int *lda, cdouble *b, int *ldb, double *tola, 
	double *tolb, double *alpha, double *beta, cdouble *
	u, int *ldu, cdouble *v, int *ldv, cdouble *q, 
	int *ldq, cdouble *work, int *ncycle, int *info);
 
/* Subroutine */ int ztgsna_(const char *job, const char *howmny, logical *select, 
	int *n, cdouble *a, int *lda, cdouble *b, int 
	*ldb, cdouble *vl, int *ldvl, cdouble *vr, int *
	ldvr, double *s, double *dif, int *mm, int *m, 
	cdouble *work, int *lwork, int *iwork, int *info);
 
/* Subroutine */ int ztgsy2_(const char *trans, int *ijob, int *m, int *
	n, cdouble *a, int *lda, cdouble *b, int *ldb, 
	cdouble *c__, int *ldc, cdouble *d__, int *ldd, 
	cdouble *e, int *lde, cdouble *f, int *ldf, 
	double *scale, double *rdsum, double *rdscal, int *
	info);
 
/* Subroutine */ int ztgsyl_(const char *trans, int *ijob, int *m, int *
	n, cdouble *a, int *lda, cdouble *b, int *ldb, 
	cdouble *c__, int *ldc, cdouble *d__, int *ldd, 
	cdouble *e, int *lde, cdouble *f, int *ldf, 
	double *scale, double *dif, cdouble *work, int *
	lwork, int *iwork, int *info);
 
/* Subroutine */ int ztpcon_(const char *norm, const char *uplo, const char *diag, int *n, 
	cdouble *ap, double *rcond, cdouble *work, double 
	*rwork, int *info);
 
/* Subroutine */ int ztprfs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *nrhs, cdouble *ap, cdouble *b, int *ldb, 
	cdouble *x, int *ldx, double *ferr, double *berr, 
	cdouble *work, double *rwork, int *info);
 
/* Subroutine */ int ztptri_(const char *uplo, const char *diag, int *n, 
	cdouble *ap, int *info);
 
/* Subroutine */ int ztptrs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *nrhs, cdouble *ap, cdouble *b, int *ldb, 
	int *info);
 
/* Subroutine */ int ztrcon_(const char *norm, const char *uplo, const char *diag, int *n, 
	cdouble *a, int *lda, double *rcond, cdouble *
	work, double *rwork, int *info);
 
/* Subroutine */ int ztrevc_(const char *side, const char *howmny, logical *select, 
	int *n, cdouble *t, int *ldt, cdouble *vl, 
	int *ldvl, cdouble *vr, int *ldvr, int *mm, int 
	*m, cdouble *work, double *rwork, int *info);
 
/* Subroutine */ int ztrexc_(const char *compq, int *n, cdouble *t, 
	int *ldt, cdouble *q, int *ldq, int *ifst, int *
	ilst, int *info);
 
/* Subroutine */ int ztrrfs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *nrhs, cdouble *a, int *lda, cdouble *b, 
	int *ldb, cdouble *x, int *ldx, double *ferr, 
	double *berr, cdouble *work, double *rwork, int *
	info);
 
/* Subroutine */ int ztrsen_(const char *job, const char *compq, logical *select, int 
	*n, cdouble *t, int *ldt, cdouble *q, int *ldq, 
	cdouble *w, int *m, double *s, double *sep, 
	cdouble *work, int *lwork, int *info);
 
/* Subroutine */ int ztrsna_(const char *job, const char *howmny, logical *select, 
	int *n, cdouble *t, int *ldt, cdouble *vl, 
	int *ldvl, cdouble *vr, int *ldvr, double *s, 
	double *sep, int *mm, int *m, cdouble *work, 
	int *ldwork, double *rwork, int *info);
 
/* Subroutine */ int ztrsyl_(const char *trana, const char *tranb, int *isgn, int 
	*m, int *n, cdouble *a, int *lda, cdouble *b, 
	int *ldb, cdouble *c__, int *ldc, double *scale, 
	int *info);
 
/* Subroutine */ int ztrti2_(const char *uplo, const char *diag, int *n, 
	cdouble *a, int *lda, int *info);
 
/* Subroutine */ int ztrtri_(const char *uplo, const char *diag, int *n, 
	cdouble *a, int *lda, int *info);
 
/* Subroutine */ int ztrtrs_(const char *uplo, const char *trans, const char *diag, int *n, 
	int *nrhs, cdouble *a, int *lda, cdouble *b, 
	int *ldb, int *info);
 
/* Subroutine */ int ztzrqf_(int *m, int *n, cdouble *a, 
	int *lda, cdouble *tau, int *info);
 
/* Subroutine */ int ztzrzf_(int *m, int *n, cdouble *a, 
	int *lda, cdouble *tau, cdouble *work, int *lwork,
	 int *info);
 
/* Subroutine */ int zung2l_(int *m, int *n, int *k, 
	cdouble *a, int *lda, cdouble *tau, cdouble *
	work, int *info);
 
/* Subroutine */ int zung2r_(int *m, int *n, int *k, 
	cdouble *a, int *lda, cdouble *tau, cdouble *
	work, int *info);
 
/* Subroutine */ int zungbr_(const char *vect, int *m, int *n, int *k, 
	cdouble *a, int *lda, cdouble *tau, cdouble *
	work, int *lwork, int *info);
 
/* Subroutine */ int zunghr_(int *n, int *ilo, int *ihi, 
	cdouble *a, int *lda, cdouble *tau, cdouble *
	work, int *lwork, int *info);
 
/* Subroutine */ int zungl2_(int *m, int *n, int *k, 
	cdouble *a, int *lda, cdouble *tau, cdouble *
	work, int *info);
 
/* Subroutine */ int zunglq_(int *m, int *n, int *k, 
	cdouble *a, int *lda, cdouble *tau, cdouble *
	work, int *lwork, int *info);
 
/* Subroutine */ int zungql_(int *m, int *n, int *k, 
	cdouble *a, int *lda, cdouble *tau, cdouble *
	work, int *lwork, int *info);
 
/* Subroutine */ int zungqr_(int *m, int *n, int *k, 
	cdouble *a, int *lda, cdouble *tau, cdouble *
	work, int *lwork, int *info);
 
/* Subroutine */ int zungr2_(int *m, int *n, int *k, 
	cdouble *a, int *lda, cdouble *tau, cdouble *
	work, int *info);
 
/* Subroutine */ int zungrq_(int *m, int *n, int *k, 
	cdouble *a, int *lda, cdouble *tau, cdouble *
	work, int *lwork, int *info);
 
/* Subroutine */ int zungtr_(const char *uplo, int *n, cdouble *a, 
	int *lda, cdouble *tau, cdouble *work, int *lwork,
	 int *info);
 
/* Subroutine */ int zunm2l_(const char *side, const char *trans, int *m, int *n, 
	int *k, cdouble *a, int *lda, cdouble *tau, 
	cdouble *c__, int *ldc, cdouble *work, int *info);
 
/* Subroutine */ int zunm2r_(const char *side, const char *trans, int *m, int *n, 
	int *k, cdouble *a, int *lda, cdouble *tau, 
	cdouble *c__, int *ldc, cdouble *work, int *info);
 
/* Subroutine */ int zunmbr_(const char *vect, const char *side, const char *trans, int *m, 
	int *n, int *k, cdouble *a, int *lda, cdouble 
	*tau, cdouble *c__, int *ldc, cdouble *work, int *
	lwork, int *info);
 
/* Subroutine */ int zunmhr_(const char *side, const char *trans, int *m, int *n, 
	int *ilo, int *ihi, cdouble *a, int *lda, 
	cdouble *tau, cdouble *c__, int *ldc, cdouble *
	work, int *lwork, int *info);
 
/* Subroutine */ int zunml2_(const char *side, const char *trans, int *m, int *n, 
	int *k, cdouble *a, int *lda, cdouble *tau, 
	cdouble *c__, int *ldc, cdouble *work, int *info);
 
/* Subroutine */ int zunmlq_(const char *side, const char *trans, int *m, int *n, 
	int *k, cdouble *a, int *lda, cdouble *tau, 
	cdouble *c__, int *ldc, cdouble *work, int *lwork,
	 int *info);
 
/* Subroutine */ int zunmql_(const char *side, const char *trans, int *m, int *n, 
	int *k, cdouble *a, int *lda, cdouble *tau, 
	cdouble *c__, int *ldc, cdouble *work, int *lwork,
	 int *info);
 
/* Subroutine */ int zunmqr_(const char *side, const char *trans, int *m, int *n, 
	int *k, cdouble *a, int *lda, cdouble *tau, 
	cdouble *c__, int *ldc, cdouble *work, int *lwork,
	 int *info);
 
/* Subroutine */ int zunmr2_(const char *side, const char *trans, int *m, int *n, 
	int *k, cdouble *a, int *lda, cdouble *tau, 
	cdouble *c__, int *ldc, cdouble *work, int *info);
 
/* Subroutine */ int zunmr3_(const char *side, const char *trans, int *m, int *n, 
	int *k, int *l, cdouble *a, int *lda, cdouble 
	*tau, cdouble *c__, int *ldc, cdouble *work, int *
	info);
 
/* Subroutine */ int zunmrq_(const char *side, const char *trans, int *m, int *n, 
	int *k, cdouble *a, int *lda, cdouble *tau, 
	cdouble *c__, int *ldc, cdouble *work, int *lwork,
	 int *info);
 
/* Subroutine */ int zunmrz_(const char *side, const char *trans, int *m, int *n, 
	int *k, int *l, cdouble *a, int *lda, cdouble 
	*tau, cdouble *c__, int *ldc, cdouble *work, int *
	lwork, int *info);
 
/* Subroutine */ int zunmtr_(const char *side, const char *uplo, const char *trans, int *m, 
	int *n, cdouble *a, int *lda, cdouble *tau, 
	cdouble *c__, int *ldc, cdouble *work, int *lwork,
	 int *info);
 
/* Subroutine */ int zupgtr_(const char *uplo, int *n, cdouble *ap, 
	cdouble *tau, cdouble *q, int *ldq, cdouble *
	work, int *info);
 
/* Subroutine */ int zupmtr_(const char *side, const char *uplo, const char *trans, int *m, 
	int *n, cdouble *ap, cdouble *tau, cdouble *c__,
	 int *ldc, cdouble *work, int *info);

void dgemm_(const char *TRANSA, const char *TRANSB, int *M, int *N, int *K,
            double *ALPHA,double *A,int *LDA, double *B, int *LDB,
            double *BETA, double *C, int *LDC);

void dgemv_(const char *TRANS, int *M, int *N, double *Alpha,
            double *A, int *lda, double *X, int *incx, double *beta,
            double *Y, int *incy);

void zgemm_(const char *TRANSA, const char *TRANSB, int *M, int *N, int *K,
            cdouble *ALPHA, cdouble *A, int *LDA, cdouble *B, int *LDB,
            cdouble *BETA, cdouble *C, int *LDC);

void zgemv_(const char *TRANS, int *M, int *N, cdouble *Alpha,
            cdouble *A, int *lda, cdouble *X, int *incx, cdouble *beta,
            cdouble *Y, int *incy);

#endif /* __CLAPACK_H */

#ifdef __cplusplus
}
#endif
