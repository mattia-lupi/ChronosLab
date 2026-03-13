#ifndef EMIN_BLAS
#define EMIN_BLAS

#ifdef MATLAB_MEX_FILE
#include <mex.h>
#include <blas.h>
#include <lapack.h>
// To be consistent with Matlab's blas/lapack
using lapack_int = ptrdiff_t;
#else // MATLAB_MEX_FILE

#include <cblas.h>
#include <lapacke.h>

inline void dgemv(const char *trans, const lapack_int *m, const lapack_int *n, const double *alpha,
                  const double *a, const lapack_int *lda, const double *x, const lapack_int *incx,
                  const double *beta, double *y, const lapack_int *incy)
{
   CBLAS_TRANSPOSE cblas_trans = (*trans == 'T' || *trans == 't') ? CblasTrans : CblasNoTrans;
   cblas_dgemv(CblasColMajor, cblas_trans, *m, *n, *alpha, a, *lda, x, *incx, *beta, y, *incy);
}

inline void dgeqrf(const lapack_int *m, const lapack_int *n, double *a,
                   const lapack_int *lda, double *tau, double *work,
                   const lapack_int *lwork, lapack_int *err)
{
   *err = LAPACKE_dgeqrf_work(LAPACK_COL_MAJOR, *m, *n, a, *lda, tau, work, *lwork);
}

inline void dorgqr(const lapack_int *m, const lapack_int *n, const lapack_int *k,
                   double *a, const lapack_int *lda, const double *tau,
                   double *work, const lapack_int *lwork, lapack_int *err)
{
   *err = LAPACKE_dorgqr_work(LAPACK_COL_MAJOR, *m, *n, *k, a, *lda, tau, work, *lwork);
}

inline void dgesvd(const char *jobu, const char *jobvt, const lapack_int *m,
                   const lapack_int *n, double *a, const lapack_int *lda,
                   double *s, double *u, const lapack_int *ldu,
                   double *vt, const lapack_int *ldvt,
                   double *work, const lapack_int *lwork, lapack_int *err)
{
   *err = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, *jobu, *jobvt, *m, *n, a, *lda, s, u, *ldu,
                              vt, *ldvt, work, *lwork);
}

inline void dtrtrs(const char *uplo, const char *trans, const char *diag,
                   const lapack_int *n, const lapack_int *nrhs,
                   const double *a, const lapack_int *lda,
                   double *b, const lapack_int *ldb, lapack_int *err)
{
   *err = LAPACKE_dtrtrs(LAPACK_COL_MAJOR, *uplo, *trans, *diag, *n, *nrhs, a, *lda, b, *ldb);
}

#endif // MATLAB_MEX_FILE
#endif // EMIN_BLAS
