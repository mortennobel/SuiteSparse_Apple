//
//  BlasWrapper.h
//  BlasWrapper
//
//  Created by Morten Nobel-Jørgensen on 26/01/16.
//  Copyright © 2016 Morten Nobel-Joergensen. All rights reserved.
//

#define BLAS_INT int
// _dgemm_, _dgemv_, _dsyrk_, _dtrsm_, _dtrsv_, _zgemm_, _zgemv_, _zherk_, _ztrsm_, _ztrsv_

void BLAS_DGEMM_WRAP (char *transa, char *transb, BLAS_INT *m, BLAS_INT *n,
                 BLAS_INT *k, double *alpha, double *A, BLAS_INT *lda, double *B,
                 BLAS_INT *ldb, double *beta, double *C, BLAS_INT *ldc) ;

void BLAS_DGEMV_WRAP (char *trans, BLAS_INT *m, BLAS_INT *n, double *alpha,
                 double *A, BLAS_INT *lda, double *X, BLAS_INT *incx, double *beta,
                 double *Y, BLAS_INT *incy) ;

void BLAS_DSYRK_WRAP (char *uplo, char *trans, BLAS_INT *n, BLAS_INT *k,
                 double *alpha, double *A, BLAS_INT *lda, double *beta, double *C,
                 BLAS_INT *ldc) ;

void BLAS_DTRSM_WRAP (char *side, char *uplo, char *transa, char *diag, BLAS_INT *m,
                 BLAS_INT *n, double *alpha, double *A, BLAS_INT *lda, double *B,
                 BLAS_INT *ldb) ;

void BLAS_DTRSV_WRAP (char *uplo, char *trans, char *diag, BLAS_INT *n, double *A,
                 BLAS_INT *lda, double *X, BLAS_INT *incx) ;

void BLAS_DGER_WRAP (BLAS_INT *m, BLAS_INT *n, double *alpha,
                     double *X, BLAS_INT *incx, double *Y, BLAS_INT *incy,
                     double *A, BLAS_INT *lda) ;

void BLAS_DSCAL_WRAP (BLAS_INT *n, double *alpha, double *Y, BLAS_INT *incy) ;


void BLAS_ZGEMM_WRAP (char *transa, char *transb, BLAS_INT *m, BLAS_INT *n,
                 BLAS_INT *k, double *alpha, double *A, BLAS_INT *lda, double *B,
                 BLAS_INT *ldb, double *beta, double *C, BLAS_INT *ldc) ;

void BLAS_ZGEMV_WRAP (char *trans, BLAS_INT *m, BLAS_INT *n, double *alpha,
                 double *A, BLAS_INT *lda, double *X, BLAS_INT *incx, double *beta,
                 double *Y, BLAS_INT *incy) ;

void BLAS_ZHERK_WRAP (char *uplo, char *trans, BLAS_INT *n, BLAS_INT *k,
                 double *alpha, double *A, BLAS_INT *lda, double *beta, double *C,
                 BLAS_INT *ldc) ;

void BLAS_ZTRSM_WRAP (char *side, char *uplo, char *transa, char *diag, BLAS_INT *m,
                 BLAS_INT *n, double *alpha, double *A, BLAS_INT *lda, double *B,
                 BLAS_INT *ldb) ;

void BLAS_ZTRSV_WRAP (char *uplo, char *trans, char *diag, BLAS_INT *n, double *A,
                 BLAS_INT *lda, double *X, BLAS_INT *incx) ;

void BLAS_ZGER_WRAP (BLAS_INT *m, BLAS_INT *n, double *alpha,
                     double *X, BLAS_INT *incx, double *Y, BLAS_INT *incy,
                     double *A, BLAS_INT *lda) ;

void BLAS_ZSCAL_WRAP (BLAS_INT *n, double *alpha, double *Y, BLAS_INT *incy) ;
