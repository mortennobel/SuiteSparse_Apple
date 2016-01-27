//
//  BlasWrapper.m
//  BlasWrapper
//
//  Created by Morten Nobel-Jørgensen on 26/01/16.
//  Copyright © 2016 Morten Nobel-Joergensen. All rights reserved.
//

#include "BlasWrapper.h"
#include <Accelerate/Accelerate.h>
#include <stdio.h>

#include <stdlib.h>

#define ORDER CblasColMajor

#define PRINT_DEBUG_X


// _dgemm_, _dgemv_, _dsyrk_, _dtrsm_, _dtrsv_, _zgemm_, _zgemv_, _zherk_, _ztrsm_, _ztrsv_

enum CBLAS_SIDE BLAS_GET_SIDE(char* X){
    return ((X[0] == 'R')?CblasRight:CblasLeft);
}
enum CBLAS_DIAG BLAS_GET_DIAG(char *X){
    return ((X[0] == 'N' ||X[0] == 'n')?CblasNonUnit:CblasUnit);
}

enum CBLAS_UPLO BLAS_GET_UPLO(char* X){
    return ((X[0] == 'L' || X[0] == 'l')?CblasLower:CblasUpper);
}

enum CBLAS_TRANSPOSE BLAS_GET_TRANSPOSE(char *X){
    return ((X[0] == 'N' || X[0] == 'n')?CblasNoTrans:((X[0] == 'C' || X[0] == 'c')?CblasConjTrans:CblasTrans));
}

void dgemm_(char *transa, char *transb, BLAS_INT *m, BLAS_INT *n,
                      BLAS_INT *k, double *alpha, double *A, BLAS_INT *lda, double *B,
            BLAS_INT *ldb, double *beta, double *C, BLAS_INT *ldc);

double sum(double *v, int length){
    double s = 0;
    for (int i=0;i<length;i++){
        s+= v[i];
    }
    return s;
}

void BLAS_DGEMM_WRAP (char *transa, char *transb, BLAS_INT *m, BLAS_INT *n,
                      BLAS_INT *k, double *alpha, double *A, BLAS_INT *lda, double *B,
                      BLAS_INT *ldb, double *beta, double *C, BLAS_INT *ldc)
{
    enum CBLAS_TRANSPOSE TransA = BLAS_GET_TRANSPOSE(transa);
    enum CBLAS_TRANSPOSE TransB = BLAS_GET_TRANSPOSE(transb);
    BLAS_INT M = *m;
    BLAS_INT N = *n;
    BLAS_INT K = *k;
    double Alpha = *alpha;
    BLAS_INT Lda = *lda;
    BLAS_INT Ldb = *ldb;
    double Beta = *beta;
    BLAS_INT Ldc = *ldc;
    
#ifdef USE_FORTRAN
    dgemm_(transa, transb, m, n,
          k, alpha, A, lda, B,
          ldb, beta, C, ldc);
#else
    cblas_dgemm(ORDER,TransA,TransB, M, N, K, Alpha, A, Lda, B, Ldb, Beta,C, Ldc);
#endif
#ifdef PRINT_DEBUG
    double s = sum(C, M*N);
    printf("%f dgemm(transa %d,transb %d, M %d, N %d, Alpha %f,A[0] %f, lda %d, B[0] %f, ldb %d, beta %f, C[0] %f, ldc %d)\n", s ,TransA,TransB, M, N,Alpha,A[0],Lda, B[0],Ldb, Beta, C[0], Ldc);
#endif
    
}


void dgemv_(char *trans, BLAS_INT *m, BLAS_INT *n, double *alpha,
           double *A, BLAS_INT *lda, double *X, BLAS_INT *incx, double *beta,
           double *Y, BLAS_INT *incy);

void BLAS_DGEMV_WRAP (char *trans, BLAS_INT *m, BLAS_INT *n, double *alpha,
                      double *A, BLAS_INT *lda, double *X, BLAS_INT *incx, double *beta,
                      double *Y, BLAS_INT *incy)
{
    enum CBLAS_TRANSPOSE Trans = BLAS_GET_TRANSPOSE(trans);
    BLAS_INT M = *m;
    BLAS_INT N = *n;
    double Alpha = *alpha;
    BLAS_INT Lda = *lda;
    BLAS_INT Incx = *incx;
    double Beta = *beta;
    BLAS_INT Incy = *incy;
#ifdef USE_FORTRAN
    dgemv_(trans, m, n, alpha,
          A, lda, X, incx, beta,
          Y, incy);
#else
    cblas_dgemv(ORDER,Trans,M, N, Alpha, A, Lda, X, Incx, Beta, Y, Incy);
#endif
#ifdef PRINT_DEBUG
    double s = sum(Y, M);
    printf("%f dgemv(trans %d, M %d, N %d, Alpha %f,A[0] %f, lda %d, X[0] %f, incx %d, beta %f, Y[0] %f, incy %d)\n", s, Trans, M, N,Alpha,A[0],Lda, X[0],Incx, Beta, Y[0], Incy);
#endif
}

void dsyrk_(char *uplo, char *trans, BLAS_INT *n, BLAS_INT *k,
           double *alpha, double *A, BLAS_INT *lda, double *beta, double *C,
           BLAS_INT *ldc);

void BLAS_DSYRK_WRAP (char *uplo, char *trans, BLAS_INT *n, BLAS_INT *k,
                      double *alpha, double *A, BLAS_INT *lda, double *beta, double *C,
                      BLAS_INT *ldc)
{
    enum CBLAS_UPLO Uplo = BLAS_GET_UPLO(uplo);
    enum CBLAS_TRANSPOSE Trans = BLAS_GET_TRANSPOSE(trans);
    BLAS_INT N = *n;
    BLAS_INT K = *k;
    double Alpha = *alpha;
    BLAS_INT Lda = *lda;
    double Beta = *beta;
    BLAS_INT Ldc = *ldc;
#ifdef USE_FORTRAN
    dsyrk_(uplo, trans, n, k,
          alpha, A, lda, beta, C,
          ldc);
#else
    cblas_dsyrk(ORDER,Uplo, Trans, N, K, Alpha, A,Lda, Beta,C,Ldc);
#endif
#ifdef PRINT_DEBUG
    double s = sum(C, N*N);
    printf("%f dsyrk(uplo %d, trans %d, N %d, K %d, Alpha %f,A[0] %f, lda %d, beta %f, C[0] %f, ldc %d)\n", s, Uplo, Trans, N,K, Alpha,A[0],Lda, Beta, C[0],Ldc);
#endif
}

void dtrsm_(char *side, char *uplo, char *transa, char *diag, BLAS_INT *m,
           BLAS_INT *n, double *alpha, double *A, BLAS_INT *lda, double *B,
           BLAS_INT *ldb);

void BLAS_DTRSM_WRAP (char *side, char *uplo, char *transa, char *diag, BLAS_INT *m,
                      BLAS_INT *n, double *alpha, double *A, BLAS_INT *lda, double *B,
                      BLAS_INT *ldb)
{
    enum CBLAS_SIDE Side = BLAS_GET_SIDE(side);
    enum CBLAS_UPLO Uplo = BLAS_GET_UPLO(uplo);
    enum CBLAS_TRANSPOSE TransA = BLAS_GET_TRANSPOSE(transa);
    enum CBLAS_DIAG Diag = BLAS_GET_DIAG(diag);
    BLAS_INT M = *m;
    BLAS_INT N = *n;
    double Alpha = *alpha;
    BLAS_INT Lda = *lda;
    BLAS_INT Ldb = *ldb;
#ifdef USE_FORTRAN
    dtrsm_(side, uplo, transa, diag, m,
          n, alpha, A, lda, B,
          ldb);
#else
    cblas_dtrsm(ORDER,Side, Uplo, TransA, Diag, M,N,Alpha,A,Lda,B,Ldb);
#endif
#ifdef PRINT_DEBUG
    double s = sum(B, M*N);
    printf("%f dtrsm(side %d, uplo %d, transa %d, diag %d, M %d, N %d, Alpha %f, A[0] %f, lda %d, B[0] %f, ldb %d)\n",s, Side, Uplo, TransA, Diag, M,N,Alpha,A[0],Lda,B[0],Ldb);
#endif
    
}

void dtrsv_(char *uplo, char *trans, char *diag, BLAS_INT *n, double *A,
           BLAS_INT *lda, double *X, BLAS_INT *incx);

void BLAS_DTRSV_WRAP (char *uplo, char *trans, char *diag, BLAS_INT *n, double *A,
                      BLAS_INT *lda, double *X, BLAS_INT *incx)
{
    enum CBLAS_UPLO Uplo = BLAS_GET_UPLO(uplo);
    enum CBLAS_TRANSPOSE TransA = BLAS_GET_TRANSPOSE(trans);
    enum CBLAS_DIAG Diag = BLAS_GET_DIAG(diag);
    BLAS_INT N = *n;
    BLAS_INT Lda = *lda;
    BLAS_INT IncX = *incx;
    
    double* tmpX = malloc(N*sizeof(double));
    double* tmpX2 = malloc(N*sizeof(double));
    for (int i=0;i<N;i++){ // todo memcpy
        tmpX[i] = X[i];
        tmpX2[i] = X[i];
    }
    dtrsv_(uplo, trans, diag, n, A,
          lda, X, incx);
    cblas_dtrsv(ORDER,Uplo, TransA, Diag, N, A, Lda, tmpX, IncX);
    bool equal = true;
    for (int i=0;i<N;i++){
        if (X[i] != tmpX[i]){
            equal = false;
        }
    }

    if (!equal){
        printf("\n\ndtrsv(uplo %d, trans %d, diag %d, N %d, A[0] %f, lda %d, X[0] %f, incx %d)\n", Uplo, TransA, Diag, N,A[0],Lda,X[0],IncX);
        
        printf("\nA: ");
        for (int i=0;i<N*N;i++){
            if (i%N == 0){
                printf("\n");
            }
            printf("%.40e,", A[i]);
        }
        printf("\n\n\nX: ");
        for (int i=0;i<N;i++){
            printf("%.40e,", tmpX2[i]);
        }
        
        printf("\n\ndtrsv: ");
        for (int i=0;i<N;i++){
            printf("%.20e ", X[i]);
        }
        printf("\n\n\ncblas_dtrsv: ");
        for (int i=0;i<N;i++){
            printf("%.20e ", tmpX[i]);
        }
        printf("\n\n\n");
    }
    
    free(tmpX);
    free(tmpX2);
#ifdef PRINT_DEBUG
    double s = sum(X, N*N);
    printf("%f dtrsv(uplo %d, trans %d, diag %d, N %d, A[0] %f, lda %d, X[0] %f, incx %d)\n", s, Uplo, TransA, Diag, N,A[0],Lda,X[0],IncX);
#endif
}

void BLAS_DGER_WRAP (BLAS_INT *m, BLAS_INT *n, double *alpha,
                     double *X, BLAS_INT *incx, double *Y, BLAS_INT *incy,
                     double *A, BLAS_INT *lda) {
    printf("BLAS_DGER not implemented");

}

void BLAS_DSCAL_WRAP (BLAS_INT *n, double *alpha, double *Y, BLAS_INT *incy) {
    printf("BLAS_DSCAL not implemented");
}


void BLAS_ZGEMM_WRAP (char *transa, char *transb, BLAS_INT *m, BLAS_INT *n,
                      BLAS_INT *k, double *alpha, double *A, BLAS_INT *lda, double *B,
                      BLAS_INT *ldb, double *beta, double *C, BLAS_INT *ldc)
{
    enum CBLAS_TRANSPOSE TransA = BLAS_GET_TRANSPOSE(transa);
    enum CBLAS_TRANSPOSE TransB = BLAS_GET_TRANSPOSE(transb);
    BLAS_INT M = *m;
    BLAS_INT N = *n;
    BLAS_INT K = *k;

    BLAS_INT Lda = *lda;
    BLAS_INT Ldb = *ldb;

    BLAS_INT Ldc = *ldc;
    cblas_zgemm(ORDER,TransA,TransB, M, N, K, alpha, A, Lda, B, Ldb, beta,C, Ldc);
    printf("ZGEMM not debugable");
}

void BLAS_ZGEMV_WRAP (char *trans, BLAS_INT *m, BLAS_INT *n, double *alpha,
                      double *A, BLAS_INT *lda, double *X, BLAS_INT *incx, double *beta,
                      double *Y, BLAS_INT *incy)
{
    enum CBLAS_TRANSPOSE Trans = BLAS_GET_TRANSPOSE(trans);
    BLAS_INT M = *m;
    BLAS_INT N = *n;

    BLAS_INT Lda = *lda;
    BLAS_INT Incx = *incx;

    BLAS_INT Incy = *incy;
    
    cblas_zgemv(ORDER,Trans,M, N, alpha, A, Lda, X, Incx, beta, Y, Incy);
    printf("ZGEMV not debugable");
}

void BLAS_ZHERK_WRAP (char *uplo, char *trans, BLAS_INT *n, BLAS_INT *k,
                      double *alpha, double *A, BLAS_INT *lda, double *beta, double *C,
                      BLAS_INT *ldc)
{
    enum CBLAS_UPLO Uplo = BLAS_GET_UPLO(uplo);
    enum CBLAS_TRANSPOSE Trans = BLAS_GET_TRANSPOSE(trans);
    BLAS_INT N = *n;
    BLAS_INT K = *k;
    BLAS_INT Lda = *lda;
    BLAS_INT Ldc = *ldc;
    cblas_zsyrk(ORDER,Uplo, Trans, N, K, alpha, A,Lda, beta,C,Ldc);
    printf("ZHERK not debugable");
}

void BLAS_ZTRSM_WRAP (char *side, char *uplo, char *transa, char *diag, BLAS_INT *m,
                      BLAS_INT *n, double *alpha, double *A, BLAS_INT *lda, double *B,
                      BLAS_INT *ldb)
{
    enum CBLAS_SIDE Side = BLAS_GET_SIDE(side);
    enum CBLAS_UPLO Uplo = BLAS_GET_UPLO(uplo);
    enum CBLAS_TRANSPOSE TransA = BLAS_GET_TRANSPOSE(transa);
    enum CBLAS_DIAG Diag = BLAS_GET_DIAG(diag);
    BLAS_INT M = *m;
    BLAS_INT N = *n;
    BLAS_INT Lda = *lda;
    BLAS_INT Ldb = *ldb;
    cblas_ztrsm(ORDER,Side, Uplo, TransA, Diag, M,N,alpha,A,Lda,B,Ldb);
    printf("ZTRSM not debugable");
}

void BLAS_ZTRSV_WRAP (char *uplo, char *trans, char *diag, BLAS_INT *n, double *A,
                      BLAS_INT *lda, double *X, BLAS_INT *incx)
{
    enum CBLAS_UPLO Uplo = BLAS_GET_UPLO(uplo);
    enum CBLAS_TRANSPOSE TransA = BLAS_GET_TRANSPOSE(trans);
    enum CBLAS_DIAG Diag = BLAS_GET_DIAG(diag);
    BLAS_INT N = *n;
    BLAS_INT Lda = *lda;
    BLAS_INT IncX = *incx;
    cblas_ztrsv(ORDER,Uplo, TransA, Diag, N, A, Lda, X, IncX);
    printf("ZTRSV not debugable");
}

void BLAS_ZGER_WRAP (BLAS_INT *m, BLAS_INT *n, double *alpha,
                     double *X, BLAS_INT *incx, double *Y, BLAS_INT *incy,
                     double *A, BLAS_INT *lda) {
    printf("BLAS_ZGER not implemented\n");
}

void BLAS_ZSCAL_WRAP (BLAS_INT *n, double *alpha, double *Y, BLAS_INT *incy) {
    printf("BLAS_ZSCAL_WRAP not implemented\n");
}
