/* ========================================================================== */
/* === Include/cholmod_blas.h =============================================== */
/* ========================================================================== */

/* -----------------------------------------------------------------------------
 * CHOLMOD/Include/cholmod_blas.h.
 * Copyright (C) 2005-2006, Univ. of Florida.  Author: Timothy A. Davis
 * CHOLMOD/Include/cholmod_blas.h is licensed under Version 2.1 of the GNU
 * Lesser General Public License.  See lesser.txt for a text of the license.
 * CHOLMOD is also available under other licenses; contact authors for details.
 * -------------------------------------------------------------------------- */

/* This does not need to be included in the user's program. */

#ifndef CHOLMOD_BLAS_H
#define CHOLMOD_BLAS_H

// #define FORCE_OSX_CBLAS_PREFIX
#define FORCE_WRAP

/* ========================================================================== */
/* === Architecture ========================================================= */
/* ========================================================================== */

#if defined (__sun) || defined (MSOL2) || defined (ARCH_SOL2)
#define CHOLMOD_SOL2
#define CHOLMOD_ARCHITECTURE "Sun Solaris"

#elif defined (__sgi) || defined (MSGI) || defined (ARCH_SGI)
#define CHOLMOD_SGI
#define CHOLMOD_ARCHITECTURE "SGI Irix"

#elif defined (__linux) || defined (MGLNX86) || defined (ARCH_GLNX86)
#define CHOLMOD_LINUX
#define CHOLMOD_ARCHITECTURE "Linux"

#elif defined (__APPLE__)
#define CHOLMOD_MAC
#define CHOLMOD_ARCHITECTURE "Mac"

#elif defined (_AIX) || defined (MIBM_RS) || defined (ARCH_IBM_RS)
#define CHOLMOD_AIX
#define CHOLMOD_ARCHITECTURE "IBM AIX"
/* recent reports from IBM AIX seem to indicate that this is not needed: */
/* #define BLAS_NO_UNDERSCORE */

#elif defined (__alpha) || defined (MALPHA) || defined (ARCH_ALPHA)
#define CHOLMOD_ALPHA
#define CHOLMOD_ARCHITECTURE "Compaq Alpha"

#elif defined (_WIN32) || defined (WIN32) || defined (_WIN64) || defined (WIN64)
#if defined (__MINGW32__) || defined (__MINGW32__)
#define CHOLMOD_MINGW
#elif defined (__CYGWIN32__) || defined (__CYGWIN32__)
#define CHOLMOD_CYGWIN
#else
#define CHOLMOD_WINDOWS
#define BLAS_NO_UNDERSCORE
#endif
#define CHOLMOD_ARCHITECTURE "Microsoft Windows"

#elif defined (__hppa) || defined (__hpux) || defined (MHPUX) || defined (ARCH_HPUX)
#define CHOLMOD_HP
#define CHOLMOD_ARCHITECTURE "HP Unix"
#define BLAS_NO_UNDERSCORE

#elif defined (__hp700) || defined (MHP700) || defined (ARCH_HP700)
#define CHOLMOD_HP
#define CHOLMOD_ARCHITECTURE "HP 700 Unix"
#define BLAS_NO_UNDERSCORE

#else
/* If the architecture is unknown, and you call the BLAS, you may need to */
/* define BLAS_BY_VALUE, BLAS_NO_UNDERSCORE, and/or BLAS_CHAR_ARG yourself. */
#define CHOLMOD_ARCHITECTURE "unknown"
#endif

/* ========================================================================== */
/* === BLAS and LAPACK names ================================================ */
/* ========================================================================== */

/* Prototypes for the various versions of the BLAS.  */

/* Determine if the 64-bit Sun Performance BLAS is to be used */
#if defined(CHOLMOD_SOL2) && !defined(NSUNPERF) && defined(BLAS64)
#define SUN64
#endif

#ifdef SUN64

#define BLAS_DTRSV dtrsv_64_
#define BLAS_DGEMV dgemv_64_
#define BLAS_DTRSM dtrsm_64_
#define BLAS_DGEMM dgemm_64_
#define BLAS_DSYRK dsyrk_64_
#define BLAS_DGER  dger_64_
#define BLAS_DSCAL dscal_64_
#define LAPACK_DPOTRF dpotrf_64_

#define BLAS_ZTRSV ztrsv_64_
#define BLAS_ZGEMV zgemv_64_
#define BLAS_ZTRSM ztrsm_64_
#define BLAS_ZGEMM zgemm_64_
#define BLAS_ZHERK zherk_64_
#define BLAS_ZGER  zgeru_64_
#define BLAS_ZSCAL zscal_64_
#define LAPACK_ZPOTRF zpotrf_64_

#elif defined(FORCE_WRAP)

#define BLAS_DTRSV BLAS_DTRSV_WRAP
#define BLAS_DGEMV BLAS_DGEMV_WRAP
#define BLAS_DTRSM BLAS_DTRSM_WRAP
#define BLAS_DGEMM BLAS_DGEMM_WRAP
#define BLAS_DSYRK BLAS_DSYRK_WRAP
#define BLAS_DGER  BLAS_DGER_WRAP
#define BLAS_DSCAL BLAS_DSCAL_WRAP
#define LAPACK_DPOTRF dpotrf_

#define BLAS_ZTRSV BLAS_ZTRSV_WRAP
#define BLAS_ZGEMV BLAS_ZGEMV_WRAP
#define BLAS_ZTRSM BLAS_ZTRSM_WRAP
#define BLAS_ZGEMM BLAS_ZGEMM_WRAP
#define BLAS_ZHERK BLAS_ZHERK_WRAP
#define BLAS_ZGER  BLAS_ZGER_WRAP
#define BLAS_ZSCAL BLAS_ZSCAL_WRAP
#define LAPACK_ZPOTRF zpotrf_

#elif defined (CHOLMOD_MAC) && defined (FORCE_OSX_CBLAS_PREFIX)

#define BLAS_DTRSV cblas_dtrsv
#define BLAS_DGEMV cblas_dgemv
#define BLAS_DTRSM cblas_dtrsm
#define BLAS_DGEMM cblas_dgemm
#define BLAS_DSYRK cblas_dsyrk
#define BLAS_DGER  cblas_dger
#define BLAS_DSCAL cblas_dscal
#define LAPACK_DPOTRF dpotrf_

#define BLAS_ZTRSV cblas_ztrsv
#define BLAS_ZGEMV cblas_zgemv
#define BLAS_ZTRSM cblas_ztrsm
#define BLAS_ZGEMM cblas_zgemm
#define BLAS_ZHERK cblas_zherk
#define BLAS_ZGER  cblas_zgeru
#define BLAS_ZSCAL cblas_zscal
#define LAPACK_ZPOTRF zpotrf_

#elif defined (BLAS_NO_UNDERSCORE)

#define BLAS_DTRSV dtrsv
#define BLAS_DGEMV dgemv
#define BLAS_DTRSM dtrsm
#define BLAS_DGEMM dgemm
#define BLAS_DSYRK dsyrk
#define BLAS_DGER  dger
#define BLAS_DSCAL dscal
#define LAPACK_DPOTRF dpotrf

#define BLAS_ZTRSV ztrsv
#define BLAS_ZGEMV zgemv
#define BLAS_ZTRSM ztrsm
#define BLAS_ZGEMM zgemm
#define BLAS_ZHERK zherk
#define BLAS_ZGER  zgeru
#define BLAS_ZSCAL zscal
#define LAPACK_ZPOTRF zpotrf

#else

#define BLAS_DTRSV dtrsv_
#define BLAS_DGEMV dgemv_
#define BLAS_DTRSM dtrsm_
#define BLAS_DGEMM dgemm_
#define BLAS_DSYRK dsyrk_
#define BLAS_DGER  dger_
#define BLAS_DSCAL dscal_
#define LAPACK_DPOTRF dpotrf_

#define BLAS_ZTRSV ztrsv_
#define BLAS_ZGEMV zgemv_
#define BLAS_ZTRSM ztrsm_
#define BLAS_ZGEMM zgemm_
#define BLAS_ZHERK zherk_
#define BLAS_ZGER  zgeru_
#define BLAS_ZSCAL zscal_
#define LAPACK_ZPOTRF zpotrf_

#endif

/* ========================================================================== */
/* === BLAS and LAPACK integer arguments ==================================== */
/* ========================================================================== */

/* Compile CHOLMOD, UMFPACK, and SPQR with -DBLAS64 if you have a BLAS that
 * uses 64-bit integers */

#if defined (LONGBLAS) || defined (BLAS64)
#define BLAS_INT SuiteSparse_long
#else
#define BLAS_INT int
#endif

/* If the BLAS integer is smaller than the basic CHOLMOD integer, then we need
 * to check for integer overflow when converting from Int to BLAS_INT.  If
 * any integer overflows, the externally-defined BLAS_OK variable is
 * set to FALSE.  BLAS_OK should be set to TRUE before calling any
 * BLAS_* macro.
 */



#define CHECK_BLAS_INT (sizeof (BLAS_INT) < sizeof (Int))
#define EQ(K,k) (((BLAS_INT) K) == ((Int) k))

#if defined (CHOLMOD_MAC) && defined (FORCE_OSX_CBLAS_PREFIX) //************-----------------------------------------------------

//#import "/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Headers/cblas.h"
//#import "/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Headers/clapack.h"

//enum X_CBLAS_ORDER     {CblasRowMajor=101, CblasColMajor=102};
//enum X_CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113,AtlasConj=114};
//enum X_CBLAS_UPLO      {CblasUpper=121, CblasLower=122};
//enum X_CBLAS_DIAG      {CblasNonUnit=131, CblasUnit=132};
//enum X_CBLAS_SIDE      {CblasLeft=141, CblasRight=142};

/* column major order */
#define BLAS_ORDER CblasColMajor
// side,uplo,transa,diag
#define BLAS_GET_SIDE(X) ((X[0] == 'R')?CblasRight:CblasLeft)
#define BLAS_GET_DIAG(X) ((X[0] == 'N' ||X[0] == 'n')?CblasNonUnit:CblasUnit)
#define BLAS_GET_UPLO(X) ((X[0] == 'L' || X[0] == 'l')?CblasLower:CblasUpper)
#define BLAS_GET_TRANSPOSE(X) ((X[0] == 'N' || X[0] == 'n')?CblasNoTrans:((X[0] == 'C' || X[0] == 'c')?CblasConjTrans:CblasTrans))

/* ========================================================================== */
/* === BLAS and LAPACK prototypes and macros ================================ */
/* ========================================================================== */

//void BLAS_DGEMV (enum X_CBLAS_ORDER order, enum X_CBLAS_TRANSPOSE trans, BLAS_INT m, BLAS_INT n, double *alpha,
//	double *A, BLAS_INT lda, double *X, BLAS_INT incx, double *beta,
//	double *Y, BLAS_INT incy) ;
//void cblas_dgemv(const enum X_CBLAS_ORDER __Order,
//        const enum X_CBLAS_TRANSPOSE __TransA, const BLAS_INT __M, const BLAS_INT __N,
//        const double __alpha, const double *__A, const BLAS_INT __lda,
//        const double *__X, const BLAS_INT __incX, const double __beta, double *__Y,
//        const BLAS_INT __incY);
#define BLAS_dgemv(trans,m,n,alpha,A,lda,X,incx,beta,Y,incy) \
{ \
    BLAS_INT M = m, N = n, LDA = lda, INCX = incx, INCY = incy ; \
    double *ALPHA = alpha, *BETA = beta; \
    if (CHECK_BLAS_INT && !(EQ (M,m) && EQ (N,n) && EQ (LDA,lda) && \
        EQ (INCX,incx) && EQ (INCY,incy))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_DGEMV (BLAS_ORDER, BLAS_GET_TRANSPOSE(trans), M, N, ALPHA[0], A, LDA, X, INCX, BETA[0], Y, INCY) ; \
    } \
}

//void BLAS_ZGEMV (enum X_CBLAS_ORDER order, enum X_CBLAS_TRANSPOSE trans, BLAS_INT m, BLAS_INT n, double *alpha,
//	double *A, BLAS_INT lda, double *X, BLAS_INT incx, double *beta,
//	double *Y, BLAS_INT incy) ;
//void cblas_zgemv(const enum X_CBLAS_ORDER __Order,
//        const enum X_CBLAS_TRANSPOSE __TransA, const BLAS_INT __M, const BLAS_INT __N,
//        const void *__alpha, const void *__A, const BLAS_INT __lda, const void *__X,
//        const BLAS_INT __incX, const void *__beta, void *__Y, const BLAS_INT __incY) ;

#define BLAS_zgemv(trans,m,n,alpha,A,lda,X,incx,beta,Y,incy) \
{ \
    BLAS_INT M = m, N = n, LDA = lda, INCX = incx, INCY = incy ; \
    if (CHECK_BLAS_INT && !(EQ (M,m) && EQ (N,n) && EQ (LDA,lda) && \
        EQ (INCX,incx) && EQ (INCY,incy))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_ZGEMV (BLAS_ORDER, BLAS_GET_TRANSPOSE(trans), M, N, alpha, A, LDA, X, INCX, beta, Y, INCY) ; \
    } \
}

//void BLAS_DTRSV (enum X_CBLAS_ORDER order, enum X_CBLAS_UPLO uplo, enum X_CBLAS_TRANSPOSE trans, enum X_CBLAS_DIAG diag, BLAS_INT n, double *A,
//	BLAS_INT lda, double *X, BLAS_INT incx) ;
//void cblas_dtrsv(const enum X_CBLAS_ORDER __Order, const enum X_CBLAS_UPLO __Uplo,
//        const enum X_CBLAS_TRANSPOSE __TransA, const enum X_CBLAS_DIAG __Diag,
//        const BLAS_INT __N, const double *__A, const BLAS_INT __lda, double *__X,
//        const BLAS_INT __incX);
        
#define BLAS_dtrsv(uplo,trans,diag,n,A,lda,X,incx) \
{ \
    BLAS_INT N = n, LDA = lda, INCX = incx ; \
    if (CHECK_BLAS_INT && !(EQ (N,n) && EQ (LDA,lda) && EQ (INCX,incx))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_DTRSV (BLAS_ORDER, BLAS_GET_UPLO(uplo), BLAS_GET_TRANSPOSE(trans), BLAS_GET_DIAG(diag), N, A, LDA, X, INCX) ; \
    } \
}

//void BLAS_ZTRSV (enum X_CBLAS_ORDER order, enum X_CBLAS_UPLO uplo, enum X_CBLAS_TRANSPOSE trans, enum X_CBLAS_DIAG diag, BLAS_INT n, double *A,
//	BLAS_INT lda, double *X, BLAS_INT incx) ;
//void cblas_ztrsv(const enum X_CBLAS_ORDER __Order, const enum X_CBLAS_UPLO __Uplo,
//        const enum X_CBLAS_TRANSPOSE __TransA, const enum X_CBLAS_DIAG __Diag,
//        const BLAS_INT __N, const void *__A, const BLAS_INT __lda, void *__X,
//        const BLAS_INT __incX);
#define BLAS_ztrsv(uplo,trans,diag,n,A,lda,X,incx) \
{ \
    BLAS_INT N = n, LDA = lda, INCX = incx ; \
    if (CHECK_BLAS_INT && !(EQ (N,n) && EQ (LDA,lda) && EQ (INCX,incx))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_ZTRSV (BLAS_ORDER, BLAS_GET_UPLO(uplo), BLAS_GET_TRANSPOSE(trans), BLAS_GET_DIAG(diag), N, A, LDA, X, INCX) ; \
    } \
}

//void BLAS_DTRSM (enum X_CBLAS_ORDER order, enum X_CBLAS_SIDE side, enum X_CBLAS_UPLO uplo, enum X_CBLAS_TRANSPOSE transa, enum X_CBLAS_DIAG diag, BLAS_INT m,
//	BLAS_INT n, double *alpha, double *A, BLAS_INT lda, double *B,
//	BLAS_INT ldb) ;
//void cblas_dtrsm(const enum X_CBLAS_ORDER __Order, const enum X_CBLAS_SIDE __Side,
//        const enum X_CBLAS_UPLO __Uplo, const enum X_CBLAS_TRANSPOSE __TransA,
//        const enum X_CBLAS_DIAG __Diag, const BLAS_INT __M, const BLAS_INT __N,
//        const double __alpha, const double *__A, const BLAS_INT __lda, double *__B,
//        const BLAS_INT __ldb);
#define BLAS_dtrsm(side,uplo,transa,diag,m,n,alpha,A,lda,B,ldb) \
{ \
    BLAS_INT M = m, N = n, LDA = lda, LDB = ldb ; \
    double* ALPHA = alpha; \
    if (CHECK_BLAS_INT && !(EQ (M,m) && EQ (N,n) && EQ (LDA,lda) && \
        EQ (LDB,ldb))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_DTRSM (BLAS_ORDER, BLAS_GET_SIDE(side), BLAS_GET_UPLO(uplo), BLAS_GET_TRANSPOSE(transa), BLAS_GET_DIAG(diag), M, N, ALPHA[0], A, LDA, B, LDB);\
    } \
}

//void BLAS_ZTRSM (enum X_CBLAS_ORDER order, enum X_CBLAS_SIDE side, enum X_CBLAS_UPLO uplo, enum X_CBLAS_TRANSPOSE transa, enum X_CBLAS_DIAG diag, BLAS_INT m,
//	BLAS_INT n, double *alpha, double *A, BLAS_INT lda, double *B,
//	BLAS_INT ldb) ;
//void cblas_ztrsm(const enum X_CBLAS_ORDER __Order, const enum X_CBLAS_SIDE __Side,
//        const enum X_CBLAS_UPLO __Uplo, const enum X_CBLAS_TRANSPOSE __TransA,
//        const enum X_CBLAS_DIAG __Diag, const BLAS_INT __M, const BLAS_INT __N,
//        const void *__alpha, const void *__A, const BLAS_INT __lda, void *__B,
//        const BLAS_INT __ldb);
#define BLAS_ztrsm(side,uplo,transa,diag,m,n,alpha,A,lda,B,ldb) \
{ \
    BLAS_INT M = m, N = n, LDA = lda, LDB = ldb ; \
    if (CHECK_BLAS_INT && !(EQ (M,m) && EQ (N,n) && EQ (LDA,lda) && \
        EQ (LDB,ldb))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_ZTRSM (BLAS_ORDER, BLAS_GET_SIDE(side), BLAS_GET_UPLO(uplo), BLAS_GET_TRANSPOSE(transa), BLAS_GET_DIAG(diag), M, N, alpha, A, LDA, B, LDB);\
    } \
}

//void BLAS_DGEMM (enum X_CBLAS_ORDER order, enum X_CBLAS_TRANSPOSE transa, enum X_CBLAS_TRANSPOSE transb, BLAS_INT m, BLAS_INT n,
//	BLAS_INT k, double *alpha, double *A, BLAS_INT lda, double *B,
//	BLAS_INT ldb, double *beta, double *C, BLAS_INT ldc) ;

//void cblas_dgemm(const enum X_CBLAS_ORDER __Order,
//        const enum X_CBLAS_TRANSPOSE __TransA,
//        const enum X_CBLAS_TRANSPOSE __TransB, const BLAS_INT __M, const BLAS_INT __N,
//        const BLAS_INT __K, const double __alpha, const double *__A,
//        const BLAS_INT __lda, const double *__B, const BLAS_INT __ldb,
//        const double __beta, double *__C, const BLAS_INT __ldc) ;
#define BLAS_dgemm(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) \
{ \
    BLAS_INT M = m, N = n, K = k, LDA = lda, LDB = ldb, LDC = ldc ; \
    double* ALPHA = alpha, *BETA = beta; \
    if (CHECK_BLAS_INT && !(EQ (M,m) && EQ (N,n) && EQ (K,k) && \
        EQ (LDA,lda) && EQ (LDB,ldb) && EQ (LDC,ldc))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_DGEMM (BLAS_ORDER, BLAS_GET_TRANSPOSE(transa), BLAS_GET_TRANSPOSE(transb), M, N, K, ALPHA[0], A, LDA, B, LDB, BETA[0], \
	    C, LDC) ; \
    } \
}

//void BLAS_ZGEMM (enum X_CBLAS_ORDER order, enum X_CBLAS_TRANSPOSE transa, enum X_CBLAS_TRANSPOSE transb, BLAS_INT m, BLAS_INT n,
//	BLAS_INT k, double *alpha, double *A, BLAS_INT lda, double *B,
//	BLAS_INT ldb, double *beta, double *C, BLAS_INT ldc) ;
//void cblas_zgemm(const enum X_CBLAS_ORDER __Order,
//        const enum X_CBLAS_TRANSPOSE __TransA,
//        const enum X_CBLAS_TRANSPOSE __TransB, const BLAS_INT __M, const BLAS_INT __N,
//        const BLAS_INT __K, const void *__alpha, const void *__A, const BLAS_INT __lda,
//        const void *__B, const BLAS_INT __ldb, const void *__beta, void *__C,
//        const BLAS_INT __ldc);
#define BLAS_zgemm(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) \
{ \
    BLAS_INT M = m, N = n, K = k, LDA = lda, LDB = ldb, LDC = ldc ; \
    if (CHECK_BLAS_INT && !(EQ (M,m) && EQ (N,n) && EQ (K,k) && \
        EQ (LDA,lda) && EQ (LDB,ldb) && EQ (LDC,ldc))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_ZGEMM (BLAS_ORDER, BLAS_GET_TRANSPOSE(transa), BLAS_GET_TRANSPOSE(transb), M, N, K, alpha, A, LDA, B, LDB, beta, \
	    C, LDC) ; \
    } \
}

//void BLAS_DSYRK (enum X_CBLAS_ORDER order, enum X_CBLAS_UPLO uplo, enum X_CBLAS_TRANSPOSE trans, BLAS_INT n, BLAS_INT k,
//	double alpha, double *A, BLAS_INT lda, double beta, double *C,
//	BLAS_INT ldc) ;
//void cblas_dsyrk(const enum X_CBLAS_ORDER __Order, const enum X_CBLAS_UPLO __Uplo,
//        const enum X_CBLAS_TRANSPOSE __Trans, const BLAS_INT __N, const BLAS_INT __K,
//        const double __alpha, const double *__A, const BLAS_INT __lda,
//        const double __beta, double *__C, const BLAS_INT __ldc) ;
#define BLAS_dsyrk(uplo,trans,n,k,alpha,A,lda,beta,C,ldc) \
{ \
    BLAS_INT N = n, K = k, LDA = lda, LDC = ldc ; \
    if (CHECK_BLAS_INT && !(EQ (N,n) && EQ (K,k) && EQ (LDA,lda) && \
        EQ (LDC,ldc))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_DSYRK (BLAS_ORDER, BLAS_GET_UPLO(uplo), BLAS_GET_TRANSPOSE(trans), N, K, alpha[0], A, LDA, beta[0], C, LDC) ; \
    } \
} \

//void BLAS_ZHERK (enum X_CBLAS_ORDER order, enum X_CBLAS_UPLO uplo, enum X_CBLAS_TRANSPOSE trans, BLAS_INT n, BLAS_INT k,
//	double alpha, double *A, BLAS_INT lda, double beta, double *C,
//	BLAS_INT ldc) ;
//void cblas_zherk(const enum X_CBLAS_ORDER __Order, const enum X_CBLAS_UPLO __Uplo,
//        const enum X_CBLAS_TRANSPOSE __Trans, const BLAS_INT __N, const BLAS_INT __K,
//        const double __alpha, const void *__A, const BLAS_INT __lda,
//        const double __beta, void *__C, const BLAS_INT __ldc) ;
#define BLAS_zherk(uplo,trans,n,k,alpha,A,lda,beta,C,ldc) \
{ \
    BLAS_INT N = n, K = k, LDA = lda, LDC = ldc ; \
    if (CHECK_BLAS_INT && !(EQ (N,n) && EQ (K,k) && EQ (LDA,lda) && \
        EQ (LDC,ldc))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_ZHERK (BLAS_ORDER, BLAS_GET_UPLO(uplo), BLAS_GET_TRANSPOSE(trans), N, K, alpha[0], A, LDA, beta[0], C, LDC) ; \
    } \
} \

//void LAPACK_DPOTRF (char *uplo, BLAS_INT *n, double *A, BLAS_INT *lda,
//	BLAS_INT *info) ;

#define LAPACK_dpotrf(uplo,n,A,lda,info) \
{ \
    BLAS_INT N = n, LDA = lda, INFO = 1 ; \
    if (CHECK_BLAS_INT && !(EQ (N,n) && EQ (LDA,lda))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	LAPACK_DPOTRF (uplo, &N, A, &LDA, &INFO) ; \
    } \
    info = INFO ; \
}

//void LAPACK_ZPOTRF (char *uplo, BLAS_INT *n, double *A, BLAS_INT *lda,
//	BLAS_INT *info) ;

#define LAPACK_zpotrf(uplo,n,A,lda,info) \
{ \
    BLAS_INT N = n, LDA = lda, INFO = 1 ; \
    if (CHECK_BLAS_INT && !(EQ (N,n) && EQ (LDA,lda))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	LAPACK_ZPOTRF (uplo, &N, A, &LDA, &INFO) ; \
    } \
    info = INFO ; \
}

/* ========================================================================== */

//void BLAS_DSCAL (BLAS_INT n, double *alpha, double *Y, BLAS_INT incy) ;
//void cblas_dscal(const BLAS_INT __N, const double __alpha, double *__X,
//        const BLAS_INT __incX);
#define BLAS_dscal(n,alpha,Y,incy) \
{ \
    BLAS_INT N = n, INCY = incy ; \
    if (CHECK_BLAS_INT && !(EQ (N,n) && EQ (INCY,incy))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_DSCAL (N, alpha, Y, INCY) ; \
    } \
}

//void BLAS_ZSCAL (BLAS_INT n, double *alpha, double *Y, BLAS_INT incy) ;
//void cblas_zscal(const BLAS_INT __N, const void *__alpha, void *__X,
//        const BLAS_INT __incX);
#define BLAS_zscal(n,alpha,Y,incy) \
{ \
    BLAS_INT N = n, INCY = incy ; \
    if (CHECK_BLAS_INT && !(EQ (N,n) && EQ (INCY,incy))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_ZSCAL (N, alpha, Y, INCY) ; \
    } \
}

//void BLAS_DGER (enum X_CBLAS_ORDER order, BLAS_INT m, BLAS_INT n, double *alpha,
//	double *X, BLAS_INT incx, double *Y, BLAS_INT incy,
//	double *A, BLAS_INT lda) ;
//void cblas_dger(const enum X_CBLAS_ORDER __Order, const BLAS_INT __M, const BLAS_INT __N,
//        const double __alpha, const double *__X, const BLAS_INT __incX,
//        const double *__Y, const BLAS_INT __incY, double *__A, const BLAS_INT __lda) ;
#define BLAS_dger(m,n,alpha,X,incx,Y,incy,A,lda) \
{ \
    BLAS_INT M = m, N = n, LDA = lda, INCX = incx, INCY = incy ; \
    double *ALPHA = alpha; \
    if (CHECK_BLAS_INT && !(EQ (M,m) && EQ (N,n) && EQ (LDA,lda) && \
          EQ (INCX,incx) && EQ (INCY,incy))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_DGER (BLAS_ORDER, M, N, ALPHA[0], X, INCX, Y, INCY, A, LDA) ; \
    } \
}

//void BLAS_ZGER (enum X_CBLAS_ORDER order, BLAS_INT m, BLAS_INT n, double *alpha,
//	double *X, BLAS_INT incx, double *Y, BLAS_INT incy,
//	double *A, BLAS_INT lda) ;
//void cblas_zgeru(const enum X_CBLAS_ORDER __Order, const BLAS_INT __M, const BLAS_INT __N,
//        const void *__alpha, const void *__X, const BLAS_INT __incX,
//        const void *__Y, const BLAS_INT __incY, void *__A, const BLAS_INT __lda) ;
        
#define BLAS_zgeru(m,n,alpha,X,incx,Y,incy,A,lda) \
{ \
    BLAS_INT M = m, N = n, LDA = lda, INCX = incx, INCY = incy ; \
    if (CHECK_BLAS_INT && !(EQ (M,m) && EQ (N,n) && EQ (LDA,lda) && \
          EQ (INCX,incx) && EQ (INCY,incy))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_ZGER (BLAS_ORDER, M, N, alpha, X, INCX, Y, INCY, A, LDA) ; \
    } \
}
#else // ------------------------------------------------------------------------------------

/* ========================================================================== */
/* === BLAS and LAPACK prototypes and macros ================================ */
/* ========================================================================== */

void BLAS_DGEMV (char *trans, BLAS_INT *m, BLAS_INT *n, double *alpha,
	double *A, BLAS_INT *lda, double *X, BLAS_INT *incx, double *beta,
	double *Y, BLAS_INT *incy) ;

#define BLAS_dgemv(trans,m,n,alpha,A,lda,X,incx,beta,Y,incy) \
{ \
    BLAS_INT M = m, N = n, LDA = lda, INCX = incx, INCY = incy ; \
    if (CHECK_BLAS_INT && !(EQ (M,m) && EQ (N,n) && EQ (LDA,lda) && \
        EQ (INCX,incx) && EQ (INCY,incy))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_DGEMV (trans, &M, &N, alpha, A, &LDA, X, &INCX, beta, Y, &INCY) ; \
    } \
}

void BLAS_ZGEMV (char *trans, BLAS_INT *m, BLAS_INT *n, double *alpha,
	double *A, BLAS_INT *lda, double *X, BLAS_INT *incx, double *beta,
	double *Y, BLAS_INT *incy) ;

#define BLAS_zgemv(trans,m,n,alpha,A,lda,X,incx,beta,Y,incy) \
{ \
    BLAS_INT M = m, N = n, LDA = lda, INCX = incx, INCY = incy ; \
    if (CHECK_BLAS_INT && !(EQ (M,m) && EQ (N,n) && EQ (LDA,lda) && \
        EQ (INCX,incx) && EQ (INCY,incy))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_ZGEMV (trans, &M, &N, alpha, A, &LDA, X, &INCX, beta, Y, &INCY) ; \
    } \
}

void BLAS_DTRSV (char *uplo, char *trans, char *diag, BLAS_INT *n, double *A,
	BLAS_INT *lda, double *X, BLAS_INT *incx) ;

#define BLAS_dtrsv(uplo,trans,diag,n,A,lda,X,incx) \
{ \
    BLAS_INT N = n, LDA = lda, INCX = incx ; \
    if (CHECK_BLAS_INT && !(EQ (N,n) && EQ (LDA,lda) && EQ (INCX,incx))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_DTRSV (uplo, trans, diag, &N, A, &LDA, X, &INCX) ; \
    } \
}

void BLAS_ZTRSV (char *uplo, char *trans, char *diag, BLAS_INT *n, double *A,
	BLAS_INT *lda, double *X, BLAS_INT *incx) ;

#define BLAS_ztrsv(uplo,trans,diag,n,A,lda,X,incx) \
{ \
    BLAS_INT N = n, LDA = lda, INCX = incx ; \
    if (CHECK_BLAS_INT && !(EQ (N,n) && EQ (LDA,lda) && EQ (INCX,incx))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_ZTRSV (uplo, trans, diag, &N, A, &LDA, X, &INCX) ; \
    } \
}

void BLAS_DTRSM (char *side, char *uplo, char *transa, char *diag, BLAS_INT *m,
	BLAS_INT *n, double *alpha, double *A, BLAS_INT *lda, double *B,
	BLAS_INT *ldb) ;

#define BLAS_dtrsm(side,uplo,transa,diag,m,n,alpha,A,lda,B,ldb) \
{ \
    BLAS_INT M = m, N = n, LDA = lda, LDB = ldb ; \
    if (CHECK_BLAS_INT && !(EQ (M,m) && EQ (N,n) && EQ (LDA,lda) && \
        EQ (LDB,ldb))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_DTRSM (side, uplo, transa, diag, &M, &N, alpha, A, &LDA, B, &LDB);\
    } \
}

void BLAS_ZTRSM (char *side, char *uplo, char *transa, char *diag, BLAS_INT *m,
	BLAS_INT *n, double *alpha, double *A, BLAS_INT *lda, double *B,
	BLAS_INT *ldb) ;

#define BLAS_ztrsm(side,uplo,transa,diag,m,n,alpha,A,lda,B,ldb) \
{ \
    BLAS_INT M = m, N = n, LDA = lda, LDB = ldb ; \
    if (CHECK_BLAS_INT && !(EQ (M,m) && EQ (N,n) && EQ (LDA,lda) && \
        EQ (LDB,ldb))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_ZTRSM (side, uplo, transa, diag, &M, &N, alpha, A, &LDA, B, &LDB);\
    } \
}

void BLAS_DGEMM (char *transa, char *transb, BLAS_INT *m, BLAS_INT *n,
	BLAS_INT *k, double *alpha, double *A, BLAS_INT *lda, double *B,
	BLAS_INT *ldb, double *beta, double *C, BLAS_INT *ldc) ;

#define BLAS_dgemm(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) \
{ \
    BLAS_INT M = m, N = n, K = k, LDA = lda, LDB = ldb, LDC = ldc ; \
    if (CHECK_BLAS_INT && !(EQ (M,m) && EQ (N,n) && EQ (K,k) && \
        EQ (LDA,lda) && EQ (LDB,ldb) && EQ (LDC,ldc))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_DGEMM (transa, transb, &M, &N, &K, alpha, A, &LDA, B, &LDB, beta, \
	    C, &LDC) ; \
    } \
}

void BLAS_ZGEMM (char *transa, char *transb, BLAS_INT *m, BLAS_INT *n,
	BLAS_INT *k, double *alpha, double *A, BLAS_INT *lda, double *B,
	BLAS_INT *ldb, double *beta, double *C, BLAS_INT *ldc) ;

#define BLAS_zgemm(transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) \
{ \
    BLAS_INT M = m, N = n, K = k, LDA = lda, LDB = ldb, LDC = ldc ; \
    if (CHECK_BLAS_INT && !(EQ (M,m) && EQ (N,n) && EQ (K,k) && \
        EQ (LDA,lda) && EQ (LDB,ldb) && EQ (LDC,ldc))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_ZGEMM (transa, transb, &M, &N, &K, alpha, A, &LDA, B, &LDB, beta, \
	    C, &LDC) ; \
    } \
}

void BLAS_DSYRK (char *uplo, char *trans, BLAS_INT *n, BLAS_INT *k,
	double *alpha, double *A, BLAS_INT *lda, double *beta, double *C,
	BLAS_INT *ldc) ;

#define BLAS_dsyrk(uplo,trans,n,k,alpha,A,lda,beta,C,ldc) \
{ \
    BLAS_INT N = n, K = k, LDA = lda, LDC = ldc ; \
    if (CHECK_BLAS_INT && !(EQ (N,n) && EQ (K,k) && EQ (LDA,lda) && \
        EQ (LDC,ldc))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_DSYRK (uplo, trans, &N, &K, alpha, A, &LDA, beta, C, &LDC) ; \
    } \
} \

void BLAS_ZHERK (char *uplo, char *trans, BLAS_INT *n, BLAS_INT *k,
	double *alpha, double *A, BLAS_INT *lda, double *beta, double *C,
	BLAS_INT *ldc) ;

#define BLAS_zherk(uplo,trans,n,k,alpha,A,lda,beta,C,ldc) \
{ \
    BLAS_INT N = n, K = k, LDA = lda, LDC = ldc ; \
    if (CHECK_BLAS_INT && !(EQ (N,n) && EQ (K,k) && EQ (LDA,lda) && \
        EQ (LDC,ldc))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_ZHERK (uplo, trans, &N, &K, alpha, A, &LDA, beta, C, &LDC) ; \
    } \
} \

void LAPACK_DPOTRF (char *uplo, BLAS_INT *n, double *A, BLAS_INT *lda,
	BLAS_INT *info) ;

#define LAPACK_dpotrf(uplo,n,A,lda,info) \
{ \
    BLAS_INT N = n, LDA = lda, INFO = 1 ; \
    if (CHECK_BLAS_INT && !(EQ (N,n) && EQ (LDA,lda))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	LAPACK_DPOTRF (uplo, &N, A, &LDA, &INFO) ; \
    } \
    info = INFO ; \
}

void LAPACK_ZPOTRF (char *uplo, BLAS_INT *n, double *A, BLAS_INT *lda,
	BLAS_INT *info) ;

#define LAPACK_zpotrf(uplo,n,A,lda,info) \
{ \
    BLAS_INT N = n, LDA = lda, INFO = 1 ; \
    if (CHECK_BLAS_INT && !(EQ (N,n) && EQ (LDA,lda))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	LAPACK_ZPOTRF (uplo, &N, A, &LDA, &INFO) ; \
    } \
    info = INFO ; \
}

/* ========================================================================== */

void BLAS_DSCAL (BLAS_INT *n, double *alpha, double *Y, BLAS_INT *incy) ;

#define BLAS_dscal(n,alpha,Y,incy) \
{ \
    BLAS_INT N = n, INCY = incy ; \
    if (CHECK_BLAS_INT && !(EQ (N,n) && EQ (INCY,incy))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_DSCAL (&N, alpha, Y, &INCY) ; \
    } \
}

void BLAS_ZSCAL (BLAS_INT *n, double *alpha, double *Y, BLAS_INT *incy) ;

#define BLAS_zscal(n,alpha,Y,incy) \
{ \
    BLAS_INT N = n, INCY = incy ; \
    if (CHECK_BLAS_INT && !(EQ (N,n) && EQ (INCY,incy))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_ZSCAL (&N, alpha, Y, &INCY) ; \
    } \
}

void BLAS_DGER (BLAS_INT *m, BLAS_INT *n, double *alpha,
	double *X, BLAS_INT *incx, double *Y, BLAS_INT *incy,
	double *A, BLAS_INT *lda) ;

#define BLAS_dger(m,n,alpha,X,incx,Y,incy,A,lda) \
{ \
    BLAS_INT M = m, N = n, LDA = lda, INCX = incx, INCY = incy ; \
    if (CHECK_BLAS_INT && !(EQ (M,m) && EQ (N,n) && EQ (LDA,lda) && \
          EQ (INCX,incx) && EQ (INCY,incy))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_DGER (&M, &N, alpha, X, &INCX, Y, &INCY, A, &LDA) ; \
    } \
}

void BLAS_ZGER (BLAS_INT *m, BLAS_INT *n, double *alpha,
	double *X, BLAS_INT *incx, double *Y, BLAS_INT *incy,
	double *A, BLAS_INT *lda) ;

#define BLAS_zgeru(m,n,alpha,X,incx,Y,incy,A,lda) \
{ \
    BLAS_INT M = m, N = n, LDA = lda, INCX = incx, INCY = incy ; \
    if (CHECK_BLAS_INT && !(EQ (M,m) && EQ (N,n) && EQ (LDA,lda) && \
          EQ (INCX,incx) && EQ (INCY,incy))) \
    { \
	BLAS_OK = FALSE ; \
    } \
    if (!CHECK_BLAS_INT || BLAS_OK) \
    { \
	BLAS_ZGER (&M, &N, alpha, X, &INCX, Y, &INCY, A, &LDA) ; \
    } \
}
#endif
#endif
