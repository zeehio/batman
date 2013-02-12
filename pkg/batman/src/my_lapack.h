// written by Dr. Jie Hao, Dr William Astle
#ifndef MY_LAPACK    
#define MY_LAPACK 

typedef int MY_LAPACK_INT;
typedef double MY_LAPACK_DBL;
typedef char MY_LAPACK_CHAR;


extern "C" {
    void dgeqp3_(MY_LAPACK_INT* m, MY_LAPACK_INT* n, MY_LAPACK_DBL* a, MY_LAPACK_INT* lda,MY_LAPACK_INT* jpvt, MY_LAPACK_DBL* tau, MY_LAPACK_DBL* work, MY_LAPACK_INT* lwork, MY_LAPACK_INT* info); 
    void dorgqr_(MY_LAPACK_INT* m, MY_LAPACK_INT* n, MY_LAPACK_INT* k, MY_LAPACK_DBL* a, MY_LAPACK_INT* lda, MY_LAPACK_DBL* tau, MY_LAPACK_DBL*  work, MY_LAPACK_INT* lwork, MY_LAPACK_INT* info);
}
#endif
