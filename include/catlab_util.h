#ifndef __CATLAB_UTIL_H__
#define __CATLAB_UTIL_H__

//#include <mkl_cblas.h>
//#include <mkl_lapacke.h>
#include <math.h>
#include <cblas.h>
#include <lapacke.h>
#include <string.h>
#include <complex.h>

#ifndef PI
#define PI             3.14159265358979323846264338327950288419716939937510
#endif//PI

#define CAT_UNKNOWN     0
#define CAT_TRUE        1
#define CAT_FALSE       -1

#define CAT_S           0x00000000
#define CAT_D           0x00000001
#define CAT_C           0x00000002
#define CAT_Z           0x00000004
#define CAT_SIZEOF_S    sizeof(float)
#define CAT_SIZEOF_D    sizeof(double)
#define CAT_SIZEOF_C    sizeof(lapack_complex_float)
#define CAT_SIZEOF_Z    sizeof(lapack_complex_double)

typedef int cat_bool;
typedef unsigned int cat_flag;

typedef float cat_s;
typedef double cat_d;
typedef lapack_complex_float cat_c;
typedef lapack_complex_double cat_z;

//#define CAT_REAL(x) (x).real
//#define CAT_IMAG(x) (x).imag
#define CAT_REAL(x) creal(x)
#define CAT_IMAG(x) cimag(x) 
#define CAT_CONJ(x) conj(x)

#define CAT_SIZEOF(x) (                                         \
                ((x) == CAT_D) ? CAT_SIZEOF_D : (               \
                    ((x) == CAT_Z) ? CAT_SIZEOF_Z : (           \
                        ((x) == CAT_S) ? CAT_SIZEOF_S : (       \
                            ((x) == CAT_C) ? CAT_SIZEOF_C : 0   \
                        )                                       \
                    )                                           \
                )                                               \
            )

#define CAT_CHECKFLAG(flag,check) ((flag & check) != 0)

#define CAT_ERR(x) printf("Error:"#x"=%d\n",x);

#define CAT_MAX(x,y) ( ((x)>=(y))?(x):(y) )
#define CAT_MIN(x,y) ( ((x)<=(y))?(x):(y) )

#ifdef CAT_DEBUG

#define CAT_ASSERT(assertion) if (!(assertion)) printf("CAT:assert:abort:"#assertion"\n");
#define CAT_LOG(x) printf(#x"\n");

#else

#define CAT_ASSERT(assertion) 
#define CAT_LOG(x) 

#endif//CAT_DEBUG



#endif//__CATLAB_UTIL_H__