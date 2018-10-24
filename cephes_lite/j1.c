/*							j1.c
 *
 *	Bessel function of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, j1();
 *
 * y = j1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order one of the argument.
 *
 * The domain is divided into the intervals [0, 8] and
 * (8, infinity). In the first interval a 24 term Chebyshev
 * expansion is used. In the second, the asymptotic
 * trigonometric representation is employed using two
 * rational functions of degree 5/5.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain      # trials      peak         rms
 *    DEC       0, 30       10000       4.0e-17     1.1e-17
 *    IEEE      0, 30       30000       2.6e-16     1.1e-16
 *
 *
 */
/*							y1.c
 *
 *	Bessel function of second kind of order one
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, y1();
 *
 * y = y1( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of the second kind of order one
 * of the argument.
 *
 * The domain is divided into the intervals [0, 8] and
 * (8, infinity). In the first interval a 25 term Chebyshev
 * expansion is used, and a call to j1() is required.
 * In the second, the asymptotic trigonometric representation
 * is employed using two rational functions of degree 5/5.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain      # trials      peak         rms
 *    DEC       0, 30       10000       8.6e-17     1.3e-17
 *    IEEE      0, 30       30000       1.0e-15     1.3e-16
 *
 * (error criterion relative when |y1| > 1).
 *
 */


/*
 * Cephes Math Library Release 2.8:  June, 2000
 * Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
 */

/*
 * #define PIO4 .78539816339744830962
 * #define THPIO4 2.35619449019234492885
 * #define SQ2OPI .79788456080286535588
 */

#include "mconf.h"

static unsigned short RP[16] = {
    0x4fc5,0x4cda,0xd23c,0xc1ca,
    0x8c6b,0x0d43,0x52ba,0x425a,
    0xf175,0xe6cc,0x8a92,0xc2d0,
    0xc482,0xa697,0x2b42,0x432a,
};
static unsigned short RQ[32] = {
    /*0x0000,0x0000,0x0000,0x3ff0,*/
    0x86e7,0x1b70,0x66b1,0x4083,
    0x01b2,0x0dd7,0x5eda,0x410f,
    0xa1b1,0xdc92,0xe954,0x4193,
    0xeac1,0x7bef,0xa13f,0x4214,
    0xffa8,0x8076,0x46fb,0x4291,
    0xf45f,0x3ecc,0x4b0a,0x4306,
    0x3f81,0xf465,0xe0bf,0x4373,
    0x2939,0x7670,0x7795,0x43d2,
};

static unsigned short PP[28] = {
    0x651b,0x4c6c,0xf92c,0x3f48,
    0xc4b6,0xa3fe,0xb948,0x3fb2,
    0x96d6,0xc215,0x08fe,0x3ff2,
    0x3a6a,0xf8b1,0x72c4,0x4014,
    0x2f16,0x8b5d,0xd91c,0x4020,
    0x81a2,0x142f,0xdbaa,0x4014,
    0x0000,0x0000,0x0000,0x3ff0,
};
static unsigned short PQ[28] = {
    0x3d69,0x1344,0xb89b,0x3f42,
    0xaa83,0x5948,0x9fdd,0x3fb1,
    0xeed6,0xb850,0xaea9,0x3ff1,
    0x51a1,0xf7d2,0x4ba2,0x4014,
    0xfd65,0xdda2,0xccb9,0x4020,
    0xb4d9,0x4762,0xd6dd,0x4014,
    0x0000,0x0000,0x0000,0x3ff0,
};

static unsigned short QP[32] = {
    0xba40,0x6b70,0x27fa,0x3faa,
    0x8fd6,0xc66d,0xedb5,0x4013,
    0x1c67,0x9acf,0xf4b9,0x4052,
    0x180d,0x47aa,0xec79,0x4076,
    0x6e50,0xb66f,0x36d9,0x4086,
    0x02d0,0xb9e8,0xabea,0x4082,
    0xbb0b,0x4c54,0x760a,0x406a,
    0x9eb5,0x4d15,0x34ff,0x4039,
};
static unsigned short QQ[28] = {
    /*0x0000,0x0000,0x0000,0x3ff0,*/
    0x5077,0x6089,0x8f30,0x4052,
    0x5f6f,0xa20e,0x81cb,0x4090,
    0xfe81,0x1bfd,0x7a69,0x40b3,
    0xd118,0xd280,0xad28,0x40c2,
    0x3d14,0xa697,0x3d0a,0x40bf,
    0x1781,0xb4bd,0x1462,0x40a6,
    0x5997,0x6ae7,0x017f,0x4075,
};

static unsigned short YP[24] = {
    0xb6c5,0x62f9,0xd2be,0x41d2,
    0x6521,0x5883,0xd72d,0xc262,
    0x0fef,0xb091,0x0954,0x42da,
    0xb083,0x37a1,0xe01a,0xc33c,
    0x66b1,0x0b73,0x79ad,0x4386,
    0x2b17,0x5dde,0x9e41,0xc3a5,
};
static unsigned short YQ[32] = {
    /*0x0000,0x0000,0x0000,0x3ff0,*/
    0x7ac2,0xa93f,0x9269,0x4082,
    0xef7f,0xbe58,0xc160,0x410c,
    0xacee,0xa9c8,0x84ef,0x4191,
    0x7b83,0x906b,0x78c3,0x4211,
    0x9316,0xfda9,0x3f5e,0x428c,
    0x1e4e,0xd71d,0xa326,0x4301,
    0xa488,0xc547,0x83e3,0x436e,
    0x747f,0x90f6,0x90f1,0x43cb,
};

static unsigned short DZ1[] = {0x822c,0x4189,0x5d2b,0x402d};
static unsigned short DZ2[] = {0xa432,0x6072,0x9bf6,0x4048};
#define Z1 (*(double *)DZ1)
#define Z2 (*(double *)DZ2)


extern double polevl ( double, void *, int );
extern double p1evl ( double, void *, int );
extern double log ( double );
extern double sin ( double );
extern double cos ( double );
extern double sqrt ( double );
double j1 ( double );

extern double TWOOPI, THPIO4, SQ2OPI;

double j1(double x)
{
    double w, z, p, q, xn;
    
    w = x;
    if( x < 0 )
        w = -x;
    
    if( w <= 5.0 )
    {
        z = x * x;
        w = polevl( z, RP, 3 ) / p1evl( z, RQ, 8 );
        w = w * x * (z - Z1) * (z - Z2);
        return( w );
    }
    
    w = 5.0/x;
    z = w * w;
    p = polevl( z, PP, 6)/polevl( z, PQ, 6 );
    q = polevl( z, QP, 7)/p1evl( z, QQ, 7 );
    xn = x - THPIO4;
    p = p * cos(xn) - w * q * sin(xn);
    return( p * SQ2OPI / sqrt(x) );
}


extern double MAXNUM;

double y1(double x)
{
    double w, z, p, q, xn;
    
    if( x <= 5.0 )
    {
        if( x <= 0.0 )
        {
            mtherr( "y1", DOMAIN );
            return( -MAXNUM );
        }
        z = x * x;
        w = x * (polevl( z, YP, 5 ) / p1evl( z, YQ, 8 ));
        w += TWOOPI * ( j1(x) * log(x)  -  1.0/x );
        return( w );
    }
    
    w = 5.0/x;
    z = w * w;
    p = polevl( z, PP, 6)/polevl( z, PQ, 6 );
    q = polevl( z, QP, 7)/p1evl( z, QQ, 7 );
    xn = x - THPIO4;
    p = p * sin(xn) + w * q * cos(xn);
    return( p * SQ2OPI / sqrt(x) );
}
