/*							j0.c
 *
 *	Bessel function of order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, j0();
 *
 * y = j0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of order zero of the argument.
 *
 * The domain is divided into the intervals [0, 5] and
 * (5, infinity). In the first interval the following rational
 * approximation is used:
 *
 *
 *        2         2
 * (w - r  ) (w - r  ) P (w) / Q (w)
 *       1         2    3       8
 *
 *            2
 * where w = x  and the two r's are zeros of the function.
 *
 * In the second interval, the Hankel asymptotic expansion
 * is employed with two rational functions of degree 6/6
 * and 7/7.
 *
 *
 *
 * ACCURACY:
 *
 *                      Absolute error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30       10000       4.4e-17     6.3e-18
 *    IEEE      0, 30       60000       4.2e-16     1.1e-16
 *
 */
/*							y0.c
 *
 *	Bessel function of the second kind, order zero
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, y0();
 *
 * y = y0( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns Bessel function of the second kind, of order
 * zero, of the argument.
 *
 * The domain is divided into the intervals [0, 5] and
 * (5, infinity). In the first interval a rational approximation
 * R(x) is employed to compute
 *   y0(x)  = R(x)  +   2 * log(x) * j0(x) / PI.
 * Thus a call to j0() is required.
 *
 * In the second interval, the Hankel asymptotic expansion
 * is employed with two rational functions of degree 6/6
 * and 7/7.
 *
 *
 *
 * ACCURACY:
 *
 *  Absolute error, when y0(x) < 1; else relative error:
 *
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30        9400       7.0e-17     7.9e-18
 *    IEEE      0, 30       30000       1.3e-15     1.6e-16
 *
 */

/*
 * Cephes Math Library Release 2.8:  June, 2000
 * Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
 */

/* Note: all coefficients satisfy the relative error criterion
 * except YP, YQ which are designed for absolute error. */

#include "mconf.h"


static unsigned short PP[28] = {
    0x6b27,0x983b,0x1d30,0x3f4a,
    0xd1cf,0xb35d,0x34b0,0x3fb5,
    0x0b98,0x4e68,0xd521,0x3ff3,
    0x0956,0xe97a,0xc9fb,0x4015,
    0x9888,0x6940,0x7e8c,0x4021,
    0x25a1,0xa594,0x3684,0x4015,
    0x0000,0x0000,0x0000,0x3ff0,
};
static unsigned short PQ[28] = {
    0x9737,0xce03,0x4a80,0x3f4e,
    0x54e3,0xab54,0xebc5,0x3fb5,
    0x069f,0xc9b3,0x0e72,0x3ff4,
    0x62bb,0xe681,0xe247,0x4015,
    0x21a1,0xea1b,0x8618,0x4021,
    0x3a19,0xed42,0x3965,0x4015,
    0x0000,0x0000,0x0000,0x3ff0,
};

static unsigned short QP[32] = {
    0x384a,0x38a5,0x4742,0xbf87,
    0x1174,0x3a32,0x853b,0xbff4,
    0x2c0c,0xf50e,0x8dcf,0xc033,
    0xe8c4,0x5a6d,0x4d2f,0xc057,
    0xe8ea,0x20ca,0x35cc,0xc066,
    0x392d,0xec17,0x627a,0xc062,
    0x18cd,0x55b2,0xb48c,0xc049,
    0xa1dd,0xd1b9,0x3358,0xc018,
};
static unsigned short QQ[28] = {
    /*0x0000,0x0000,0x0000,0x3ff0,*/
    0x25ac,0x413c,0x1457,0x4050,
    0x9c7f,0xb175,0xc370,0x408a,
    0x8cb5,0xbd74,0x54cd,0x40ae,
    0xd63e,0xbdef,0x4877,0x40bc,
    0x3b11,0x1d73,0x2aba,0x40b7,
    0x9e82,0xc731,0x1c2f,0x40a0,
    0x0a54,0x0628,0x402f,0x406e,
};

static unsigned short YP[32] = {
    0x898f,0xe896,0x7437,0x40ce,
    0x8896,0x32e4,0xf81f,0xc16b,
    0x4cdd,0xf028,0x3f78,0x41f4,
    0xbd2b,0xe1d6,0x957b,0xc26c,
    0xac2d,0x3cc3,0xea72,0x42d3,
    0xcc02,0xd1d8,0xa121,0xc328,
    0x4003,0x660b,0xa94b,0x4363,
    0x367b,0x5906,0x6d4b,0xc350,
};
static unsigned short YQ[28] = {
    /*0x0000,0x0000,0x0000,0x3ff0,*/
    0xfcb6,0x576d,0x4522,0x4090,
    0xbc0c,0xa907,0x1b76,0x4123,
    0xd101,0x5164,0x0763,0x41b0,
    0x64bc,0x2b86,0x1ddb,0x4234,
    0x828e,0xc57e,0x75fc,0x42b2,
    0x596d,0xdfeb,0x8910,0x4326,
    0xb5d0,0xbcf9,0xd25f,0x438b,
};

static unsigned short R1[] = {0x2bbb,0x8046,0x21fb,0x4017};
#define DR1 *(double *)R1
static unsigned short R2[] = {0xdd6f,0xa621,0x78a4,0x403e};
#define DR2 *(double *)R2



static unsigned short RP[16] = {
    0x8325,0xad1c,0xdc53,0xc1f1,
    0x990d,0xc772,0x7751,0x427c,
    0x00f7,0xe0d9,0x5614,0xc2ec,
    0x5fb4,0x69ff,0x3ef8,0x4341,
};
static unsigned short RQ[32] = {
    /*0x0000,0x0000,0x0000,0x3ff0,*/
    0xb78c,0xa696,0x3902,0x407f,
    0x1a67,0x36a2,0x36cb,0x4105,
    0x0634,0x2eac,0x1934,0x4187,
    0x4914,0x0944,0xd5b0,0x4204,
    0x2e46,0x7218,0xbeb3,0x427e,
    0x48e9,0x8c97,0xa6a2,0x42f1,
    0x2e9c,0x7e7b,0x4141,0x435c,
    0x62cc,0xc7b6,0xbe34,0x43b7,
};

extern double polevl ( double, void *, int );
extern double p1evl ( double, void *, int );
extern double log ( double );
extern double sin ( double );
extern double cos ( double );
extern double sqrt ( double );
double j0 ( double );

extern double TWOOPI, SQ2OPI, PIO4;

double j0(double x)

{
    double w, z, p, q, xn;
    
    if( x < 0 )
        x = -x;
    
    if( x <= 5.0 )
    {
        z = x * x;
        if( x < 1.0e-5 )
            return( 1.0 - z/4.0 );
        
        p = (z - DR1) * (z - DR2);
        p = p * polevl( z, RP, 3)/p1evl( z, RQ, 8 );
        return( p );
    }
    
    w = 5.0/x;
    q = 25.0/(x*x);
    p = polevl( q, PP, 6)/polevl( q, PQ, 6 );
    q = polevl( q, QP, 7)/p1evl( q, QQ, 7 );
    xn = x - PIO4;
    p = p * cos(xn) - w * q * sin(xn);
    return( p * SQ2OPI / sqrt(x) );
}

/*							y0() 2	*/
/* Bessel function of second kind, order zero	*/

/* Rational approximation coefficients YP[], YQ[] are used here.
 * The function computed is  y0(x)  -  2 * log(x) * j0(x) / PI,
 * whose value at x = 0 is  2 * ( log(0.5) + EUL ) / PI
 * = 0.073804295108687225.
 */

/*
 * #define PIO4 .78539816339744830962
 * #define SQ2OPI .79788456080286535588
 */
extern double MAXNUM;

double y0(double x)

{
    double w, z, p, q, xn;
    
    if( x <= 5.0 )
    {
        if( x <= 0.0 )
        {
            mtherr( "y0", DOMAIN );
            return( -MAXNUM );
        }
        z = x * x;
        w = polevl( z, YP, 7) / p1evl( z, YQ, 7 );
        w += TWOOPI * log(x) * j0(x);
        return( w );
    }
    
    w = 5.0/x;
    z = 25.0 / (x * x);
    p = polevl( z, PP, 6)/polevl( z, PQ, 6 );
    q = polevl( z, QP, 7)/p1evl( z, QQ, 7 );
    xn = x - PIO4;
    p = p * sin(xn) + w * q * cos(xn);
    return( p * SQ2OPI / sqrt(x) );
}
