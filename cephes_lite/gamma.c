/*							gamma.c
 *
 *	Gamma function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, gamma();
 * extern int sgngam;
 *
 * y = gamma( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns gamma function of the argument.  The result is
 * correctly signed, and the sign (+1 or -1) is also
 * returned in a global (extern) variable named sgngam.
 * This variable is also filled in by the logarithmic gamma
 * function lgam().
 *
 * Arguments |x| <= 34 are reduced by recurrence and the function
 * approximated by a rational function of degree 6/7 in the
 * interval (2,3).  Large arguments are handled by Stirling's
 * formula. Large negative arguments are made positive using
 * a reflection formula.
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC      -34, 34      10000       1.3e-16     2.5e-17
 *    IEEE    -170,-33      20000       2.3e-15     3.3e-16
 *    IEEE     -33,  33     20000       9.4e-16     2.2e-16
 *    IEEE      33, 171.6   20000       2.3e-15     3.2e-16
 *
 * Error for arguments outside the test range will be larger
 * owing to error amplification by the exponential function.
 *
 */
/*							lgam()
 *
 *	Natural logarithm of gamma function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, lgam();
 * extern int sgngam;
 *
 * y = lgam( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the base e (2.718...) logarithm of the absolute
 * value of the gamma function of the argument.
 * The sign (+1 or -1) of the gamma function is returned in a
 * global (extern) variable named sgngam.
 *
 * For arguments greater than 13, the logarithm of the gamma
 * function is approximated by the logarithmic version of
 * Stirling's formula using a polynomial approximation of
 * degree 4. Arguments between -33 and +33 are reduced by
 * recurrence to the interval [2,3] of a rational approximation.
 * The cosecant reflection formula is employed for arguments
 * less than -33.
 *
 * Arguments greater than MAXLGM return MAXNUM and an error
 * message.  MAXLGM = 2.035093e36 for DEC
 * arithmetic or 2.556348e305 for IEEE arithmetic.
 *
 *
 *
 * ACCURACY:
 *
 *
 * arithmetic      domain        # trials     peak         rms
 *    DEC     0, 3                  7000     5.2e-17     1.3e-17
 *    DEC     2.718, 2.035e36       5000     3.9e-17     9.9e-18
 *    IEEE    0, 3                 28000     5.4e-16     1.1e-16
 *    IEEE    2.718, 2.556e305     40000     3.5e-16     8.3e-17
 * The error criterion was relative when the function magnitude
 * was greater than one but absolute when it was less than one.
 *
 * The following test used the relative error criterion, though
 * at certain points the relative error could be much higher than
 * indicated.
 *    IEEE    -200, -4             10000     4.8e-16     1.3e-16
 *
 */

/*							gamma.c	*/
/*	gamma function	*/

/*
 * Cephes Math Library Release 2.8:  June, 2000
 * Copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
 */


#include "mconf.h"



static unsigned short P[] = {
    0x2153,0x3998,0xfcb8,0x3f24,
    0xbfab,0xe686,0x84e3,0x3f53,
    0x14b0,0xe9db,0x57cd,0x3f85,
    0x23d3,0x18c4,0x63d9,0x3fa8,
    0x7d31,0xdcae,0x8da9,0x3fca,
    0xe312,0x3993,0xa137,0x3fdf,
    0x0000,0x0000,0x0000,0x3ff0
};
static unsigned short Q[] = {
    0xd3af,0x8400,0x487a,0xbef8,
    0x2573,0x2915,0xae8a,0x3f41,
    0xb44a,0xe750,0x40e4,0xbf72,
    0xb117,0x5b1b,0x31ed,0x3f88,
    0xde67,0xe33f,0x5779,0x3fa2,
    0x87c2,0x9d42,0x071a,0xbfce,
    0x3c51,0xc9cd,0x4944,0x3fb2,
    0x0000,0x0000,0x0000,0x3ff0
};
#define MAXGAM 171.624376956302725
static unsigned short LPI[4] = {
    0xa1bd,0x48e7,0x50d0,0x3ff2,
};
#define LOGPI *(double *)LPI



/* Stirling's formula for the gamma function */


static unsigned short STIR[20] = {
    0x7293,0x592d,0xcc72,0x3f49,
    0x1d7c,0x27e6,0x166b,0xbf2e,
    0x4fd7,0x07d4,0xf726,0xbf65,
    0xc5fd,0x1b98,0x71c7,0x3f6c,
    0x5986,0x5555,0x5555,0x3fb5,
};
#define MAXSTIR 143.01608
static unsigned short SQT[4] = {
    0x2706,0x1ff6,0x0d93,0x4004,
};
#define SQTPI *(double *)SQT



int sgngam = 0;
extern int sgngam;
extern double MAXLOG, MAXNUM, PI;

extern double pow ( double, double );
extern double log ( double );
extern double exp ( double );
extern double sin ( double );
extern double polevl ( double, void *, int );
extern double p1evl ( double, void *, int );
extern double floor ( double );
extern double fabs ( double );
extern int isnan ( double );
extern int isfinite ( double );
static double stirf ( double );
double lgam ( double );




/* Gamma function computed by Stirling's formula.
 * The polynomial STIR is valid for 33 <= x <= 172.
 */
static double stirf(double x)

{
    double y, w, v;
    
    w = 1.0/x;
    w = 1.0 + w * polevl( w, STIR, 4 );
    y = exp(x);
    if( x > MAXSTIR )
    { /* Avoid overflow in pow() */
        v = pow( x, 0.5 * x - 0.25 );
        y = v * (v / y);
    }
    else
    {
        y = pow( x, x - 0.5 ) / y;
    }
    y = SQTPI * y * w;
    return( y );
}



double gamma(double x)
{
    double p, q, z;
    int i;
    
    sgngam = 1;
    
    
    q = fabs(x);
    
    if( q > 33.0 )
    {
        if( x < 0.0 )
        {
            p = floor(q);
            if( p == q )
            {
                
                goto goverf;
                
            }
            i = p;
            if( (i & 1) == 0 )
                sgngam = -1;
            z = q - p;
            if( z > 0.5 )
            {
                p += 1.0;
                z = q - p;
            }
            z = q * sin( PI * z );
            if( z == 0.0 )
            {
                
                goverf:
                    mtherr( "gamma", OVERFLOW );
                    return( sgngam * MAXNUM);
                    
            }
            z = fabs(z);
            z = PI/(z * stirf(q) );
        }
        else
        {
            z = stirf(x);
        }
        return( sgngam * z );
    }
    
    z = 1.0;
    while( x >= 3.0 )
    {
        x -= 1.0;
        z *= x;
    }
    
    while( x < 0.0 )
    {
        if( x > -1.E-9 )
            goto small;
        z /= x;
        x += 1.0;
    }
    
    while( x < 2.0 )
    {
        if( x < 1.e-9 )
            goto small;
        z /= x;
        x += 1.0;
    }
    
    if( x == 2.0 )
        return(z);
    
    x -= 2.0;
    p = polevl( x, P, 6 );
    q = polevl( x, Q, 7 );
    return( z * p / q );
    
    small:
        if( x == 0.0 )
        {
            
            mtherr( "gamma", SING );
            return( MAXNUM );
            
        }
        else
            return( z/((1.0 + 0.5772156649015329 * x) * x) );
}



/* A[]: Stirling's formula expansion of log gamma
 * B[], C[]: log gamma function between 2 and 3
 */


static unsigned short A[] = {
    0x6661,0x2733,0x9850,0x3f4a,
    0xe943,0xb580,0x7fbd,0xbf43,
    0x5ebb,0x20dc,0x019f,0x3f4a,
    0xa5a1,0x16b0,0xc16c,0xbf66,
    0x554b,0x5555,0x5555,0x3fb5
};
static unsigned short B[] = {
    0x6761,0x8ff3,0x8901,0xc095,
    0xb93e,0x355b,0xf234,0xc0e2,
    0x89e5,0xf890,0x3d73,0xc114,
    0xdb51,0xf994,0xbc82,0xc131,
    0xf20b,0x0219,0x4589,0xc13a,
    0x055e,0x5418,0x0c67,0xc12a
};
static unsigned short C[] = {
    /*0x0000,0x0000,0x0000,0x3ff0,*/
    0x12b2,0x1cf3,0xfd0d,0xc075,
    0xd757,0x7b89,0xaa0d,0xc0d0,
    0x4c9b,0xb974,0xeb84,0xc10a,
    0x0043,0x7195,0x6286,0xc131,
    0xf34c,0x892f,0x5255,0xc143,
    0xe14a,0x6a11,0xce4b,0xc13e
};
/* log( sqrt( 2*pi ) ) */
static unsigned short LS2P[] = {
    0xbeb5,0xc864,0x67f1,0x3fed
};
#define LS2PI *(double *)LS2P
#define MAXLGM 2.556348e305


/* Logarithm of gamma function */


double lgam(double x)
{
    double p, q, u, w, z;
    int i;
    
    sgngam = 1;
    
    if( x < -34.0 )
    {
        q = -x;
        w = lgam(q); /* note this modifies sgngam! */
        p = floor(q);
        if( p == q )
        {
            lgsing:
                
                goto loverf;
                
        }
        i = p;
        if( (i & 1) == 0 )
            sgngam = -1;
        else
            sgngam = 1;
        z = q - p;
        if( z > 0.5 )
        {
            p += 1.0;
            z = p - q;
        }
        z = q * sin( PI * z );
        if( z == 0.0 )
            goto lgsing;
        /*	z = log(PI) - log( z ) - w;*/
        z = LOGPI - log( z ) - w;
        return( z );
    }
    
    if( x < 13.0 )
    {
        z = 1.0;
        p = 0.0;
        u = x;
        while( u >= 3.0 )
        {
            p -= 1.0;
            u = x + p;
            z *= u;
        }
        while( u < 2.0 )
        {
            if( u == 0.0 )
                goto lgsing;
            z /= u;
            p += 1.0;
            u = x + p;
        }
        if( z < 0.0 )
        {
            sgngam = -1;
            z = -z;
        }
        else
            sgngam = 1;
        if( u == 2.0 )
            return( log(z) );
        p -= 2.0;
        x = x + p;
        p = x * polevl( x, B, 5 ) / p1evl( x, C, 6);
        return( log(z) + p );
    }
    
    if( x > MAXLGM )
    {
        
        loverf:
            mtherr( "lgam", OVERFLOW );
            return( sgngam * MAXNUM );
            
    }
    
    q = ( x - 0.5 ) * log(x) - x + LS2PI;
    if( x > 1.0e8 )
        return( q );
    
    p = 1.0/(x*x);
    if( x >= 1000.0 )
        q += ((   7.9365079365079365079365e-4 * p
                - 2.7777777777777777777778e-3) *p
                + 0.0833333333333333333333) / x;
    else
        q += polevl( p, A, 4 ) / x;
    return( q );
}
