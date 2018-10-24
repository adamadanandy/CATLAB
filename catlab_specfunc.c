#include "catlab_specfunc.h"
#include "cephes.h"

extern void zbesj_(void*, void*, const double*, const int*, const int*, void*, void*, void*, void*);
extern void zbesy_(void*, void*, void*, void*, void*, void*, void*, void*, void*, void*, void*);



cat_z cat_z_besselj1(cat_z z)
{
    static double zr,zi,cyr,cyi;
    const double nu=1.0;
    const int n=1;
    const int kode=1;
    static int ierr,nz;
    zr = CAT_REAL(z);
    zi = CAT_IMAG(z);
    zbesj_(&zr,&zi,&nu,&kode,&n,&cyr,&cyi,&nz,&ierr);
    if (ierr != 0)
    {
        //printf("zr=%15.9e,zi=%15.9e,ierr=%d\n",zr,zi,ierr);
        return INFINITY+INFINITY*I;
    }
    
    return cyr+cyi*I;
}

cat_z cat_z_besselj0(cat_z z)
{
    static double zr,zi,cyr,cyi;
    const double nu=0.0;
    const int n=1;
    const int kode=1;
    static int ierr,nz;
    zr = CAT_REAL(z);
    zi = CAT_IMAG(z);
    zbesj_(&zr,&zi,&nu,&kode,&n,&cyr,&cyi,&nz,&ierr);
    return cyr+cyi*I;
}

cat_d cat_d_besselj1(cat_d d)
{
    return j1(d);
}

cat_d cat_d_besselj0(cat_d d)
{
    return j0(d);
}

cat_z cat_z_besselj(double nu, cat_z z)
/*return J_nu(z)*/
{
    static double zr, zi, cyr, cyi, cwrkr, cwrki;
    static int n, kode, nz, ierr;
    cat_z rtval, zache;
    n = nu;
    if ((double)n == nu && CAT_IMAG(z) == 0)
    {
        if (n == 0)
        {
            zr = CAT_REAL(z);
            rtval = j0(zr);
        }
        else if (n == 1)
        {
            zr = CAT_REAL(z);
            rtval = j1(zr);
        }
        else if (n == -1)
        {
            zr = CAT_REAL(z);
            rtval = j1(zr);
            rtval = -rtval;
        }
        else
        {
            zr = CAT_REAL(z);
            rtval = jn(n,zr);
            if (n%2)
            {
                rtval = -rtval;
            }
        }
    }
    else
    {
        n = 1;
        zr = CAT_REAL(z);
        zi = CAT_IMAG(z);
        kode = 1;
        if (nu >= 0)
        {
            zbesj_(&zr,&zi,&nu,&kode,&n,&cyr,&cyi,&nz,&ierr);
            rtval = cyr+cyi*I;
        }
        else
        {
            nu = -nu;
            zbesj_(&zr,&zi,&nu,&kode,&n,&cyr,&cyi,&nz,&ierr);
            rtval = cyr+cyi*I;
            zbesy_(&zr,&zi,&nu,&kode,&n,&cyr,&cyi,&nz,&cwrkr,&cwrki,&ierr);
            zache = cyr+cyi*I;
            rtval = rtval*cos(PI*nu)-zache*sin(PI*nu);
        }
    }
    return rtval;
}

double cat_d_besselj(double nu, double z)
{
    static double  zr, zi, cyr, cyi, cwrkr, cwrki;
    static int n, kode, nz, ierr;
    double rtval, dache;
    n = nu;
    if ((double)n == nu)
    {
        if (n == 0)
        {
            rtval = j0(z);
        }
        else if (n == 1)
        {
            rtval = j1(z);
        }
        else if (n == -1)
        {
            rtval = j1(z);
            rtval = -rtval;
        }
        else
        {
            zr = CAT_REAL(z);
            rtval = jn(n,zr);
        }
    }
    else
    {
        n = 1;
        kode = 1;
        zr = z;
        zi = 0;
        if (nu >= 0)
        {
            zbesj_(&zr,&zi,&nu,&kode,&n,&cyr,&cyi,&nz,&ierr);
            rtval = cyr;
        }
        else
        {
            nu = -nu;
            zbesj_(&zr,&zi,&nu,&kode,&n,&cyr,&cyi,&nz,&ierr);
            rtval = cyr;
            zbesy_(&zr,&zi,&nu,&kode,&n,&cyr,&cyi,&nz,&cwrkr,&cwrki,&ierr);
            dache = cyr;
            rtval = rtval*cos(PI*nu)-dache*sin(PI*nu);
        }
    }
    return rtval;
}
