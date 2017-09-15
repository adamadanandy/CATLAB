#include <catlab_specfunc.h>

extern void zbesj_(void*, void*, void*, void*, void*, void*, void*, void*, void*);
extern void zbesy_(void*, void*, void*, void*, void*, void*, void*, void*, void*, void*, void*);
extern double dbesj0_(void*);
extern double dbesj1_(void*);
extern double dbesjn_(void*,void*);

cat_z cat_z_besselj(double nu, cat_z z)
/*return J_nu(z)*/
{
    static double zr, zi, cyr, cyi, cwrkr, cwrki;
    static int n, kode, nz, ierr;
    cat_z rtval, zache;
    n = nu;
    if ((double)n == nu && CAT_IMAG(z) == 0) {
        if (n == 0) {
            zr = CAT_REAL(z);
            rtval = dbesj0_(&zr);
        }
        else if (n == 1) {
            zr = CAT_REAL(z);
            rtval = dbesj1_(&zr);
        }
        else if (n == -1) {
            zr = CAT_REAL(z);
            rtval = dbesj1_(&zr);
            rtval = -rtval;
        }
        else {
            zr = CAT_REAL(z);
            rtval = dbesjn_(&n,&zr);
            if (n%2) {
                rtval = -rtval;
            }
        }
    }
    else {
        n = 1;
        zr = CAT_REAL(z);
        zi = CAT_IMAG(z);
        kode = 1;
        if (nu >= 0) {
            zbesj_(&zr,&zi,&nu,&kode,&n,&cyr,&cyi,&nz,&ierr);
            rtval = cyr+cyi*I;
        }
        else {
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
    if ((double)n == nu && CAT_IMAG(z) == 0) {
        if (n == 0) {
            rtval = dbesj0_(&z);
        }
        else if (n == 1) {
            rtval = dbesj1_(&z);
        }
        else if (n == -1) {
            rtval = dbesj1_(&z);
            rtval = -rtval;
        }
        else {
            zr = CAT_REAL(z);
            rtval = dbesjn_(&n,&zr);
        }
    }
    else {
        n = 1;
        kode = 1;
        zr = z;
        zi = 0;
        if (nu >= 0) {
            zbesj_(&zr,&zi,&nu,&kode,&n,&cyr,&cyi,&nz,&ierr);
            rtval = cyr;
        }
        else {
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