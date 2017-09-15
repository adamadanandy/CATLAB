#include <stdio.h>
#include <catlab.h>
#include <complex.h>

int main(int argc, char** argv)
{
    int itri,itrj;
    double cache1,cache2;
    double _Complex zache1,zache2;
    printf("----------------------Testing BesselJ---------------------\n");
    printf("J_%g(%g)=%g\n",0.0,1.0,cat_d_besselj(0.0,1.0));
    printf("J_%g(%g)=%g\n",1.0,1.0,cat_d_besselj(1.0,1.0));
    printf("J_%g(%g)=%g\n",2.0,1.0,cat_d_besselj(2.0,1.0));
    printf("J_%g(%g)=%g\n",3.0,1.0,cat_d_besselj(3.0,1.0));
    printf("J_%g(%g)=%g\n",0.1,1.0,cat_d_besselj(0.1,1.0));
    printf("J_%g(%g)=%g\n",0.0,1.0,CAT_REAL(cat_z_besselj(0.0,1.0)));
    printf("J_%g(%g)=%g\n",1.0,1.0,CAT_REAL(cat_z_besselj(1.0,1.0)));
    printf("J_%g(%g)=%g\n",2.0,1.0,CAT_REAL(cat_z_besselj(2.0,1.0)));
    printf("J_%g(%g)=%g\n",3.0,1.0,CAT_REAL(cat_z_besselj(3.0,1.0)));
    printf("J_%g(%g)=%g\n",0.1,1.0,CAT_REAL(cat_z_besselj(0.1,1.0)));
    printf("J_%g(%g)=%g\n",-1.0,1.0,cat_d_besselj(-1.0,1.0));
    printf("J_%g(%g)=%g\n",-2.0,1.0,cat_d_besselj(-2.0,1.0));
    printf("J_%g(%g)=%g\n",-3.0,1.0,cat_d_besselj(-3.0,1.0));
    printf("J_%g(%g)=%g\n",-0.1,1.0,cat_d_besselj(-0.1,1.0));
    printf("J_%g(%g)=%g\n",-1.0,1.0,CAT_REAL(cat_z_besselj(-1.0,1.0)));
    printf("J_%g(%g)=%g\n",-2.0,1.0,CAT_REAL(cat_z_besselj(-2.0,1.0)));
    printf("J_%g(%g)=%g\n",-3.0,1.0,CAT_REAL(cat_z_besselj(-3.0,1.0)));
    printf("J_%g(%g)=%g\n",-0.1,1.0,CAT_REAL(cat_z_besselj(-0.1,1.0)));
    
    zache1 = 1.0+1.0*I;
    zache2 = cat_z_besselj(0.0,zache1);
    printf("J_%g((%g,%g))=(%g,%g)\n",0.0,1.0,1.0,CAT_REAL(zache2),CAT_IMAG(zache2));
    zache1 = 1.0-1.0*I;
    zache2 = cat_z_besselj(0.0,zache1);
    printf("J_%g((%g,%g))=(%g,%g)\n",0.0,1.0,-1.0,CAT_REAL(zache2),CAT_IMAG(zache2));
    zache1 = 1.0+1.0*I;
    zache2 = cat_z_besselj(1.0,zache1);
    printf("J_%g((%g,%g))=(%g,%g)\n",1.0,1.0,1.0,CAT_REAL(zache2),CAT_IMAG(zache2));
    zache1 = 1.0-1.0*I;
    zache2 = cat_z_besselj(1.0,zache1);
    printf("J_%g((%g,%g))=(%g,%g)\n",1.0,1.0,-1.0,CAT_REAL(zache2),CAT_IMAG(zache2));
    zache1 = 1.0+1.0*I;
    zache2 = cat_z_besselj(3.0,zache1);
    printf("J_%g((%g,%g))=(%g,%g)\n",3.0,1.0,1.0,CAT_REAL(zache2),CAT_IMAG(zache2));
    zache1 = 1.0-1.0*I;
    zache2 = cat_z_besselj(3.0,zache1);
    printf("J_%g((%g,%g))=(%g,%g)\n",3.0,1.0,-1.0,CAT_REAL(zache2),CAT_IMAG(zache2));
    zache1 = 1.0+1.0*I;
    zache2 = cat_z_besselj(0.1,zache1);
    printf("J_%g((%g,%g))=(%g,%g)\n",0.1,1.0,1.0,CAT_REAL(zache2),CAT_IMAG(zache2));
    zache1 = 1.0-1.0*I;
    zache2 = cat_z_besselj(0.1,zache1);
    printf("J_%g((%g,%g))=(%g,%g)\n",0.1,1.0,-1.0,CAT_REAL(zache2),CAT_IMAG(zache2));
    
    zache1 = 1.0+1.0*I;
    zache2 = cat_z_besselj(-1.0,zache1);
    printf("J_%g((%g,%g))=(%g,%g)\n",-1.0,1.0,1.0,CAT_REAL(zache2),CAT_IMAG(zache2));
    zache1 = 1.0-1.0*I;
    zache2 = cat_z_besselj(-1.0,zache1);
    printf("J_%g((%g,%g))=(%g,%g)\n",-1.0,1.0,-1.0,CAT_REAL(zache2),CAT_IMAG(zache2));
    zache1 = 1.0+1.0*I;
    zache2 = cat_z_besselj(-3.0,zache1);
    printf("J_%g((%g,%g))=(%g,%g)\n",-3.0,1.0,1.0,CAT_REAL(zache2),CAT_IMAG(zache2));
    zache1 = 1.0-1.0*I;
    zache2 = cat_z_besselj(-3.0,zache1);
    printf("J_%g((%g,%g))=(%g,%g)\n",-3.0,1.0,-1.0,CAT_REAL(zache2),CAT_IMAG(zache2));
    zache1 = 1.0+1.0*I;
    zache2 = cat_z_besselj(-0.1,zache1);
    printf("J_%g((%g,%g))=(%g,%g)\n",-0.1,1.0,1.0,CAT_REAL(zache2),CAT_IMAG(zache2));
    zache1 = 1.0-1.0*I;
    zache2 = cat_z_besselj(-0.1,zache1);
    printf("J_%g((%g,%g))=(%g,%g)\n",-0.1,1.0,-1.0,CAT_REAL(zache2),CAT_IMAG(zache2));
}