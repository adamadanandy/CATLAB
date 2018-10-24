#ifndef __CATLAB_SPECFUNC_H__
#define __CATLAB_SPECFUNC_H__

#include <catlab_util.h>
#include <math.h>
#include <string.h>
#include <complex.h>

#ifndef PI
#define PI             3.14159265358979323846264338327950288419716939937510
#endif//PI

cat_z cat_z_besselj(double nu, cat_z z);
double cat_d_besselj(double nu, double z);
cat_z cat_z_besselj1(cat_z z);
cat_z cat_z_besselj0(cat_z z);
cat_d cat_d_besselj1(cat_d d);
cat_d cat_d_besselj0(cat_d d);

#endif//__CATLAB_SPECFUNC_H__