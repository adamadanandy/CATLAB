#ifndef __CATLAB_EIG_H__
#define __CATLAB_EIG_H__

#include <catlab_util.h>
#include <catlab_matrix.h>
#include <float.h>

typedef struct _tagCatzEigsWorkspace {
    int             isreuse;
    int             ndim;
    int             nev;
    int             ncv;
    int             lworkl;
    int             maxit;
    double          tol;
    int             iparam[11];
    int             ipntr[14];

    double _Complex *resid;
    double _Complex *v;
    double _Complex *workd;
    double _Complex *workl;
    double          *rwork;
    
    int            *selectarray;
    double _Complex *d;
    double _Complex *z;
    double _Complex *workev;
    
    int             bPrepared;

} catz_eigs_workspace, *ptr_catz_eigs_workspace;

cat_bool cat_eig(ptr_cat_mat mata, ptr_cat_zmat d, ptr_cat_zmat v);
cat_bool cat_zmat_eig(ptr_cat_zmat mata, ptr_cat_zmat d, ptr_cat_zmat v);
cat_bool cat_dmat_eig(ptr_cat_dmat mata, ptr_cat_zmat d, ptr_cat_zmat v);
cat_bool cat_zmat_geig(ptr_cat_zmat mata, ptr_cat_zmat matb,
            ptr_cat_zmat d, ptr_cat_zmat v);
cat_bool cat_dmat_geig(ptr_cat_dmat mata, ptr_cat_dmat matb,
            ptr_cat_zmat d, ptr_cat_zmat v);


cat_bool    cat_zmat_eigs(ptr_cat_zmat mata, ptr_cat_zmat matd,
                ptr_cat_zmat matv, cat_z target, int neig,
                ptr_catz_eigs_workspace* pzw);
void        cat_zmat_eigs_workspace_dealloc(ptr_catz_eigs_workspace zw);
ptr_catz_eigs_workspace cat_zmat_eigs_workspace_alloc(int n, int neig);
#endif//__CATLAB_EIG_H__