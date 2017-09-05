/*  <title title="eig">
 *      <section title="double precision eigenvalue problem">
 *          <subsection title="real symmetric A">
 *              <routine>
 *                  <lapack subroutine="DSYEV">
 *              </routine>
 *          </subsection>
 *          <subsection title="real nonsymmetric A">
 *              <list>
 *                  <item title="with preliminary balance step">
 *                      <routine>
 *                          <lapack subroutine="DGEEV">
 *                      </routine>
 *                  </item>
 *                  <item title="d = eig(A,'nobalance')">
 *                      <routine>
 *                          <lapack	subroutine="DGEHRD">
 *                          <lapack	subroutine="DHSEQR">
 *                      </routine>
 *                  </item>
 *                  <item title="[V,D] = eig(A,'nobalance')">
 *                      <routine>
 *                          <lapack	subroutine="DGEHRD">
 *                          <lapack	subroutine="DORGHR">
 *                          <lapack	subroutine="DHSEQR">
 *                          <lapack	subroutine="DTREVC">
 *                      </routine>
 *              </list>
 *          </subsection>
 *          <subsection title="Hermitian A">
 *              <routine>
 *                  <lapack	subroutine="ZHEEV">
 *              </routine>
 *          </subsection>
 *          <subsection title="Non-Hermitian A">
 *              <list>
 *                  <item title="With preliminary balance step">
 *                      <routine>
 *                          <lapack subroutine="ZGEEV">
 *                      </routine>
 *                  </item>
 *                  <item title="d = eig(A,'nobalance')">
 *                      <routine>
 *                          <lapack	subroutine="ZGEHRD">
 *                          <lapack subroutine="ZHSEQR">
 *                      </routine>
 *                  </item>
 *                  <item title="[V,D] = eig(A,'nobalance')">
 *                      <routine>
 *                          <lapack	subroutine="ZGEHRD">
 *                          <lapack	subroutine="ZORGHR">
 *                          <lapack	subroutine="ZHSEQR">
 *                          <lapack	subroutine="ZTREVC">
 *                      </routine>
 *                  </item>
 *              </list>
 *          </subsection>
 *      </section>
 *      <section title="double precision general eigenvalue problem">
 *          <subsection title="Real symmetric A, symmetric positive definite B">
 *              <routine>
 *                  <lapack subroutine="DSYGV">
 *              </routine>
 *          </subsection>
 *          <subsection title="eig(A,B,'qz') for real A,B">
 *              <routine>
 *                  <lapack subroutine="DGGEV">
 *              </routine>
 *          </subsection>
 *          <subsection title="Real nonsymmetric A, real general B">
 *              <routine>
 *                  <lapack subroutine="DGGEV">
 *              </routine>
 *          </subsection>
 *          <subsection title="Complex Hermitian A, Hermitian positive definite B">
 *              <routine>
 *                  <lapack subroutine="ZHEGV">
 *              </routine>
 *          </subsection>
 *          <subsection title="eig(A,B,'qz') for complex A or B">
 *              <routine>
 *                  <lapack subroutine="ZGGEV">
 *              </routine>
 *          </subsection>
 *          <subsection title="Complex non-Hermitian A, complex B">
 *              <routine>
 *                  <lapack subroutine="ZGGEV">
 *              </routine>
 *          </subsection>
 *      </section>
 *  </title>
 */

#include <catlab_eig.h>
#include <string.h>

/*
 * MATLAB: eig
 *
 * CATLAB: cat_eig, cat_geig
 *
 * MATLAB: eigs
 * 
 * CATLAB: cat_eigs, cat_geigs
 */




cat_bool cat_zmat_eig(ptr_cat_zmat matb, ptr_cat_zmat d, ptr_cat_zmat v)
{
    int info,n;
    cat_z *w,*vl,*vr;
    
    CAT_ASSERT(matb->bPrepared == CAT_TRUE);
    
    ptr_cat_zmat mata = (ptr_cat_zmat)cat_matrixDuplicate((ptr_cat_mat)matb);
    n = mata->shape[0];
    w = malloc(sizeof(cat_z)*n);
    vl = malloc(sizeof(cat_z)*1*n);
    vr = malloc(sizeof(cat_z)*n*n);
    info = LAPACKE_zgeev(LAPACK_COL_MAJOR, 'N', 'V',
                mata->shape[0], mata->data, mata->shape[0],
                w, vl, 1, vr, n);
    if (info != 0) {
        printf("info == %d\n",info);
        free(w);
        free(vl);
        free(vr);
        cat_matrixDestructor((cat_mat**)&mata);
        return CAT_FALSE;
    }
    else {
        d->data = w;
        d->ldata = n;
        d->lshape = 1;
        d->shape = malloc(sizeof(int));
        d->shape[0] = n;
        d->bPrepared = CAT_TRUE;
        d->datatype = CAT_Z;
        d->matrixtype = CAT_GEM;
        v->data = vr;
        v->ldata = n*n;
        v->lshape = 2;
        v->shape = malloc(sizeof(int)*2);
        v->shape[0] = n;
        v->shape[1] = n;
        v->bPrepared = CAT_TRUE;
        v->datatype = CAT_Z;
        v->matrixtype = CAT_GEM;
        free(vl);
        cat_matrixDestructor((cat_mat**)&mata);
        return CAT_TRUE;
    }
}

cat_bool cat_dmat_eig(ptr_cat_dmat matb, ptr_cat_zmat d, ptr_cat_zmat v)
{
    int info, n;
    int itr, itrj;
    cat_d *wr, *wi, *vl, *vr;
    double _Complex *wc, *wv;
    
    CAT_ASSERT(matb->bPrepared == CAT_TRUE);
    
    ptr_cat_dmat mata = (ptr_cat_dmat)cat_matrixDuplicate((ptr_cat_mat)matb);
    n = mata->shape[0];
    wr = malloc(sizeof(cat_d)*n);
    wi = malloc(sizeof(cat_d)*n);
    vl = malloc(sizeof(cat_d)*1*n);
    vr = malloc(sizeof(cat_d)*n*n);
    info = LAPACKE_dgeev(LAPACK_COL_MAJOR, 'N', 'V',
                mata->shape[0], mata->data, mata->shape[0],
                wr, wi, vl, 1, vr, n);
    if (info != 0) {
        printf("info == %d\n",info);
        free(wr);
        free(wi);
        free(vl);
        free(vr);
        cat_matrixDestructor((cat_mat**)&mata);
        return CAT_FALSE;
    }
    else {
        
        d->ldata = n;
        d->lshape = 1;
        d->shape = malloc(sizeof(int));
        d->shape[0] = n;
        d->bPrepared = CAT_TRUE;
        d->matrixtype = CAT_GEM;
        d->datatype = CAT_Z;

        wc = malloc(sizeof(cat_z)*n);
        wv = malloc(sizeof(cat_z)*n*n);
        
        for (itr = 0; itr < n; itr++) {
            wc[itr] = wr[itr] + wi[itr]* _Complex_I;
            if (wi[itr] != 0.0) {
                wc[itr+1] = wr[itr+1] + wi[itr+1]* _Complex_I;
                for (itrj = 0; itrj < n; itrj++) {
                    wv[itrj+itr*n]      = vr[itrj+itr*n] + vr[itrj+(itr+1)*n]*_Complex_I;
                    wv[itrj+(itr+1)*n]  = vr[itrj+itr*n] - vr[itrj+(itr+1)*n]*_Complex_I;
                }
                itr++;
            }
            else {
                for (itrj = 0; itrj < n; itrj++) {
                    wv[itrj+itr*n] = vr[itrj+itr*n];
                }
            }
        }
        
        d->data = (cat_z*)wc;
        
        free(wr);
        free(wi);
        free(vr);
        
        v->data = (cat_z*)wv;
        v->ldata = n*n;
        v->lshape = 2;
        v->shape = malloc(sizeof(int)*2);
        v->shape[0] = n;
        v->shape[1] = n;
        v->bPrepared = CAT_TRUE;
        v->datatype = CAT_Z;
        v->matrixtype = CAT_GEM;
        free(vl);
        cat_matrixDestructor((cat_mat**)&mata);
        return CAT_TRUE;
    }
}

cat_bool cat_eig(ptr_cat_mat mata, ptr_cat_zmat d, ptr_cat_zmat v)
{
    CAT_ASSERT(mata->bPrepared == CAT_TRUE);
    
    if (CAT_CHECKFLAG(mata->datatype,CAT_S)) {
        CAT_LOG(TODO);
    }
    else if (CAT_CHECKFLAG(mata->datatype,CAT_C)) {
        CAT_LOG(TODO);
    }
    else if (CAT_CHECKFLAG(mata->datatype,CAT_D)) {
        if (CAT_CHECKFLAG(mata->matrixtype,CAT_GEM)) {
            return cat_dmat_eig((ptr_cat_dmat)mata, d, v);
        }
    }
    else if (CAT_CHECKFLAG(mata->datatype,CAT_Z)) {
        if (CAT_CHECKFLAG(mata->matrixtype,CAT_GEM)) {
            return cat_zmat_eig((ptr_cat_zmat)mata, d, v);
        }
    }
    return CAT_FALSE;
}

cat_bool cat_dmat_geig(ptr_cat_dmat mata, ptr_cat_dmat matb,
            ptr_cat_zmat d, ptr_cat_zmat v)
{
    int n, itr, info, itrj;
    cat_d *alphar, *alphai, *beta, *vl, *vr;
    double _Complex *datad, *datav;
    CAT_ASSERT(mata->bPrepared == CAT_TRUE);
    CAT_ASSERT(matb->bPrepared == CAT_TRUE);
    ptr_cat_dmat a = (ptr_cat_dmat)cat_matrixDuplicate((ptr_cat_mat)mata);
    ptr_cat_dmat b = (ptr_cat_dmat)cat_matrixDuplicate((ptr_cat_mat)matb);
    n = mata->shape[0];
    alphar = malloc(sizeof(cat_d)*n);
    alphai = malloc(sizeof(cat_d)*n);
    beta = malloc(sizeof(cat_d)*n);
    vl = malloc(sizeof(cat_d)*1*n);
    vr = malloc(sizeof(cat_d)*n*n);
    info = LAPACKE_dggev(LAPACK_COL_MAJOR, 'N', 'V', n, a->data, n,
         	        b->data, n, alphar, alphai, beta, vl, 1, vr, n );
    if (info != 0) {
        CAT_ERR(info);
        free(alphar);
        free(alphai);
        free(beta);
        free(vl);
        free(vr);
        cat_matrixDestructor((cat_mat**)&a);
        cat_matrixDestructor((cat_mat**)&b);
        return CAT_FALSE;
    }
    else {
        datad = malloc(sizeof(double _Complex)*n);
        datav = malloc(sizeof(double _Complex)*n*n);
        
        d->ldata = n;
        d->lshape = 1;
        d->shape = malloc(sizeof(int));
        d->shape[0] = n;
        d->bPrepared = CAT_TRUE;
        d->matrixtype = CAT_GEM;
        d->datatype = CAT_Z;
        
        for (itr = 0; itr < n; itr++) {
            datad[itr] = (alphar[itr]+alphai[itr]*_Complex_I)/beta[itr];
            if (alphai[itr] != 0.0) {
                datad[itr+1] = (alphar[itr+1]+alphai[itr+1]*_Complex_I)/beta[itr+1];
                for (itrj = 0; itrj < n; itrj++) {
                    datav[itrj+itr*n]       = vr[itrj+itr*n] + vr[itrj+(itr+1)*n]*_Complex_I;
                    datav[itrj+(itr+1)*n]   = vr[itrj+itr*n] - vr[itrj+(itr+1)*n]*_Complex_I;
                }
                itr++;
            }
            else {
                for (itrj = 0; itrj < n; itrj++) {
                    datav[itrj+itr*n] = vr[itrj+itr*n];
                }
            }
        }
        d->data = (cat_z*)datad;

        v->data = (cat_z*)datav;
        v->ldata = n*n;
        v->lshape = 2;
        v->shape = malloc(sizeof(int)*2);
        v->shape[0] = n;
        v->shape[1] = n;
        v->bPrepared = CAT_TRUE;
        v->datatype = CAT_Z;
        v->matrixtype = CAT_GEM;
        
        free(alphar);
        free(alphai);
        free(beta);
        free(vl);
        free(vr);
        
        cat_matrixDestructor((cat_mat**)&a);
        cat_matrixDestructor((cat_mat**)&b);
        return CAT_TRUE;
    }
}

cat_bool cat_zmat_geig(ptr_cat_zmat mata, ptr_cat_zmat matb,
            ptr_cat_zmat d, ptr_cat_zmat v)
{
    int n, itr, info;
    cat_z *alpha, *beta, *vl, *vr;
    double _Complex *pz1, *pz2, *pz3;
    CAT_ASSERT(mata->bPrepared == CAT_TRUE);
    CAT_ASSERT(matb->bPrepared == CAT_TRUE);
    ptr_cat_zmat a = (ptr_cat_zmat)cat_matrixDuplicate((ptr_cat_mat)mata);
    ptr_cat_zmat b = (ptr_cat_zmat)cat_matrixDuplicate((ptr_cat_mat)matb);
    n = mata->shape[0];
    alpha = malloc(sizeof(cat_z)*n);
    beta = malloc(sizeof(cat_z)*n);
    vl = malloc(sizeof(cat_z)*1*n);
    vr = malloc(sizeof(cat_z)*n*n);
    info = LAPACKE_zggev(LAPACK_COL_MAJOR, 'N', 'V', n, a->data, n,
         	        b->data, n, alpha, beta, vl, 1, vr, n );
    if (info != 0) {
        CAT_ERR(info);
        free(alpha);
        free(beta);
        free(vl);
        free(vr);
        cat_matrixDestructor((cat_mat**)&a);
        cat_matrixDestructor((cat_mat**)&b);
        return CAT_FALSE;
    }
    else {
        d->ldata = n;
        d->lshape = 1;
        d->shape = malloc(sizeof(int));
        d->shape[0] = n;
        d->bPrepared = CAT_TRUE;
        d->matrixtype = CAT_GEM;
        d->datatype = CAT_Z;
        
        pz1 = (double _Complex*)vl;
        pz2 = (double _Complex*)alpha;
        pz3 = (double _Complex*)beta;
        for (itr = 0; itr < n; itr++) {
            pz1[itr] = pz2[itr]/ pz3[itr];
        }
        d->data = vl;

        v->data = vr;
        v->ldata = n*n;
        v->lshape = 2;
        v->shape = malloc(sizeof(int)*2);
        v->shape[0] = n;
        v->shape[1] = n;
        v->bPrepared = CAT_TRUE;
        v->datatype = CAT_Z;
        v->matrixtype = CAT_GEM;
        
        free(alpha);
        free(beta);
        
        cat_matrixDestructor((cat_mat**)&a);
        cat_matrixDestructor((cat_mat**)&b);
        return CAT_TRUE;
    }
}

cat_bool cat_geig(ptr_cat_mat mata, ptr_cat_mat matb,
            ptr_cat_mat d, ptr_cat_mat v)
{

    CAT_ASSERT(mata->bPrepared == CAT_TRUE);
    CAT_ASSERT(matb->bPrepared == CAT_TRUE);
    if (CAT_CHECKFLAG(mata->datatype,CAT_S)) {
        CAT_LOG(TODO);
    }
    else if (CAT_CHECKFLAG(mata->datatype,CAT_C)) {
        CAT_LOG(TODO);
    }
    else if (CAT_CHECKFLAG(mata->datatype,CAT_D)) {
        if (CAT_CHECKFLAG(mata->matrixtype,CAT_GEM)) {
            return cat_dmat_geig((ptr_cat_dmat)mata, (ptr_cat_dmat)matb,
                (ptr_cat_zmat)d, (ptr_cat_zmat)v);
        }
    }
    else if (CAT_CHECKFLAG(mata->datatype,CAT_Z)) {
        if (CAT_CHECKFLAG(mata->matrixtype,CAT_GEM)) {
            return cat_zmat_geig((ptr_cat_zmat)mata, (ptr_cat_zmat)matb,
                (ptr_cat_zmat)d, (ptr_cat_zmat)v);
        }
    }
    return CAT_FALSE;
}

/***********************************************************************/
extern void znaupd_(void*, void*, void*, void*, void*, void*, void*, 
                    void*, void*, void*, void*, void*, void*, void*, 
                    void*, void*, void*);

extern void dnaupd_(void*, void*, void*, void*, void*, void*, void*,
                    void*, void*, void*, void*, void*, void*, void*,
                    void*, void*);

extern void zneupd_(void*, void*, void*, void*, void*, void*, void*, void*,
                    void*, void*, void*, void*, void*, void*, void*, void*,
                    void*, void*, void*, void*, void*, void*, void*, void*);

extern void dneupd_(void*, void*, void*, void*, void*, void*, void*, void*,
                    void*, void*, void*, void*, void*, void*, void*, void*,
                    void*, void*, void*, void*, void*, void*, void*, void*,
                    void*);
/***********************************************************************/

ptr_catz_eigs_workspace cat_zmat_eigs_workspace_alloc(int n, int neig)
{
    ptr_catz_eigs_workspace zw;
    zw = malloc(sizeof(catz_eigs_workspace));
    memset(zw,0,sizeof(catz_eigs_workspace));
    zw->tol = DBL_EPSILON;
    zw->nev = neig;
    zw->ndim = n;
    zw->ncv = CAT_MIN(n,CAT_MAX(2*neig,20));
    zw->lworkl = (3*zw->ncv+5)*zw->ncv;
    zw->maxit = CAT_MAX(300,2*zw->ndim/neig);
    zw->resid = malloc(sizeof(double _Complex)*n);
    zw->v     = malloc(sizeof(double _Complex)*n*zw->ncv);
    zw->workd = malloc(sizeof(double _Complex)*3*n);
    zw->rwork = malloc(sizeof(double)*zw->ncv);
    zw->workl = malloc(sizeof(double _Complex)*zw->lworkl);
    
    zw->selectarray = malloc(sizeof(int)*zw->ncv);
    zw->d           = malloc(sizeof(double _Complex)*(zw->nev+1));
    zw->z           = malloc(sizeof(double _Complex)*n*zw->nev);
    zw->workev      = malloc(sizeof(double _Complex)*2*zw->ncv);
    
    zw->iparam[0] = 1;
    zw->iparam[2] = zw->maxit;
    zw->iparam[6] = 3;
    
    zw->bPrepared = CAT_TRUE;
    return zw;
}

void cat_zmat_eigs_workspace_dealloc(ptr_catz_eigs_workspace zw)
{
    if (zw == NULL) {
        return;
    }
    if (zw->bPrepared == CAT_TRUE) {
        free(zw->resid);
        free(zw->v);
        free(zw->workd);
        free(zw->rwork);
        free(zw->selectarray);
        free(zw->d);
        free(zw->z);
        free(zw->workev);
    }
    free(zw);
}

cat_bool cat_zmat_eigs(ptr_cat_zmat mata, ptr_cat_zmat matd,
                ptr_cat_zmat matv, cat_z target, int neig,
                ptr_catz_eigs_workspace* pzw)
{
    static char *bmat = "I";
    static char *which = "LM";
    static char *howmny = "A";
    int ido = 0;
    int ndim;
    int info;
    int itri,itrj;
    char rvec;
    double _Complex sigma;
    ptr_cat_zmat cpmat;
    ptr_catz_eigs_workspace zw;
    ptr_cat_workspace divwp = cat_mat_csc_mldivide_alloc();
    
    CAT_ASSERT(mata->lshape == 2 && mata->shape[0] == mata->shape[1]);
    CAT_ASSERT(mata->datatype == CAT_Z);
    
    ndim = mata->shape[0];
    
    if (pzw != NULL && *pzw != NULL) {
        /*the workspace is already given*/
        zw = *pzw;
    }
    else {
        zw = cat_zmat_eigs_workspace_alloc(ndim, neig);
    }
    for (itri = 0; itri < ndim; itri++) {
        zw->resid[itri] = drand48()+drand48()*I;
    }
    /*
     * In CATLAB, the CSC matrix must contain all the diagonal elements
     * the CSC matrix will be convert from GBM & COO & DIA matrix, which
     * the diagonal elements will be kept by CATLAB automatically.
     */
    cpmat = (ptr_cat_zmat)cat_matrixDuplicate((ptr_cat_mat)mata);
    if (cpmat->matrixtype != CAT_CSC) {
        cat_zmat_matrixGB2CC(cpmat);
    }
    //cpmat = mata-target*eye;
    for (itri = 0; itri < ndim; itri++) {
        for (itrj = cpmat->ptrind1[itri]; itrj < cpmat->ptrind1[itri+1]; itrj++) {
            if (itri == cpmat->ptrind2[itrj]) {
                cpmat->data[itrj] = cpmat->data[itrj] - target;
                break;
            }
        }
    }
    //prepare the solver of [A-lambda*I]
    cat_matz_csc_mldivide(cpmat, NULL, NULL, divwp);
    divwp->isreuse = CAT_TRUE;
    
    while (ido != 99) {
        znaupd_(&ido, bmat, &(zw->ndim), which, &(zw->nev), &(zw->tol),
                zw->resid, &(zw->ncv), zw->v, &(zw->ndim), zw->iparam,
                zw->ipntr, zw->workd, zw->workl, &(zw->lworkl),
                zw->rwork, &info);
        
/* MODE == 3 :> OP := inv[A-lambda*I] :> M := I */
        if (ido == -1) {
/* Y = OP * X :: X := WORKD(IPNTR(1)) :: Y := WORKD(IPNTR(2))*/
            cat_matz_csc_mldivide(cpmat, zw->workd+zw->ipntr[0]-1,
                zw->workd+zw->ipntr[1]-1, divwp);
        }
        else if (ido == 1) {
/* Y = OP * X :: X := WORKD(IPNTR(1)) :: Y := WORKD(IPNTR(2)) :: B * X ::>> WORKD(IPNTR(3)) not need to be recomputed*/
            cat_matz_csc_mldivide(cpmat, zw->workd+zw->ipntr[0]-1,
                zw->workd+zw->ipntr[1]-1, divwp);
        }
        else if (ido == 2) {
/* Y = M * X :: X := WORKD(IPNTR(1)) :: Y := WORKD(IPNTR(2))      */
            memcpy(zw->workd+zw->ipntr[1]-1,zw->workd+zw->ipntr[0]-1,
                sizeof(double _Complex)*ndim);
        }
        else if (ido == 99) {
/*CONVERGED*/
        }
        else {
            printf("Diverge\n");
            return CAT_FALSE;
        }
    }
    if (matv == NULL) {
        rvec = 0;
    }
    else {
        rvec = 1;
    }
    zneupd_(&rvec, howmny, zw->selectarray, zw->d, zw->z,
            &(zw->ndim), &sigma, zw->workev,
            bmat, &ndim, which, &(zw->nev), &(zw->tol),
            zw->resid, &(zw->ncv), zw->v, &(zw->ndim), zw->iparam,
            zw->ipntr, zw->workd, zw->workl, &(zw->lworkl),
            zw->rwork, &info);
    
    matd->datatype = CAT_Z;
    matd->matrixtype = CAT_GEM;
    matd->bPrepared = CAT_TRUE;
    matd->lshape = 1;
    matd->shape = malloc(sizeof(int));
    matd->shape[0] = neig;
    matd->ldata = neig+1;
    matd->data = zw->d;
    zw->d = NULL;
    for (itri = 0; itri < neig; itri++) {
        matd->data[itri] = matd->data[itri]+target;
    }
    
    matv->datatype = CAT_Z;
    matv->matrixtype = CAT_GEM;
    matv->bPrepared = CAT_TRUE;
    if (neig == 1) {
        matv->lshape = 1;
        matv->shape = malloc(sizeof(int));
        matv->shape[0] = ndim;
        matv->ldata = ndim;
        
    }
    else {
        matv->lshape = 2;
        matv->shape = malloc(sizeof(int)*2);
        matv->shape[0] = ndim;
        matv->shape[1] = neig;
        matv->ldata = ndim*neig;
    }
    matv->data = zw->z;
    zw->z = NULL;
    
    if (pzw == NULL) {
        /*User doesn't need the workspace*/
        cat_zmat_eigs_workspace_dealloc(zw);
    }
    else {
        *pzw = zw;
    }
    cat_matrixDestructor((ptr_cat_mat*)&cpmat);
    return CAT_TRUE;
}

cat_bool cat_zmat_geigs(ptr_cat_zmat mata, ptr_cat_zmat matb,
                ptr_cat_zmat matd, ptr_cat_zmat matv,
                cat_z target, int neig, ptr_catz_eigs_workspace* pzw)
{
    static char *bmat = "G";
    static char *which = "LM";
    static char *howmny = "A";
    int ido = 0;
    int ndim;
    int info;
    int itri;
    double _Complex *cachedata;
    char rvec;
    double scaleB,invscaleB;
    double _Complex sigma;
    ptr_cat_zmat cpmat1,cpmat2;
    ptr_catz_eigs_workspace zw;
    ptr_cat_workspace divwp = cat_mat_csc_mldivide_alloc();
    
    CAT_ASSERT(mata->lshape == 2 && mata->shape[0] == mata->shape[1]);
    CAT_ASSERT(matb->lshape == 2 && matb->shape[0] == matb->shape[1]);
    CAT_ASSERT(mata->datatype == CAT_Z);
    CAT_ASSERT(matb->datatype == CAT_Z);
    CAT_ASSERT(mata->shape[0] == matb->shape[0]);
    
    ndim = mata->shape[0];
    
    if (pzw != NULL && *pzw != NULL) {
        /*the workspace is already given*/
        zw = *pzw;
    }
    else {
        zw = cat_zmat_eigs_workspace_alloc(ndim, neig);
    }
    for (itri = 0; itri < ndim; itri++) {
        zw->resid[itri] = drand48()+drand48()*I;
    }
    /*
     * In CATLAB, the CSC matrix must contain all the diagonal elements
     * the CSC matrix will be convert from GBM & COO & DIA matrix, which
     * the diagonal elements will be kept by CATLAB automatically.
     */
    sigma = target;
    cpmat1 = (ptr_cat_zmat)cat_matrixDuplicate((ptr_cat_mat)mata);
    cpmat2 = (ptr_cat_zmat)cat_matrixDuplicate((ptr_cat_mat)matb);
    if (cpmat1->matrixtype != CAT_CSC) {
        cat_zmat_matrixGB2CC(cpmat1);
    }
    if (cpmat2->matrixtype != CAT_CSC) {
        cat_zmat_matrixGB2CC(cpmat2);
    }
    //cpmat1 = cpmat1-sigma*cpmat2
    cat_zmat_AXPBY(-sigma, cpmat2, 1.0, cpmat1);
    //scaleB = norm2(B=>vector)/sqrt(ndim)
    scaleB = cblas_dznrm2(cpmat2->ptrind1[ndim],cpmat2->data,1)/sqrt(ndim);
    //scaleB = 2^(floor(log2(scaleB+1)))
    scaleB = exp2((int)(log2(scaleB+1)));
    printf("scaleB=%g\n",scaleB);
    invscaleB = 1.0/scaleB;
    //cpmat2 = cpmat2/scaleB
    cblas_zscal(cpmat2->ptrind1[ndim],&invscaleB, cpmat2->data, 1);
    //sigma = sigma*scaleB
    sigma = sigma*scaleB;
    

    //prepare the solver of [A-sigma*B]
    cat_matz_csc_mldivide(cpmat1, NULL, NULL, divwp);
    divwp->isreuse = CAT_TRUE;
    
    cachedata = malloc(sizeof(double _Complex)*zw->ndim);
    
    while (ido != 99) {
        znaupd_(&ido, bmat, &(zw->ndim), which, &(zw->nev), &(zw->tol),
                zw->resid, &(zw->ncv), zw->v, &(zw->ndim), zw->iparam,
                zw->ipntr, zw->workd, zw->workl, &(zw->lworkl),
                zw->rwork, &info);
        
/* MODE == 3 :> OP := inv[A-lambda*I] :> M := I */
        if (ido == -1) {
/* Y = OP * X :: X := WORKD(IPNTR(1)) :: Y := WORKD(IPNTR(2))*/
            cat_zcscgemv(CblasNoTrans, cpmat2->shape[1], 
                cpmat2->data, cpmat2->ptrind1, cpmat2->ptrind2,
                zw->workd+zw->ipntr[0]-1, cachedata);
            cat_matz_csc_mldivide(cpmat1, cachedata,
                zw->workd+zw->ipntr[1]-1, divwp);
        }
        else if (ido == 1) {
/* Y = OP * X :: X := WORKD(IPNTR(1)) :: Y := WORKD(IPNTR(2)) :: B * X ::>> WORKD(IPNTR(3)) not need to be recomputed*/
            cat_zcscgemv(CblasNoTrans, cpmat2->shape[1], 
                cpmat2->data, cpmat2->ptrind1, cpmat2->ptrind2,
                zw->workd+zw->ipntr[0]-1, cachedata);
            cat_matz_csc_mldivide(cpmat1, cachedata,
                zw->workd+zw->ipntr[1]-1, divwp);
        }
        else if (ido == 2) {
/* Y = M * X :: X := WORKD(IPNTR(1)) :: Y := WORKD(IPNTR(2))      */
            memcpy(zw->workd+zw->ipntr[1]-1,zw->workd+zw->ipntr[0]-1,
                sizeof(double _Complex)*ndim);
        }
        else if (ido == 99) {
/*CONVERGED*/
        }
        else {
            printf("Diverge\n");
            return CAT_FALSE;
        }
    }
    if (matv == NULL) {
        rvec = 0;
    }
    else {
        rvec = 1;
    }
    zneupd_(&rvec, howmny, zw->selectarray, zw->d, zw->z,
            &(zw->ndim), &sigma, zw->workev,
            bmat, &ndim, which, &(zw->nev), &(zw->tol),
            zw->resid, &(zw->ncv), zw->v, &(zw->ndim), zw->iparam,
            zw->ipntr, zw->workd, zw->workl, &(zw->lworkl),
            zw->rwork, &info);
    
    matd->datatype = CAT_Z;
    matd->matrixtype = CAT_GEM;
    matd->bPrepared = CAT_TRUE;
    matd->lshape = 1;
    matd->shape = malloc(sizeof(int));
    matd->shape[0] = neig;
    matd->ldata = neig+1;
    matd->data = zw->d;
    zw->d = NULL;
    for (itri = 0; itri < neig; itri++) {
        matd->data[itri] = matd->data[itri]/scaleB;
    }
    
    matv->datatype = CAT_Z;
    matv->matrixtype = CAT_GEM;
    matv->bPrepared = CAT_TRUE;
    if (neig == 1) {
        matv->lshape = 1;
        matv->shape = malloc(sizeof(int));
        matv->shape[0] = ndim;
        matv->ldata = ndim;
        
    }
    else {
        matv->lshape = 2;
        matv->shape = malloc(sizeof(int)*2);
        matv->shape[0] = ndim;
        matv->shape[1] = neig;
        matv->ldata = ndim*neig;
    }
    matv->data = zw->z;
    zw->z = NULL;
    
    if (pzw == NULL) {
        /*User doesn't need the workspace*/
        cat_zmat_eigs_workspace_dealloc(zw);
    }
    else {
        *pzw = zw;
    }
    cat_matrixDestructor((ptr_cat_mat*)&cpmat1);
    return CAT_TRUE;
}

/***********************************************************************/
/*
ptr_catd_eigs_workspace cat_dmat_eigs_workspace_alloc(int n, int neig)
{
    ptr_catd_eigs_workspace zw;
    zw = malloc(sizeof(catd_eigs_workspace));
    memset(zw,0,sizeof(catd_eigs_workspace));
    zw->tol = DBL_EPSILON;
    zw->nev = neig;
    zw->ndim = n;
    zw->ncv = CAT_MIN(n,CAT_MAX(2*neig,20));
    zw->lworkl = (3*zw->ncv+5)*zw->ncv;
    zw->maxit = CAT_MAX(300,2*zw->ndim/neig);
    zw->resid = malloc(sizeof(double _Complex)*n);
    zw->v     = malloc(sizeof(double _Complex)*n*zw->ncv);
    zw->workd = malloc(sizeof(double _Complex)*3*n);
    zw->rwork = malloc(sizeof(double)*zw->ncv);
    zw->workl = malloc(sizeof(double _Complex)*zw->lworkl);
    
    zw->selectarray = malloc(sizeof(int)*zw->ncv);
    zw->d           = malloc(sizeof(double _Complex)*(zw->nev+1));
    zw->z           = malloc(sizeof(double _Complex)*n*zw->nev);
    zw->workev      = malloc(sizeof(double _Complex)*2*zw->ncv);
    
    zw->iparam[0] = 1;
    zw->iparam[2] = zw->maxit;
    zw->iparam[6] = 3;
    
    zw->bPrepared = CAT_TRUE;
    return zw;
}

void cat_zmat_eigs_workspace_dealloc(ptr_catz_eigs_workspace zw)
{
    if (zw == NULL) {
        return;
    }
    if (zw->bPrepared == CAT_TRUE) {
        free(zw->resid);
        free(zw->v);
        free(zw->workd);
        free(zw->rwork);
        free(zw->selectarray);
        free(zw->d);
        free(zw->z);
        free(zw->workev);
    }
    free(zw);
}
*/


