#include <catlab_matrix.h>

ptr_cat_mat cat_EmptyMatConstructor()
{
    ptr_cat_mat rtmat;
    rtmat = malloc(sizeof(cat_mat));
    memset(rtmat,0,sizeof(cat_mat));
    return rtmat;
}

void cat_zvec_add(ptr_cat_zmat vec, double _Complex addnum)
{
    int itr;
    CAT_ASSERT(vec->lshape == 1);
    for (itr = 0; itr < vec->shape[0]; itr++)
    {
        vec->data[itr] += addnum;
    }
}
extern int izamin_(void*, void*, void*);
void cat_zvec_absmin_index(ptr_cat_zmat vec, int* ptrindex)
{
    static int incx = 1;
    CAT_ASSERT(vec->lshape == 1);
    *ptrindex = izamin_(vec->shape, vec->data, &incx);
    *ptrindex -= 1;
}
// Y = aX+bY
void cat_zmat_gem_AXPBY(double _Complex a, ptr_cat_zmat matx,
                double _Complex b, ptr_cat_zmat maty)
{
    cat_zmat_SCALE(maty,b);
    cblas_zaxpy(matx->ldata,&a,matx->data,1,maty->data,1);
}

void cat_zmat_gbm_AXPBY(double _Complex a, ptr_cat_zmat matx,
                double _Complex b, ptr_cat_zmat maty)
{
    int maxkl, maxku;
    cat_zmat_SCALE(maty,b);
    if (maty->kl != matx->kl || maty->ku != matx->ku) {
        maxkl = CAT_MAX(matx->kl,maty->kl);
        maxku = CAT_MAX(matx->ku,maty->ku);
        cat_zmat_gbm_Adjust(matx,maxku,maxkl);
        cat_zmat_gbm_Adjust(maty,maxku,maxkl);
    }
    cblas_zaxpy(matx->ldata,&a,matx->data,1,maty->data,1);
}

void cat_zmat_gbm_Adjust(ptr_cat_zmat pmat, int newku, int newkl)
{
    double _Complex *ptrdata;
    int lcol;
    int itri,itrj,ofst;
    if (pmat->ku == newku && pmat->kl == newkl) {
        return;
    }
    if (pmat->ku > newku || pmat->kl > newkl) {
        printf("ILLEGAL\n");
        return;
    }
    lcol = pmat->shape[1];
    ptrdata = malloc(sizeof(double _Complex)*(newku+newkl+1)*lcol);
    memset(ptrdata,0,sizeof(double _Complex)*(newku+newkl+1)*lcol);
    for (itrj = 0; itrj < lcol; itrj++) {
        ofst = newku - pmat->ku;
        for (itri = 0; itri < pmat->ku+pmat->kl+1; itri++) {
            ptrdata[itri+ofst+itrj*(newku+newkl+1)] = pmat->data[itri+itrj*(pmat->ku+pmat->kl+1)];
        }
    }
    
    free(pmat->data);
    pmat->data = ptrdata;
    pmat->ldata = (newku+newkl+1)*lcol;
    pmat->ku = newku;
    pmat->kl = newkl;
}

void cat_zmat_SCALE(ptr_cat_zmat pmat, double _Complex alpha)
{
    cblas_zscal(pmat->ldata,&alpha,pmat->data,1);
}

void cat_zmat_csc_AXPBY(double _Complex a, ptr_cat_zmat matx,
    double _Complex b, ptr_cat_zmat maty) {
    CAT_ASSERT(matx->datatype == CAT_Z && maty->datatype == CAT_Z);
    CAT_ASSERT(matx->lshape == 2 && maty->lshape == 2);
    CAT_ASSERT(matx->shape[0] == maty->shape[0]);
    CAT_ASSERT(matx->shape[1] == maty->shape[1]);
    cat_zcscaxpby(matx->shape[1], &a, matx->data, matx->ptrind1, matx->ptrind2,
        &b, &(maty->data), &(maty->ptrind1), &(maty->ptrind2));
    maty->ldata = maty->ptrind1[10];
}

/* matrix should be stored in column priority (fortran like) form*/
ptr_cat_zmat cat_zmat_matrixCombineCC2GE(int nrow, int ncol, ptr_cat_zmat *ptrmat)
{

    int itr, itri, itrj, cache1, offsetr, offsetc, matrow, matcol;
    int nmat, itrmatr, itrmatc;
    int totalrow, totalcol;
    int *matrows, *matcols, *rows, *cols;
    int matshape[2];
    ptr_cat_zmat rtmat, cachemat;


    nmat = nrow*ncol;
    matrows = malloc(sizeof(int)*nmat);
    matcols = malloc(sizeof(int)*nmat);
    rows = malloc(sizeof(int)*nrow);
    cols = malloc(sizeof(int)*ncol);
    
    for (itr = 0; itr<nmat; itr++)
    {
        if (ptrmat[itr] == NULL)
        {
            matrows[itr] = 0;
            matcols[itr] = 0;
        }
        else
        {
            matrows[itr] = ptrmat[itr]->shape[0];
            matcols[itr] = ptrmat[itr]->shape[1];
        }
    }


    totalrow = 0;
    totalcol = 0;
    
    for (itrmatc = 0; itrmatc<ncol; itrmatc++)
    {
        for (itrmatr = 0; itrmatr<nrow; itrmatr++)
        {
            cache1 = matcols[itrmatc*nrow+itrmatr];
            if (cache1 != 0)
            {
                break;
            }
        }
        
        for (itrmatr = 0; itrmatr<nrow; itrmatr++)
        {
            if (matcols[itrmatc*nrow+itrmatr] == 0)
            {
                matcols[itrmatc*nrow+itrmatr] = cache1;
            }
            else if (matcols[itrmatc*nrow+itrmatr] != cache1)
            {
                printf("Warning!! Matrix Combination don't match!\n");
                return NULL;
            }
        }
    }

    for (itrmatr = 0; itrmatr<nrow; itrmatr++)
    {
        itr = 0;
        
        for (itrmatc = 0; itrmatc<ncol; itrmatc++)
        {
            cache1 = matrows[itrmatr+itrmatc*nrow];
            if (cache1 != 0)
            {
                break;
            }
        }
        
        for (itrmatc = 0; itrmatc<ncol; itrmatc++)
        {
            if (matrows[itrmatc*nrow+itrmatr] == 0)
            {
                matrows[itrmatc*nrow+itrmatr] = cache1;
            }
            else if (matrows[itrmatc*nrow+itrmatr] != cache1)
            {
                printf("Warning!! Matrix Combination don't match!\n");
                return NULL;
            }
        }
    }

    rows[0] = 0;
    cols[0] = 0;
    for (itrmatr = 0; itrmatr<nrow; itrmatr++)
    {
        totalrow += matrows[itrmatr];
        if (itrmatr != nrow-1)
        {
            rows[itrmatr+1] = rows[itrmatr] + matrows[itrmatr];
        }
    }
    
    for (itrmatc = 0; itrmatc<ncol; itrmatc++)
    {
        totalcol += matcols[itrmatc*nrow];
        if (itrmatc != ncol-1)
        {
            cols[itrmatc+1] = cols[itrmatc] + matcols[itrmatc*nrow];
        }
    }

    /*rows,cols stores the initial index of each matrix*/
    matshape[0] = totalrow;
    matshape[1] = totalcol;
    rtmat = (ptr_cat_zmat) cat_GMatConstructor(CAT_Z, 2, matshape);
    if (rtmat == NULL)
    {
        printf("Not enough Memory!\n");
        return NULL;
    }
    
    /*copy data*/

    for (itrmatc = 0; itrmatc<ncol; itrmatc++)
    {
        for (itrmatr = 0; itrmatr<nrow; itrmatr++)
        {
            cachemat = ptrmat[itrmatc*nrow+itrmatr];
            
            if (cachemat == NULL)
            {
                continue;
            }
            offsetr = rows[itrmatr];
            offsetc = cols[itrmatc];
            matrow = cachemat->shape[0];
            matcol = cachemat->shape[1];
            //printf("\n%d %d %d %d %d %d\n",itrmatr,itrmatc,offsetr,offsetc,matrow,matcol);
            /*
             * the small matrix <cachemat> has shape [matrow, matcol]
             * located at [offsetr:offsetr+matrow-1, offsetc:offsetc+matcol-1]
             */

            for (itri = 0; itri < matcol; itri++)
            {
                for (itrj = cachemat->ptrind1[itri]; itrj < cachemat->ptrind1[itri+1]; itrj++)
                {
                    // data row index in small matrix: cachemat->ptrind2[itrj]
                    // data col index in small matrix: itri
                    // data row index in total matrix: cachemat->ptrind2[itrj] + offsetr
                    // data col index in total matrix: itri + offsetc
                    rtmat->data[(itri + offsetc)*totalrow + cachemat->ptrind2[itrj] + offsetr]
                     = cachemat->data[itrj];
                }
            }
        }
    }

    return rtmat;
}

void cat_zmat_matrixGE2CC(ptr_cat_zmat tmat)
{
    double _Complex *newdata;
    int itri,itrj,itr;
    int m,n;
    CAT_ASSERT(tmat->datatype == CAT_Z);
    m = tmat->shape[0];
    n = tmat->shape[1];
    newdata = malloc(sizeof(double _Complex)*m*n);
    memset(newdata,0,sizeof(double _Complex)*m*n);

    // for (itri = 0; itri < n; itri++) {
    //     for (itrj = tmat->ptrind1[itri]; itrj < tmat->ptrind1[itri+1]; itrj++) {
    //         newdata[itri*m+tmat->ptrind2[itrj]] = tmat->data[itrj];
    //     }
    // }
    // free(tmat->data);
    // free(tmat->ptrind1);
    // free(tmat->ptrind2);

    tmat->ptrind1 = malloc(sizeof(int)*(m+1));
    tmat->ptrind2 = malloc(sizeof(int)*m*n);
    itr = 0;
    for (itri = 0; itri < m; itri++)
    {
        tmat->ptrind1[itri] = itr;
        for (itrj = 0; itrj < n; itrj++)
        {
            tmat->ptrind2[itr] = itrj;
            itr++;
        }
    }
    tmat->ptrind1[m] = itr;
    // tmat->data = newdata;
    tmat->matrixtype = CAT_CSC;
    tmat->ldata = n*m;
}

void cat_zmat_matrixCC2GE(ptr_cat_zmat tmat)
{
    double _Complex *newdata;
    int itri,itrj;
    int m,n;
    CAT_ASSERT(tmat->datatype == CAT_Z);
    m = tmat->shape[0];
    n = tmat->shape[1];
    newdata = malloc(sizeof(double _Complex)*m*n);
    memset(newdata,0,sizeof(double _Complex)*m*n);

    for (itri = 0; itri < n; itri++) {
        for (itrj = tmat->ptrind1[itri]; itrj < tmat->ptrind1[itri+1]; itrj++) {
            newdata[itri*m+tmat->ptrind2[itrj]] = tmat->data[itrj];
        }
    }
    free(tmat->data);
    free(tmat->ptrind1);
    free(tmat->ptrind2);
    tmat->ptrind1 = NULL;
    tmat->ptrind2 = NULL;
    tmat->data = newdata;
    tmat->matrixtype = CAT_GEM;
    tmat->ldata = n*m;
}

ptr_cat_mat cat_CSCMatConstructor(cat_flag matdatatype, int row, int col, 
            int* ai, int* aj, void* av) {
    ptr_cat_mat rtmat;
    rtmat = malloc(sizeof(cat_mat));
    memset(rtmat,0,sizeof(cat_mat));
    rtmat->lshape = 2;
    rtmat->shape = malloc(sizeof(int)*2);
    rtmat->shape[0] = row;
    rtmat->shape[1] = col;
    rtmat->ldata = ai[col];
    rtmat->ptrind1 = malloc(sizeof(int)*(col+1));
    rtmat->ptrind2 = malloc(sizeof(int)*rtmat->ldata);
    rtmat->datatype = matdatatype;
    rtmat->matrixtype = CAT_CSC;
    rtmat->data = malloc(CAT_SIZEOF(matdatatype)*rtmat->ldata);
    memcpy(rtmat->ptrind1,ai,sizeof(int)*(col+1));
    memcpy(rtmat->ptrind2,aj,sizeof(int)*rtmat->ldata);
    memcpy(rtmat->data,av,CAT_SIZEOF(matdatatype)*rtmat->ldata);
    return rtmat;
}

/*sparse blas like routine: y=ax+by*/
void cat_zcscaxpby(int m, const void* alpha,
    const cat_z *x, const int *ix, const int *jx, const void* beta,
    cat_z **y, int **iy, int** jy) {
    
    int *newi,*newj;
    double _Complex *newv;
    int *io,*jo;
    double _Complex *o;
    int itri,itrj,itrk,itrl;
    double _Complex a,b;
    
    io = *iy;
    jo = *jy;
    o  = *y;
    a = *((double _Complex*)alpha);
    b = *((double _Complex*)beta);
    newi = malloc(sizeof(int)*(m+1));
    memset(newi,0,sizeof(int)*(m+1));
    newi[0] = 0;
    itri = 0;
    newj = malloc(sizeof(int)*(ix[m]+io[m]));
    memset(newj,0,sizeof(int)*(ix[m]+io[m]));
    itrj = 0;
    newv = malloc(sizeof(double _Complex)*(ix[m]+io[m]));
    memset(newv,0,sizeof(double _Complex)*(ix[m]+io[m]));
    itrk = 0;
    itrl = 0;
    for (itri = 0; itri < m; itri++) {
        while (itrk < ix[itri+1] ||
                itrl < io[itri+1]) {
            if (itrk < ix[itri+1] && itrl < io[itri+1] && jx[itrk] == jo[itrl]) {
                newj[itrj] = jx[itrk];
                newv[itrj] = a*x[itrk]+b*o[itrl];
                itrk++;
                itrl++;
                itrj++;
            }
            else if (itrl == io[itri+1] ||
                (itrk < ix[itri+1] && jx[itrk] < jo[itrl])) {
                newj[itrj] = jx[itrk];
                newv[itrj] = a*x[itrk];
                itrk++;
                itrj++;
            }
            else {
                newj[itrj] = jo[itrl];
                newv[itrj] = b*o[itrl];
                itrl++;
                itrj++;
            }
        }
        newi[itri+1] = itrj;
    }
    itrl = newi[m];
    if (itrl < (ix[m]+io[m])) {
        newj = realloc(newj,sizeof(int)*itrl);
        newv = realloc(newv,sizeof(double _Complex)*itrl);
    }
    free(*y);
    free(*iy);
    free(*jy);
    *y = newv;
    *iy = newi;
    *jy = newj;
}
/*sparse blas routine*/
void cat_zcscgemv(CBLAS_TRANSPOSE TransA, int m,
    const cat_z *a, const int *ia, const int *ja,
    const cat_z *x, cat_z* y)
{
    int itri,itrj;
    if (TransA == CblasNoTrans) {
        memset(y,0,sizeof(cat_z)*m);
        for (itri = 0; itri < m; itri++) {
            for (itrj = ia[itri]; itrj < ia[itri+1]; itrj++) {
                y[ja[itrj]] += a[itrj]*x[itri];
            }
        }
    }
    else if (TransA == CblasTrans) {
        memset(y,0,sizeof(cat_z)*m);
        for (itri = 0; itri < m; itri++) {
            for (itrj = ia[itri]; itrj < ia[itri+1]; itrj++) {
                y[itri] += a[itrj]*x[ja[itrj]];
            }
        }
    }
    else if (TransA == CblasConjTrans) {
        memset(y,0,sizeof(cat_z)*m);
        for (itri = 0; itri < m; itri++) {
            for (itrj = ia[itri]; itrj < ia[itri+1]; itrj++) {
                y[itri] += CAT_CONJ(a[itrj])*x[ja[itrj]];
            }
        }
    }
}
/*catlab routine*/
void cat_matz_csc_mldivide(ptr_cat_zmat matA,
        cat_z* ivec, cat_z* ovec, ptr_cat_workspace pworkspace)
{
    int row,col;
    
    row = matA->shape[0];
    col = matA->shape[1];
    pworkspace->umftype = CAT_Z;
    CAT_ASSERT(matA->datatype == CAT_Z && matA->matrixtype == CAT_CSC);
    if (pworkspace->isreuse != CAT_TRUE && pworkspace->Symbolic != NULL) {
        umfpack_zi_free_symbolic(&(pworkspace->Symbolic));
        pworkspace->Symbolic = NULL;
    }
    if (pworkspace->Symbolic == NULL) {
        umfpack_zi_symbolic(row, col, matA->ptrind1, matA->ptrind2,
            (double*)matA->data, NULL, &(pworkspace->Symbolic),
            pworkspace->Control, pworkspace->Info);
    }
    if (pworkspace->isreuse != CAT_TRUE && pworkspace->Numeric != NULL) {
        umfpack_zi_free_numeric(&(pworkspace->Numeric));
        pworkspace->Numeric = NULL;
    }
    if (pworkspace->Numeric == NULL) {
        umfpack_zi_numeric(matA->ptrind1, matA->ptrind2, (double*)matA->data,
            NULL, pworkspace->Symbolic, &(pworkspace->Numeric),
            pworkspace->Control, pworkspace->Info);
    }
    //Just prepare a solver, don't solve.
    if (ivec == NULL) {
        return;
    }
    umfpack_zi_solve(UMFPACK_A, matA->ptrind1, matA->ptrind2,
        (double*)matA->data, NULL, (double*)ovec, NULL,
        (double*)ivec, NULL, pworkspace->Numeric,
        pworkspace->Control, pworkspace->Info);
}

void cat_matd_csc_mldivide(ptr_cat_zmat matA,
        cat_z* ivec, cat_z* ovec, ptr_cat_workspace pworkspace)
{
    int row,col;
    
    row = matA->shape[0];
    col = matA->shape[1];
    pworkspace->umftype = CAT_D;
    CAT_ASSERT(matA->datatype == CAT_D && matA->matrixtype == CAT_CSC);
    if (pworkspace->isreuse != CAT_TRUE && pworkspace->Symbolic != NULL) {
        umfpack_di_free_symbolic(&(pworkspace->Symbolic));
        pworkspace->Symbolic = NULL;
    }
    if (pworkspace->Symbolic == NULL) {
        umfpack_di_symbolic(row, col, matA->ptrind1, matA->ptrind2,
            (double*)matA->data, &(pworkspace->Symbolic),
            pworkspace->Control, pworkspace->Info);
    }
    if (pworkspace->isreuse != CAT_TRUE && pworkspace->Numeric != NULL) {
        umfpack_di_free_numeric(&(pworkspace->Numeric));
        pworkspace->Numeric = NULL;
    }
    if (pworkspace->Numeric == NULL) {
        umfpack_di_numeric(matA->ptrind1, matA->ptrind2,
            (double*)matA->data, pworkspace->Symbolic,
            &(pworkspace->Numeric), pworkspace->Control,
            pworkspace->Info);
    }
    //Just prepare a solver, don't solve.
    if (ivec == NULL) {
        return;
    }
    umfpack_di_solve(UMFPACK_A, matA->ptrind1, matA->ptrind2,
        (double*)matA->data, (double*)ovec, (double*)ivec,
        pworkspace->Numeric, pworkspace->Control, pworkspace->Info);
}

ptr_cat_workspace cat_mat_csc_mldivide_alloc()
{
    ptr_cat_workspace rtspace;
    rtspace = malloc(sizeof(cat_workspace));
    memset(rtspace,0,sizeof(cat_workspace));
    return rtspace;
}

void cat_mat_csc_mldivide_dealloc(ptr_cat_workspace psp)
{
    if (psp->Symbolic != NULL) {
        if (psp->umftype == CAT_UMFPACK_DI) {
            umfpack_di_free_symbolic(&(psp->Symbolic));
        }
        else if (psp->umftype == CAT_UMFPACK_ZI) {
            umfpack_zi_free_symbolic(&(psp->Symbolic));
        }
    }
    if (psp->Numeric != NULL) {
        if (psp->umftype == CAT_UMFPACK_DI) {
            umfpack_di_free_numeric(&(psp->Numeric));
        }
        else if (psp->umftype == CAT_UMFPACK_ZI) {
            umfpack_zi_free_numeric(&(psp->Numeric));
        }
    }
    free(psp);
}

ptr_cat_mat cat_matrixDuplicate(ptr_cat_mat smat)
{
    ptr_cat_mat tmat;
    CAT_ASSERT(smat->bPrepared == CAT_TRUE);
    tmat                = malloc(sizeof(cat_mat));
    memset(tmat,0,sizeof(cat_mat));
    tmat->lshape        = smat->lshape;
    tmat->ldata         = smat->ldata;
    tmat->datatype      = smat->datatype;
    tmat->bPrepared     = smat->bPrepared;
    tmat->matrixtype    = smat->matrixtype;
    tmat->shape         = malloc(sizeof(int)*tmat->lshape);
    tmat->ku            = smat->ku;
    tmat->kl            = smat->kl;
    memcpy(tmat->shape,smat->shape,sizeof(int)*tmat->lshape);
    tmat->data          = malloc(CAT_SIZEOF(tmat->datatype)*tmat->ldata);
    memcpy(tmat->data,smat->data,CAT_SIZEOF(tmat->datatype)*tmat->ldata);
    if (smat->ptrind1 != NULL) {
        if (smat->matrixtype == CAT_CSC) {
            //CSC matrix short ind array
            tmat->ptrind1 = malloc(sizeof(int)*(tmat->shape[0]+1));
            memcpy(tmat->ptrind1,smat->ptrind1,sizeof(int)*(tmat->shape[0]+1));
        }
        else {
            printf("Not Supported\n");
        }
    }
    if (smat->ptrind2 != NULL) {
        if (smat->matrixtype == CAT_CSC) {
            tmat->ptrind2 = malloc(sizeof(int)*tmat->ldata);
            memcpy(tmat->ptrind2,smat->ptrind2,sizeof(int)*tmat->ldata);
        }
        else {
            printf("Not Supported\n");
        }
    }
    return tmat;
}

ptr_cat_mat cat_GMatConstructor(cat_flag matdatatype, int dimension, int* shape)
{
    ptr_cat_mat tmat;
    int itri,cachen;
    tmat = malloc(sizeof(cat_mat));
    memset(tmat,0,sizeof(cat_mat));
    tmat->matrixtype = CAT_GEM;
    tmat->datatype = matdatatype;
    tmat->lshape = dimension;
    tmat->shape = malloc(sizeof(int)*dimension);
    memcpy(tmat->shape,shape,sizeof(int)*dimension);
    cachen = 1;
    for (itri = 0; itri < dimension; itri++) {
        cachen *= shape[itri];
    }
    tmat->ldata = cachen;
    tmat->data = malloc(CAT_SIZEOF(tmat->datatype)*cachen);
    memset(tmat->data,0,CAT_SIZEOF(tmat->datatype)*cachen);
    tmat->bPrepared = CAT_TRUE;
    return tmat;
}

ptr_cat_mat cat_BMatConstructor(cat_flag matdatatype, int rol, int col, int ku, int kl)
{
    ptr_cat_mat tmat;
    int * newshape;
    tmat = malloc(sizeof(cat_mat));
    memset(tmat,0,sizeof(cat_mat));
    newshape = malloc(sizeof(int)*2);
    tmat->matrixtype = CAT_GBM;
    tmat->datatype = matdatatype;
    tmat->lshape = 2;
    tmat->shape = newshape;
    tmat->shape[0] = rol;
    tmat->shape[1] = col;
    tmat->ku = ku;
    tmat->kl = kl;
    tmat->ldata = col*(ku+kl+1);
    tmat->data = malloc(CAT_SIZEOF(matdatatype)*tmat->ldata);
    memset(tmat->data,0,CAT_SIZEOF(tmat->datatype)*tmat->ldata);
    tmat->bPrepared = CAT_TRUE;
    return tmat;
}

ptr_cat_zmat cat_matz_BMat_Eye(int n)
{
    ptr_cat_zmat tmat;
    int * newshape;
    int itr;
    tmat = malloc(sizeof(cat_mat));
    memset(tmat,0,sizeof(cat_mat));
    newshape = malloc(sizeof(int)*2);
    tmat->matrixtype = CAT_GBM;
    tmat->datatype = CAT_Z;
    tmat->lshape = 2;
    tmat->shape = newshape;
    tmat->shape[0] = n;
    tmat->shape[1] = n;
    tmat->ku = 0;
    tmat->kl = 0;
    tmat->ldata = n*(tmat->ku+tmat->kl+1);
    tmat->data = malloc(CAT_SIZEOF_Z*tmat->ldata);
    for (itr = 0; itr < tmat->ldata; itr++)
    {
        tmat->data[itr] = 1.0;
    }

    tmat->bPrepared = CAT_TRUE;
    return tmat;
}

void cat_matrixDestructor(cat_mat **pptrmat)
{
    /*TODO:
     * 100% not complete!
     * need rewrite for specific matrix type
     */
    if ((*pptrmat) != NULL) {
        return;
    }
    if ((*pptrmat)->shape != NULL) {
        free((*pptrmat)->shape);
    }
    if ((*pptrmat)->data != NULL) {
        free((*pptrmat)->data);
    }
    if ((*pptrmat)->matrixtype == CAT_CSC) {
        free((*pptrmat)->ptrind1);
        free((*pptrmat)->ptrind2);
    }
    free(*pptrmat);
    *pptrmat = NULL;    
}

void cat_matrixD2Z(ptr_cat_dmat tmat) 
{
    double _Complex *ptrzdata;
    int itr;
    ptrzdata = malloc(sizeof(double _Complex)*tmat->ldata);
    for (itr = 0; itr < tmat->ldata; itr++) {
        ptrzdata[itr] = tmat->data[itr];
    }
    free(tmat->data);
    tmat->data = (double*)ptrzdata;
    tmat->datatype = CAT_Z;
}

void cat_matrixZ2DReal(ptr_cat_zmat tmat)
{
    double *ptrzdata;
    int itr;
    ptrzdata = malloc(sizeof(double)*tmat->ldata);
    for (itr = 0; itr < tmat->ldata; itr++) {
        ptrzdata[itr] = CAT_REAL(tmat->data[itr]);
    }
    free(tmat->data);
    tmat->data = (cat_z*)ptrzdata;
    tmat->datatype = CAT_D;
}

void cat_matrixZ2DImag(ptr_cat_zmat tmat)
{
    double *ptrzdata;
    int itr;
    ptrzdata = malloc(sizeof(double)*tmat->ldata);
    for (itr = 0; itr < tmat->ldata; itr++) {
        ptrzdata[itr] = CAT_IMAG(tmat->data[itr]);
    }
    free(tmat->data);
    tmat->data = (cat_z*)ptrzdata;
    tmat->datatype = CAT_D;
}

void cat_zmat_matrixGB2CC(ptr_cat_zmat tmat)
{
    //double _Complex *newdata;
    int *newind, *newicc;
    int col,row,ku,kl,datacol;
    int itri,itrj,itrdata;
    CAT_ASSERT(tmat->lshape == 2 && tmat->datatype == CAT_Z && CAT_CHECKFLAG(tmat->matrixtype,CAT_GBM));
    row = tmat->shape[0];
    col = tmat->shape[1];
    ku = tmat->ku;
    kl = tmat->kl;
    datacol = ku+kl+1;
    newind = malloc(sizeof(int)*(col+1));
    newicc = malloc(sizeof(int)*tmat->ldata);
    itrdata = 0;
    newind[0] = 0;
    for (itrj = 0; itrj < col; itrj++) {
        for (itri = CAT_MAX(0,ku-itrj); itri < CAT_MIN(datacol,row+ku-itrj); itri++) {
            tmat->data[itrdata] = tmat->data[itri+itrj*datacol];
            newicc[itrdata] = itrj + itri - ku;
            itrdata++;
        }
        newind[itrj+1] = itrdata;
    }

    //free(tmat->data); tmat->data = newdata;
    tmat->ptrind1 = newind;
    tmat->ptrind2 = newicc;
    tmat->matrixtype = CAT_CSC;
    //tmat->ldata = newind[col];
}

void cat_zmat_matrixGB2GE(ptr_cat_zmat tmat)
{
    double _Complex *newdata;
    int itri,itrj;
    int ku,kl;
    int m,n;
    CAT_ASSERT(tmat->datatype == CAT_Z);
    m = tmat->shape[0];
    n = tmat->shape[1];
    newdata = malloc(sizeof(double _Complex)*m*n);
    memset(newdata,0,sizeof(double _Complex)*m*n);
    ku = tmat->ku;
    kl = tmat->kl;
    for (itrj = 0; itrj < n; itrj++) {
        for (itri = CAT_MAX(0,itrj-ku); itri < CAT_MIN(m,itrj+kl+1); itri++) {
            newdata[itri+itrj*m] = tmat->data[itri-itrj+ku+itrj*(ku+kl+1)];
        }
    }
    free(tmat->data);
    tmat->data = newdata;
    tmat->matrixtype = CAT_GEM;
    tmat->ldata = n*m;
}

void cat_dmat_matrixGB2GE(ptr_cat_dmat tmat)
{
    double *newdata;
    int itri,itrj;
    int ku,kl;
    int m,n;
    CAT_ASSERT(tmat->datatype == CAT_D);
    m = tmat->shape[0];
    n = tmat->shape[1];
    newdata = malloc(sizeof(double)*m*n);
    memset(newdata,0,sizeof(double)*m*n);
    ku = tmat->ku;
    kl = tmat->kl;
    for (itrj = 0; itrj < n; itrj++) {
        for (itri = CAT_MAX(0,itrj-ku); itri < CAT_MIN(m,itrj+kl+1); itri++) {
            newdata[itri+itrj*m] = tmat->data[itri-itrj+ku+itrj*(ku+kl+1)];
        }
    }
    free(tmat->data);
    tmat->data = newdata;
    tmat->matrixtype = CAT_GEM;
    tmat->ldata = n*m;
}

/***********************************************************************/

ptr_cat_dmat cat_GmatRealZ(ptr_cat_zmat tmat)
{
    ptr_cat_dmat rtmat;
    int itr;
    if (CAT_CHECKFLAG(tmat->matrixtype,CAT_GEM)) {
        rtmat = (ptr_cat_dmat)cat_GMatConstructor(CAT_D, tmat->lshape, tmat->shape);
        for (itr = 0; itr < tmat->ldata; itr++) {
            rtmat->data[itr] = CAT_REAL(tmat->data[itr]);
        }
        return rtmat;
    }
    return NULL;
}

ptr_cat_dmat cat_GmatImagZ(ptr_cat_zmat tmat)
{
    ptr_cat_dmat rtmat;
    int itr;
    if (CAT_CHECKFLAG(tmat->matrixtype,CAT_GEM)) {
        rtmat = (ptr_cat_dmat)cat_GMatConstructor(CAT_D, tmat->lshape, tmat->shape);
        for (itr = 0; itr < tmat->ldata; itr++) {
            rtmat->data[itr] = CAT_IMAG(tmat->data[itr]);
        }
        return rtmat;
    }
    return NULL;
}

void cat_zmat_GMatMultVec(ptr_cat_zmat pmat, cat_z* ivec, cat_z* ovec) 
{
    static double _Complex one = 1.0;
    static double _Complex zero= 0.0;
    CAT_ASSERT(pmat->bPrepared==CAT_TRUE && pmat->datatype==CAT_Z && CAT_CHECKFLAG(mata->matrixtype,CAT_GEM));
    cblas_zgemv(CblasColMajor, CblasNoTrans, pmat->shape[0], pmat->shape[1],
        &one, pmat->data, pmat->shape[0], ivec, 1, &zero, ovec, 1);
}

void cat_dmat_GMatMultVec(ptr_cat_dmat pmat, cat_d* ivec, cat_d* ovec) 
{
    CAT_ASSERT(pmat->bPrepared==CAT_TRUE && pmat->datatype==CAT_D && CAT_CHECKFLAG(mata->matrixtype,CAT_GEM));
    cblas_dgemv(CblasColMajor, CblasNoTrans, pmat->shape[0], pmat->shape[1],
        1.0, pmat->data, pmat->shape[0], ivec, 1, 0.0, ovec, 1);
}

void cat_mat_GMatMultVec(ptr_cat_mat pmat, void* ivec, void* ovec) 
{
    if (pmat->datatype == CAT_D) {
        cat_dmat_GMatMultVec((ptr_cat_dmat)pmat,ivec,ovec);
    }
    else if (pmat->datatype == CAT_Z) {
        cat_zmat_GMatMultVec((ptr_cat_zmat)pmat,ivec,ovec);
    }
    else {
        printf("Not supported\n");
    }
}

void cat_zmat_BMatMultVec(ptr_cat_zmat pmat, cat_z* ivec, cat_z* ovec)
{
    static double _Complex one = 1.0;
    static double _Complex zero= 0.0;
    CAT_ASSERT(pmat->bPrepared==CAT_TRUE && pmat->datatype==CAT_Z && CAT_CHECKFLAG(mata->matrixtype,CAT_GBM));
    cblas_zgbmv(CblasColMajor, CblasNoTrans, pmat->shape[0], pmat->shape[1],
        pmat->kl, pmat->ku, &one, pmat->data, pmat->kl+pmat->ku+1, ivec, 1, &zero, ovec, 1);
}

void cat_dmat_BMatMultVec(ptr_cat_dmat pmat, double* ivec, double* ovec)
{
    CAT_ASSERT(pmat->bPrepared==CAT_TRUE && pmat->datatype==CAT_D && CAT_CHECKFLAG(mata->matrixtype,CAT_GBM));
    cblas_dgbmv(CblasColMajor, CblasNoTrans, pmat->shape[0], pmat->shape[1],
        pmat->kl, pmat->ku, 1.0, pmat->data, pmat->kl+pmat->ku+1, ivec, 1, 0.0, ovec, 1);
}

void cat_mat_BMatMultVec(ptr_cat_mat pmat, void* ivec, void* ovec)
{
    if (pmat->datatype == CAT_D) {
        cat_dmat_BMatMultVec((ptr_cat_dmat)pmat,ivec,ovec);
    }
    else if (pmat->datatype == CAT_Z) {
        cat_zmat_BMatMultVec((ptr_cat_zmat)pmat,ivec,ovec);
    }
    else {
        printf("Not supported\n");
    }
}

void cat_mat_MatMultVec(ptr_cat_mat pmat, void* ivec, void* ovec)
{
    if (pmat->matrixtype == CAT_GEM) {
        cat_mat_GMatMultVec(pmat,ivec,ovec);
    }
    else if (pmat->matrixtype == CAT_GBM) {
        cat_mat_BMatMultVec(pmat,ivec,ovec);
    }
    else {
        printf("Not supported\n");
    }
}

/***********************************************************************/

void cat_matrix_fprint_s_1d(FILE *fp, ptr_cat_smat smat)
{
    int itri;
    for (itri = 0; itri < smat->shape[0]; itri++) {
        fprintf(fp, "%15.9e\n",smat->data[itri]);
    }
}

void cat_matrix_fprint_d_1d(FILE *fp, ptr_cat_dmat smat)
{
    int itri;
    for (itri = 0; itri < smat->shape[0]; itri++) {
        fprintf(fp, "%15.9e\n",smat->data[itri]);
    }
}

void cat_matrix_fprint_c_1d(FILE *fp, ptr_cat_cmat smat)
{
    int itri;
    for (itri = 0; itri < smat->shape[0]; itri++) {
        fprintf(fp, " %15.9e %15.9e\n",CAT_REAL(smat->data[itri]),CAT_IMAG(smat->data[itri]));
    }
}

void cat_matrix_fprint_z_1d(FILE *fp, ptr_cat_zmat smat)
{
    int itri;
    for (itri = 0; itri < smat->shape[0]; itri++) {
        fprintf(fp, " %15.9e %15.9e\n",CAT_REAL(smat->data[itri]),CAT_IMAG(smat->data[itri]));
    }
}

void cat_matrix_fprint_s_2d(FILE *fp, ptr_cat_smat smat)
{
    int itri,itrj;
    int nrow,ncol;
    nrow = smat->shape[0];
    ncol = smat->shape[1];
    for (itri = 0; itri < nrow; itri++) {
        for (itrj = 0; itrj < ncol; itrj++) {
            fprintf(fp, "%15.9e ",smat->data[itri+itrj*nrow]);
        }
        fprintf(fp, "\n");
    }
}

void cat_matrix_fprint_d_2d(FILE *fp, ptr_cat_dmat smat)
{
    int itri,itrj;
    int nrow,ncol;
    nrow = smat->shape[0];
    ncol = smat->shape[1];
    for (itri = 0; itri < nrow; itri++) {
        for (itrj = 0; itrj < ncol; itrj++) {
            fprintf(fp, "%15.9e ",smat->data[itri+itrj*nrow]);
        }
        fprintf(fp, "\n");
    }
}

void cat_matrix_fprint_c_2d(FILE *fp, ptr_cat_cmat smat)
{
    int itri,itrj;
    int nrow,ncol;
    nrow = smat->shape[0];
    ncol = smat->shape[1];
    for (itri = 0; itri < nrow; itri++) {
        for (itrj = 0; itrj < ncol; itrj++) {
            fprintf(fp, " %15.9e %15.9e",CAT_REAL(smat->data[itri+itrj*nrow]),CAT_IMAG(smat->data[itri+itrj*nrow]));
        }
        fprintf(fp, "\n");
    }
}

void cat_matrix_fprint_z_2d(FILE *fp, ptr_cat_zmat smat)
{
    int itri,itrj;
    int nrow,ncol;
    nrow = smat->shape[0];
    ncol = smat->shape[1];
    for (itri = 0; itri < nrow; itri++) {
        for (itrj = 0; itrj < ncol; itrj++) {
            //printf("(%15.9e,%15.9e) ",CAT_REAL(smat->data[itri+itrj*nrow]),CAT_IMAG(smat->data[itri+itrj*nrow]));
            fprintf(fp, " %15.9e %15.9e",CAT_REAL(smat->data[itri+itrj*nrow]),CAT_IMAG(smat->data[itri+itrj*nrow]));
        }
        fprintf(fp, "\n");
    }
}

void cat_matrix_fprint(FILE *fp, ptr_cat_mat smat)
{
    CAT_ASSERT(smat->bPrepared == CAT_TRUE && smat->lshape > 0);
    if (smat->lshape == 1) {
        if (smat->datatype == CAT_S) {
            cat_matrix_fprint_s_1d(fp, (ptr_cat_smat)smat);
        }
        else if (smat->datatype == CAT_D) {
            cat_matrix_fprint_d_1d(fp, (ptr_cat_dmat)smat);
        }
        else if (smat->datatype == CAT_C) {
            cat_matrix_fprint_c_1d(fp, (ptr_cat_cmat)smat);
        }
        else if (smat->datatype == CAT_Z) {
            cat_matrix_fprint_z_1d(fp, (ptr_cat_zmat)smat);
        }
    }
    else if (smat->lshape == 2) {
        if (smat->datatype == CAT_S) {
            cat_matrix_fprint_s_2d(fp, (ptr_cat_smat)smat);
        }
        else if (smat->datatype == CAT_D) {
            cat_matrix_fprint_d_2d(fp, (ptr_cat_dmat)smat);
        }
        else if (smat->datatype == CAT_C) {
            cat_matrix_fprint_c_2d(fp, (ptr_cat_cmat)smat);
        }
        else if (smat->datatype == CAT_Z) {
            cat_matrix_fprint_z_2d(fp, (ptr_cat_zmat)smat);
        }
    }
    else {
        fprintf(fp, "Not supported, larger than 2 dimensions\n");
    }
}


void cat_matrix_fwrite_z_2d(FILE *fp, ptr_cat_zmat smat)
{
    fwrite(smat->data, sizeof(double _Complex), smat->shape[0]*smat->shape[1], fp);
}

void cat_matrix_fwrite_z_1d(FILE *fp, ptr_cat_zmat smat)
{
    fwrite(smat->data, sizeof(double _Complex), smat->shape[0], fp);
}

void cat_matrix_filenamewrite_z_2d(char* filename, ptr_cat_zmat smat)
{
    FILE *fp = fopen(filename,"wb");
    cat_matrix_fwrite_z_2d(fp,smat);
    fclose(fp);
}

void cat_matrix_filenamewrite_z_1d(char* filename, ptr_cat_zmat smat)
{
    FILE *fp = fopen(filename,"wb");
    cat_matrix_fwrite_z_1d(fp,smat);
    fclose(fp);
}

