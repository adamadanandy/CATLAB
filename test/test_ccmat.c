#include <catlab.h>
#include <time.h>

int main(int argc, char** argv)
{
    ptr_cat_zmat pmat, qmat, rmat, smat, tmat, omat;
    int itri,itrj;
    clock_t t1,t2,t3,t4;
    double _Complex ivec[10], ovec[10];
    int ai[11],aj[10];
    double _Complex av[10];
    double _Complex alpha = 1.0+1.0*I;
    double _Complex beta = 1.0-1.0*I;
    ptr_cat_workspace pcw = cat_mat_csc_mldivide_alloc();
    pmat = (ptr_cat_zmat)cat_BMatConstructor(CAT_Z, 10, 10, 2, 1);
    
    for (itri = 0; itri < 10; itri++) {
        pmat->data[itri*4+0] = itri;
        pmat->data[itri*4+1] = itri*I;
        pmat->data[itri*4+2] = 1.0;
        pmat->data[itri*4+3] = 9-itri+1*I;
        ivec[itri] = itri*(1.0+I);
        ai[itri] = itri;
        aj[itri] = 9-itri;
        av[itri] = itri*(1.0+I);
    }
    ai[10] = 10;
    smat = (ptr_cat_zmat)cat_CSCMatConstructor(CAT_Z, 10, 10, ai, aj, av);
    tmat = (ptr_cat_zmat)cat_matrixDuplicate((ptr_cat_mat)smat);
    omat = (ptr_cat_zmat)cat_matrixDuplicate((ptr_cat_mat)smat);
    cat_zmat_matrixCC2GE(tmat);
    qmat = (ptr_cat_zmat)cat_matrixDuplicate((ptr_cat_mat)pmat);
    rmat = (ptr_cat_zmat)cat_matrixDuplicate((ptr_cat_mat)pmat);
    cat_zmat_matrixGB2GE(qmat);
    cat_zmat_matrixGB2CC(rmat);
    printf("-------------------Test convert GBM2GEM-------------------\n");
    CAT_MATRIX_PRINT(qmat);
    printf("-------------------Test mldivide reused-------------------\n");
    t1 = clock();
    for (itri = 0; itri < 10; itri++) {
        ivec[itri] = itri*(1.0+I);
    }
    for (itrj = 0; itrj < 500000; itrj++) {
        cat_matz_csc_mldivide(rmat, ivec, ovec, pcw);
        for (itri = 0; itri < 10; itri++) {
            ivec[itri] = ovec[itri]*1.62481;
        }
    }
    t2 = clock();
    printf("Not Reused Elapse: %g\n",(double)(t2-t1)/CLOCKS_PER_SEC);
    t3 = clock();
    pcw->isreuse = CAT_TRUE;
    for (itri = 0; itri < 10; itri++) {
        ivec[itri] = itri*(1.0+I);
    }
    for (itrj = 0; itrj < 500000; itrj++) {
        cat_matz_csc_mldivide(rmat, ivec, ovec, pcw);
        for (itri = 0; itri < 10; itri++) {
            ivec[itri] = ovec[itri]*1.62481;
        }
    }
    t4 = clock();
    printf("With Reused Elapse: %g\n",(double)(t4-t3)/CLOCKS_PER_SEC);
    for (itri = 0; itri < 10; itri++) {
        ivec[itri] = itri*(1.0+I);
    }
    cat_zcscgemv(CblasNoTrans, rmat->shape[0], rmat->data,
        rmat->ptrind1, rmat->ptrind2, ivec, ovec);
    printf("---------------------Test multiply------------------------\n");
    for (itri = 0; itri < 10; itri++) {
        printf("(%g,%g)\n",CAT_REAL(ivec[itri]),CAT_IMAG(ivec[itri]));
    }
    printf("----------------------------------------------------------\n");
    for (itri = 0; itri < 10; itri++) {
        printf("(%g,%g)\n",CAT_REAL(ovec[itri]),CAT_IMAG(ovec[itri]));
    }
    printf("------------------------Test axpby------------------------\n");
    CAT_MATRIX_PRINT(qmat);
    printf("----------------------------------------------------------\n");
    CAT_MATRIX_PRINT(tmat);
    printf("----------------------------------------------------------\n");
//    cat_zcscaxpby(10, &alpha, rmat->data, rmat->ptrind1, rmat->ptrind2,
//        &beta, &(omat->data), &(omat->ptrind1), &(omat->ptrind2));
//    omat->ldata = omat->ptrind1[10];
    cat_zmat_csc_AXPBY(alpha, rmat, beta, omat);
    cat_zmat_matrixCC2GE(omat);
    CAT_MATRIX_PRINT(omat);
}