#include <catlab.h>
#include <time.h>

int main(int argc, char** argv)
{
    ptr_cat_zmat pmat, qmat, rmat;
    int itri,itrj;
    clock_t t1,t2,t3,t4;
    double _Complex ivec[10], ovec[10];
    ptr_cat_workspace pcw = cat_matz_csc_mldivide_alloc();
    pmat = (ptr_cat_zmat)cat_BMatConstructor(CAT_Z, 10, 10, 2, 1);
    for (itri = 0; itri < 10; itri++) {
        pmat->data[itri*4+0] = itri;
        pmat->data[itri*4+1] = itri*I;
        pmat->data[itri*4+2] = 1.0;
        pmat->data[itri*4+3] = 9-itri+1*I;
        ivec[itri] = itri*(1.0+I);
    }
    qmat = (ptr_cat_zmat)cat_matrixDuplicate((ptr_cat_mat)pmat);
    rmat = (ptr_cat_zmat)cat_matrixDuplicate((ptr_cat_mat)pmat);
    cat_zmat_matrixGB2GE(qmat);
    cat_zmat_matrixGB2CC(rmat);
    printf("----------------------------------------------------------\n");
    CAT_MATRIX_PRINT(qmat);
    printf("----------------------------------------------------------\n");
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
}