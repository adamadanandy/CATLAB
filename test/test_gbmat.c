#include <catlab.h>

int main(int argc, char** argv)
{
    ptr_cat_zmat pmat,qmat,rmat,smat;
    int itri;
    double _Complex ivec[10], ovec[10];
    pmat = (ptr_cat_zmat)cat_BMatConstructor(CAT_Z, 10, 10, 2, 3);
    smat = (ptr_cat_zmat)cat_BMatConstructor(CAT_Z, 10, 10, 1, 4);
    for (itri = 0; itri < 10; itri++) {
        pmat->data[itri*6+0] = itri;
        pmat->data[itri*6+1] = itri*I;
        pmat->data[itri*6+2] = 1.0;
        pmat->data[itri*6+3] = 9-itri+1*I;
        pmat->data[itri*6+4] = 1+(9-itri)*I;
        pmat->data[itri*6+5] = itri+itri*I;
        smat->data[itri*6+0] = itri;
        smat->data[itri*6+1] = itri*I;
        smat->data[itri*6+2] = 1.0;
        smat->data[itri*6+3] = 9-itri+1*I;
        smat->data[itri*6+4] = 1+(9-itri)*I;
        smat->data[itri*6+5] = itri+itri*I;
        ivec[itri] = itri*I;
        ovec[itri] = 0.0;
    }
    qmat = (ptr_cat_zmat)cat_matrixDuplicate((ptr_cat_mat)pmat);
    rmat = (ptr_cat_zmat)cat_matrixDuplicate((ptr_cat_mat)smat);
    printf("-------------Test Mat Mult Vec------------\n");
    CAT_VECZ_PRINT(ivec,itri,10);
    printf("------------------------------------------\n");
    cat_mat_MatMultVec((ptr_cat_mat)pmat, ivec, ovec);
    CAT_VECZ_PRINT(ovec,itri,10);
    printf("------------------------------------------\n");
    cat_zmat_matrixGB2GE(pmat);
    CAT_MATRIX_PRINT(pmat);
    printf("--------Check With GE Mat Mult Vec--------\n");
    for (itri = 0; itri < 10; itri++) ovec[itri] = 0.0;
    cat_mat_MatMultVec((ptr_cat_mat)pmat, ivec, ovec);
    CAT_VECZ_PRINT(ovec,itri,10);
    printf("------------Check GB Mat AXPBY------------\n");
    CAT_MATRIX_PRINT(pmat);
    printf("------------------------------------------\n");
    cat_zmat_matrixGB2GE(smat);
    CAT_MATRIX_PRINT(smat);
    printf("------------------------------------------\n");
    cat_zmat_gbm_AXPBY(1.0,qmat,1.0,rmat);
    cat_zmat_matrixGB2GE(rmat);
    CAT_MATRIX_PRINT(rmat);
}