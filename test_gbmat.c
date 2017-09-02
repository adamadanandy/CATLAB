#include <catlab.h>

int main(int argc, char** argv)
{
    ptr_cat_zmat mymat;
    int itri,itrj;
    double _Complex ivec[10], ovec[10];
    mymat = (ptr_cat_zmat)cat_BMatConstructor(CAT_Z, 10, 10, 2, 3);
    for (itri = 0; itri < 10; itri++) {
        mymat->data[itri*6+0] = itri;
        mymat->data[itri*6+1] = itri*I;
        mymat->data[itri*6+2] = 1.0;
        mymat->data[itri*6+3] = 9-itri+1*I;
        mymat->data[itri*6+4] = 1+(9-itri)*I;
        mymat->data[itri*6+5] = itri+itri*I;
        ivec[itri] = itri*I;
        ovec[itri] = 0.0;
    }
    CAT_VECZ_PRINT(ivec,itri,10);
    printf("--------------------------------\n");
    cat_mat_MatMultVec((ptr_cat_mat)mymat, ivec, ovec);
    CAT_VECZ_PRINT(ovec,itri,10);
    printf("--------------------------------\n");
    cat_zmat_matrixGB2GE(mymat);
    CAT_MATRIX_PRINT(mymat);
    printf("--------------------------------\n");
    for (itri = 0; itri < 10; itri++) ovec[itri] = 0.0;
    cat_mat_MatMultVec((ptr_cat_mat)mymat, ivec, ovec);
    CAT_VECZ_PRINT(ovec,itri,10);
}