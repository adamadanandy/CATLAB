#include <stdio.h>
#include <catlab.h>
#include <complex.h>

int main(int argc, char** argv)
{
    ptr_cat_zmat gbmat, gemat, ccmat;
    cat_zmat dmat, vmat;
    int itri,itrj;
    
    double _Complex ivec[10], ovec[10];
    ptr_cat_workspace pcw = cat_matz_csc_mldivide_alloc();
    gbmat = (ptr_cat_zmat)cat_BMatConstructor(CAT_Z, 10, 10, 2, 1);
    for (itri = 0; itri < 10; itri++) {
        gbmat->data[itri*4+0] = itri;
        gbmat->data[itri*4+1] = itri*I;
        gbmat->data[itri*4+2] = 1.0;
        gbmat->data[itri*4+3] = 9-itri+1*I;
        ivec[itri] = itri*(1.0+I);
    }
    gemat = (ptr_cat_zmat)cat_matrixDuplicate((ptr_cat_mat)gbmat);
    ccmat = (ptr_cat_zmat)cat_matrixDuplicate((ptr_cat_mat)gbmat);
    cat_zmat_matrixGB2GE(gemat);
    cat_zmat_matrixGB2CC(ccmat);
    printf("----------------------------------------------------------\n");
    CAT_MATRIX_PRINT(gemat);
    printf("----------------------------------------------------------\n");
    cat_zmat_eigs(gbmat, &dmat, &vmat, 1.0+1.0*I, 3, NULL);
    CAT_MATRIX_PRINT(&dmat);
    printf("----------------------------------------------------------\n");
}