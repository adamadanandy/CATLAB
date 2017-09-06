#include <stdio.h>
#include <catlab.h>
#include <complex.h>

int main(int argc, char** argv)
{
    ptr_cat_zmat gbmat1, gemat1, ccmat1;
    ptr_cat_zmat gbmat2, gemat2, ccmat2;
    cat_zmat dmat, vmat;
    int itri;
    
    gbmat1 = (ptr_cat_zmat)cat_BMatConstructor(CAT_Z, 10, 10, 2, 1);
    gbmat2 = (ptr_cat_zmat)cat_BMatConstructor(CAT_Z, 10, 10, 2, 1);
    for (itri = 0; itri < 10; itri++) {
        gbmat1->data[itri*4+0] = itri;
        gbmat1->data[itri*4+1] = itri*I;
        gbmat1->data[itri*4+2] = 1.0;
        gbmat1->data[itri*4+3] = 9-itri+1*I;
        
        gbmat2->data[itri*4+0] = 11-itri;
        gbmat2->data[itri*4+1] = (10-itri)*I;
        gbmat2->data[itri*4+2] = -1.0;
        gbmat2->data[itri*4+3] = itri+1+1*I;
    }
    gemat1 = (ptr_cat_zmat)cat_matrixDuplicate((ptr_cat_mat)gbmat1);
    ccmat1 = (ptr_cat_zmat)cat_matrixDuplicate((ptr_cat_mat)gbmat1);
    cat_zmat_matrixGB2GE(gemat1);
    cat_zmat_matrixGB2CC(ccmat1);
    gemat2 = (ptr_cat_zmat)cat_matrixDuplicate((ptr_cat_mat)gbmat2);
    ccmat2 = (ptr_cat_zmat)cat_matrixDuplicate((ptr_cat_mat)gbmat2);
    cat_zmat_matrixGB2GE(gemat2);
    cat_zmat_matrixGB2CC(ccmat2);
    printf("----------------------------------------------------------\n");
    CAT_MATRIX_PRINT(gemat1);
    printf("----------------------------------------------------------\n");
    CAT_MATRIX_PRINT(gemat2);
    printf("----------------------------------------------------------\n");
    cat_zmat_geigs(gbmat1, gbmat2, &dmat, &vmat, 1.0+1.0*I, 3, NULL);
    CAT_MATRIX_PRINT(&dmat);
    printf("----------------------------------------------------------\n");
}