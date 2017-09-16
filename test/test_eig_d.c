#define CAT_DEBUG
#include <stdio.h>
#include <catlab.h>

int main(int argc, char** argv)
{
    ptr_cat_dmat pmat;
    int matshape[2] = {10,10};
    int itri,itrj;
    ptr_cat_zmat dmat,vmat;
    dmat = (ptr_cat_zmat)cat_EmptyMatConstructor();
    vmat = (ptr_cat_zmat)cat_EmptyMatConstructor();
    pmat = (ptr_cat_dmat)cat_GMatConstructor(CAT_D,2,matshape);
    for (itri = 0; itri < 10; itri++) {
        for (itrj = 0; itrj < 10; itrj++) {
            if (itri == itrj) {
                pmat->data[itri+itrj*10] = 2.0;
            }
            else if (itri == itrj-1) {
                pmat->data[itri+itrj*10] = 1.0;
            }
            else if (itri == itrj+1) {
                pmat->data[itri+itrj*10] = -itri;
            }
            else {
                pmat->data[itri+itrj*10] = 0.0;
            }
        }
    }
    cat_dmat_eig(pmat,dmat,vmat);
    CAT_MATRIX_PRINT(pmat);
    CAT_MATRIX_PRINT(dmat);
    CAT_MATRIX_PRINT(vmat);
}
