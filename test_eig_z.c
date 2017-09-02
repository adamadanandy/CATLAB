#define CAT_DEBUG
#include <stdio.h>
#include <catlab.h>

int main(int argc, char** argv)
{
    ptr_cat_zmat pmat;
    int matshape[2] = {10,10};
    int itri,itrj;
    ptr_cat_zmat dmat,vmat;
    dmat = malloc(sizeof(cat_zmat));
    vmat = malloc(sizeof(cat_zmat));

    pmat = (ptr_cat_zmat)cat_GMatConstructor(CAT_Z,2, matshape);
    for (itri = 0; itri < 10; itri++) {
        for (itrj = 0; itrj < 10; itrj++) {
            if (itri == itrj) {
                pmat->data[itri+itrj*10] = 1.0;
            }
            else if (itri == itrj-1) {
                pmat->data[itri+itrj*10] = itri + I;
            }
            else if (itri == itrj+1) {
                pmat->data[itri+itrj*10] = itri - 2.0*I;
            }
        }
    }
    cat_zmat_eig(pmat,dmat,vmat);
    CAT_MATRIX_PRINT(pmat);
    CAT_MATRIX_PRINT(dmat);
    CAT_MATRIX_PRINT(vmat);
}
