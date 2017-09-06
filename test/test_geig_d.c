#define CAT_DEBUG
#include <stdio.h>
#include <catlab.h>

int main(int argc, char** argv)
{
    ptr_cat_dmat pmat, qmat;
    int matshape[2] = {10,10};
    int itri,itrj;
    ptr_cat_zmat dmat, vmat;
    dmat = malloc(sizeof(cat_zmat));
    vmat = malloc(sizeof(cat_dmat));
    pmat = (ptr_cat_dmat)cat_GMatConstructor(CAT_D,2,matshape);
    qmat = (ptr_cat_dmat)cat_GMatConstructor(CAT_D,2,matshape);
    for (itri = 0; itri < 10; itri++) {
        for (itrj = 0; itrj < 10; itrj++) {
            if (itri == itrj) {
                pmat->data[itri+itrj*10] = 2.0;
                qmat->data[itri+itrj*10] = -1.0;
            }
            else if (itri == itrj-1) {
                pmat->data[itri+itrj*10] = 1.0;
                qmat->data[itri+itrj*10] = 1.0;
            }
            else if (itri == itrj+1) {
                pmat->data[itri+itrj*10] = -itri;
                qmat->data[itri+itrj*10] = 1.0;
            }
            else {
                pmat->data[itri+itrj*10] = 0.0;
                qmat->data[itri+itrj*10] = 0.0;
            }
        }
    }
    cat_dmat_geig(pmat, qmat, dmat, vmat);
    CAT_MATRIX_PRINT(pmat);
    CAT_LOG(------------------------------------)
    CAT_MATRIX_PRINT(qmat);
    CAT_LOG(------------------------------------)
    CAT_MATRIX_PRINT(dmat);
    CAT_LOG(------------------------------------)
    CAT_MATRIX_PRINT(vmat);
}
