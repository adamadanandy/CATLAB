#include <stdio.h>
#include <catlab.h>

int main(int argc, char** argv)
{
    ptr_cat_zmat testzmat;
    ptr_cat_dmat testdmat;
    double _Complex izvec[10],ozvec[10];
    double idvec[10],odvec[10];
    int shape10x10[2] = {10,10};
    int itri,itrj;
    
    testzmat = (ptr_cat_zmat)cat_GMatConstructor(CAT_Z, 2, shape10x10);
    
    for (itri = 0; itri < 10; itri++) {
        for (itrj = 0; itrj < 10; itrj++) {
            if (itri - itrj == 1 || itri-itrj == -1) {
                testzmat->data[itri+itrj*10] = itri+itrj*I;
            }
            else if (itri == itrj) {
                testzmat->data[itri+itrj*10] = 1.0;
            }
            else {
                testzmat->data[itri+itrj*10] = 0.0;
            }
        }
        izvec[itri] = 1.0+itri*I;
        idvec[itri] = itri;
    }
    CAT_MATRIX_PRINT(testzmat);
    printf("---------------------------------\n");
    for (itri = 0; itri < 10; itri++) {
        printf("(%g,%g)\n",CAT_REAL(izvec[itri]),CAT_IMAG(izvec[itri]));
    }
    printf("---------------------------------\n");
    cat_zmat_GMatMultVec(testzmat, izvec, ozvec);
    for (itri = 0; itri < 10; itri++) {
        printf("(%g,%g)\n",CAT_REAL(ozvec[itri]),CAT_IMAG(ozvec[itri]));
    }
    printf("---------------------------------\n");
    testdmat = cat_GmatRealZ(testzmat);
    CAT_MATRIX_PRINT(testdmat);
    printf("---------------------------------\n");
    for (itri = 0; itri < 10; itri++) {
        printf("%g\n",idvec[itri]);
    }
    printf("---------------------------------\n");
    cat_dmat_GMatMultVec(testdmat,idvec,odvec);
    for (itri = 0; itri < 10; itri++) {
        printf("%g\n",odvec[itri]);
    }
}