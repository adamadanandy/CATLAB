gcc -I/home/yyli/staticlib/cblas/include/ -I/home/yyli/staticlib/lapacke/include/ -I./ -c catlab_eig.c
gcc test_eig.o catlab_eig.o /home/yyli/staticlib/lapacke/liblapacke.a /home/yyli/staticlib/lapack/liblapack.a /home/yyli/staticlib/cblas/libcblas.a /home/yyli/staticlib/blas/librefblas.a -lgfortran -lm -o test_eig
