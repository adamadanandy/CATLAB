default: catlab_eig.c catlab_matrix.c
	gcc -c catlab_matrix.c -Wall -I./ -I/home/yyli/include/
	gcc -c catlab_eig.c -Wall -I./ -I/home/yyli/include/
	ar rv libcatlab.a *.o
	ranlib libcatlab.a

test: test_eig_d test_eig_z test_geig_d test_geig_z test_matmul test_gbmat test_ccmat test_eigs_z

test_eig_d: test_eig_d.c
	gcc -I. -I/home/yyli/include -c test_eig_d.c
	gcc test_eig_d.o -L./ -lcatlab -L/home/yyli/lib -lopenblas -llapack -lumfpack -larpack -lgfortran -lm -pthread -o test_eig_d

test_eig_z: test_eig_z.c
	gcc -I. -I/home/yyli/include -c test_eig_z.c
	gcc test_eig_z.o -L./ -lcatlab -L/home/yyli/lib -lopenblas -llapack -lumfpack -larpack -lgfortran -lm -pthread -o test_eig_z

test_geig_z: test_geig_z.c
	gcc -I. -I/home/yyli/include -c test_geig_z.c
	gcc test_geig_z.o -L./ -lcatlab -L/home/yyli/lib -lopenblas -llapack -lumfpack -larpack -lgfortran -lm -pthread -o test_geig_z

test_geig_d: test_geig_d.c
	gcc -I. -I/home/yyli/include -c test_geig_d.c
	gcc test_geig_d.o -L./ -lcatlab -L/home/yyli/lib -lopenblas -llapack -lumfpack -larpack -lgfortran -lm -pthread -o test_geig_d

test_matmul: test_matmul.c
	gcc -I. -I/home/yyli/include -c test_matmul.c
	gcc test_matmul.o -L./ -lcatlab -L/home/yyli/lib -lopenblas -llapack -lumfpack -larpack -lgfortran -lm -pthread -o test_matmul

test_gbmat: test_gbmat.c
	gcc -I. -I/home/yyli/include -c test_gbmat.c
	gcc test_gbmat.o -L./ -lcatlab -L/home/yyli/lib -lopenblas -llapack -lumfpack -larpack -lgfortran -lm -pthread -o test_gbmat

test_ccmat: test_ccmat.c
	gcc -I. -I/home/yyli/include -c test_ccmat.c
	gcc test_ccmat.o -L./ -lcatlab -L/home/yyli/lib -lopenblas -llapack -lumfpack -larpack -lgfortran -lm -pthread -o test_ccmat

test_eigs_z: test_eigs_z.c
	gcc -I. -I/home/yyli/include -c test_eigs_z.c
	gcc test_eigs_z.o -L./ -lcatlab -L/home/yyli/lib -lopenblas -llapack -lumfpack -larpack -lgfortran -lm -pthread -o test_eigs_z
clean:
	-rm *.o
	-rm *.a
	-rm test_eig_d test_eig_z test_geig_d test_geig_z test_matmul test_gbmat test_ccmat test_eigs_z
