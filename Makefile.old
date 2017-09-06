COPT = -Wall -g
CINCLUDE = -I. -I/home/yyli/include
LOPT = -lcatlab -lopenblas -llapack -lumfpack -larpack -lgfortran -lm -pthread
LINCLUDE = -L. -L/home/yyli/lib

default: catlab_eig.c catlab_matrix.c
	gcc -c catlab_matrix.c ${COPT} ${CINCLUDE}
	gcc -c catlab_eig.c ${COPT} ${CINCLUDE}
	ar rv libcatlab.a *.o
	ranlib libcatlab.a

test: test_eig_d test_eig_z test_geig_d test_geig_z test_matmul test_gbmat test_ccmat test_eigs_z

test_eig_d: test_eig_d.c
	gcc ${COPT} ${CINCLUDE} -c test_eig_d.c
	gcc test_eig_d.o ${LINCLUDE} ${LOPT} -o test_eig_d

test_eig_z: test_eig_z.c
	gcc ${COPT} ${CINCLUDE} -c test_eig_z.c
	gcc test_eig_z.o ${LINCLUDE} ${LOPT} -o test_eig_z

test_geig_z: test_geig_z.c
	gcc ${COPT} ${CINCLUDE} -c test_geig_z.c
	gcc test_geig_z.o ${LINCLUDE} ${LOPT} -o test_geig_z

test_geig_d: test_geig_d.c
	gcc ${COPT} ${CINCLUDE} -c test_geig_d.c
	gcc test_geig_d.o ${LINCLUDE} ${LOPT} -o test_geig_d

test_matmul: test_matmul.c
	gcc ${COPT} ${CINCLUDE} -c test_matmul.c
	gcc test_matmul.o ${LINCLUDE} ${LOPT} -o test_matmul

test_gbmat: test_gbmat.c
	gcc ${COPT} ${CINCLUDE} -c test_gbmat.c
	gcc test_gbmat.o ${LINCLUDE} ${LOPT} -o test_gbmat

test_ccmat: test_ccmat.c
	gcc ${COPT} ${CINCLUDE} -c test_ccmat.c
	gcc test_ccmat.o ${LINCLUDE} ${LOPT} -o test_ccmat

test_eigs_z: test_eigs_z.c
	gcc ${COPT} ${CINCLUDE} -c test_eigs_z.c
	gcc test_eigs_z.o ${LINCLUDE} ${LOPT} -o test_eigs_z

test_geigs_z: test_geigs_z.c
	gcc ${COPT} ${CINCLUDE} -c test_geigs_z.c
	gcc test_geigs_z.o ${LINCLUDE} ${LOPT} -o test_geigs_z

clean:
	-rm *.o
	-rm *.a
	-rm test_eig_d test_eig_z test_geig_d test_geig_z test_matmul test_gbmat test_ccmat test_eigs_z test_geigs_z
