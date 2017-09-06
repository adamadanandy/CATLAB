include ../make.inc

EXE = test_eig_d test_eig_z test_geig_d test_geig_z test_matmul test_gbmat test_ccmat test_eigs_z

test: $(EXE)

test_eig_d: test_eig_d.o

test_eig_z: test_eig_z.o

test_geig_d: test_geig_d.o

test_geig_z: test_geig_z.o

test_matmul: test_matmul.o

test_gbmat: test_gbmat.o

test_ccmat: test_ccmat.o

test_eigs_z: test_eigs_z.o

clean:
	-rm *.o
	-rm $(EXE)