include ../make.inc

OBJS = catlab_eig.o catlab_matrix.o catlab_specfunc.o
SRCS = catlab_eig.c catlab_matrix.c catlab_specfunc.c

default: $(OBJS)
	ar rv ../libcatlab.a $(OBJS)
	ranlib ../libcatlab.a

clean:
	-rm *.o
	-rm ../libcatlab.a