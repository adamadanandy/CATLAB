include ../make.inc

OBJS = catlab_eig.o catlab_matrix.o
SRCS = catlab_eig.c catlab_matrix.c

default: $(OBJS)
	ar rv ../libcatlab.a $(OBJS)
	ranlib ../libcatlab.a

clean:
	-rm *.o
	-rm ../libcatlab.a