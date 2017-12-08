


CBLASLIB = $(HOME)/lib
LAPACKLIB = $(HOME)/lib
UMFPACKLIB = $(HOME)/lib
ARPACKLIB = $(HOME)/lib

CBLASINC   = /usr/local/intel/Compiler/11.1/ifort/mkl/include 
LAPACKEINC = $(HOME)/include
UMFPACKINC = $(HOME)/include

INCLUDES = -I$(CBLASINC) -I$(LAPACKEINC) -I$(UMFPACKINC) -I../include

CC = icc
CFLAGS = -g 
LDDIR = -L.. -L$(CBLASLIB) -L$(LAPACKLIB) -L$(UMFPACKLIB) -L$(ARPACKLIB) -L/usr/local/intel/Compiler/11.1/ifort/lib/intel64
LDFLAGS = -lcatlab -llapacke -larpack -lumfpack -lamos -lmachcon -lifcore -lifport -lm -lpthread -mkl

.SUFFIXES:
.SUFFIXES: .c .o

%.o: %.c
	${CC} -c ${CFLAGS} ${INCLUDES} $<

%: %.o
	${CC} $< ${LDDIR} ${LDFLAGS}  -o $@
