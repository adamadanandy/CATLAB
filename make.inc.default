


CBLASLIB = $(HOME)/lib
LAPACKLIB = $(HOME)/lib
UMFPACKLIB = $(HOME)/lib
ARPACKLIB = $(HOME)/lib

CBLASINC   = $(HOME)/include
LAPACKEINC = $(HOME)/include
UMFPACKINC = $(HOME)/include

INCLUDES = -I$(CBLASINC) -I$(LAPACKEINC) -I$(UMFPACKINC) -I../include

CC = gcc
CFLAGS = -Wall -g
LDDIR = -L.. -L$(CBLASLIB) -L$(LAPACKLIB) -L$(UMFPACKLIB) -L$(ARPACKLIB)
LDFLAGS = -lcatlab -lopenblas -llapack -larpack -lumfpack -lmachcon -lamos -lgfortran -lm -pthread

.SUFFIXES:
.SUFFIXES: .c .o

%.o: %.c
	${CC} -c ${CFLAGS} ${INCLUDES} $<

%: %.o
	${CC} $< ${LDDIR} ${LDFLAGS}  -o $@