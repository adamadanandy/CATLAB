all: lib test

default: lib

lib: src
	cd src; make
test: lib
	cd test; make

clean:
	cd src; make clean
	cd test; make clean