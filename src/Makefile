CFLAGS = -I../build/include
LDFLAGS = -lblas -llapack -L../build/ -lexpm_static

expm:
	MAKEFLAGS='CFLAGS=${CFLAGS}' R CMD SHLIB *.c $(LDFLAGS)

.PHONY: expm

clean:
	rm -rf *.so *.o
