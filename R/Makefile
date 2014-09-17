CFLAGS = -I../build/include
LDFLAGS = -lblas -llapack -L../build/ -lexpm_static

expm:
	MAKEFLAGS='CFLAGS=${CFLAGS}' R CMD SHLIB expm_wrap.c $(LDFLAGS)

.PHONY: expm

clean:
	rm -rf *.so *.o
