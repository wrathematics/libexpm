CFLAGS = -I../build/include
LDFLAGS = -lblas -llapack -L../build/ -lexpm_static

expm:
	R CMD COMPILE expm_wrap.c CFLAGS=$(CFLAGS); \
	R CMD SHLIB expm_wrap.o $(LDFLAGS)

.PHONY: expm

clean:
	rm -rf *.so *.o
