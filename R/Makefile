LDFLAGS = -lblas -llapack -L../build -lexpm_static

expm:
	R CMD SHLIB expm_wrap.c $(LDFLAGS)

.PHONY: expm

clean:
	rm -rf expm.so expm.o
