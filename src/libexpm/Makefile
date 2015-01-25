all: documentation libexpm R

documentation:
	doxygen libexpm.conf

libexpm:
	mkdir -p build;
	cd build; \
	mkdir include; \
	cmake ..; \
	make install

R: libexpm
	cd R; \
	make

clean:
	rm -rf ./build ./docs/html/; \
	cd R; make clean
