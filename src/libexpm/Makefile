all: documentation libexpm

documentation:
	doxygen libexpm.conf

libexpm:
	mkdir -p build;
	cd build; \
	mkdir include; \
	cmake ..; \
	make install

clean:
	rm -rf ./build ./docs/html/; \
	cd R; make clean
