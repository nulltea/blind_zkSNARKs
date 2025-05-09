.PHONY: all clean

all: _fhe_example_params_cffi.o

_fhe_example_params_cffi.o: fhe_example_params.h ../_lazer_cffi.o
	python3 ../params_cffi_build.py fhe_example_params.h ../..

fhe_example_params.h: fhe_example_params.py
	cd ../../scripts && sage lin-codegen.sage ../python/fhe_example/fhe_example_params.py > ../python/fhe_example/fhe_example_params.h

../_lazer_cffi.o:
	cd .. && make

clean:
	rm -rf __pycache__
	rm -f *.c *.o *.so
