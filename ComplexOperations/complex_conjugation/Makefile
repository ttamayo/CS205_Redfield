OS := $(shell uname)
OPTIONS:= 

ifeq ($(OS),Darwin)
	OPTIONS += -framework OpenCL
else
	OPTIONS += -l OpenCL
endif

main: mat_complex_conjugate.c
	gcc -Wall -g mat_complex_conjugate.c -o mat_complex_conjugate $(OPTIONS)

clean:
	rm -rf mat_complex_conjugate
