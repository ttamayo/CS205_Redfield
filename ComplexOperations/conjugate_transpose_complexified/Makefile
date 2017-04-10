OS := $(shell uname)
OPTIONS:= 

ifeq ($(OS),Darwin)
	OPTIONS += -framework OpenCL
else
	OPTIONS += -l OpenCL
endif

main: mat_conjugate_transpose_complexified.c
	gcc -Wall -g mat_conjugate_transpose_complexified.c -o mat_conjugate_transpose_complexified $(OPTIONS)

clean:
	rm -rf mat_conjugate_transpose_complexified
