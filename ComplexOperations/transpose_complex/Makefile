OS := $(shell uname)
OPTIONS:= 

ifeq ($(OS),Darwin)
	OPTIONS += -framework OpenCL
else
	OPTIONS += -l OpenCL
endif

main: matrix_transpose.c
	gcc -Wall -g matrix_transpose.c -o matrix_transpose $(OPTIONS)

clean:
	rm -rf matrix_transpose
