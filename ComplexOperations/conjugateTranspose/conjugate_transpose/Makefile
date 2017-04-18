OS := $(shell uname)
OPTIONS:= 

ifeq ($(OS),Darwin)
	OPTIONS += -framework OpenCL
else
	OPTIONS += -l OpenCL
endif

main: mat_conjugate_transpose.c
	gcc -Wall -g mat_conjugate_transpose.c -o mat_conjugate_transpose $(OPTIONS)

clean:
	rm -rf mat_conjugate_transpose
