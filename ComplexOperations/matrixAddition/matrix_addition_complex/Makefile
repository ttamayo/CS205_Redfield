OS := $(shell uname)
OPTIONS:= 

ifeq ($(OS),Darwin)
	OPTIONS += -framework OpenCL
else
	OPTIONS += -l OpenCL
endif

main: mat_add_complex.c
	gcc -Wall -g mat_add_complex.c -o mat_add_complex $(OPTIONS)

clean:
	rm -rf mat_add_complex
