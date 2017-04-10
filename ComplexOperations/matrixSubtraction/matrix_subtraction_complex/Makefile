OS := $(shell uname)
OPTIONS:= 

ifeq ($(OS),Darwin)
	OPTIONS += -framework OpenCL
else
	OPTIONS += -l OpenCL
endif

main: mat_sub_complex.c
	gcc -Wall -g mat_sub_complex.c -o mat_sub_complex $(OPTIONS)

clean:
	rm -rf mat_sub_complex
