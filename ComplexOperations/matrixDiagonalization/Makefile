OS := $(shell uname)
OPTIONS:= 

ifeq ($(OS),Darwin)
	OPTIONS += -framework OpenCL
else
	OPTIONS += -l OpenCL
endif

main: mat_diagonalize.c
	gcc -Wall -g mat_diagonalize.c -o mat_diagonalize $(OPTIONS) -llapack -lm

clean:
	rm -rf matAdd
