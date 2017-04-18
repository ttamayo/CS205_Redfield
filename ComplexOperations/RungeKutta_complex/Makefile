OS := $(shell uname)
OPTIONS:= 

ifeq ($(OS),Darwin)
	OPTIONS += -framework OpenCL
else
	OPTIONS += -l OpenCL
endif

main: RungeKutta_complex.c
	gcc -Wall -g RungeKutta_complex.c -o RungeKutta_complex $(OPTIONS)

clean:
	rm -rf RungeKutta_complex
