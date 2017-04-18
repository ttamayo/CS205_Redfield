OS := $(shell uname)
OPTIONS:= 

ifeq ($(OS),Darwin)
	OPTIONS += -framework OpenCL
else
	OPTIONS += -l OpenCL
endif

main: mat_add_complexified.c
	gcc -Wall -g mat_add_complexified.c -o mat_add_complexified $(OPTIONS)

clean:
	rm -rf mat_add_complexified
