OS := $(shell uname)
OPTIONS:= 

ifeq ($(OS),Darwin)
	OPTIONS += -framework OpenCL
else
	OPTIONS += -l OpenCL
endif

main: mat_sub_complexified.c
	gcc -Wall -g mat_sub_complexified.c -o mat_sub_complexified $(OPTIONS)

clean:
	rm -rf mat_sub_complexified
