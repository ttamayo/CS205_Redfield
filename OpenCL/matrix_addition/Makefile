OS := $(shell uname)
OPTIONS:= 

ifeq ($(OS),Darwin)
	OPTIONS += -framework OpenCL
else
	OPTIONS += -l OpenCL
endif

main: main.c
	gcc -Wall -g main.c -o matAdd $(OPTIONS)

clean:
	rm -rf matAdd
