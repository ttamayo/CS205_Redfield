OS := $(shell uname)
OPTIONS:= 

ifeq ($(OS),Darwin)
	OPTIONS += -framework OpenCL
else
	OPTIONS += -l OpenCL
endif

main: complexTypeTranslate.c
	gcc -Wall -g complexTypeTranslate.c -o complexTypeTranslate $(OPTIONS)

clean:
	rm -rf complexTypeTranslate
