CFLAGS   = -c -O2 -openmp
COMPILER = icc
PRO         = *
OBJECTS     = ${PRO}.o

${PRO}.x : $(OBJECTS)
	$(COMPILER) -o ${PRO}.x $(OBJECTS) -lblas -llapack -lm -lrt -openmp

%.o : %.c
	$(COMPILER) $(CFLAGS) $(<F)

clean :
	rm -fr *.o *.x
