#=======================================================================
# Make file for omp_diag.f90
#=======================================================================
CCFLAGS   = -c -O3
CCOMPILER = icc
PRO         = cgeev
OBJECTS     = cgeev.o

LIB = -mkl -lblas -llapack 



${PRO}.x : $(OBJECTS)
	$(CCOMPILER) -o ${PRO}.x $(OBJECTS) $(LIB)

redfield.o : 
	$(CCOMPILER) $(CCFLAGS) cgeev.c $(LIB)

clean :
	rm -fr *.o *.x
