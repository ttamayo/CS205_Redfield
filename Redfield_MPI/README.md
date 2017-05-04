Routines for propagating an open quantum system according to the secular Redfield equations
in nodes with a gpu

Copile with mpi and open acc:
```
	mpicc  *.c -lblas -llapack -lm -lpgftnrtl -lrt -Minfo=accel  -O2 -fast -acc  
```	

	Maybe:

```
	pgcc -acc *.c -lblas -llapack -lm -lpgftnrtl -lrt -Minfo=accel -ta=tesla,cc35 -Mmpi=mpich1
```

Environmental variables for running on odyssey

```
	export PGI=/n/seasfs03/IACS/cs205/pgi
	export PATH=$PGI/linux86-64/16.10/bin:$PATH
	export MANPATH=$MANPATH:$PGI/linux86-64/16.10/man
	export LM_LICENSE_FILE=$LM_LICENSE_FILE:$PGI/license.dat
```

Modules for running in odyssey
```
	module load pgi/14.10-fasrc01 openmpi/1.10.0-fasrc01 netcdf/4.1.3-fasrc05
```
