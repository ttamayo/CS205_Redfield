Compilation on Oddisey:

### Loading

module load pgi/14.10-fasrc01 openmpi/1.10.0-fasrc01 netcdf/4.1.3-fasrc05
module swap openmpi mvapich2/2.2-fasrc01
 
### Setting some env variables and licence

export PGI=/n/seasfs03/IACS/cs205/pgi
export PATH=$PGI/linux86-64/16.10/bin:$PATH
export MANPATH=$MANPATH:$PGI/linux86-64/16.10/man
export LM_LICENSE_FILE=$LM_LICENSE_FILE:$PGI/license.dat
export MPIDIR=/usr/lib64/slurm

### Compiling
#### Last 
mpicc  *.c -lblas -llapack -lm -lpgftnrtl -lrt -Minfo=accel  -O2 -fast -acc  


 mpicc  *.c -lblas -llapack -lm -lpgftnrtl -lrt -Minfo=accel -ta=tesla,cc35 -O2 -fast 

#### Maybe?
pgcc -acc *.c -lblas -llapack -lm -lpgftnrtl -lrt -Minfo=accel -ta=tesla,cc35 -Mmpi=mpich1
