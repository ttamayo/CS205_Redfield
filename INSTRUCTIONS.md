
MAIN COMPONENTS/DELIVERABLES
----------------------------

1. `./Python_FirstSteps` ... Our first implementation of the Redfield code. Was able to identify bottlenecks and confirm problem scaling.

2. `./ComplexOperations` ... Various matrix operations written in C for optimizing and testing. Copied parts of these operations into final Redfield code.

3. `./RedfieldPropagation_OpenACC` ... Our main/final Redfield code. Please follow the steps below to test the code. 

4. `./MatrixMultiply_OpenMP` ... Advanced feature: Matrix-matrix multiplication implemented using OpenMP to optimize the bottleneck portion of the Redfield code ("Lindblad term" from README.md).

5. `./LindbladTerm_OpenMP` ... Advanced feature (continued): A single iteration of the Lindblad term of the Redfield equations implemented in serial and parallel. 

6. `./Redfield_MPI` ... A work version of an hybrid implementation of MPI + OpenACC.


INSTRUCTIONS TO RUN MAIN/FINAL REDFIELD CODE
--------------------------------------------

The GPU accelerated Redfield propagation is located in the `./RedfieldPropagation_OpenACC` directory. Adapt the global variables `NSITES` and `number_of_steps` in `main.c` to set a particular problem size and a number of integration steps. The code can be compiled for serial execution on a single CPU using the command 

-    gcc *.c -llapack -lm

and for parallel execution on a GPU using the command

-    pgcc -acc *.c -lblas  -llapack  -lm -lpgftnrtl -lrt  -ta=nvidia:cc35,nocache -Msafeptr -Minfo=all

Please note that you must have set the following environment variables to use the `pgcc` compiler on odyssey:

```
export PGI=/n/seasfs03/IACS/cs205/pgi
export PATH=$PGI/linux86-64/16.10/bin:$PATH
export MANPATH=$MANPATH:$PGI/linux86-64/16.10/man
export LM_LICENSE_FILE=$LM_LICENSE_FILE:$PGI/license.dat
```
Run an interactive job on odyssey, for instance on the `holyseasgpu` partition with the command 

```
srun -p holyseasgpu -N 1 -n 1 --gres=gpu:1  --pty --mem 4000 -t 0-12:00 /bin/bash
```

to execute the serial or the parallel version of the code. Redirect the output of the executable to a text file of your choice to store the population dynamics trajectories for your specified system. 
