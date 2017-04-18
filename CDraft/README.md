# Redfield serial version
This folder cointains a naive version of Redfield equiations
in C.

## Oddisey
module load legacy/0.0.1-fasrc01
module load hpc/intel-mkl-10.0.4.023
module load intel/17.0.2-fasrc01

### Interactive session

* Compilation
1. Run configure file
2. make clean
3. make

* Running
srun -n 1 --mem 1000 ./cgeev.x > output.txt 2> errors.txt
