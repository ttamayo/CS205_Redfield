make clean
make
srun -n 1 --mem 1000 ./redfield.x > output.txt 2> errors.txt
