#PBS -S /bin/bash

# set name of job
#PBS -N ben_reiss_ex_3

# copy enviroment variables
#PBS -V

# ressources
#PBS -l nodes=4:ppn=1
#PBS -l walltime=00:10:00

# write error and standard output in one file
#PBS -j oe

# change to working directory
cd $PBS_O_WORKDIR

#program
make clean
make

mpirun -np 4 -machinefile $PBS_NODEFILE caseq 14 1
