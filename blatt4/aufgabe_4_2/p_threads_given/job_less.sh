#PBS -S /bin/bash

# set name of job
#PBS -N kubicek_reissaus_ex_2

# copy enviroment variables
#PBS -V

# ressources
#PBS -l nodes=1:ppn=8
#PBS -l walltime=02:59:00

# write error and standard output in one file
#PBS -j oe

# change to working directory
cd $PBS_O_WORKDIR

#program
make clean
make

samples=10000000
./pi 1 $samples 
./pi 4 $samples
