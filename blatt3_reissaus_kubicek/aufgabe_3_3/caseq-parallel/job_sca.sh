#PBS -S /bin/bash

# set name of job
#PBS -N reiss_kub_caseq_par

# copy enviroment variables
#PBS -V

# ressources
#PBS -l nodes=3:ppn=8
#PBS -l walltime=00:59:00

# write error and standard output in one file
#PBS -j oe

# change to working directory
cd $PBS_O_WORKDIR

#program
make -f Makefile_Scalasca clean
make -f Makefile_Scalasca

/usr/bin/time -p scalasca -analyze mpiexec -np 24 -machinefile $PBS_NODEFILE caseq 10000 20000
