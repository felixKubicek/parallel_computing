#PBS -S /bin/bash

# set name of job
#PBS -N kubicek_reissaus_ex_2

# copy enviroment variables
#PBS -V

# ressources
#PBS -l nodes=1:ppn=8
#PBS -l walltime=00:10:00

# write error and standard output in one file
#PBS -j oe

# change to working directory
cd $PBS_O_WORKDIR

#program
make clean
make

samples=1000000
echo ""
echo "Executing with 1,2,4,8,16 threads."
/usr/bin/time -p ./pi 1 $samples 
/usr/bin/time -p ./pi 2 $samples 
/usr/bin/time -p ./pi 4 $samples 
/usr/bin/time -p ./pi 8 $samples 
/usr/bin/time -p ./pi 16 $samples 
echo ""
echo "results for 8 threads 4 times:"
/usr/bin/time -p ./pi 8 $samples
/usr/bin/time -p ./pi 8 $samples
/usr/bin/time -p ./pi 8 $samples
/usr/bin/time -p ./pi 8 $samples

