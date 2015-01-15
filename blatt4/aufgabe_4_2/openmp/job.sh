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

samples=1000000000
echo ""
echo "Executing with 1,2,4,8,16 threads."

for t in 1 2 4 8 16
do
  echo "threads $t"
  for x in {1..5}
  do
    ./pi $t $samples; 
  done
done
