#!/bin/bash
#$ -pe omp 4
# set default value for n; override with qsub -v at runtime
#$ -v n=100

# loading the version of matlab of interest
functionHandle=$1
echo "Running $1 on scc-lite"
currentPath=$(pwd) 
cd /usr3/graduate/nithin/MATLAB/
matlab -nodisplay -nosplash -r "run_script($1); exit">/dev/null 2>&1 &
#nohup matlab -nodisplay -nosplash -r "run_script($1); quit">"/usr3/graduate/nithin/git-repos/nohup.txt"&
cd "$currentPath"
echo "Now you can run something else"