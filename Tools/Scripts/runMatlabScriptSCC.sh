#!/bin/bash
functionHandle=$1
echo "Running $1 on scc-lite"
currentPath=$(pwd) 
cd /usr3/graduate/$USER/MATLAB/
mkdir -p /projectnb/semetergrp/$USER/local_cluster_object
nohup matlab -nodisplay -nosplash -r "run_script($1); quit">/dev/null 2>&1 &
#nohup matlab -nodisplay -nosplash -r "run_script($1); quit">"/usr3/graduate/nithin/git-repos/nohup.txt"&
rm -rf /projectnb/semetergrp/$USER/local_cluster_object
cd "$currentPath"
echo "Now you can run something else"