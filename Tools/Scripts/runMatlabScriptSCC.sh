#!/bin/bash
functionHandle=$1
echo "Running $1"
currentPath=$(pwd) 
cd /usr3/graduate/nithin/MATLAB/
nohup matlab -nodisplay -nosplash -r "run_script($1); quit">/dev/null 2>&1 &
cd $currentPath
#nohup /usr/local/MATLAB/R2018b/bin/matlab -nodisplay -nosplash -r "run_script($1); quit">"/home/nithin/Documents/git-repos/nohup.txt" &
echo "Run complete"