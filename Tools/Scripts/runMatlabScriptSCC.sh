#!/bin/bash
functionHandle=$1
echo "Running $1 on scc-lite"
currentPath=$(pwd) 
cd /usr3/graduate/nithin/MATLAB/
nohup matlab -nodisplay -nosplash -r "run_script($1); quit">/dev/null 2>&1 &
#nohup matlab -nodisplay -nosplash -r "run_script($1); quit">"/usr3/graduate/nithin/git-repos/nohup.txt"&
cd "$currentPath"
echo "Now you can run something else"