#!/bin/bash
functionHandle=$1
echo "Running $1"
nohup /usr/local/MATLAB/R2018b/bin/matlab -nodisplay -nosplash -r "run_script($1); quit">/dev/null 2>&1 &
#nohup /usr/local/MATLAB/R2018b/bin/matlab -nodisplay -nosplash -r "run_script($1); quit">"/home/nithin/Documents/git-repos/nohup.txt" &
echo "Run complete"