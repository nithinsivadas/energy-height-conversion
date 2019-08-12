#!/bin/bash
currentPath=$(pwd)
cd $gitRootDir
git fetch 
git pull
git add --all :/
git commit -m 'backup from scc-lite'
git push origin master
cd $currentPath
echo "Backed up energy-height-conversion to git"