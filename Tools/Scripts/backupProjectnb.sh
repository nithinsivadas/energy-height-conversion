#!/bin/bash
# rsyncs all files in $dataDir with auroraOptiplex except for those excluded in 'exclude-data.txt' 
rsync -auv --info=progress2 --exclude-from='/usr3/graduate/nithin/git-repos/energy-height-conversion/Tools/Scripts/exclude-data.txt' /projectnb/semetergrp/nithin/Data/ nithin@128.197.173.50:/media/nithin/Elements/Nithin/Data/projectnb
