# MaginputmodelONE.f is a slightly modified version of original file ../MagParameterProgram
# It was modified to run non-interactively.  The differences between the original and this version
# are shown in MaginputmodelONE.f.diff.  To run, execute on the command line
# echo "60" | MaginputmodelONE
# echo "5" | MaginputmodelONE
# echo "1" | MaginputmodelONE
# To test, execute
# make testhour
# make test5min
# make test1min

This program is to calculate W and G parameters using either minute data, 5-minute data or hourly data.
Input: multiple year data in a single file.

Compiler: G77

Under the current directory (in which the Magmodelinput.f is located), there are 3 directories 

directory 1min: kpdst.lst (hourly resolution) and omni_min.asc (1 min resolution)
directory 5min: kpdst.lst (hourly resolution) and omni_5min.asc (5 min resolution)
   NOTE:kpdst.lst and 1min/omni_min.asc or 5min/omni_5min.asc should begin and end at the same time,
     For instance, if omni_min.asc starts at 2000/001/00/00(mm), and ends at 2000/366/23/59, 
     the kpdst.lst file must include data from 2000/001/00(yyyy/ddd/hh) to 2000/366/23.
directory hour: omni2_hour.dat (doesn't need separate file for Kp and Dst values because Kp and 
   Dst are already in this file)

omni_min.asc and omni_5min.asc are downloaded from:
ftp://nssdcftp.gsfc.nasa.gov/spacecraft_data/omni/high_res_omni/
Then the yearly files with year YYYY (omni_minYYYY.asc or omni_5minYYYY.asc) are combined into a 
  single file (omni_min.asc or omni_5min.asc, respectively)

omni2_hour.dat are downloaded from:
ftp://nssdcftp.gsfc.nasa.gov/spacecraft_data/omni/
Then the yearly files with year YYYY (omni2_xxxx.dat) are combined into a single file omni2_hour.dat


kpdst.lst files are downloaded from:
http://omniweb.gsfc.nasa.gov/form/dx1.html
Select "Create file"
Select "Hourly averaged"
Start with YYYY0101
 stop with XXXX1231         
Select "Kp*10 index" and "Dst index"
Submit

Then you will see:

"Location:  nssdcftp.gsfc.nasa.gov
 Directory: /staging/ftpweb

 Filename          Description         Size (bytes)
 * omni2_*****.lst ASCII data file     183960
 * omni2_*****.fmt Format description  303           "
Download omni2_*****.lst and change its name to kpdst.lst

Compile the program MagmodelinputONE.f with:

  g77 -o MagmodelinputONE.x MagmodelinputONE.f

When you run the program, there are two input choices.  First, you select the number of minutes 
  per data point (1 for 1 min resolution data, 5 for 5 min resolution data, or 60 for hourly data).
Secondly, you choose an initialization option for the initial W parameters, 1 for starting from zero 
  values or 2 for specifying the W values at a time interval immediately before the first time 
  interval of the file being analyzed
If initialization option 2 is chosen, you also need to enter the preceding W values

After running the program, the output file will be in the directory with the input data (1min, 5min,   or hour).  The output file will be called WGparametersmin.d (1 min data), 
  WGparameters5min.d (5 min data), or WGhour.d (hourly data).

For a test, there are sample input files (Year 2000) and output files in each directory. WG*_0.d is the output file with the starting option 1 (all the W values start from 0), and WG*_3.d is the output file with the starting option 2 using starting values (W_1,W_2,W_3,W_4,W_5,W_6) = (3.,3.,3.,3.,3.,3.).  If you run the program, you should be able to reproduce these results.  
