## GeoDataMATLAB
![alt text](https://raw.github.com/jswoboda/GeoDataMATLAB/master/logo/logo1.png "GeoDataMATLAB")
#Overview
This is the repository for the MATLAB version of GeoData to plot and analyze data from geophysics sources such as radar and optical systems.

#Installation
To install first clone repository:

	$ git clone https://github.com/jswoboda/GeoDataMATLAB.git
	
The user should then start MATLAB	 and go to the folder containing the code. 

The software package can be installed in MATLAB by running setup.m, which will add the tools to the MATLAB path. The user can specify they’re developing the toolbox further by adding the string “develop” as the second argument. This will create a directory Test that will be added within the main GeoData directory. This folder will not be added to the path to allow the user to test new functions. There is also an option to create a new path file that can be saved where ever the user wants. It is suggested that this is saved in the MATLAB folder in the user's Documents directory.

~~~matlab
cd GeoData 
setUpPath('permanent','develop'); % perminately saves path and creates test directory.
setUpPath('permanent','develop',' ~/Documents/MATLAB'); % perminately saves path and creates test directory.

~~~

#Software Structure

The code is broken up into three main directories: Main, Utilities and Reading. Main holds code related to the GeoData class such as the class def file and any other functions that directly use or impact the class. The folder Utilities holds code related to functions that would be used to support the class such as coordinate transforms. Lastly the Reading directory is to be used to store functions that will be used to read in data from new data types.

#Style Guide
This style guide will cover conventions and elements specific to this codebase. For more general tips on proper MATLAB style guidelines see The Elements of MATLAB Style by Richard K. Johnson.


The code to read in data will be within functions and output the class data variables in an order shown in the code. These read functions will be placed in the Reading folder and be within a specific file. The names of the functions will start with read_ and then followed by a descriptive name.

Code to impact the class will be placed in the class def file. For example if it is desired to interpolate the data to a different set of coordinates this should be implemented as a method in the class. The code for doing this can be written outside of the class def file if needed but this should only be done to keep the code neat.

The properties names will be all lower case. While all function names will be lower case with _ to separate words. The classes will be have capitalized words with no spaces. Directories will be capitalized with _ to separate words.

If the user would like to create test code please do this in the Test folder. Also this code is not be uploaded to the main code base on GitHub. 

