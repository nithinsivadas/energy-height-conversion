NOTE:	This file is adapted from a similar introductory text 
	composed by Mr Clive G. Page (cgp@star.le.ac.uk) for his 
	compact packaging of the EMX-based G77 for DOS
	(available from ftp://ftp.star.le.ac.uk/pub/fortran/ )
----------------------------------------------------------------------------
	Last update of this document:  August 1999
----------------------------------------------------------------------------


This file contains a brief introduction to g77, the GNU Fortran compiler,
here packaged for use on PCs running Windows 95/98 or Windows NT.

----------------------------------------------------------------------------
1. What is g77?
----------------------------------------------------------------------------

G77 is a Fortran compiler from the Free Software Foundation (GNU). The g77
compiler is designed to be portable, and there are implementations  on a
wide variety of systems including Linux, other Unixes, OS/2, MS-DOS and
Windows 95/98 and Windows NT.

The g77 compiler supports the whole of the (now obsolete) Fortran77
Standard plus many commonly-used extensions.  

Typographical note: in this introduction Fortran keywords and code are
shown in UPPER CASE merely to make them stand out; g77 is equally happy
with Fortran in lower-case.

----------------------------------------------------------------------------
2. How to download.
----------------------------------------------------------------------------

The compact G77 fop Win32 (Windows 95/NT) package is available from 
  http://www.geocities.com/Athens/Olympus/5564

There are three zip files:

  g77exe.zip    [1,574Kb]  The main executables
  g77lib.zip    [  208Kb]  The libraries
  g77doc.zip    [  301Kb]  Documentation and other stuff.


This port is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
more details.  The GNU Public licence can be obtained by writing to
the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.  
The GNU Public Licence requires the corresponding source code to be 
available. Alternatively it can be obtained from a number of Internet sites 
including: 	
http://www.gnu.org/copyleft/copyleft.html


----------------------------------------------------------------------------
3. How to install.
----------------------------------------------------------------------------

You will need at least 6.5 MB of free space on your hard disc (11 Mb
if you install the optional SLATEC math library and PSPLOT scientific
graphics library.)  If you start from your root directory the files 
will unzip into subdirectories \G77, \G77\BIN, \G77\LIB and \G77\DOC 
as necessary.  You are free to put them somewhere else but if so you 
will need to update the file G77SETUP.BAT to reflect the changes.  



The simplest commands are:

  CD \
  UNZIP G77EXE.ZIP
  UNZIP G77LIB.ZIP	
  UNZIP G77DOC.ZIP

NOTE:  Using the free 32-bit "unzip" utility that's available 
       from:
         http://www.geocities.com/Athens/Olympus/5564/g77.htm
       ensures that unzipping the above zipfiles will preserve
       the directory structure of the archives _and_ the
       long filenames that some of the library files have.	
       Do NOT use a 16-bit unzip utility (e.g., an older,
       DOS version, of PKUNZIP) as then the long filenames 
       will not be preserved for sure and the compiler will 
       complain that it cannot find some of the library files.

To use the system you need to execute the very short batch file 
G77SETUP.BAT to set up the G77 environment.  It adds the folder 
that holds the executable files, \G77\BIN, to your path and defines 
an environment variable LIBRARY_PATH that points to the library 
folder \G77\LIB (this variable tells the compiler where to find the 
library files.)  The default locations are assumed to be in the
folder C:\G77 . If a user prefers to install G77 in another drive
(e.g., in F:\G77) then the file G77SETUP.BAT must be corrected.
The file G77SETUP.BAT _must_ be run _everytime_ a user opens 
a "console" window (aka, "DOS" window) in order to run G77.  
Alternatively a user might prefer to add the contents of the 
G77SETUP.BAT in his/her C:\AUTOEXEC.BAT file, in which case 
the G77 environment will be automatically visible in all "DOS" 
windows everytime the PC is booted.


----------------------------------------------------------------------------
4. How to use the compiler.
----------------------------------------------------------------------------

The command G77 will compile and link your program and produce an
executable that can only be run under Windows 95 or Windows NT,
not under DOS (nor in a DOS box of Windows 3.1x)  There are various 
command options (see below) which need to be specified before the 
name of the program file(s).  A simple demonstration program called 
MYTEST1.F is included in the package.  To compile this just use:
  G77 MYTEST1.F

This generates a file A.EXE, which can be executed simply using
  A

The compiler by default names the produced executable file as "A.EXE".
If one wishes to override this default behavior, then one can specify
the desired name of the produced executable, using the "-o" compiler
option.  E.g., to produce an executable with the name MYTEST1.EXE
the above command should be:
  G77 MYTEST1.F -o MYTEST1.EXE
 
It is good practice to get the compiler to warn you of unused and
undeclared variables, the switch to do this is "-Wall" which can be used
like this:
  G77 -Wall MYTEST1.F

Note that the switches are case-sensitive: -WALL will not work.

The g77 compiler has a large number of other command switches - a few of
the most useful are shown here:

-c               Compile-only: produces .OBJ files.
-ffree-form      Selects free-format source code
-fpedantic       Warns of non-portable/non-standard code.
-fno-automatic   Static storage for all variables, like universal SAVE
-fno-backslash   Interprets "\" as a normal character in strings
-fvxt            Use VAX Fortran interpretation of certain syntax
-g               Produces debugging information.
-Idirectory      Specifies directory to search for INCLUDE files
-O               Optimise code generation
-Wimplicit       Warns of any names with no explicit data type
-Wuninitialised  Warns of some cases of unset variables (if -O also set).
-Wall            Warns of both of above cases.

Any number of source-code files can be used on the command-line, and
wild-cards (*) are supported in file-names.  The command-line can also
contain compiled object modules (.O) and library files (called "archives"
in Unix-speak) which have .A as the extension.  File-names are, of course,
not case-sensitive in under Windows 95 or NT.


NOTE:	How To Capture Errors Produced By The Compiler In a File:
        ---------------------------------------------------------
	The command shell of Windows 95/98 does not provide for re-direction
        of the standard error output.  So, if a user wishes to capture the
	output of the compiler in a file, then a possible solution is by
        using a utility like "etime" (available as 
	http://www.geocities.com/Athens/Olympus/5564/etime.zip ) 
	for this purpose.  In this case, running the compiler as:
		ETIME -2ERRFILE G77 MYFILE.F
	will capture all the errors in the file ERRFILE.

	Under Windows NT the command shell supports Unix-like redirection
	of the standard error output, so a user only needs to enter:
		G77 MYFILE.F 2>ERRFILE
	in order to capture the errors in the file ERRFILE.

----------------------------------------------------------------------------
5. Extensions to g77 compatible with Fortran90:
----------------------------------------------------------------------------

The g77 compiler conforms fully only to the Fortran77 Standard.  Although
this is technically obsolete, the vast majority of Fortran programs in the
world were written (more or less) to this Standard.  The new standard,
Fortran90, has many more advanced features.  The g77 compiler already
supports some features of Fortran90, the most notable of these being:

IMPLICIT NONE (to flag non-explicit data types in a program unit).
INCLUDE 'filename' 
Automatic arrays (a form of dynamic storage within subprograms).
Free-format input (if the switch -ffree-form is used) 
Multiple statements on a line (with ";" as the separator).
Symbolic names can be up to 31-characters long, can contain lower-case
  letters (equivalent to upper-case) and underscores.
End-of-line comments may be used starting with ! (but ! must not be
  in column 6 with fixed-format source-code).
Relational operators > >= < <= == /= can be used instead of .GT. etc.
Character constants can use "double" or 'single' quotes.
Program unit names permitted on END, e.g. END SUBROUTINE MYSUB.
DO without labels and END DO are permitted, also indefinite DO (but a
  conditional EXIT or STOP is obviously needed in such loops).
DO WHILE(logical expression) with END DO is permitted.
EXIT and CYCLE are allowed in DO loops.
SELECT CASE structure is supported but only with integer/logical selectors.
Construct names are allowed with IF/DO/CYCLE/EXIT/SELECT CASE.
Zero-length strings are valid.
Substrings of character constants are permitted.
Character intrinsic functions ACHAR, IACHAR, and LEN_TRIM are provided.
Bit-wise integer functions BTEST, IAND, IBCLR, IBITS, IBSET, IEOR, IOR,
  ISHFT, ISHFT, MVBITS, NOT (the MIL-STD 1753 intrinsics) are provided.
OPEN with STATUS='REPLACE' is supported.
NAMELIST input/output is also supported.
Type declarations may use KIND values (but this is of limited use because
  kind-selection functions are not yet provided).

----------------------------------------------------------------------------
6. Other g77 extensions (NOT compatible with Fortran90)
----------------------------------------------------------------------------

Many extensions to the official Fortran77 Standard were introduced by
companies which produced Fortran compilers for sale, but not all of these 
were incorporated into Fortran90.  You may find that existing "Fortran77"
code makes use of some of these non-standard features.  Fortunately g77
supports some of the more common extensions, especially those of VAX
Fortran.  The most important ones are listed below.  They make it possible
to use "legacy" code with a minimum of alteration, but are NOT recommended
if you are writing new code.

Data types BYTE, INTEGER*1, INTEGER*2, INTEGER*4 (default), INTEGER*8,
 REAL*4 (default), REAL*8, DOUBLE COMPLEX, COMPLEX*8 (default), COMPLEX*16,
 LOGICAL*1, LOGICAL*2, LOGICAL*4 (default) are supported.
DATA statements may be intermixed with specifications statements.
Arguments of procedure calls may use %VAL, %REF, %LOC, %DESCR functions.
Additional intrinsic functions: LOC, AND, OR, LSHIFT, RSHIFT, XOR,
 CDABS, CDCOS, CDEXP, CDLOG, CDSIN, CDSQRT, DCMPLX, DCONJG, DFLOAT,
 DIMAG, DREAL, IMAG, ZABS, ZCOS, ZEXP, ZLOG, ZSIN, and ZSQRT supported.
Any number of continuations lines may be used.
Symbolic names may include "$" if -fdollar-ok switch is specified..
Integer constants may be specified in other number bases: e.g. B'01', 
 O'01234567', X'01EF', Z'FFFF' etc.; in addition "typeless" constants  may
 be given in a similar form but with the letter following the string of 
 digits.
FORMAT specifications may include $ to suppress the final carriage-return.
Debug lines starting "D" or "d" are treated as comments.

----------------------------------------------------------------------------
7. Notes on Input/Output 
----------------------------------------------------------------------------

I/O unit numbers must be in the range 0 to 99, with 0 and 6 pre-connected to
  the screen, 5 to the keyboard (but best to use UNIT=* for both).
Unformatted direct-access files have bytes as units of record length.
Output to the screen does not count as "printing" in Fortran terms, so
  the first character of each record is never removed for carriage-control.
Output to unit N is sent to file "FORT.N" if no file was opened for it.


----------------------------------------------------------------------------
8. The Fortran Library 
----------------------------------------------------------------------------

G77 supports, even under Windows 95/98/NT, most of the routines commonly used 
to access system services on Unix.  These include: ABORT, GETARG, IARGC, 
EXIT, CTIME, DATE, DTIME, ETIME, FDATE, GMTIME, LTIME, IDATE, ITIME, MCLOCK,
SECNDS, SECOND, SYSTEM, TIME,  ERF, DERF, ERFC, DERFC, IRAND, RAND,
SRAND, ACCESS, CHDIR, CHMOD, GETCWD, FSTAT, SSTAT, LSTAT, RENAME, UMASK,
UNLINK, GERROR, IERRNO, PERROR, GETENV, FGETC, GFET, FPUTC, FPUT, FLUSH,
FNUM, FSEEK, FTELL, ISATTY, TTYNAM, SLEEP.

Note that I have not checked all of these.  HOSTNM and GETLOG are  not
supported, and SYSTEM (which executes an operating system command line)
seems to work in a DOS-box under Windows 95/98/NT

For details of the calling sequences see the g77 documentation.


----------------------------------------------------------------------------
9. How to use object libraries
----------------------------------------------------------------------------

A version of the standard Unix library manager, unhelpfully called AR, 
is provided.  It manages what it calls "archives" which are just 
collections of object modules.  These libraries have ".A" as their 
filename extension and the prefix "LIB" to their name.  E.g., the library 
with the "basic" name SLATEC will be stored in the file LIBSLATEC.A  
To create a library you need to compile the source files with the -c switch 
to produce one or more object files (.O).  These can be put in a library 
like this:

  AR -r  LIBMYLIB.A  FILE1.O FILE2.O ...etc  (note the prefix "LIB")
or, even:
  AR -r  LIBMYLIB.A  *.O

To list the contents of an archive use
  AR -t LIBMYLIB.A

More details regarding the usage of the librarian utility "AR" are
given in the text file \G77\AR.DOC.
The generated libraries are best placed in the same directory where
all the G77 libraries were placed during installation (\G77\LIB\). 
This directory is pointed at by the environment variable LIB_PATH
as set by the setup file G77SETUP.BAT   Then, the libraries can be 
referenced in a G77 command line by the "-l" option followed (without 
a space) with the library _basic_ name without the extension ".A", 
_nor_ the prefix "LIB".  E.g., a library file LIBMYLIB.A that exists 
in \G77\LIB\ has the _basic_ library name "MYLIB".  Then if one 
specifies -lMYLIB on the G77 command-line, it will be searched during 
the linking phase  and link only those modules which are called 
directly or indirectly, e.g.,  
  G77 MYFILE.F -o MYFILE.EXE -lMYLIB
Alternatively, one can specify the _full_ pathname of a library file 
on the command line.  E.g., to link with the modules contained with 
the library file LIBXYZ.A which exists in the directory E:\QWERTY, 
one can write:
  G77 MYFILE.F -o MYFILE.EXE E:\QWERTY\LIBXYZ.A
If the library file is present in the current directory then one
can write:
  G77 MYFILE.F -o MYFILE.EXE LIBXYZ.A
In general, libraries that are not to be used temporarily are best 
placed in \G77\LIB and referenced as explained above with the "-l" 
option, using only the basic library name.
	

	

----------------------------------------------------------------------------
10. Debugging
----------------------------------------------------------------------------

[to be written]

----------------------------------------------------------------------------
11. Scientific Graphics using Fortran
----------------------------------------------------------------------------
I recommend the PSPLOT package written by Mr Kevin Kohler
(http://www.nova.edu/cwis/oceanography/psplot.html)  
The ready-to-use binary form of the library is available from:
  http://www.geocities.com/Athens/Olympus/5564/
The PSPLOT library is free to use for non-commercial purposes.




----------------------------------------------------------------------------
12. Source Code and Further Information
----------------------------------------------------------------------------

All of this software has been released by the various authors under the
terms of the GNU Public Licence.  The licence requires source code to be 
provided or easily available.  These sources are widely available on the 
Internet.


The G77 documentation is in the archive g77doc.zip, available from:   
http://www.geocities.com/Athens/Olympus/5564/
It is strongly recommended that the user installs these G77 documentation
files for future reference (these are just ordinary .htm files and can be 
viewed using any web browser software as Netscape or Microsoft Internet 
Explorer.)
-------------------------------------------------------------------------

K Kourakis
kkou@geocities.com
August 1999