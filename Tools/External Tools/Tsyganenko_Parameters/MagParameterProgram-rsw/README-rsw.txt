# ./MaginputmodelONE.f and ./MaginputmodelONE_corrected.f in this
# directory are slightly modified version (by rsw) of the original
# files in ../MagParameterProgram. (They were modified to run
# non-interactively.)  The differences between the original and this
# version are shown in MaginputmodelONE.f.diff.

# On 4/20/2009, Richard Denton provided
# ./MagmodelinputONE_corrected.f, a slightly modified version of
# ./MagmodelinputONE.f that only affects the iG2 variable for hourly
# data.  The difference bewteen the new and old files are shown in
# ./MagmodelinputONE_corrected.f.diff

# To run, execute on the command line
#
# make all

# To test, execute
# make testhour
# make testhour_corrected
# make test5min
# make test1min

# make test2000 will show the differences between the result of
# Richard's run on omni2000, and my run.  The differences are probably
# due to the fact that the (input) omni2000 data file changed between runs.


