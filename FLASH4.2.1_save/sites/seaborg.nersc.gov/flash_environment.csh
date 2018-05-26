# Environment settings for University of Chicago Flash Center operations

umask 07

setenv FLASH_GROUP_ID incite13
setenv MP_COREFILE_FORMAT 1   # Generate lightweight core file format,
                              # to file rather than to stdout
			      
setenv FLASH_SITE_NAME seaborg.nersc.gov
setenv FLASH_PROJECT_HOME /project/projectdirs/incite13/project/seaborg
setenv FLASH_SOURCE_TREE $FLASH_PROJECT_HOME/Flash3_wd_def

set path =  ($path $FLASH_PROJECT_HOME/bin)

module load python
