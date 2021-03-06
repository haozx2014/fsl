# FSL configuration file
#  - to be sourced by the user, typically in .bashrc or equivalent
#  - note that the user should set

# Written by Mark Jenkinson, FMRIB Analysis Group, University of Oxford
# Modified for Debian by Michael Hanke <michael.hanke@gmail.com>

# clean out previous fsl PATH components: DO NOT EDIT THE NEXT TWO LINES
PATH=$( echo $PATH | tr ":" "\n" | grep  -v "/usr/lib/fsl/" | tr -s "\n" ":" | sed 's/:$//')
LD_LIBRARY_PATH=$( echo $LD_LIBRARY_PATH | tr ":" "\n" | grep  -v "/usr/lib/fsl/" | tr -s "\n" ":" | sed 's/:$//')

#### Set up standard FSL user environment variables ####

# Debian has a fixed FSLDIR
FSLDIR=/usr/share/fsl/5.0

# Possum is installed in the same directory
POSSUMDIR=$FSLDIR

# add the fsl binary path to the search path
PATH=/usr/lib/fsl/5.0:$PATH

# The following variable selects the default output image type
# Legal values are:
# NIFTI, NIFTI_PAIR, NIFTI_GZ, NIFTI_PAIR_GZ
# This would typically be overwritten in ${HOME}/.fsl/fsl.sh if the user
# wished to write files with a different format
FSLOUTPUTTYPE=NIFTI_GZ

# Comment out the definition of FSLMULTIFILEQUIT to enable
#  FSL programs to soldier on after detecting multiple image
#  files with the same basename ( e.g. epi.hdr and epi.nii )
FSLMULTIFILEQUIT=TRUE

# The following variables specify paths for programs and can be changed
# or replaced by different programs, by default set sensible Debian-defaults
FSLTCLSH=/usr/bin/tclsh
FSLWISH=/usr/bin/wish
FSLBROWSER=/etc/alternatives/x-www-browser



# The following variables are used for running code in parallel across
#  several machines ( i.e. for FDT )
# for a cluster engine setup see below
FSLLOCKDIR=
FSLMACHINELIST=
FSLREMOTECALL=

# If set, tell FSL to use Sun Gridengine to submit jobs instead of running them
# directly on the machine. If unset, no attempt will be made to utilize
# gridengine, even if it is running. By default SGE is not utilized.
#FSLPARALLEL=1

# Mail setup for gridengine jobs. See man qsub (-m option) for all possible
# settings. By default no email is sent.
FSLCLUSTER_MAILOPTS="n"

# default queue for job submissions
#FSLCLUSTER_DEFAULT_QUEUE="all.q"



###################################################
####    DO NOT ADD ANYTHING BELOW THIS LINE    ####
###################################################

export FSLDIR POSSUMDIR PATH FSLMULTIFILEQUIT FSLOUTPUTTYPE FSLTCLSH \
       FSLWISH FSLBROWSER FSLLOCKDIR FSLMACHINELIST FSLREMOTECALL


# Configure the linker search path for Debian FSLs internal shared libraries
LD_LIBRARY_PATH=/usr/lib/fsl/5.0${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
export LD_LIBRARY_PATH


# load user configuration
if [ -f "${HOME}/.fslconf/fsl.sh" ] ; then
  . "${HOME}/.fslconf/fsl.sh" ;
fi
