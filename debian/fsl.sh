# FSL configuration file 
#  - to be sourced by the user, typically in .bashrc or equivalent
#  - note that the user should set 

# Written by Mark Jenkinson, FMRIB Analysis Group, University of Oxford
# Modified for Debian by Michael Hanke <michael.hanke@gmail.com>

#### Set up standard FSL user environment variables ####

# Debian has a fixed FSLDIR
FSLDIR=/usr/share/fsl
export FSLDIR

# add the fsl binary path to the search path
PATH=$PATH:/usr/lib/fsl
export PATH

# The following variable selects the default output image type
# Legal values are:
# NIFTI, NIFTI_PAIR, NIFTI_GZ, NIFTI_PAIR_GZ
# This would typically be overwritten in ${HOME}/.fsl/fsl.sh if the user
# wished to write files with a different format
FSLOUTPUTTYPE=NIFTI_GZ
export FSLOUTPUTTYPE

# Comment out the definition of FSLMULTIFILEQUIT to enable 
#  FSL programs to soldier on after detecting multiple image
#  files with the same basename ( e.g. epi.hdr and epi.nii )
FSLMULTIFILEQUIT=TRUE ; export FSLMULTIFILEQUIT

# The following variables specify paths for programs and can be changed
#  or replaced by different programs ( e.g. FSLDISPLAY=open   for MacOSX)

FSLTCLSH=/usr/bin/tclsh
FSLWISH=/usr/bin/wish

FSLBROWSER=/etc/alternatives/x-www-browser

export FSLTCLSH FSLWISH FSLBROWSER


# The following variables are used for running code in parallel across
#  several machines ( i.e. for FDT )

FSLLOCKDIR=
FSLMACHINELIST=
FSLREMOTECALL=

export FSLLOCKDIR FSLMACHINELIST FSLREMOTECALL


###################################################
####    DO NOT ADD ANYTHING BELOW THIS LINE    ####
###################################################

# Configure the linker search path for Debian FSLs internal shared libraries
LD_LIBRARY_PATH=/usr/lib/fsl${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}
export LD_LIBRARY_PATH


# load user configuration
if [ -f "${HOME}/.fslconf/fsl.sh" ] ; then
  . "${HOME}/.fslconf/fsl.sh" ;
fi
