# Assumes that ${FSLDIR} is set as an environment variable

# use Debian specific build configuration
include ${FSLDIR}/debian/systemvars.mk
include ${FSLDIR}/debian/externallibs.mk
include ${FSLCONFDIR}/common/vars.mk
include ${FSLCONFDIR}/common/rules.mk
