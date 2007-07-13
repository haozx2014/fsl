# $Id: externallibs.mk,v 1.2 2004/03/11 13:51:20 mark Exp $

# External Library and Include Paths

# debian settings
DEB_LIBPATH=/usr/lib
DEB_INCPATH=/usr/include

FSLEXTLIB=${FSLDIR}/extras/lib
FSLEXTINC=${FSLDIR}/extras/include

# CEPHES library
LIB_CEPHES = ${FSLEXTLIB}
INC_CEPHES = ${FSLEXTINC}/cephes

# FFTW library
LIB_FFTW = ${DEB_LIBPATH}
INC_FFTW = ${DEB_INCPATH}

# GD library
LIB_GD = ${DEB_LIBPATH}
INC_GD = ${DEB_INCPATH}

# GDC library
LIB_GDC = ${DEB_LIBPATH}
INC_GDC = ${DEB_INCPATH}

# PNG library
LIB_PNG = ${DEB_LIBPATH}
INC_PNG = ${DEB_INCPATH}/libpng

# PROB library
LIB_PROB = ${FSLEXTLIB}
INC_PROB = ${FSLEXTINC}/libprob

# CPROB library
LIB_CPROB = ${FSLEXTLIB}
INC_CPROB = ${FSLEXTINC}/libcprob

# NEWMAT library
LIB_NEWMAT = ${DEB_LIBPATH}
INC_NEWMAT = ${DEB_INCPATH}/newmat

# FFTW3 library
LIB_FFTW3 = ${DEB_LIBPATH}
INC_FFTW3 = ${DEB_INCPATH}

# BOOST library
#BOOSTDIR = /usr/local/boost
#LIB_BOOST = ${BOOSTDIR}
INC_BOOST = ${DEB_INCPATH}/boost
