# $Id: Makefile,v 1.19 2007/06/14 14:50:21 duncan Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = first

USRINCFLAGS = -I${INC_NEWMAT}

USRLDFLAGS = -L${LIB_NEWMAT} 


LIBS= -lshapeModel  -lmeshclass -lnewimage -lprob -lmiscmaths -lfslio -lniftiio -lznz -lnewmat  -lutils -lz 
SWOBJS= first.o
XFILES= first boundaryCorr first_utils

SCRIPTS = run_first run_first_sh_prior run_first_all first_flirt first_bcorr

all: ${XFILES}
	@if [ ! -d ${DESTDIR}/etc/standard ] ; then ${MKDIR} -p ${DESTDIR}/etc/standard ; ${CHMOD} -R g+w ${DESTDIR}/etc ; fi
	${CP} -r data/first_models_317_bin ${DESTDIR}/etc/standard 
	${CP} data/*.nii.gz ${DESTDIR}/etc/standard

first: ${SWOBJS} 
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${SWOBJS} ${LIBS}

boundaryCorr: boundaryCorr.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ boundaryCorr.o ${LIBS}

first_utils: first_utils.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ first_utils.o ${LIBS}
