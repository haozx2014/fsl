# $Id: Makefile,v 1.21 2007/07/21 01:04:53 brian Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = first

USRINCFLAGS = -I${INC_NEWMAT}

USRLDFLAGS = -L${LIB_NEWMAT} 


LIBS= -lshapeModel  -lmeshclass -lnewimage -lprob -lmiscmaths -lfslio -lniftiio -lznz -lnewmat  -lutils -lz 
SWOBJS= first.o
XFILES= first first_utils

SCRIPTS = run_first run_first_intref run_first_sh_prior run_first_all first_flirt concat_bvars

all: ${XFILES}
	@if [ ! -d ${DESTDIR}/etc/standard ] ; then ${MKDIR} -p ${DESTDIR}/etc/standard ; ${CHMOD} -R g+w ${DESTDIR}/etc ; fi
	${CP} -r data/first_models_317_bin ${DESTDIR}/etc/standard 
	${CP} data/*.nii.gz ${DESTDIR}/etc/standard

first: ${SWOBJS} 
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${SWOBJS} ${LIBS}

first_utils: first_utils.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ first_utils.o ${LIBS}
