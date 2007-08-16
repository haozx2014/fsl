# $Id: Makefile,v 1.24 2007/08/12 22:32:40 brian Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = first

USRINCFLAGS = -I${INC_NEWMAT}

USRLDFLAGS = -L${LIB_NEWMAT} 


LIBS= -lshapeModel  -lmeshclass -lnewimage -lprob -lmiscmaths -lfslio -lniftiio -lznz -lnewmat  -lutils -lz 
SWOBJS= first.o
XFILES= first first_utils

SCRIPTS = run_first run_first_all first_flirt concat_bvars first_roi_slicesdir

all: ${XFILES}

first: ${SWOBJS} 
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${SWOBJS} ${LIBS}

first_utils: first_utils.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ first_utils.o ${LIBS}
