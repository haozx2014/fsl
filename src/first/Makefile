# $Id: Makefile,v 1.26 2009/03/20 10:17:02 brian Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = first

USRINCFLAGS = -I${INC_NEWMAT}

USRLDFLAGS = -L${LIB_NEWMAT} 


LIBS= -lshapeModel -lfirst_lib -lfslvtkio -lmeshclass -lnewimage -lprob -lmiscmaths -lfslio -lniftiio -lznz -lnewmat  -lutils -lz 

XFILES= first first_utils

SCRIPTS = run_first run_first_all first_flirt concat_bvars first_roi_slicesdir

all: ${XFILES}

first: first.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ first.o ${LIBS}

first_utils: first_utils.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ first_utils.o ${LIBS}
