# $Id: Makefile,v 1.5 2008/11/19 19:09:31 brian Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = fslvtkio

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB} -I${HOME}/fsldev/lib
USRLDFLAGS = -L${LIB_NEWMAT}  -L${LIB_ZLIB} 

LIBS= -lmeshclass -lnewimage -lmiscmaths -lfslio -lniftiio -lznz -lnewmat -lutils -lz

SHAPEOBJS=fslvtkio.o
 
MODELOBJS=fslvtkio.o


all: fslvtkio.a 

fslvtkio.a: ${SHAPEOBJS}
	${AR} -r libfslvtkio.a ${SHAPEOBJS}




