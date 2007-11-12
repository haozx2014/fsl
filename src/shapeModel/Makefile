# $Id: Makefile,v 1.4 2006/11/06 14:14:23 brian Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = shapeModel

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB} 
USRLDFLAGS = -L${LIB_NEWMAT}  -L${LIB_ZLIB} 

LIBS= -lmeshclass -lnewimage -lmiscmaths -lfslio -lniftiio -lznz -lnewmat -lutils -lz

SHAPEOBJS=mshape.o shapeModel.o
 
MODELOBJS=mshape.o shapeModel.o

DBGFLAGS=-g

all: shapeModel

shapeModel: ${SHAPEOBJS}
	${CXX} -shared -o libshapeModel.so ${SHAPEOBJS} ${LDFLAGS} -fPIC -lmeshclass -lnewimage -lnewmat




