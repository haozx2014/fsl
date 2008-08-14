# $Id: Makefile,v 1.5 2008/08/10 09:21:34 brian Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = shapeModel

SHAPEOBJS=shapeModel.o
 

all: shapeModel

shapeModel: ${SHAPEOBJS}
	${CXX} -shared -o libshapeModel.so ${SHAPEOBJS} ${LDFLAGS} -fPIC -lmeshclass -lnewimage -lnewmat




