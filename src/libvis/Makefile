# A Makefile for miscplot

include ${FSLCONFDIR}/default.mk

OPTFLAG = -O3

PROJNAME = libvis

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_GD} -I${INC_GDC} -I${INC_PNG} -I${INC_ZLIB} -I${INC_BOOST} -fPIC
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_GD} -L${LIB_GDC} -L${LIB_PNG} -L${LIB_ZLIB}

LIBS = -lnewimage -lmiscplot -lmiscpic -lmiscmaths -lfslio -lniftiio -lznz -lnewmat -lprob -lm  -lgdc -lgd -lpng -lz

PLOTOBJS =  miscplot.o  

PICOBJS  =  miscpic.o

TOBJS = test.o  

TESTXFILES = test

all: miscplot miscpic

miscplot: ${PLOTOBJS}
	${CXX} -shared -o libmiscplot.so ${PLOTOBJS} ${LDFLAGS} -lnewmat -lgd -lgdc -lmiscmaths -lfslio -lprob -lniftiio

miscpic: ${PICOBJS}
	${CXX} -shared -o libmiscpic.so ${PICOBJS} ${LDFLAGS} -lnewimage -lgd -lfslio -lmiscmaths -lnewmat -lniftiio

test: ${TOBJS}
	$(CXX) ${CXXFLAGS} ${LDFLAGS} -o $@ ${TOBJS}
