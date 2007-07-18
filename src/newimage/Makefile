# A Makefile for the image class

include ${FSLCONFDIR}/default.mk

PROJNAME = newimage

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB}

LIBS = -lmiscmaths -lfslio -lniftiio -lutils -lnewmat -lznz -lm -lz

OBJS = lazy.o newimage.o generalio.o newimagefns.o complexvolume.o imfft.o warpfns.o
NOBJS = ${OBJS} costfns.o

TOBJS = ${NOBJS} testprog.o  
SOBJS = ${NOBJS} skeletonapp.o  

TESTXFILES = testprog skeletonapp
XFILES = 

all: libnewimage

libnewimage: ${NOBJS}
	${CXX} -shared -o libnewimage.so ${NOBJS} ${LDFLAGS} -fPIC -lnewmat -lfslio -lmiscmaths -lniftiio

skeletonapp: ${SOBJS}
	$(CXX)  ${CXXFLAGS} ${LDFLAGS} -o $@ ${SOBJS} ${LIBS}

testprog: ${TOBJS}
	$(CXX)  ${CXXFLAGS} ${LDFLAGS} -o $@ ${TOBJS} ${LIBS}

