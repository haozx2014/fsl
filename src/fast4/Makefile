# A Makefile for the fast_two project

include ${FSLCONFDIR}/default.mk

PROJNAME = fast4

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB}

LIBS = -lnewimage -lmiscmaths -lprob -lfslio -lniftiio -lutils -lnewmat -lznz -lm -lz

 
SOBJS = ${NOBJS} fast_two.o mriseg_two.o multi_mriseg_two.o



XFILES = fast4
all: ${XFILES}


fast4: ${SOBJS}
	$(CXX)  ${CXXFLAGS} ${LDFLAGS} -o $@ ${SOBJS} ${LIBS}


