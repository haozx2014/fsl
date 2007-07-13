include ${FSLCONFDIR}/default.mk

PROJNAME = flame

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_ZLIB}

LIBS = -lutils -lnewimage -lmiscmaths -lm -lnewmat -lfslio -lniftiio -lznz -lprob -lz

XFILES = flame

OBJS = design.o mcmc_mh.o gsoptions.o gsmanager.o

#OPTFLAGS =

all:	${XFILES}

flame:${OBJS} flame.o
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} flame.o ${LIBS}

testgam:${OBJS} testgam.o
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} testgam.o ${LIBS}

