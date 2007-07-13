include ${FSLCONFDIR}/default.mk

PROJNAME = mcflirt

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_ZLIB}

MOBJS = mcflirt.o Globaloptions.o Log.o
VOBJS = vmapper.o Globaloptions.o
DOBJS = displacement.o Log.o
ROBJS = reform.o
R2OBJS = reform2.o
TOBJS = testprog.o Log.o
NOBJS = nullmotion.o Log.o
N2OBJS = nulldata.o
RMASK = resize_mask.o

LIBS = -lnewimage -lmiscmaths -lutils -lfslio -lniftiio -lznz -lnewmat -lprob -lm -lz

XFILES = mcflirt vmapper
TESTXFILES = testprog displacement reform reform2 nullmotion nulldata resize_mask

#OPTFLAGS = 

all:	${XFILES}

test:	${TESTXFILES}

mcflirt:${MOBJS}
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${MOBJS} ${LIBS}

vmapper:	${VOBJS}
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${VOBJS} ${LIBS}

displacement:	${DOBJS}
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${DOBJS} ${LIBS}

reform: ${ROBJS}
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${ROBJS} ${LIBS}

reform2: ${R2OBJS}
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${R2OBJS} ${LIBS}

testprog: ${TOBJS}
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${TOBJS} ${LIBS}

nullmotion:  ${NOBJS}
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${NOBJS} ${LIBS}

nulldata:  ${N2OBJS}
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${N2OBJS} ${LIBS}

resize_mask:  ${RMASK}
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@  ${RMASK} ${LIBS}
