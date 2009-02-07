include ${FSLCONFDIR}/default.mk

PROJNAME = bint

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_ZLIB}

LIBS = -lutils -lnewimage -lmiscmaths -lm -lnewmat -lfslio -lniftiio -lznz -lprob -lz

XFILES = mean

OBJS = model.o lsmcmcmanager.o lslaplacemanager.o bintoptions.o meanoptions.o

#OPTFLAGS =

all:	${XFILES} libbint

libbint: ${OBJS}
	${CXX} -shared -o libbint.so ${OBJS} ${LDFLAGS} -fPIC -lutils -lnewmat -lnewimage -lmiscmaths -lprob

mean:${OBJS} mean.o
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} mean.o ${LIBS}


