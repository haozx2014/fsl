include ${FSLCONFDIR}/default.mk

PROJNAME = fslio

USRINCFLAGS = -I${INC_ZLIB}

OBJS=fslio.o

SCRIPTS = remove_ext fsloutputtype imtest imglob imcp imln imrm immv

all: libfslio

libfslio: ${OBJS}
	${CXX} -shared -o libfslio.so $(OBJS) -fPIC ${LDFLAGS} -lniftiio -lznz -lz

TESTOBJS = exttest.o
TESTLIBS = -lfslio -lniftiio -lznz -lz

test: all exttest

exttest: ${TESTOBJS}
	${CC} ${LDFLAGS} -o exttest ${TESTOBJS} ${TESTLIBS}

