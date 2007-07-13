include ${FSLCONFDIR}/default.mk

PROJNAME = istats

OBJS = istats.o statsfunc.o

LIBS = -lavwio -lm

XFILES = istats
HFILES =

RUNTCLS = Istats

all: ${XFILES}

istats: ${OBJS}
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} ${LIBS}
