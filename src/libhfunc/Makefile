include ${FSLCONFDIR}/default.mk

PROJNAME = libhfunc

AFILES = libhfunc.so
HFILES = libhfunc.h dbh.h
XFILES = 

OBJS = libhfunc.o

all: libhfunc

libhfunc: ${OBJS}
	${CXX} -shared -o libhfunc.so ${OBJS} ${LDFLAGS} -fPIC
