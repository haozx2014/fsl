include ${FSLCONFDIR}/default.mk

PROJNAME = cprob

OBJS = bdtr.o btdtr.o chdtr.o drand.o expx2.o fdtr.o gamma.o gdtr.o \
igam.o igami.o incbet.o incbi.o mtherr.o nbdtr.o ndtr.o ndtri.o pdtr.o \
stdtr.o unity.o polevl.o const.o xmath.o

USRCFLAGS = -std=gnu99

all: libprob.a

dftest: all dftest.o
	$(CC) $(CFLAGS) -L. -o dftest dftest.o -lprob

libprob.a: $(OBJS) $(INCS)
	${CXX} -shared -o libprob.so $(OBJS) -fPIC ${LDFLAGS}
	rm -f libcprob.so
	ln -s libprob.so libcprob.so
