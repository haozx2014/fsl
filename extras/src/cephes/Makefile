include ${FSLCONFDIR}/default.mk

PROJNAME = cephes

OBJS = zetac.o j0.o i0.o chbevl.o j1.o i1.o

#OBJS = airy.o beta.o chbevl.o clog.o cmplx.o \
#dawsn.o ei.o ellie.o ellik.o ellpe.o ellpj.o ellpk.o \
#exp10.o exp2.o expn.o fac.o fresnl.o hyp2f1.o hyperg.o i0.o i1.o iv.o \
#jv.o k0.o k1.o kn.o kolmogorov.o log2.o lrand.o planck.o \
#polmisc.o polylog.o polyn.o powi.o psi.o rgamma.o round.o \
#shichi.o sici.o sindg.o spence.o struve.o tandg.o zeta.o zetac.o

# Non-POSIX functions that nevertheless exist in math.h on mace (and others)
# cbrt.o isnan.o j0.o j1.o jn.o yn.o 

LIBS = -lm -lcprob

MOBJS = mtst.o libcephes.a
COBJS = ctest.o libcephes.a

USRINCFLAGS = -I${INC_CPROB} -I../cprob

XFILES = 
AFILES = libcephes.a
HFILES = libcephes.h
TESTXFILES = mtst ctest

all: libcephes.a

libcephes.a: ${OBJS}
	${AR} -r libcephes.a ${OBJS}

mtst: all ${MOBJS}
	${CC} ${CFLAGS} ${LDFLAGS} -o $@ ${MOBJS} ${LIBS} 

ctest: all ${COBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${COBJS} ${LIBS} 

