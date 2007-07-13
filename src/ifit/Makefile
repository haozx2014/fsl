include ${FSLCONFDIR}/default.mk

PROJNAME = ifit

OBJS = covsrt.o gaussj.o mrqcof.o mrqmin.o nrerror.o matrix.o \
        vector.o powell_lim.o linmin.o brent.o f1dim.o mnbrak.o

##Remove these comments to compile on cayenne
#CC = cc
#ANSI_FLAGS = -std1 
#USRCFLAGS = -ieee
#OPTFLAGS = 
#DEPENDFLAGS = -M
##

USRINCFLAGS = -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_ZLIB}

LIBS = -lfslio -lniftiio -lznz -lm -lz

XFILES = linfit lmfit minfit
HFILES =

all: linfit lmfit minfit mfilt

linfit: linfit.o ${OBJS}
	${CC}  ${CFLAGS} ${LDFLAGS} -o $@ linfit.o ${OBJS} ${LIBS}

lmfit: lmfit.o ${OBJS}
	${CC}  ${CFLAGS} ${LDFLAGS} -o $@ lmfit.o ${OBJS} ${LIBS}

minfit: minfit.o ${OBJS}
	${CC}  ${CFLAGS} ${LDFLAGS} -o $@ minfit.o ${OBJS} ${LIBS}

mfilt: mfilt.o
	${CC}  ${CFLAGS} ${LDFLAGS} -o $@ mfilt.o ${LIBS}
