# $Id: Makefile,v 1.32 2006/11/08 09:08:49 steve Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = fugue

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB}

LIBS=-lnewimage -lmiscmaths -lprob -lfslio -lniftiio -lznz -lnewmat -lutils -lz

PRELUDEOBJS=prelude.o unwarpfns.o
FUGUEOBJS=fugue.o unwarpfns.o
APPLYOBJS=applywarp.o
CONVERTOBJS=convertwarp.o
VOBJS=variation.o unwarpfns.o
IOBJS=invwarptest.o unwarpfns.o
SOBJS=sigloss.o
MOBJS=macfnirt.o

XFILES=prelude fugue applywarp convertwarp sigloss
TESTXFILES=variation invwarptest macflirt
SCRIPTS = extractfidparams dof2warp

all: ${XFILES}

prelude: ${PRELUDEOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${PRELUDEOBJS} ${LIBS} 

fugue:  ${FUGUEOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${FUGUEOBJS} ${LIBS} 

applywarp:  ${APPLYOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${APPLYOBJS} ${LIBS} 

convertwarp:  ${CONVERTOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${CONVERTOBJS} ${LIBS} 

variation: ${VOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${VOBJS} ${LIBS} 

invwarptest: ${IOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${IOBJS} ${LIBS} 

sigloss: ${SOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${SOBJS} ${LIBS} 

macfnirt: ${MOBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ ${MOBJS} ${LIBS} 
