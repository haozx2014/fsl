include ${FSLCONFDIR}/default.mk

PROJNAME = utils

USRINCFLAGS = -I${INC_NEWMAT}
USRLIBFLAGS = -L${LIB_NEWMAT}

USRCXXFLAGS = 

OptionsTestTarget = opttst
OptionsTestObjects = opttst.o
OptionsTestLibrarys = -lutils

LibraryTarget = libutils.a
LibraryObjects = matches.o functions.o usage.o check.o parse.o log.o time_tracer.o

all: ${LibraryTarget}
	rm -rf ${FSLDIR}/etc/standard
	(cd ${FSLDIR}/etc; ln -s ../data/standard .)
	mkdir -p ${FSLDIR}/lib
	rm -rf ${FSLDIR}/lib/atlases
	(cd ${FSLDIR}/lib; ln -s ../data/atlases .)

${OptionsTestTarget}: ${LibraryTarget} ${OptionsTestObjects}
	${CXX} ${LDFLAGS} -o $@  ${OptionsTestObjects} \
		${OptionsTestLibrarys}

${LibraryTarget}: ${LibraryObjects}
	${AR} -r $@ ${LibraryObjects}
	${RANLIB} ${LibraryTarget}

XFILES = 
AFILES = ${LibraryTarget}

