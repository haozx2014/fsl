# $Id: Makefile,v 1.16 2012/07/25 15:56:33 cowboy Exp $
include ${FSLCONFDIR}/default.mk

PROJNAME = fslsurface

LINUX_LIBS=-lGL -lGLU 

ifeq ($(FSLMACHTYPE),apple-darwin8-gcc4.0)
        ARCHFLAGS =  -arch i386 -isysroot /Developer/SDKs/MacOSX10.4u.sdk -I/usr/X11R6/include/
        ARCHLDFLAGS = -Wl,-search_paths_first -arch i386 -isysroot /Developer/SDKs/MacOSX10.4u.sdk -L/Developer/SDKs/MacOSX10.4u.sdk/usr/X11R6/lib/ -framework OpenGL -framework GLUT
	LINUX_LIBS=""
endif 

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_ZLIB} -I${FSLDIR}/include/niftiio 

USRLDFLAGS = -L${LIB_NEWMAT}  -L${INC_PROB} -L${LIB_ZLIB} -L./ 


LIBS=-lgiftiio -lexpat -lfirst_lib -lmeshclass
FMLIBS= -lshapeModel -lfslvtkio -lmeshclass -lnewimage -lmiscmaths  -lfslio -lniftiio -lznz -lnewmat -lutils -lprob -lz

APPLY_ZLIB = -DHAVE_ZLIB


SHAPEOBJS=fslsurface.o 
MODELOBJS=fslsurface.o

BRIANLIBFILES=fslsurfacegl.a
LIBFILES=fslsurface.a fslsurface_backcompat.a
XFILES=fslsurfacemaths

all:  ${LIBFILES} ${XFILES}

fslsurface.a: fslsurfacefns.o fslsurfaceio.o fslsurface.o 
	${AR} -r libfslsurface.a  fslsurfacefns.o fslsurfaceio.o fslsurface.o
fslsurfacegl.a: fslsurfacegl.o fslsurface.a
	${AR} -r libfslsurfacegl.a fslsurfacegl.o
fslsurface_backcompat.a: fslsurface_first.o fslsurface_dataconv.o fslsurface.a
	${AR} -r libfslsurface_backcompat.a fslsurface_first.o fslsurface_dataconv.o
fslsurfacemaths: fslsurface.a fslsurfacemaths.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslsurfacemaths.o fslsurfaceio.o fslsurfacefns.o fslsurface_first.o fslsurface.o ${LIBS} ${FMLIBS} ${LINUX_LIBS}
marching_cubes: fslsurface.a marching_cubes.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ marching_cubes.o fslsurfaceio.o fslsurfacefns.o ${LIBS} ${FMLIBS} ${LINUX_LIBS}
