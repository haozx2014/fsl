
# for SHELL, do not change the type of shell - only use Bourne or BASH
SHELL = /bin/sh

# System dependent commands (NB: the first two are the most platform dependent)

INSTALL = install -p
RANLIB = ranlib

RM = /bin/rm
CP = /bin/cp
CHMOD = /bin/chmod
MKDIR = /bin/mkdir
TCLSH = /usr/bin/tclsh

# Compiler dependent variables

CC = gcc
CXX = g++
CSTATICFLAGS = 
CXXSTATICFLAGS = 

ARCHFLAGS = ${DEB_ARCH_OPT_FLAGS}

DEPENDFLAGS = -MM

OPTFLAGS =  -O3 -fexpensive-optimizations ${ARCHFLAGS}
MACHDBGFLAGS =
GNU_ANSI_FLAGS = -Wall -ansi -pedantic
SGI_ANSI_FLAGS = -ansi -fullwarn
ANSI_FLAGS = ${GNU_ANSI_FLAGS} -I/usr/include/nifti -I/usr/include/newmat \
             -I${FSLDIR}/extras/include -I${FSLDIR}/include/libprob \
             -DHAVE_LIBFREETYPE -fPIC
#ARCHLDFLAGS = -Wl,--no-undefined -Wl,--as-needed -Wl,--rpath ${FSLDIR}/lib -Wl,--rpath /usr/lib/fsl
ARCHLDFLAGS = -Wl,--no-undefined -Wl,--as-needed
