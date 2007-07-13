include ${FSLCONFDIR}/default.mk

PROJNAME = load_varian

OBJS = main.o process.o hfunc.o phasefunc.o procfunc.o read_varian.o image_ft.o scaleslice.o pssreorder.o segment.o byteorder.o reslice.o kspace_mask.o pocs.o hermit.o margosian.o image_fftw.o tep.o dropouts.o fermi.o entghost.o navigator.o resamp.o concat_info.o fdf.o

USRINCFLAGS = -I${INC_FFTW3} -I${INC_NEWMAT} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_FFTW3} -L${LIB_NEWMAT} -L${LIB_ZLIB}

LIBS = -lfftw3 -lmiscmaths -lfslio -lniftiio -lznz -lnewmat -lm -lz

XFILES = load_varian_exe set_hdr_info get_orient get_dim
SCRIPTS = load_varian
TFILES = *.tcl *.scale
HFILES = 

RUNTCLS = LoadVarian MRinfo

all: ${XFILES}

load_varian_exe: ${OBJS}
	${CXX} ${LDFLAGS} -o $@ ${OBJS} ${LIBS}

set_hdr_info: set_hdr_info.o
	${CXX} ${LDFLAGS} -o $@ set_hdr_info.o hfunc.o ${LIBS}

get_orient: get_orient.o
	${CXX} ${LDFLAGS} -o $@ get_orient.o hfunc.o ${LIBS}

get_dim: get_dim.o byteorder.o
	${CXX} ${LDFLAGS} -o $@ get_dim.o byteorder.o ${LIBS}

avw2fdf: avw2fdf.o
	${CXX} ${LDFLAGS} -o $@ avw2fdf.o -lfslio -lznz -lniftiio -lm -lz
