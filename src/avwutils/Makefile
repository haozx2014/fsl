
include ${FSLCONFDIR}/default.mk

PROJNAME = avwutils

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB}

LIBS = -lnewimage -lmiscmaths -lprob -lfslio -lnewmat -lutils -lniftiio -lznz -lm -lz

IOLIBS = -lfslio -lniftiio -lznz -lm -lz
SSLIBS = -lss_32R $(IOLIBS)

XFILES = fslorient fslstats fslmerge \
         fslpspec fslroi fslnvols fsl2ascii  \
         fslsplit fslcc fslinterleave \
         fslhd fslcpgeom fslcreatehd fslmaths \
         fslcomplex fslfft fslmeants fslcheck fslslice fslswapdim

SCRIPTS = fslval fslchpixdim fslanimate fslsize sliceanimate fslinfo fsledithd fslchfiletype avw2fsl

FSCRIPTS = fmriborient fslswapdim2

all: ${XFILES}

fslcheck: fslcheck.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o fslcheck fslcheck.cc ${LIBS}

avwconv: avwconv.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o avwconv avwconv.cc ${LIBS}

fslstats: fslstats.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslstats.o ${LIBS} 

fslfft: fslfft.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o fslfft fslfft.cc ${LIBS}

fslmeants: fslmeants.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o fslmeants fslmeants.cc ${LIBS}

avwfixfloat: avwfixfloat.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o avwfixfloat avwfixfloat.cc ${LIBS}

fslcomplex: fslcomplex.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o fslcomplex fslcomplex.cc ${LIBS}

avwnvols: avwnvols.c
	${CC} ${CFLAGS} ${LDFLAGS} -o avwnvols avwnvols.c ${LIBS}

avwsplit: avwsplit.c
	${CC} ${CFLAGS} ${LDFLAGS} -o avwsplit avwsplit.c ${LIBS}

avwhd: avwhd.c
	$(CC) $(CFLAGS) -o avwhd avwhd.c $(LDFLAGS) ${LIBS}

avwcreatehd: avwcreatehd.c
	$(CC) $(CFLAGS) -o avwcreatehd avwcreatehd.c $(LDFLAGS) ${SSLIBS}

avwcpgeom: avwcpgeom.c
	$(CC) $(CFLAGS) -o avwcpgeom avwcpgeom.c $(LDFLAGS) ${SSLIBS}

avwstats: avwstats.c
	$(CC) $(CFLAGS) -DFDT="float" -o avwstats avwstats.c $(LDFLAGS) ${SSLIBS}

avwroi_8UI: avwroi.c
	for FDT in 8UI 16SI 32SI 32R 64R ; do $(DATATYPES) ; \
	$(CC) $(CFLAGS) -DFDT="$$STR" -o avwroi_$$FDT avwroi.c $(LDFLAGS) -lss_$$FDT $(IOLIBS); \
	done

avwmaths_8UI: avwmaths.c
	for FDT in 8UI 16SI 32SI 32R 64R ; do $(DATATYPES) ; \
	$(CC) $(CFLAGS) -DFDT="$$STR" -o avwmaths_$$FDT avwmaths.c $(LDFLAGS) -lss_$$FDT $(IOLIBS); \
	done

avwcc_8UI: avwcc.c
	for FDT in 8UI 16SI 32SI 32R 64R ; do $(DATATYPES) ; \
	$(CC) $(CFLAGS) -DFDT="$$STR" -o avwcc_$$FDT avwcc.c $(LDFLAGS) -lss_$$FDT $(IOLIBS); \
	done

avwinterleave_8UI: avwinterleave.c
	for FDT in 8UI 16SI 32SI 32R 64R ; do $(DATATYPES) ; \
	$(CC) $(CFLAGS) -DFDT="$$STR" -o avwinterleave_$$FDT avwinterleave.c $(LDFLAGS) -lss_$$FDT $(IOLIBS); \
	done

avw2ascii_8UI: avw2ascii.c
	for FDT in 8UI 16SI 32SI 32R 64R ; do $(DATATYPES) ; \
	$(CC) $(CFLAGS) -DFDT="$$STR" -o avw2ascii_$$FDT avw2ascii.c $(LDFLAGS) -lss_$$FDT $(IOLIBS); \
	done

fslswapdim: fslswapdim.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslswapdim.o ${LIBS} 

fslorient: fslorient.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslorient.o ${LIBS} 

fslslice: fslslice.cc
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o fslslice fslslice.cc ${LIBS}

fslpspec: fslpspec.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslpspec.o ${LIBS} 


#start of matts c++ rewrites
fslroi: fslroi.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslroi.o ${LIBS}

fslnvols: fslnvols.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslnvols.o ${LIBS}

fsl2ascii: fsl2ascii.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fsl2ascii.o ${LIBS}

fslsplit: fslsplit.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslsplit.o ${LIBS}

fslmerge: fslmerge.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslmerge.o ${LIBS}

fslcc: fslcc.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslcc.o ${LIBS}

fslinterleave: fslinterleave.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslinterleave.o ${LIBS}

fslhd: fslhd.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslhd.o ${LIBS}

fslcpgeom: fslcpgeom.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslcpgeom.o ${LIBS}

fslcreatehd: fslcreatehd.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslcreatehd.o ${LIBS}

fslmaths: fslmaths.o
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ fslmaths.o ${LIBS}

