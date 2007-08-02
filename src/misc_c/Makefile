include $(FSLCONFDIR)/default.mk

PROJNAME = misc_c

USRINCFLAGS = -I${INC_PROB} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_PROB} -L${LIB_ZLIB}

IOLIBS = -lfslio -lniftiio -lznz
#LIBS = -lss_8UI $(IOLIBS) -lprob -lm -lz

XFILES = whirlgif hist2prob morph_kernel ptoz runavw tmpnam ztop add4bytes

SCRIPTS = correctheader

all: $(XFILES)

hist2prob: hist2prob.c
	$(CC) $(CFLAGS) -o hist2prob hist2prob.c -L${LIB_ZLIB} -lm -lz

morph_kernel: morph_kernel.c
	$(CC) -o morph_kernel morph_kernel.c

ptoz: ptoz.cc
	$(CXX) $(CFLAGS) -o ptoz ptoz.cc $(LDFLAGS) -lprob -lm -lz

ztop: ztop.cc
	$(CXX) $(CFLAGS) -o ztop ztop.cc $(LDFLAGS) -lprob -lm -lz

runavw: runavw.c
	$(CC) $(CFLAGS) -o runavw runavw.c $(LDFLAGS) -lss_8UI $(IOLIBS) -lm -lz

sphere: sphere.c
	$(CC) $(CFLAGS) -o sphere sphere.c $(LDFLAGS) -lss_8UI $(IOLIBS) -lm -lz

tmpnam: tmpnam.c
	$(CC) $(CFLAGS) -o tmpnam tmpnam.c

whirlgif: gifdecod.c gifencod.c whirlgif.c
	$(CC) $(CFLAGS) -o whirlgif whirlgif.c gifdecod.c gifencod.c

add4bytes: add4bytes.cc
	$(CXX) $(CFLAGS) -o add4bytes add4bytes.cc 
