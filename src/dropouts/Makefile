#

include $(FSLCONFDIR)/default.mk

PROJNAME = dropouts

USRINCFLAGS = -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_ZLIB}

LIBS = -lss_$$FDT -lfslio -lniftiio -lznz -lm -lz

XFILES = dropouts_8UI dropouts_8SI dropouts_16UI dropouts_16SI dropouts_32UI dropouts_32SI dropouts_32R

RUNAVWS = dropouts

SCRIPTS = oddvolumes

all:	dropouts_8UI

dropouts_8UI:	dropouts.c
	for FDT in 8UI 8SI 16UI 16SI 32UI 32SI 32R ; do $(DATATYPES) ; \
        $(CC) $(CFLAGS) -DFDT="$$STR" -o dropouts_$$FDT dropouts.c $(LDFLAGS) $(LIBS) ; \
	done
