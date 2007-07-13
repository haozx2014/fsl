include $(FSLCONFDIR)/default.mk

USRINCFLAGS = -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_ZLIB}

PROJNAME = bet

LIBS = -lss_$$FDT -lfslio -lniftiio -lznz -lm -lz

XFILES = old_bet_8UI old_bet_16SI old_bet_32SI old_bet_32R

RUNAVWS = old_bet

SCRIPTS = old_betfunc old_betpremask old_betsmallz old_bet_robust

all:	old_bet_8UI

old_bet_8UI: bet.c
	for FDT in 8UI 16SI 32SI 32R ; do $(DATATYPES) ; \
        $(CC) $(CFLAGS) -DFDT="$$STR" -o old_bet_$$FDT bet.c $(LDFLAGS) $(LIBS) ; \
	done
