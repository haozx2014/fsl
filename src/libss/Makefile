#

include $(FSLCONFDIR)/default.mk

PROJNAME = libss

all: libss_8UI.a

libss_8UI.a: libss_8UI.o libavw_8UI.o libtessa_8UI.o libpicio_8UI.o
	for FDT in 8UI 8SI 16UI 16SI 32UI 32SI 32R 64R ; do $(DATATYPES) ; \
	$(AR) -r libss_$$FDT.a libss_$$FDT.o libavw_$$FDT.o libtessa_$$FDT.o libpicio_$$FDT.o ; \
	done

libss_8UI.o: libss.c libss.h libssbase.h 
	for FDT in 8UI 8SI 16UI 16SI 32UI 32SI 32R 64R ; do $(DATATYPES) ; \
	$(CC) $(CFLAGS) -I ${INC_ZLIB} -DFDT="$$STR" -c -o libss_$$FDT.o libss.c ; \
	done

libavw_8UI.o: libavw.c libavw.h libss.h libssbase.h 
	for FDT in 8UI 8SI 16UI 16SI 32UI 32SI 32R 64R ; do $(DATATYPES) ; \
	$(CC) $(CFLAGS) -I ${INC_ZLIB} -DFDT="$$STR" -c -o libavw_$$FDT.o libavw.c ; \
	done

libtessa_8UI.o:	libtessa.c libtessa.h libavw.h libss.h libssbase.h 
	for FDT in 8UI 8SI 16UI 16SI 32UI 32SI 32R 64R ; do $(DATATYPES) ; \
	$(CC) $(CFLAGS) -I ${INC_ZLIB} -DFDT="$$STR" -c -o libtessa_$$FDT.o libtessa.c ; \
	done

libpicio_8UI.o:	libpicio.c libpicio.h libss.h libssbase.h 
	for FDT in 8UI 8SI 16UI 16SI 32UI 32SI 32R 64R ; do $(DATATYPES) ; \
	$(CC) $(CFLAGS) -I ${INC_ZLIB} -DFDT="$$STR" -c -o libpicio_$$FDT.o libpicio.c ; \
	done
