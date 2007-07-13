#

include ${FSLCONFDIR}/default.mk

PROJNAME = misc_tcl

SCRIPTS = runfsl

TCLFILES = *.tcl *xbm *.sh *.ppm *.ico fslwebbrowser

all:
	@if [ ! -d ${DESTDIR}/bin ] ; then \
	  ${MKDIR} ${DESTDIR}/bin ; \
	  ${CHMOD} g+w ${DESTDIR}/bin ; \
	fi
	${INSTALL} -m 0775 Runtcl ${DESTDIR}/bin ; \
        ${CP} -r BWidget ${DESTDIR}/tcl/ ; \
	${INSTALL} -m 0775 fsl.bat FSL.lnk ${DESTDIR} ; \
	cd ${DESTDIR}/bin ; \
	${RM} -f fsl fsl_gui ; \
	ln -s Runtcl fsl ; \
	ln -s Runtcl fsl_gui 
