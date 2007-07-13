include $(FSLCONFDIR)/default.mk

PROJNAME = feat5

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_ZLIB} -I${INC_GD} -I${INC_GDC} -I${INC_PNG} 
USRLDFLAGS =  -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_ZLIB} -L${LIB_GD} -L${LIB_GDC} -L${LIB_PNG}


#USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_GD} -I${INC_GDC} -I${INC_PNG} -I${INC_ZLIB}
#USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_GD} -L${LIB_GDC} -L${LIB_PNG} -L${LIB_ZLIB}
#LIBS = -lnewimage -lmiscplot -lmiscpic -lmiscmaths -lfslio -lniftiio -lznz -lnewmat -lprob -lm  -lgdc -lgd -lpng -lz


IOLIBS = -lfslio -lniftiio -lznz
LIBS = -lutils -lnewimage -lmiscplot -lmiscpic -lmiscmaths -lprob -lnewmat -lss_32R $(IOLIBS) -lm -lz  -lgdc -lgd -lpng

XFILES = feat_model tsplot prewhiten

RUNTCLS = Feat Renderhighres Featquery Featwatcher Glm

SCRIPTS = easythresh feat mccutup renderhighres featregapply featquery fslerrorreport updatefeatreg estnoise perfusion_subtract

all: feat_model tsplot prewhiten defaultfiles

defaultfiles:
	sed s/VVV/`grep "set fmri(version)" featlib.tcl | grep "\." | awk '{print $$3}'`/g doc/ORIGtop.html > doc/top.html
	@if [ ! -d ${DESTDIR}/etc ] ; then ${MKDIR} ${DESTDIR}/etc ; ${CHMOD} g+w ${DESTDIR}/etc ; fi
	${CP} -rf default_flobs.flobs ${DESTDIR}/etc
	${CP} -f featconf doc/feat.tcl
	@if [ ! -d ${DESTDIR}/etc/fslconf ] ; then ${MKDIR} ${DESTDIR}/etc/fslconf ; ${CHMOD} g+w ${DESTDIR}/etc/fslconf ; fi
	${CP} -f featconf ${DESTDIR}/etc/fslconf/feat.tcl

featlib.o: featlib.cc featlib.h
	$(CXX) $(CXXFLAGS) -DFDT="float" -c -o featlib.o featlib.cc

tsplot: tsplot.cc featlib.o featlib.h
	$(CXX) $(CXXFLAGS) -DFDT="float" -o tsplot featlib.o tsplot.cc $(LDFLAGS) $(LIBS)

prewhiten: prewhiten.cc featlib.o featlib.h
	$(CXX) $(CXXFLAGS) -DFDT="float" -o prewhiten featlib.o prewhiten.cc $(LDFLAGS) $(LIBS)

feat_model: feat_model.cc featlib.o
	$(CXX) $(CXXFLAGS) -DFDT="float" -o feat_model featlib.o feat_model.cc $(LDFLAGS) $(LIBS)
