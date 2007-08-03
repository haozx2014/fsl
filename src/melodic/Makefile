# A Makefile for MELODIC

include ${FSLCONFDIR}/default.mk

OPTFLAGS = -O3 -Wno-deprecated 
OPTFLAGS_alphaev6-dec-osf 5.0-gcc2.95.2 = -O3 -mieee -mfp-trap-mode=sui

PROJNAME = melodic

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_GD} -I${INC_GDC} -I${INC_PNG} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB} -L${LIB_GD} -L${LIB_GDC} -L${LIB_PNG} -L${LIB_ZLIB}

LIBS = -lutils -lnewimage -lmiscplot -lmiscpic -lmiscmaths -lfslio -lniftiio -lznz -lnewmat -lprob -lm  -lgdc -lgd -lpng -lz

TEST_OBJS = melhlprfns.o test.o  

GGMIX_OBJS = ggmix.o

FSL_GLM_OBJS = melhlprfns.o fsl_glm.o

FSL_REGFILT_OBJS = melhlprfns.o fsl_regfilt.o

MELODIC_OBJS =  meloptions.o melhlprfns.o melgmix.o meldata.o melpca.o melica.o melreport.o melodic.o 

TESTXFILES = test

XFILES = fsl_glm fsl_regfilt melodic

RUNTCLS = Melodic

SCRIPTS = melodicreport

all: ggmix fsl_regfilt fsl_glm melodic

ggmix: ${GGMIX_OBJS}
	${CXX} -shared -o libggmix.so ${GGMIX_OBJS} ${LDFLAGS} -fPIC -lnewmat -lmiscmaths -lnewimage

melodic: ${MELODIC_OBJS}
	$(CXX) ${CXXFLAGS} ${LDFLAGS} -o $@ ${MELODIC_OBJS} ${LIBS}

test: ${TEST_OBJS}
	$(CXX) ${CXXFLAGS} ${LDFLAGS} -o $@ ${TEST_OBJS} ${LIBS}

fsl_glm: ${FSL_GLM_OBJS}
	$(CXX) ${CXXFLAGS} ${LDFLAGS} -o $@ ${FSL_GLM_OBJS} ${LIBS}

fsl_regfilt: ${FSL_REGFILT_OBJS}
	$(CXX) ${CXXFLAGS} ${LDFLAGS} -o $@ ${FSL_REGFILT_OBJS} ${LIBS}
