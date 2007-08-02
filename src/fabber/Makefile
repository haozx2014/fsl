include ${FSLCONFDIR}/default.mk

PROJNAME = fabber

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB}

LIBS = -lutils -lmiscmaths -lnewmat -lfslio -lniftiio -lznz -lz 

XFILES = fabber
OBJS = fabber.o fwdmodel_q2tips.o inference_spatialvb.o dataset.o inference_vb.o noisemodel.o noisemodel_white.o fwdmodel_quipss2.o fwdmodel.o fwdmodel_simple.o fwdmodel_linear.o noisemodel_ar.o inference.o dist_mvn.o easylog.o easyoptions.o fwdmodel_asl_eagle.o fwdmodel_asl_grase.o fwdmodel_asl_buxton.o

# For debugging:
#OPTFLAGS = -ggdb
OPTFLAGS =

all:	${XFILES}

fabber: ${OBJS} 
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} ${LIBS}

# DO NOT DELETE
