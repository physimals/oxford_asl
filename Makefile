include ${FSLCONFDIR}/default.mk

PROJNAME = oxford_asl

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB}


LIBS = -lutils -lnewimage -lmiscmaths -lm -lnewmat -lfslio -lniftiio -lznz -lz

XFILES = asl_file
SCRIPTS = oxford_asl asl_preproc asl_calib asl_reg quasil
RUNTCLS = Asl

OBJS = readoptions.o asl_functions.o

# Pass Git revision details
GIT_SHA1:=$(shell git describe --dirty)
GIT_DATE:=$(shell git log -1 --format=%ad --date=local)
CXXFLAGS += -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_DATE="\"${GIT_DATE}\""

# Always rebuild scripts
.PHONY: FORCE
FORCE:

all:	${XFILES} ${SCRIPTS}

asl_file: ${OBJS} asl_file.o 
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} asl_file.o ${LIBS}

$(SCRIPTS): %: %.in FORCE
	sed -e "s/\$${GIT_SHA1}/${GIT_SHA1}/" -e "s/\$${GIT_DATE}/${GIT_DATE}/" $< >$@

clean:
	rm ${SCRIPTS}

