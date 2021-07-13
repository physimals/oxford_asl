include ${FSLCONFDIR}/default.mk

PROJNAME = oxford_asl

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB}

FSLVERSION= $(shell cat ${FSLDIR}/etc/fslversion | head -c 1)
ifeq ($(FSLVERSION), 5) 
  NIFTILIB = -lfslio -lniftiio 
  MATLIB = -lnewmat
else 
  NIFTILIB = -lNewNifti
  UNAME := $(shell uname -s)
  ifeq ($(UNAME), Linux)
    MATLIB = -lopenblas
  endif
endif

LIBS = -lutils -lnewimage -lmiscmaths -lprob ${MATLIB} ${NIFTILIB} -lznz -lz

XFILES = asl_file
SCRIPTS = oxford_asl asl_calib asl_reg quasil toast oxford_asl_roi_stats.py
PYMODULES = python/oxford_asl/__init__.py python/oxford_asl/_version.py
PYGUI = python/oxford_asl/gui/*.py python/oxford_asl/gui/*.png
VERSIONED = oxford_asl asl_calib quasil asl_reg toast 

OBJS = readoptions.o asl_functions.o

# Pass Git revision details
GIT_SHA1:=$(shell git describe --dirty)
GIT_DATE:=$(shell git log -1 --format=%ad --date=local)
CXXFLAGS += -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_DATE="\"${GIT_DATE}\""

# Always rebuild scripts
.PHONY: FORCE

all:	${XFILES} ${VERSIONED}

asl_file: ${OBJS} asl_file.o 
	${CXX}  ${CXXFLAGS} ${LDFLAGS} -o $@ ${OBJS} asl_file.o ${LIBS}

$(VERSIONED): %: %.in FORCE
	sed -e "s/@GIT_SHA1@/${GIT_SHA1}/" -e "s/@GIT_DATE@/${GIT_DATE}/" $< >$@
	chmod a+x $@

postinstallscript: $(PYMODULES) $(PYGUI)
	mkdir -p $(DESTDIR)/python/oxford_asl/gui ; \
	cp $(PYMODULES) $(DESTDIR)/python/oxford_asl/ ; \
	cp $(PYGUI) $(DESTDIR)/python/oxford_asl/gui ; \
        cp asl_gui_fsl $(DESTDIR)/bin/asl_gui ; \
	cd ..

clean:
	rm -f ${VERSIONED} asl_file *.o

FORCE:
	python python/setup.py build
