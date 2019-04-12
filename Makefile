include ${FSLCONFDIR}/default.mk

PROJNAME = oxford_asl

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB}

FSLVERSION= $(shell cat ${FSLDIR}/etc/fslversion | head -c 1)
ifeq ($(FSLVERSION), 5) 
  NIFTILIB = -lfslio -lniftiio 
  MATLIB = -lnewmat
else 
  UNAME := $(shell uname -s)
  ifeq ($(UNAME), Linux)
    MATLIB = -lopenblas
  endif
  NIFTILIB = -lNewNifti
endif

LIBS = -lutils -lnewimage -lmiscmaths -lprob ${MATLIB} ${NIFTILIB} -lznz -lz

XFILES = asl_file
SCRIPTS = oxford_asl asl_calib asl_reg quasil 
PYMODULES = python/asl/__init__.py python/asl/fslhelpers.py python/asl/reg.py python/asl/fslwrap.py python/asl/image.py
PYGUI = python/asl/gui/*.py python/asl/gui/banner.png
VERSIONED = oxford_asl asl_calib quasil asl_reg python/asl/__init__.py

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
	mkdir -p $(DESTDIR)/python/asl/gui ; \
	cp $(PYMODULES) $(DESTDIR)/python/asl/ ; \
	cp $(PYGUI) $(DESTDIR)/python/asl/gui ; \
        cp asl_gui_fsl $(DESTDIR)/bin/asl_gui ; \
	cd ..

clean:
	rm -f ${VERSIONED} asl_file *.o

FORCE:

