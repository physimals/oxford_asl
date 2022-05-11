include ${FSLCONFDIR}/default.mk

PROJNAME = oxford_asl

# The FSL build system changed
# substantially in FSL 6.0.6
# FSL >= 6.0.6
ifeq (${FSL_GE_606}, true)
  LIBS = -lfsl-newimage -lfsl-miscmaths -lfsl-cprob -lfsl-utils \
         -lfsl-NewNifti -lfsl-znz
# FSL <= 6.0.5
else
  ifeq ($(shell uname -s), Linux)
	MATLIB := -lopenblas
  endif
  USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB} \
                -I${FSLDIR}/extras/include/armawrap
  USRLDFLAGS  = -L${LIB_NEWMAT} -L${LIB_ZLIB} -L${LIB_CPROB} \
                -lutils -lnewimage -lmiscmaths -lcprob       \
                ${MATLIB} -lNewNifti -lznz -lz
endif

XFILES    = asl_file
SCRIPTS   = oxford_asl asl_calib asl_reg quasil toast oxford_asl_roi_stats.py oxford_asl_hadamard_decode.py
PYMODULES = python/oxford_asl/*.py
PYGUI     = python/oxford_asl/gui/*.py python/oxford_asl/gui/*.png
VERSIONED = oxford_asl asl_calib quasil asl_reg toast

OBJS = readoptions.o asl_functions.o

# Pass Git revision details
GIT_SHA1 := $(shell git describe --dirty)
GIT_DATE := $(shell git log -1 --format=%ad --date=local)
CXXFLAGS += -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_DATE="\"${GIT_DATE}\""

# Always rebuild scripts
.PHONY: FORCE

all: ${XFILES} ${VERSIONED}

asl_file: ${OBJS} asl_file.o
	${CXX} ${CXXFLAGS} -o $@ $^ ${LDFLAGS}

$(VERSIONED): %: %.in FORCE
	sed -e "s/@GIT_SHA1@/${GIT_SHA1}/" -e "s/@GIT_DATE@/${GIT_DATE}/" $< >$@
	chmod a+x $@

# Only called in FSL <= 605
postinstallscript: $(PYMODULES) $(PYGUI)
	mkdir -p $(DESTDIR)/python/oxford_asl/gui
	cp $(PYMODULES) $(DESTDIR)/python/oxford_asl/
	cp $(PYGUI) $(DESTDIR)/python/oxford_asl/gui/
	cp wrappers/* $(DESTDIR)/bin/

# Only called in FSL >= 606
# setup.py -V is called to force creation
# of _version.py before installation
pyinstall:
	fslpython python/setup.py -V
	fslpython -m pip install --no-deps -vv ./python/

clean:
	rm -f ${VERSIONED} asl_file *.o
