include ${FSLCONFDIR}/default.mk

PROJNAME = oxford_asl

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_ZLIB}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_ZLIB}

LIBS = -lutils -lnewimage -lmiscmaths -lm -lnewmat -lfslio -lniftiio -lznz -lz

XFILES = asl_file
SCRIPTS = oxford_asl asl_calib asl_reg quasil asl_gui
PYMODULES = python/asl/__init__.py python/asl/fslhelpers.py python/asl/reg.py 
PYGUI = python/asl/gui/*.py
RUNTCLS = Asl
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
	mkdir -p $(FSLDEVDIR)/python/asl/gui ; \
	cp $(PYMODULES) $(FSLDEVDIR)/python/asl/ ; \
	cp $(PYGUI) $(FSLDEVDIR)/python/asl/gui ; \
	cd ..

clean:
	rm -f ${VERSIONED} asl_file *.o

FORCE:

