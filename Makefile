###########################
# Makefile for GLISSANDO 2
###########################

# preprocessor options 

# description of the preprocessor parameters:
# _nnwp_=0     - use the hard-sphere NN wounding profile
# _nnwp_=1     - use Gaussian NN wounding profile
# _nnwp_=2     - use NN wounding profile given by gamma approximation 
# _files_=1    - read the nuclear distributions from external files, =0 - generate randomly
# _profile_=1  - generate the nucleon profile and NN correlation data, =0 - do not
# _weight_=1   - generate data for the NN collision profiles and for the weight (RDS) ditributions, 0 - do not 
# _rapidity_=1 - generate the data for the rapidity distributions, =0 - do not
# _evout_=1    - generate the full event data, =0 - do not


# default:
# PREPROCESS = -D_nnwp_=1 -D_files_=0 -D_profile_=0 -D_weight_=0 -D_rapidity_=0 -D_evout_=0

PREPROCESS   := -D_nnwp_=1 -D_files_=0 -D_profile_=1 -D_weight_=1 -D_rapidity_=0 -D_evout_=0

# version information from file glissando2.cxx
VERSION = $(shell grep "\#define _VER_"  $(DIR_CXX)glissando2.cxx | sed 's|[\#,_,2_,a-z,A-Z, ]*||')

## directory structure
DIR_MAIN    = ./
DIR_ADDONS  = $(DIR_MAIN)addons/
DIR_BUILD   = $(DIR_MAIN)build/
DIR_DOXYGEN = $(DIR_BUILD)doxygen/
DIR_H       = $(DIR_BUILD)include/
DIR_OBJ     = $(DIR_BUILD)obj/
DIR_CXX     = $(DIR_BUILD)src/
DIR_DOC     = $(DIR_MAIN)$(shell grep "OUTPUT_DIRECTORY *=" Doxyfile | sed 's|[A-Z,_, \t,=,./]*||')
DIR_DOCHTML = $(DIR_DOC)html/
DIR_DOCTEX  = $(DIR_DOC)latex/
DIR_MACRO   = $(DIR_MAIN)macro/
DIR_INPUT   = $(DIR_MAIN)input/
DIR_OUTPUT  = $(DIR_MAIN)output/

# search paths
vpath %.h   $(DIR_H)
vpath %.cxx $(DIR_CXX)
vpath %.o   $(DIR_OBJ)

# distribution file lists
F_INCLUDE   = $(DIR_H)*.h
F_SOURCE    = $(DIR_CXX)*.cxx
F_DOXYGEN   = $(DIR_MAIN)Doxyfile
F_MACRO     = $(DIR_MACRO)*.C
F_INPUT     = $(DIR_INPUT)*.dat
F_BASH      = $(DIR_MAIN)*.sh
F_ADDONS    = $(DIR_ADDONS)*.cxx $(DIR_ADDONS)*.mk
F_PACK       = glissando2$(VERSION).tar.gz 
#F_PACK      = glissando2.101.tar.gz 
F_README    = $(DIR_MAIN)README

# compilation file lists
BIN_GLISS  = glissando2
HSRC_GLISS = glissando2.cxx
SRC_GLISS  = $(HSRC_GLISS:%=$(DIR_CXX)%) $(BIN_GLISS:%=$(DIR_CXX)%.cxx)
OBJ_GLISS  = $(SRC_GLISS:$(DIR_CXX)%.cxx=$(DIR_OBJ)%.o)

# compilation
CXX         = g++
LD          = g++
CXXFLAGS    = -O3 -g -w -Wno-deprecated -I $(DIR_H) $(PREPROCESS) `root-config --cflags`
LFLAGS      = `root-config --libs` -lm -lgcc

#################################################################################
# RULES                                                                         #












#################################################################################

all: $(BIN_GLISS:%=$(DIR_OBJ)%) 
	cp $^ $(DIR_MAIN)
	mkdir -p $(DIR_OUTPUT)
	echo "Type \"./glissando2\" to run the code"

$(DIR_OBJ)glissando2: $(OBJ_GLISS)
	echo "Linking:   $@ ($(LD))"
	$(LD)  $^ -o $@ $(LFLAGS)

$(DIR_OBJ)%.o: %.cxx $(F_INCLUDE)
	@[ -d $(DIR_OBJ) ] || mkdir -p $(DIR_OBJ)
	echo "Compiling: $< ($(CXX))"
	$(CXX) $(CXXFLAGS) -c $< -o $@

doc:
	doxygen
	rm -f $(DIR_DOXYGEN)file_*.dox
	echo "<!DOCTYPE HTML PUBLIC '-//W3C//DTD HTML 4.0 Transitional//EN'><HTML><HEAD><META http-equiv='REFRESH' content='0;url=$(DIR_DOCHTML)index.html'></HEAD><BODY></BODY></HTML>" > manual.html
	$(MAKE) -C $(DIR_DOCTEX)
	ln -s $(DIR_DOCTEX)refman.pdf manual.pdf
	echo
	echo "GLISSANDO 2 reference manual generated with Doxygen"
	echo "- manual.html"
	echo "- manual.pdf"

package: $(F_INCLUDE) $(F_SOURCE) $(F_MACRO) $(F_DOXYGEN) $(F_INPUT) $(F_BASH) $(F_ADDONS) $(F_README) doc/latex/refman.pdf one.dat Makefile
	echo "$(VERSION)" > version
	tar zcvf $(F_PACK) $^ version
	echo "Package '$(F_PACK)' created"

clean:
	rm -rf $(DIR_OBJ)
	rm -f $(DIR_OBJ)$(BIN_GLISS) $(DIR_MAIN)$(BIN_GLISS)
	rm -f err 
	rm -f version
	echo "*.o and binary files removed"

cleandoc:
	rm -rf $(DIR_DOC)
	rm -f $(DIR_MAIN)manual.html $(DIR_MAIN)manual.pdf
	echo "Doxygen documentation removed"

cleanoutput:
	rm -f $(DIR_OUTPUT)*
	echo "Output files removed"

.SILENT :

