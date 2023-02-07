#COMPILER (PARALLEL)
FC=mpif90
#gfortran
#PRECOMPILATION FLAG (leave blank for serial code)
FPP=


#EXE=gz_generate_phi_traces

#EXE=gz_1b_nR
#EXE=gz_1b_bethe_nR

#EXE=gz_2b_bethe
#EXE=gz_2b_cubic_hyb
#EXE=gz_2b_cubic_hyb_nR
#EXE=gz_2b_bethe_pair_hopping

#EXE=gz_2b_bethe_nR
#EXE=gz_2b_bethe_pair_hopping_nR
#EXE=gz_bethe_janus   !+---> temporary out of order <---+!

#EXE=gz_1b_attractiveU


#EXE=gz_1b_bethe_sc_nRQ
#EXE=gz_1b_cubic_sc_nRQ
##EXE=gz_1b_cb_sc_nRQ
#EXE=gz_1b_bethe_sc_temp
#EXE=gz_1b_eom


#EXE=gz_neq_1b_bethe
#EXE=gz_neq_1b_bethe_sc
#EXE=gz_neq_1b_cubic_sc
#EXE=gz_neq_1b_cb_sc
#EXE=gz_neq_1b_bethe_sc_tdlgr


#EXE=gz_2b_bethe_sc
#EXE=gz_2b_bethe_sc_pair_hopping
#EXE=gz_2b_bethe_sc_nRQ
#EXE=gz_2b_bethe_sc_pair_hopping_nRQ

#EXE=gz_3b_bethe
#EXE=gz_3b_bethe_nR
#EXE=gz_neq_3b_bethe

#EXE=gz_3b_bethe_sc
#EXE=gz_3b_bethe_sc_su2nlep
#EXE=gz_3b_bethe_sc_sweep
#EXE=gz_3b_bethe_sweep
#EXE=gz_3b_bethe_sc_nRQ
#EXE=gz_3b_bethe_sc_su2
#EXE=gz_3b_bethe_sc_su2isoTz


#EXE=gz_neq_sc
#EXE=gz_neq_3b_bethe_sc
#EXE=gz_neq_3b_bethe_sc_su2
#EXE=gz_neq_3b_bethe_sc_linearXY
#EXE=gz_neq_3b_bethe_sc_su2_read
#EXE=gz_neq_3b_bethe_u1su2
#EXE=gz_neq_3b_bethe_sc_su2O1

#EXE=gz_GF_pp_sc

#EXE=gz_diss_sc_bethe
EXE=bcs_diss_solitons


DIR=drivers
DIREXE=$(HOME)/.bin



#REVISION SOFTWARE GIT:
BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
ifeq ($(BRANCH),master)
BRANCH=
endif
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc
#
OBJS=RK_VIDE.o MATRIX_SPARSE.o AMOEBA.o GZ_VARS_INPUT.o GZ_VARS_GLOBAL.o ELECTRIC_FIELD.o  GZ_AUX_FUNX.o GZ_neqAUX_FUNX.o GZ_LOCAL_FOCK_SPACE.o GZ_VARIATIONAL_BASIS.o GZ_LOCAL_HAMILTONIAN.o GZ_EFFECTIVE_HOPPINGS.o GZ_ENERGY.o GZ_OPTIMIZE.o GZ_DYNAMICS.o GZ_GREENS_FUNCTIONS.o
#

GLOB_INC:=$(shell pkg-config --cflags dmft_tools scifor)
GLOB_LIB:=$(shell pkg-config --libs dmft_tools scifor)  
#
#
LIBDIR=$(HOME)/opt
GALLIBDIR  = $(LIBDIR)/galahad/objects/pc64.lnx.gfo/double
GALLIBMOD  = $(LIBDIR)/galahad/modules/pc64.lnx.gfo/double
GALLIBS = -lgalahad -lgalahad_hsl -lgalahad_metis4

GLOB_INC+=-I$(GALLIBMOD) -L$(GALLIBDIR)
GLOB_LIB+=$(GALLIBS)

#

#INCARGS =-I$(LIBDIR)/old_libs/SciFortran/gnu/include -L$(LIBDIR)/SciFortran/gnu/lib 
#INCARGS+=-I$(LIBDIR)/old_libs/DMFTtools/gnu/include -L$(LIBDIR)/DMFTtools/gnu/lib 
#INCARGS+=-I$(GALLIBMOD) -L$(GALLIBDIR)
#FFLAG += -ffree-line-length-none -cpp $(INCARGS)

#FFLAG+=-O0 -p -g -Wall -fbounds-check -fbacktrace -Wuninitialized
#ARGS=-I$(LIBDIR)/old_libs/SciFortran/gnu/include  -L$(LIBDIR)/old_libs/SciFortran/gnu/lib  -lscifor -lfftpack -lminpack  -llapack -lblas -larpack

#ARGS=-L$(GALLIBDIR) $(GALLIBS1) $(GALLIBS2) 
#ARGS+=-L$(LIBDIR)/old_libs/DMFTtools/gnu/lib  -ldmftt
#ARGS+=-L$(LIBDIR)/old_libs/SciFortran/gnu/lib  -lscifor -lfftpack -lminpack  -llapack -lblas -larpack


FFLAG = -O2 -ffree-line-length-none
DFLAG = -O0 -p -g -fimplicit-none -Wsurprising -Wuninitialized -fbounds-check  -Waliasing -Wall -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none
OFLAG = -O3 -ffast-math -march=native -funroll-loops -ffree-line-length-none
FPPSERIAL =-cpp -D_
FPPMPI =-cpp -D_MPI	




define colorecho	
	@tput setaf $2
	@tput bold
	@echo $1
	@tput sgr0
endef

.SUFFIXES: .f90

.f90.o:
	$(FC) $(FLAG) $(GLOB_INC) -c $<


all: FLAG:=${FFLAG} ${FPPMPI}
all:	$(OBJS)
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ", 6)
	@echo ""
	$(FC) $(FLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE) ${GLOB_INC} ${GLOB_LIB}
	@echo "Done"
	$(call colorecho,"created $(EXE) in  $(DIREXE)", 1)


debug: FLAG:=${DFLAG} ${FPPMPI}
debug: $(OBJS)
	@echo ""
	$(call colorecho,"compiling $(EXE).f90 ", 6)
	@echo ""
	$(FC) $(FLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE) ${GLOB_INC} ${GLOB_LIB}
	@echo "Done"
	$(call colorecho,"created $(EXE) in  $(DIREXE)", 1)


clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc

