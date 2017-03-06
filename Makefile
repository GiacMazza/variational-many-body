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
#EXE=gz_1b_cb_sc_nRQ
#EXE=gz_1b_bethe_sc
#EXE=gz_1b_eom
#EXE=gz_imt_1b_bethe
EXE=gz_imt_Hqp

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
#EXE=gz_neq_3b_bethe_sc_su2_read
#EXE=gz_neq_3b_bethe_u1su2
#EXE=gz_neq_3b_bethe_sc_su2O1

#EXE=gz_GF_pp_sc


DIR=drivers
DIREXE=$(HOME)/.project_bin

.SUFFIXES: .f90

#REVISION SOFTWARE GIT:
BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc
#
OBJS=RK_VIDE.o MATRIX_SPARSE.o AMOEBA.o GZ_VARS_INPUT.o GZ_VARS_GLOBAL.o ELECTRIC_FIELD.o  GZ_AUX_FUNX.o GZ_neqAUX_FUNX.o GZ_LOCAL_FOCK_SPACE.o GZ_VARIATIONAL_BASIS.o GZ_LOCAL_HAMILTONIAN.o GZ_EFFECTIVE_HOPPINGS.o GZ_ENERGY.o GZ_OPTIMIZE.o GZ_DYNAMICS.o GZ_imtDYNAMICS.o #GZ_GREENS_FUNCTIONS.o
#



#FFLAG +=-fpp -D_$(FPP) ONLY WITH mpif90
LIBDIR=$(HOME)/opt_local
#LIBDIR=/opt/

#GALLIBDIR  = $(LIBDIR)/galahad/objects/mac64.osx.gfo/double
#GALLIBMOD  = $(LIBDIR)/galahad/modules/mac64.osx.gfo/double
GALLIBDIR  = $(LIBDIR)/galahad/objects/pc64.lnx.gfo/double
GALLIBMOD  = $(LIBDIR)/galahad/modules/pc64.lnx.gfo/double


GALLIBS1   = -lgalahad -lgalahad_hsl 
GALLIBS2   = -lgalahad_metis 

MKLARGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm


INCARGS =-I$(LIBDIR)/SciFortran/gnu/include -L$(LIBDIR)/SciFortran/gnu/lib 
INCARGS+=-I$(LIBDIR)/DMFTtools/gnu/include -L$(LIBDIR)/DMFTtools/gnu/lib 




INCARGS+=-I$(GALLIBDIR) -L$(GALLIBDIR)
FFLAG += -ffree-line-length-none -cpp $(INCARGS)

FFLAG+=-O0 -p -g -Wall -fbounds-check -fbacktrace -Wuninitialized

ARGS= -L$(GALLIBDIR) $(GALLIBS1) $(GALLIBS2) -I$(GALLIBMOD) -ldmftt -lscifor  -lfftpack -lminpack  -llapack -lblas -larpack #-lparpack    

all:compile


lib: ed_solver

compile: version $(OBJS)
	@echo " !+------------------------------------------------- "
	@echo " ..................... compile ..................... "
	@echo " !+------------------------------------------------- "
	$(FC) $(FFLAG) $(OBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/$(EXE)_$(BRANCH) $(ARGS)
	@echo " !+--------------------------------------+! "
	@echo " .................. done .................. "
	@echo " !+--------------------------------------+! "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/$(EXE)_$(BRANCH)

ed_solver:
	@make -C ED_SOLVER/

.f90.o:	
	$(FC) $(FFLAG)  -c $< 

completion:
	sf_lib_completion.sh $(DIR)/$(EXE).f90
	@echo "run: . .bash_completion.d/$(EXE) to add completion for $(EXE) in this shell"

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc

all_clean: clean
	@make -C ED_SOLVER/ clean

version:
	@echo $(VER)

