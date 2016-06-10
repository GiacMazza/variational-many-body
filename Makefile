#COMPILER (PARALLEL)
FC=gfortran
#PRECOMPILATION FLAG (leave blank for serial code)
FPP=


EXE=gz_generate_phi_traces

#EXE=gz_1b_nR

#EXE=gz_2b_bethe
#EXE=gz_2b_cubic_hyb
#EXE=gz_2b_cubic_hyb_nR
#EXE=gz_2b_bethe_pair_hopping

#EXE=gz_2b_bethe_nR
#EXE=gz_2b_bethe_pair_hopping_nR
#EXE=gz_bethe_janus   !+---> temporary out of order <---+!

#EXE=gz_1b_attractiveU_self

#EXE=gz_2b_bethe_sc
#EXE=gz_2b_bethe_sc_pair_hopping
#EXE=gz_2b_bethe_sc_nRQ
#EXE=gz_2b_bethe_sc_pair_hopping_nRQ

#EXE=gz_3b_bethe
#EXE=gz_3b_bethe_sc
#EXE=gz_3b_bethe_sc_sweep
#EXE=gz_3b_bethe_sweep
#EXE=gz_3b_bethe_sc_nRQ
#EXE=gz_3b_bethe_sc_su2


#EXE=gz_neq_sc
EXE=gz_neq_3b_bethe_sc



DIR=drivers
DIREXE=$(HOME)/.project_bin

.SUFFIXES: .f90

#REVISION SOFTWARE GIT:
BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

#OBJS=MATRIX_SPARSE.o ARPACK_LANCZOS.o AMOEBA.o GZ_VARS_INPUT.o GZ_VARS_GLOBAL.o  GZ_AUX_FUNX.o GZ_LOCAL_FOCK_SPACE.o GZ_VARIATIONAL_BASIS.o GZ_EFFECTIVE_HOPPINGS.o GZ_ENERGY_MINIMIZATION.o GZ_OPTIMIZED_ENERGY.o

OBJS=RK_VIDE.o MATRIX_SPARSE.o AMOEBA.o GZ_VARS_INPUT.o GZ_VARS_GLOBAL.o  GZ_AUX_FUNX.o GZ_neqAUX_FUNX.o GZ_LOCAL_FOCK_SPACE.o GZ_VARIATIONAL_BASIS.o GZ_LOCAL_HAMILTONIAN.o GZ_EFFECTIVE_HOPPINGS.o GZ_ENERGY.o GZ_OPTIMIZE.o GZ_DYNAMICS.o


GALLIBDIR  = $(HOME)/opt_local/galahad/objects/pc64.lnx.gfo/double
#GALLIBDIR  = /opt/galahad/objects/pc64.lnx.gfo/double

GALLIBS1   = -lgalahad -lgalahad_hsl 
GALLIBS2   = -lgalahad_metis 

MKLARGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm

#FFLAG +=-fpp -D_$(FPP) ONLY WITH mpif90
LIBDIR=$(HOME)/opt_local
#LIBDIR=/opt

INCARGS =-I$(LIBDIR)/SciFortran/gnu/include -L$(LIBDIR)/SciFortran/gnu/lib 
INCARGS+=-I$(LIBDIR)/DMFTtools/gnu/include -L$(LIBDIR)/DMFTtools/gnu/lib 
INCARGS+=-I$(LIBDIR)/galahad/objects/pc64.lnx.gfo/double -L$(LIBDIR)/galahad/objects/pc64.lnx.gfo/double
FFLAG += -ffree-line-length-none -cpp $(INCARGS)

ARGS= -L$(GALLIBDIR) $(GALLIBS1) $(GALLIBS2) -I$(LIBDIR)/galahad/modules/pc64.lnx.gfo/double -ldmftt -lscifor  -lfftpack -lminpack  -llapack -lblas -larpack -lparpack    

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

