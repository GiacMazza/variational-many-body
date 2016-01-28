#COMPILER (PARALLEL)
FC=gfortran
#PRECOMPILATION FLAG (leave blank for serial code)
FPP=

#EXE=GZ_MB
#EXE=gz_2band_minN
#EXE=gz_2b_bethe
EXE=gz_test_symm


DIR=drivers
DIREXE=$(HOME)/.project_bin

.SUFFIXES: .f90

#REVISION SOFTWARE GIT:
BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

#OBJS=MATRIX_SPARSE.o ARPACK_LANCZOS.o AMOEBA.o GZ_VARS_INPUT.o GZ_VARS_GLOBAL.o  GZ_AUX_FUNX.o GZ_LOCAL_FOCK_SPACE.o GZ_VARIATIONAL_BASIS.o GZ_EFFECTIVE_HOPPINGS.o GZ_ENERGY_MINIMIZATION.o GZ_OPTIMIZED_ENERGY.o

OBJS=AMOEBA.o GZ_VARS_INPUT.o GZ_VARS_GLOBAL.o  GZ_AUX_FUNX.o GZ_LOCAL_FOCK_SPACE.o GZ_VARIATIONAL_BASIS.o GZ_EFFECTIVE_HOPPINGS.o GZ_ENERGY_MINIMIZATION.o GZ_OPTIMIZED_ENERGY.o


#GALLIBDIR  = /home/mazza/opt_local/galahad/objects/pc64.lnx.gfo/double
GALLIBDIR  = /opt/galahad/objects/pc64.lnx.gfo/double

GALLIBS1   = -lgalahad -lgalahad_hsl 
GALLIBS2   = -lgalahad_metis 

MKLARGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm

#FFLAG +=-fpp -D_$(FPP) ONLY WITH mpif90
#LIBDIR=/home/mazza/opt_local
LIBDIR=/opt/

INCARGS=-I$(LIBDIR)/scifor/gnu/include -L$(LIBDIR)/scifor/gnu/lib -I$(LIBDIR)/galahad/objects/pc64.lnx.gfo/double -L$(LIBDIR)/galahad/objects/pc64.lnx.gfo/double
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

