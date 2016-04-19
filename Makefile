include make.inc

OBJS=MATRIX_SPARSE.o AMOEBA.o LANCELOT_WRAP.o GZ_VARS_INPUT.o GZ_VARS_GLOBAL.o  GZ_AUX_FUNX.o GZ_LOCAL_FOCK_SPACE.o GZ_VARIATIONAL_BASIS.o GZ_EFFECTIVE_HOPPINGS.o GZ_ENERGY.o GZ_OPTIMIZE.o SOLVE_GZ.o


all:compile

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

.f90.o:	
	$(FC) $(FFLAG)  -c $< 

completion:
	sf_lib_completion.sh $(DIR)/$(EXE).f90
	@echo "run: . .bash_completion.d/$(EXE) to add completion for $(EXE) in this shell"

clean: 
	@echo "Cleaning:"
	@rm -f *.mod *.o *~ revision.inc

version:
	@echo $(VER)

