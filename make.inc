##$ COMPILER: suppprted compilers are ifort, gnu >v4.7
##$ use mpif90 for parallel compiler
FC=gfortran

##$ PLATFORM: supported platform are intel, gnu
##$ if using a parallel compiler (mpif90) check
##$ platform first using the cmd: $mpif90 -show
PLAT=gnu


##$ PRECOMPILATION FLAGS
##$ Leave blank for SERIAL CODE
##$
FPP=

##$ LOCATION OF THE scifor and dmft_tools DIRECTORIES
##$ is placed in different locations edit INCARGS here below
LIBDIR=$(HOME)/opt


##$ CHOOSE THE MODEL BY SELECTING THE PROGRAM DRIVER
##$ if needed add yours to the list:

#EXE=GZ_MB
#EXE=gz_2band_minN
EXE=gz_2b_bethe
#EXE=gz_test_symm
#EXE=gz_optimize_VS_nR
#EXE=gz_optimize_janus

##$ SET THE LOCATION OF YOU PROGRAM DRIVER (default is ./drivers)
DIR=drivers

##$ SET THE LOCATION WHERE TO PLACE THE EXECUTABLE (default is $HOME/.bin)
DIREXE=$(HOME)/.bin


ifeq ($(PLAT),intel)
FFLAG = -O2 -ftz -fpp -D_$(FPP)
DFLAG = -p -O0 -g -fpe0 -warn -warn errors  -debug extended -traceback -check all,noarg_temp_created -fpp -D_$(FPP)
OFLAG = -O3 -ftz -prof-use -prof-gen -ip -ipo -fpp -D_$(FPP)
endif

ifeq ($(PLAT),gnu)
INCARGS=-I$(LIBDIR)/scifor/gnu/include -L$(LIBDIR)/scifor/gnu/lib
INCARGS+=-I$(LIBDIR)/dmft_tools/gnu/include -L$(LIBDIR)/dmft_tools/gnu/lib
INCARGS+=-I$(GALAHAD)/modules/pc64.lnx.gfo/double -L$(GALAHAD)objects/pc64.lnx.gfo/double
FFLAG = -O2 -funroll-all-loops -ffree-line-length-none -cpp -D_$(FPP) $(INCARGS) 
DFLAG = -O0 -p -g -Wall -fimplicit-none -ffpe-trap=zero,overflow,underflow  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -pedantic  -fbacktrace -fcheck=bounds -ffree-line-length-none -cpp -D_$(FPP) $(INCARGS)
OFLAG = -O3 -ffast-math -march=native -funroll-all-loops -fno-protect-parens -flto -ffree-line-length-none -cpp -D_$(FPP) $(INCARGS 
endif




##$ CHOOSE LINKING OPTIONS:
##$ 
##$ If you intend to use mkl:
##$ 
#MKLARGS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
#ARGS= -lgalahad -lgalahad_hsl -lgalahad_metis -ldmftt -lscifor $(MKLARGS) -larpack -lparpack 
##$ 
##$ ELSE:
##$ 
ARGS=  -lgalahad -lgalahad_hsl -lgalahad_metis -ldmftt -lscifor -lfftpack -lminpack -llapack -lblas -larpack -lparpack 




##$ REVISION SOFTWARE VARIABLES
##$ 
REV=$(shell git rev-parse HEAD)
BRANCH=_$(shell git rev-parse --abbrev-ref HEAD)
VER = 'character(len=41),parameter :: revision = "$(REV)"' > revision.inc

ifeq ($(BRANCH),_master)
BRANCH=
endif


##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90
