MODULE GZ_VARS_INPUT
  USE DMFT_PARSE_INPUT
  implicit none

  public :: read_input

  integer              :: Norb                ! number of orbitals
  real(8),dimension(3) :: Uloc   !local interactions
  real(8)              :: Ust
  real(8)              ::  Jh !hund's coupling
  real(8)              :: Jsf,Jph
  real(8)              :: xmu
  real(8),parameter    :: beta=500.d0
  
  !# Minimization flags  #!
  integer          :: lancelot_verbose
  logical          :: amoeba_verbose
  logical          :: GZmin_verbose
  character(len=4) :: min_method
  character(len=6) :: lgr_method
  integer          :: wf_symmetry
  !integer :: okkkkk
  real(8)   :: amoeba_min_tol
  real(8) :: err_self
  real(8) :: Rseed 
  real(8) :: Rmix
  integer :: Niter_self
  



contains

  subroutine read_input(INPUTunit)
    character(len=*) :: INPUTunit

    call parse_input_variable(Norb,"Norb",INPUTunit,default=1)
    call parse_input_variable(Uloc,"ULOC",INPUTunit,default=[2.d0,0.d0,0.d0])
    call parse_input_variable(ust,"UST",INPUTunit,default=0.d0)
    call parse_input_variable(Jh,"JH",INPUTunit,default=0.d0)
    call parse_input_variable(Jsf,"JSF",INPUTunit,default=0.d0)
    call parse_input_variable(Jph,"JPH",INPUTunit,default=0.d0)
    call parse_input_variable(xmu,"XMU",INPUTunit,default=0.d0)
    call parse_input_variable(wf_symmetry,"SYMMETRY",INPUTunit,default=0, &
         comment="0:O1cxSU2s;1:O1cxSU2sx(SU2c/O3c)")
    call parse_input_variable(lancelot_verbose,"LANCELOT_VERBOSE",INPUTunit,default=1)
    call parse_input_variable(amoeba_verbose,"AMOEBA_VERBOSE",INPUTunit,default=.false.)
    call parse_input_variable(GZmin_verbose,"GZMIN_VERBOSE",INPUTunit,default=.false.)
    call parse_input_variable(min_method,"MIN_METHOD",INPUTunit,default='nlep')
    call parse_input_variable(lgr_method,"LGR_METHOD",INPUTunit,default='amoeba')
    call parse_input_variable(Rseed,"RSEED",INPUTunit,default=1.d0)
    call parse_input_variable(Rmix,"RMIX",INPUTunit,default=1.d0)
    call parse_input_variable(Niter_self,"NITER_SELF",INPUTunit,default=100)
    call parse_input_variable(err_self,"ERR_SELF",INPUTunit,default=1.d-10)
    call parse_input_variable(amoeba_min_tol,"AMOEBA_MIN_TOL",INPUTunit,default=1.d-12)
    
    !call parse_input_variable(new_input,"NEW",INPUTunit,default=1.d-12)



  end subroutine read_input








END MODULE GZ_VARS_INPUT
