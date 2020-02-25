MODULE GZ_VARS_INPUT
  USE SF_PARSE_INPUT
  implicit none

  public :: read_input

  integer              :: Norb                ! number of orbitals  
  real(8),dimension(3) :: Uloc   !local interactions
  real(8)              :: Ust
  real(8)              ::  Jh !hund's coupling
  real(8)              :: Jsf,Jph
  real(8)              :: xmu
  real(8)              :: beta
  integer,parameter    :: lw=512
  real(8)              :: wini,wfin
  logical              :: slater_store
  real(8)              :: k_dens_diss
  real(8)              :: k_qp_diss,beta_diss
  real(8)              :: a_nhh
  integer              :: Nsite
  
  !# Minimization flags  #!
  integer              :: lancelot_verbose
  logical              :: amoeba_verbose
  logical              :: GZmin_verbose
  character(len=4)     :: min_method
  character(len=6)     :: lgr_method
  integer              :: wf_symmetry
  logical              :: gz_superc
  real(8)              :: amoeba_min_tol
  real(8)              :: err_self
  real(8)              :: Rseed 
  real(8)              :: Rmix
  integer              :: Niter_self

  !# DYNAMICS #!
  integer              :: Nt
  real(8)              :: tstart
  real(8)              :: tstep
  logical              :: GZneq_verbose
  !# Electric field #!
  character(len=16)    :: field_type    !choose the profile of the electric field
  real(8)              :: Dpulse
  real(8)              :: Efield        !Electric field strength
  real(8),dimension(3) :: Evect         !Electric field vectors as input
  real(8)              :: Ton,Toff      !turn on/off time, t0 also center of the pulse
  integer              :: Ncycles       !Number of cycles in pulsed light packet
  real(8)              :: omega0        !parameter for the Oscilatting field and Pulsed light
  real(8)              :: E1            !Electric field strenght for the AC+DC case (tune to resonate)
  

contains

  subroutine read_input(INPUTunit)
    character(len=*) :: INPUTunit

    call parse_input_variable(Norb,"Norb",INPUTunit,default=1)
    call parse_input_variable(Uloc,"ULOC",INPUTunit,default=[0.d0,0.d0,0.d0])
    call parse_input_variable(ust,"UST",INPUTunit,default=0.d0)
    call parse_input_variable(Jh,"JH",INPUTunit,default=0.d0)
    call parse_input_variable(Jsf,"JSF",INPUTunit,default=0.d0)
    call parse_input_variable(Jph,"JPH",INPUTunit,default=0.d0)
    call parse_input_variable(xmu,"XMU",INPUTunit,default=0.d0)
    call parse_input_variable(wf_symmetry,"SYMMETRY",INPUTunit,default=0, &
         comment="0:O1cxSU2s;1:O1cxSU2sx(SU2c/O3c)")
    call parse_input_variable(gz_superc,"GZ_SUPERC",INPUTunit,default=.false.)
    call parse_input_variable(lancelot_verbose,"LANCELOT_VERBOSE",INPUTunit,default=1)
    call parse_input_variable(amoeba_verbose,"AMOEBA_VERBOSE",INPUTunit,default=.false.)
    call parse_input_variable(GZmin_verbose,"GZMIN_VERBOSE",INPUTunit,default=.false.)
    call parse_input_variable(min_method,"MIN_METHOD",INPUTunit,default='nlep')
    call parse_input_variable(lgr_method,"LGR_METHOD",INPUTunit,default='f_zero')
    call parse_input_variable(Rseed,"RSEED",INPUTunit,default=1.d0)
    call parse_input_variable(Rmix,"RMIX",INPUTunit,default=1.d0)
    call parse_input_variable(Niter_self,"NITER_SELF",INPUTunit,default=100)
    call parse_input_variable(err_self,"ERR_SELF",INPUTunit,default=1.d-10)
    call parse_input_variable(amoeba_min_tol,"AMOEBA_MIN_TOL",INPUTunit,default=1.d-12)
    call parse_input_variable(wini,"WINI",INPUTunit,default=-10.d0)
    call parse_input_variable(wfin,"WFIN",INPUTunit,default=10.d0)
    call parse_input_variable(slater_store,"SLATER_STORE",INPUTunit,default=.false.)
    call parse_input_variable(beta,"BETA",INPUTunit,default=1000.d0)           !+- fictitious temperature -+!
    call parse_input_variable(k_qp_diss,"K_DISS",INPUTunit,default=0.d0)       !+- friction dissipation -+!
    call parse_input_variable(a_nhh,"A_NHH",INPUTunit,default=1.d0)       !+- lindblat[a_nhh=1.0] -> non-hermitean [a_nhh=0.0]  -+!
    call parse_input_variable(beta_diss,"BETA_DISS",INPUTunit,default=100.d0)  !+- friction "temperature"; ie fermi-function distribution of the bath -+!

    call parse_input_variable(Nsite,"Nsite",INPUTunit,default=2)

    !
    call parse_input_variable(Nt,"NT",INPUTunit,default=100)
    call parse_input_variable(tstart,"TSTART",INPUTunit,default=0.d0)
    call parse_input_variable(tstep,"TSTEP",INPUTunit,default=1.d-2)
    call parse_input_variable(GZneq_verbose,"GZNEQ_VERBOSE",INPUTunit,default=.false.)
    !
    !ELECTRIC FIELD VARIABLES
    call parse_input_variable(field_type,"FIELD_TYPE",INPUTunit,default ='pulse',comment="profile type of the electric field ")
    call parse_input_variable(Efield,"EFIELD",INPUTunit,default=0d0,comment="electric field strength")
    call parse_input_variable(Evect,"EVECT",INPUTunit,default=[1d0,0d0,0d0],comment="electric field direction (normalized)")
    call parse_input_variable(ton,"TON",INPUTunit,default=0d0,comment="turn on time or center of the pulse")
    call parse_input_variable(toff,"TOFF",INPUTunit,default=10000d0,comment="turn off time")
    call parse_input_variable(Dpulse,"DPULSE",INPUTunit,default=10.d0,comment="time width of the field pulse")
    call parse_input_variable(omega0,"OMEGA0",INPUTunit,default=acos(-1d0) , comment="parameter for the Oscilatting field and Pulsed light")
    !
  end subroutine read_input








END MODULE GZ_VARS_INPUT
