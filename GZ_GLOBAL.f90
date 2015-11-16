MODULE GZ_VARS_GLOBAL
  implicit none

  !+--------------------+!
  !+- LOCAL FOCK SPACE -+!
  !+--------------------+!
  integer,public                            :: Norb                ! number of orbitals
  integer,public                            :: state_dim           ! dimension of a single Fock state  |(\up,\dn)_1,...,(\up,\dn)_Norb> ===> state_dim = 2*Norb  
  integer,public                            :: nFock               ! dimension of the local Fock space = 2^{2*Norb} 
  integer,public                            :: nFock_indep         ! dimension of the independent Fock space ----> this is related to the constrained minimization w/ Galahad 
  
  integer,dimension(:,:),allocatable,public :: index               ! ispin,iorb  to  istate=1,state_dim
  integer,dimension(:),allocatable,public   :: fock_indep          ! store the fock independent states
  integer,dimension(:),allocatable,public   :: full2indep_fock     ! map fock state to the corresponding independent one
  integer,dimension(:,:),allocatable,public :: indep2full_fock     ! map independent fock states to the full fock space

  !+---------------------------+!
  !+- GUTZWILLER WAVEFUNCTION -+!
  !+---------------------------+!

  !# Operators in the space (nFock X nFock) #!
  real(8),public,allocatable                :: CC(:,:,:),CA(:,:,:)          ! creation annhilation 
  real(8),public,allocatable                :: local_hamiltonian(:,:)       ! Hamiltonian of the atomic problem
  real(8),public,allocatable                :: dens_dens_interaction(:,:)   ! Density-density interaction for the atomic problem

  real(8),public,allocatable                :: UHubbard(:,:)       ! Hubbard interaction
  real(8),public,allocatable                :: docc(:,:,:)         ! Double occupancy (Norb,:,:)
  real(8),public,allocatable                :: dens(:,:,:)         ! local density
  
  !# Gutzwiller Matrices NOTE: for the moment only the DIAGONAL CASE is considered #!
  real(8),public,allocatable                :: phi_traces_basis_Rhop(:,:,:,:) ! state_dim X state_dim matrices of dimension (nFock X nFock) 
  real(8),public,allocatable                :: phi_traces_basis_dens(:,:,:,:) ! state_dim X state_dim matrices of dimension (nFock X nFock) 
  real(8),public,allocatable                :: phi_traces_basis_Hloc(:,:)   ! Single matrix of dimension (nFock X nFock) 
  real(8),public,allocatable                :: GZproj_vect(:)               ! dimension nFock
  
  !# Gutzwiller renormalization factors #!
  complex(8),public,allocatable             :: Rhop_c(:,:)         ! renormalization factors real
  real(8),public,allocatable                :: Rhop_r(:,:)         ! reormalization factors imag  
  real(8),public,allocatable                :: Rhop_diag(:)        ! reormalization factors diagonal imag    
  complex(8),public,allocatable             :: Phi_c(:,:)          ! gz_projectors -full- real
  real(8),public,allocatable                :: Phi_r(:,:)          ! gz_projectors -full- imag
  real(8),public,allocatable                :: phi_gz(:)           ! gz_projectors diagonal

  !# Slater Determinant #!
  real(8),public,dimension(:,:,:),allocatable :: slater_matrix_elements   ! (state_dim,state_dim,Lk) aka <c^+_{k\alpha} c_{k\beta}> \alpha,\beta=1,state_dim
  real(8),public,dimension(:,:),allocatable   :: slater_ground_state_deriv   ! (state_dim,state_dim) aka d<H*>/R_{\alpha,\beta} \alpha,\beta=1,state_dim
  
  !# Observables #!
  real(8),public,allocatable                :: gz_docc(:)         ! Double occupancy (1:Norb)
  real(8),public,allocatable                :: gz_dens(:)         ! local density    (1:Norb*Nspin)

  !# Minimization flags  #!
  logical :: fix_density_minimization
  integer :: lancelot_verbose
  logical :: amoeba_verbose
  logical :: GZmin_verbose

  real(8) :: Rseed,err_self
  integer :: Niter_self
  integer :: opt_energy_unit,GZmin_unit

  !+-------------------------+!
  !+- MODEL DETAILS DETAILS -+!
  !+-------------------------+!
  integer                          :: Nx,Lk
  real(8),dimension(:),allocatable :: epsik,hybik
  real(8),dimension(:),allocatable :: wtk
  real(8)                          :: Wband,Vhyb
  real(8)                          :: xmu
  real(8),dimension(:),allocatable :: Eloc
  real(8),dimension(:),allocatable :: atomic_energy_levels
  real(8),parameter                :: beta=500.d0
  real(8)                          :: U,Cfield,Nread
  

END MODULE GZ_VARS_GLOBAL
