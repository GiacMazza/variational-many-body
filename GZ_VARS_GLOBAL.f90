MODULE GZ_VARS_GLOBAL
  USE GZ_VARS_INPUT
  implicit none
  
  !+--------------------+!
  !+- LOCAL FOCK SPACE -+!
  !+--------------------+!
  integer                            :: state_dim           ! dimension of a single Fock state  |(\up,\dn)_1,...,(\up,\dn)_Norb> ===> state_dim = 2*Norb  
  integer                            :: nFock               ! dimension of the local Fock space = 2^{2*Norb} 
  integer                            :: nFock_indep         ! dimension of the independent Fock space ----> this is related to the constrained minimization w/ Galahad   
  !
  integer,dimension(:,:),allocatable :: index               ! ispin,iorb  to  istate=1,state_dim
  integer,dimension(:),allocatable   :: fock_indep          ! store the fock independent states
  integer,dimension(:),allocatable   :: full2indep_fock     ! map fock state to the corresponding independent one
  integer,dimension(:,:),allocatable :: indep2full_fock     ! map independent fock states to the full fock space

  !+---------------------------+!
  !+- GUTZWILLER WAVEFUNCTION -+!
  !+---------------------------+!
  !# Operators in the space (nFock X nFock) #!
  real(8),allocatable                :: CC(:,:,:),CA(:,:,:)          ! creation annhilation 
  real(8),allocatable                :: local_hamiltonian(:,:)       ! Hamiltonian of the atomic problem
  real(8),allocatable                :: local_hamiltonian_free(:,:)       ! free Hamiltonian of the atomic problem
  real(8),allocatable                :: dens_dens_interaction(:,:)   ! Density-density interaction for the atomic problem

  real(8),allocatable                :: UHubbard(:,:)       ! Hubbard interaction
  
  !# Gutzwiller Matrices NOTE: for the moment only the DIAGONAL CASE is considered #!
  real(8),allocatable                :: phi_traces_basis_Rhop(:,:,:,:) ! state_dim X state_dim matrices of dimension (nFock X nFock) 
  real(8),allocatable                :: phi_traces_basis_dens(:,:,:,:) ! state_dim X state_dim matrices of dimension (nFock X nFock) 
  real(8),allocatable                :: phi_traces_basis_Hloc(:,:)     ! Single matrix of dimension (nFock X nFock) 
  real(8),allocatable                :: phi_traces_basis_free_Hloc(:,:)     ! Single matrix of dimension (nFock X nFock) 
  !  real(8),allocatable                :: GZproj_vect(:)                 ! dimension nFock
  
  !# Gutzwiller renormalization factors #!
  complex(8),allocatable             :: Rhop_c(:,:)         ! renormalization factors real

  real(8),allocatable                :: Rhop_diag(:)        ! reormalization factors diagonal imag    
  complex(8),allocatable             :: Phi_c(:,:)          ! gz_projectors -full- real
  real(8),allocatable                :: Phi_r(:,:)          ! gz_projectors -full- imag
  real(8),allocatable                :: phi_gz(:)           ! gz_projectors diagonal
  
  real(8)                            :: GZ_opt_energy,GZ_opt_kinetic,GZ_opt_Eloc
  real(8),allocatable                :: GZ_opt_Rhop(:,:)         ! reormalization factors 
  real(8),allocatable                :: GZ_opt_projector_diag(:)
  

  !# Slater Determinant #!
  real(8),dimension(:,:,:),allocatable :: slater_matrix_elements   ! (state_dim,state_dim,Lk) aka <c^+_{k\alpha} c_{k\beta}> \alpha,\beta=1,state_dim
  real(8),dimension(:,:),allocatable   :: slater_ground_state_deriv   ! (state_dim,state_dim) aka d<H*>/R_{\alpha,\beta} \alpha,\beta=1,state_dim
  
  !# Observables #!
  real(8),allocatable                :: gz_docc(:)           ! Double occupancy (1:Norb)
  real(8),allocatable                :: gz_dens_dens_orb(:,:)  ! density_density different orbitals (1:Norb,1:Norb)
  real(8),allocatable                :: gz_dens(:)           ! local density    (1:state_dim)



  real(8),allocatable                :: docc(:,:,:)              ! Double occupancy (Norb,:,:)
  real(8),allocatable                :: dens(:,:,:)              ! local density    (1:state_dim,:,:)
  real(8),allocatable                :: dens_dens_orb(:,:,:,:)     ! local density    (Norb,Norb,:,:)
  real(8),allocatable                :: spin_flip(:,:,:,:)     ! spin_flip    (Norb,Norb,:,:)
  real(8),allocatable                :: pair_hopping(:,:,:,:)     ! pair_hopping    (Norb,Norb,:,:)

  integer :: opt_energy_unit,opt_rhop_unit
  integer :: GZmin_unit,GZmin_unit_,opt_GZ_unit

  !+-------------------------+!
  !+- MODEL DETAILS DETAILS -+!
  !+-------------------------+!
  integer                          :: Nx,Lk
  real(8),dimension(:),allocatable :: epsik,hybik
  real(8),dimension(:),allocatable :: wtk
  real(8)                          :: Wband,Vhyb
  real(8),dimension(:),allocatable :: Eloc
  real(8),dimension(:),allocatable :: atomic_energy_levels
  real(8)                          :: Cfield

  real(8),dimension(:,:,:),allocatable :: Hk_tb ! tight binding input matrix

  

END MODULE GZ_VARS_GLOBAL
