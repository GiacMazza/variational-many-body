MODULE GZ_VARS_GLOBAL
  USE GZ_VARS_INPUT
  USE SF_CONSTANTS
  implicit none

  !+- RELEVANT DIMENSIONS -+!
  integer                            :: nFock               ! dimension of the local Fock space = 2^{2*Norb}   
  integer                            :: Ns                  ! number of local energy levels = 2*Norb
  integer                            :: state_dim           ! dimension of a single Fock state  |(\up,\dn)_1,...,(\up,\dn)_Norb> ===> Ns = 2*Norb  
  !
  integer                            :: Nphi                ! dimension of the matrix basis for the GZprojectors

  !+- phi_basis
  complex(8),dimension(:,:,:),allocatable :: phi_basis,phi_basis_dag


  !<DEBUG
  real(8) :: e0test
  !DEBUG>


  
  !+--------------------+!
  !+- LOCAL FOCK SPACE -+!
  !+--------------------+!


  !integer                            :: nFock_indep         ! dimension of the independent Fock space ----> this is related to the constrained minimization w/ Galahad   
  !
  integer,dimension(:,:),allocatable :: index               ! ispin,iorb  to  istate=1,Ns
  !integer,dimension(:),allocatable   :: fock_indep          ! store the fock independent states
  !integer,dimension(:),allocatable   :: full2indep_fock     ! map fock state to the corresponding independent one
  !integer,dimension(:,:),allocatable :: indep2full_fock     ! map independent fock states to the full fock space

  !+---------------------------+!
  !+- GUTZWILLER WAVEFUNCTION -+!
  !+---------------------------+!
  !# Operators in the space (nFock X nFock) #!
  real(8),allocatable                :: CC(:,:,:),CA(:,:,:)          ! creation annhilation 
  real(8),allocatable                :: local_hamiltonian(:,:)       ! Hamiltonian of the atomic problem
  real(8),allocatable                :: local_hamiltonian_free(:,:)       ! free Hamiltonian of the atomic problem





  
  !# Gutzwiller Matrices #!
  ! hopping renormalization
  complex(8),allocatable                :: phi_traces_basis_Rhop(:,:,:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  
  ! density constraints
  complex(8),allocatable                :: phi_traces_basis_dens(:,:,:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  !real(8),allocatable                :: phi_traces_basis_Cdens(:,:,:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi)
  
  ! local operators
  complex(8),allocatable                :: phi_traces_basis_Hloc(:,:)     ! Single matrix of dimension (nPhi X nPhi) 
  complex(8),allocatable                :: phi_traces_basis_free_Hloc(:,:)     ! Single matrix of dimension (nPhi X nPhi) 
  complex(8),allocatable                :: phi_traces_basis_local_dens(:,:,:,:)  !Ns x Ns matrices (c+a cb)
  complex(8),allocatable                            :: phi_traces_basis_dens_dens_orb(:,:,:,:)
  complex(8),allocatable                            :: phi_traces_basis_docc_orb(:,:,:)
  complex(8),allocatable                            :: phi_traces_basis_spin_flip(:,:,:,:)
  complex(8),allocatable                            :: phi_traces_basis_pair_hopping(:,:,:,:)
  

  
  
  ! real(8),allocatable                :: gz_docc(:)           ! Double occupancy (1:Norb)
  ! real(8),allocatable                :: gz_dens_dens_orb(:,:)  ! density_density different orbitals (1:Norb,1:Norb)
  ! real(8),allocatable                :: gz_dens(:)           ! local density    (1:Ns)



  ! real(8),allocatable                :: docc(:,:,:)              ! Double occupancy (Norb,:,:)
  ! real(8),allocatable                :: dens(:,:,:)              ! local density    (1:Ns,:,:)
  ! real(8),allocatable                :: dens_dens_orb(:,:,:,:)     ! local density    (Norb,Norb,:,:)
  ! real(8),allocatable                :: spin_flip(:,:,:,:)     ! spin_flip    (Norb,Norb,:,:)
  ! real(8),allocatable                :: pair_hopping(:,:,:,:)     ! pair_hopping    (Norb,Norb,:,:)
  




  
  !# Gutzwiller renormalization factors #!
  !complex(8),allocatable             :: Rhop_c(:,:)         ! renormalization factors real

  !real(8),allocatable                :: Rhop_diag(:)        ! reormalization factors diagonal imag    
  !complex(8),allocatable             :: Phi_c(:,:)          ! gz_projectors -full- real
  !real(8),allocatable                :: Phi_r(:,:)          ! gz_projectors -full- imag
  !real(8),allocatable                :: phi_gz(:)           ! gz_projectors diagonal
  
  real(8)                            :: GZ_opt_energy,GZ_opt_kinetic,GZ_opt_Eloc
  complex(8),allocatable                :: GZ_opt_Rhop(:,:)         ! reormalization factors  
  complex(8),allocatable                :: GZ_opt_projector_diag(:)  !+- CHANGE NAME TO THIS GUY
  

  !# Slater Determinant #!
  complex(8),dimension(:,:,:),allocatable :: slater_matrix_elements   ! (Ns,Ns,Lk) aka <c^+_{k\alpha} c_{k\beta}> \alpha,\beta=1,Ns
  complex(8),dimension(:,:),allocatable   :: slater_ground_state_deriv   ! (Ns,Ns) aka d<H*>/R_{\alpha,\beta} \alpha,\beta=1,Ns
  
  !# Observables #!
  real(8),allocatable                :: gz_docc(:)           ! Double occupancy (1:Norb)
  real(8),allocatable                :: gz_dens_dens_orb(:,:)  ! density_density different orbitals (1:Norb,1:Norb)
  real(8),allocatable                :: gz_dens(:)           ! local density    (1:Ns)



  real(8),allocatable                :: op_docc(:,:,:)              ! Double occupancy (Norb,:,:)
  real(8),allocatable                :: op_dens(:,:,:)              ! local density    (1:Ns,:,:)
  real(8),allocatable                :: op_local_dens(:,:,:,:)              ! local density matrix    (1:Ns,1:Ns,:,:)
  real(8),allocatable                :: op_dens_dens_orb(:,:,:,:)     ! local density    (Norb,Norb,:,:)
  real(8),allocatable                :: op_spin_flip(:,:,:,:)     ! spin_flip    (Norb,Norb,:,:)
  real(8),allocatable                :: op_pair_hopping(:,:,:,:)     ! pair_hopping    (Norb,Norb,:,:)

  integer :: opt_energy_unit,opt_rhop_unit
  integer :: GZmin_unit,GZmin_unit_,opt_GZ_unit

  !+-------------------------+!
  !+- MODEL DETAILS DETAILS -+!
  !+-------------------------+!
  integer :: Lk
  real(8),dimension(:),allocatable :: wtk
  real(8),dimension(:),allocatable :: atomic_energy_levels 
  real(8)                          :: Cfield !+- to be defined in the driver
  real(8),dimension(:,:,:),allocatable :: Hk_tb ! tight binding input matrix
  !
  !integer                          :: Nx,Lk
  ! real(8),dimension(:),allocatable :: epsik,hybik
  
  ! real(8)                          :: Wband,Vhyb
  ! real(8),dimension(:),allocatable :: Eloc


END MODULE GZ_VARS_GLOBAL
