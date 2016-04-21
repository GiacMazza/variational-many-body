MODULE GZ_VARS_GLOBAL
  USE GZ_VARS_INPUT
  USE SF_CONSTANTS
  implicit none

  !+- RELEVANT DIMENSIONS -+!
  integer                            :: nFock               ! dimension of the local Fock space = 2^{2*Norb}   
  integer                            :: Ns                  ! number of local energy levels = 2*Norb
  integer                              :: state_dim           ! dimension of a single Fock state  |(\up,\dn)_1,...,(\up,\dn)_Norb> ===> Ns = 2*Norb  

  !+- VARIATIONAL DENSITY MATRIX DETAILS -+!
  !integer                              :: Nvdm  !number of independent entries in the variational density matrix in the natural basis 1< Nvdm <= Ns
  integer,dimension(:),allocatable     :: vdm_map
  !integer                              :: Nvdm_c!number of independent constraints for the density matrix 1< Nvdm <= Ns*Ns  
  integer,dimension(:,:),allocatable   :: vdm_c_map

  integer :: Nopt_diag,Nopt_odiag

  integer :: Nopt_lgr,Nopt_normal,Nopt_anomalous
  integer,dimension(:,:),allocatable :: opt_map
  integer,dimension(:,:),allocatable :: opt_map_anomalous


  integer :: NRhop_opt,NQhop_opt,Nvdm_NC_opt,Nvdm_NCoff_opt,Nvdm_AC_opt
  !
  integer                            :: Nphi                ! dimension of the matrix basis for the GZprojectors

  !+- phi_basis
  complex(8),dimension(:,:,:),allocatable :: phi_basis,phi_basis_dag



  logical :: optimization_flag


  !+--------------------+!
  !+- LOCAL FOCK SPACE -+!
  !+--------------------+!
  integer,dimension(:,:),allocatable :: index               ! ispin,iorb  to  istate=1,Ns

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
  complex(8),allocatable                :: phi_traces_basis_Rhop_hc(:,:,:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  complex(8),allocatable                :: phi_traces_basis_Qhop(:,:,:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  complex(8),allocatable                :: phi_traces_basis_Qhop_hc(:,:,:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 

  ! density constraints
  complex(8),allocatable                :: phi_traces_basis_dens(:,:,:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  complex(8),allocatable                :: phi_traces_basis_dens_hc(:,:,:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  complex(8),allocatable                :: phi_traces_basis_dens_anomalous(:,:,:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  complex(8),allocatable                :: phi_traces_basis_dens_anomalous_hc(:,:,:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 

  ! local operators
  complex(8),allocatable                :: phi_traces_basis_Hloc(:,:)     ! Single matrix of dimension (nPhi X nPhi) 
  complex(8),allocatable                :: phi_traces_basis_free_Hloc(:,:)     ! Single matrix of dimension (nPhi X nPhi) 
  complex(8),allocatable                :: phi_traces_basis_local_dens(:,:,:,:)  !Ns x Ns matrices (c+a cb)
  complex(8),allocatable                            :: phi_traces_basis_dens_dens_orb(:,:,:,:)
  complex(8),allocatable                            :: phi_traces_basis_docc_orb(:,:,:)
  complex(8),allocatable                            :: phi_traces_basis_spin_flip(:,:,:,:)
  complex(8),allocatable                            :: phi_traces_basis_pair_hopping(:,:,:,:)
  complex(8),allocatable                            :: phi_traces_basis_sc_order(:,:,:,:)

  complex(8),allocatable                            :: phi_traces_basis_spin2(:,:)
  complex(8),allocatable                            :: phi_traces_basis_spinZ(:,:)
  complex(8),allocatable                            :: phi_traces_basis_isospin2(:,:)
  complex(8),allocatable                            :: phi_traces_basis_isospinZ(:,:)

  !<init_lgr
  real(8),dimension(:),allocatable :: lgr_init_slater,lgr_init_gzproj



  !+- OPTIMIZED QUANTITIES -+!
  complex(8),allocatable                    :: GZ_vector(:)  !+- CHANGE NAME TO THIS GUY
  !# Slater Determinant #!
  complex(8),dimension(:,:,:),allocatable   :: GZ_opt_slater   ! (Ns,Ns,Lk) aka <c^+_{k\alpha} c_{k\beta}> \alpha,\beta=1,Ns
  complex(8),dimension(:,:,:,:),allocatable :: GZ_opt_slater_superc   ! (Ns,Ns,Lk) aka <c^+_{k\alpha} c_{k\beta}> \alpha,\beta=1,Ns
  !
  real(8)                                   :: GZ_opt_energy,GZ_opt_kinetic,GZ_opt_Eloc
  !
  complex(8),allocatable                    :: GZ_opt_Rhop(:,:)         
  complex(8),allocatable                    :: GZ_opt_Qhop(:,:)         
  !
  complex(8),allocatable                       :: GZ_opt_vdm(:,:)           
  complex(8),dimension(:,:,:),allocatable   :: GZ_Opt_vdm_superc

  complex(8),dimension(:,:),allocatable        :: GZ_opt_slater_lgr(:,:)
  complex(8),dimension(:,:,:),allocatable   :: GZ_opt_slater_lgr_superc(:,:,:)

  
  complex(8),dimension(:,:),allocatable   :: slater_ground_state_deriv   ! (Ns,Ns) aka d<H*>/R_{\alpha,\beta} \alpha,\beta=1,Ns

  !# Observables #!
  real(8),allocatable                :: gz_docc(:)           ! Double occupancy (1:Norb)
  real(8),allocatable                :: gz_dens_dens_orb(:,:)  ! density_density different orbitals (1:Norb,1:Norb)
  real(8),allocatable                :: gz_dens(:)           ! local density    (1:Ns)
  real(8),allocatable                :: gz_dens_matrix(:,:)           ! local density    (1:Ns)

  complex(8),allocatable             :: gz_sc_order(:,:)
  real(8)                :: gz_spin2,gz_spinZ,gz_isospin2,gz_isospinZ




  real(8),allocatable                :: op_docc(:,:,:)              ! Double occupancy (Norb,:,:)
  real(8),allocatable                :: op_dens(:,:,:)              ! local density    (1:Ns,:,:)
  real(8),allocatable                :: op_local_dens(:,:,:,:)              ! local density matrix    (1:Ns,1:Ns,:,:)
  real(8),allocatable                :: op_local_dens_anomalous(:,:,:,:)              ! local density matrix    (1:Ns,1:Ns,:,:)

  real(8),allocatable                :: op_dens_dens_orb(:,:,:,:)     ! local density    (Norb,Norb,:,:)
  real(8),allocatable                :: op_spin_flip(:,:,:,:)     ! spin_flip    (Norb,Norb,:,:)
  real(8),allocatable                :: op_pair_hopping(:,:,:,:)     ! pair_hopping    (Norb,Norb,:,:)
  real(8),allocatable                :: op_sc_order(:,:,:,:)     ! SC-order parameter    (Norb,Norb,:,:)
  
  complex(8),allocatable             :: op_spin2(:,:)
  complex(8),allocatable             :: op_spinZ(:,:)
  complex(8),allocatable             :: op_isospin2(:,:)
  complex(8),allocatable             :: op_isospinZ(:,:)


  integer :: opt_energy_unit,opt_rhop_unit,opt_qhop_unit
  integer :: GZmin_unit,GZmin_unit_,opt_GZ_unit

  !+-------------------------+!
  !+- MODEL DETAILS DETAILS -+!
  !+-------------------------+!
  real(8)                          :: Wband
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



  !+- TESTING POINTER FUCNTIONS -+!
  abstract interface
     subroutine opt_stride_vec2mat(vec,mat)
       complex(8),dimension(:) :: vec
       complex(8),dimension(:,:) :: mat
     end subroutine opt_stride_vec2mat
  end interface
  abstract interface
     subroutine opt_stride_mat2vec(mat,vec)
       complex(8),dimension(:) :: vec
       complex(8),dimension(:,:) :: mat
     end subroutine opt_stride_mat2vec
  end interface
  !
  procedure (opt_stride_vec2mat),pointer :: Rhop_stride_v2m => null ()
  procedure (opt_stride_vec2mat),pointer :: Qhop_stride_v2m => null ()  
  procedure (opt_stride_vec2mat),pointer :: vdm_NC_stride_v2m => null ()
  procedure (opt_stride_vec2mat),pointer :: vdm_NCoff_stride_v2m => null ()
  procedure (opt_stride_vec2mat),pointer :: vdm_AC_stride_v2m => null ()
  !
  !
  procedure (opt_stride_mat2vec),pointer :: Rhop_stride_m2v => null ()
  procedure (opt_stride_mat2vec),pointer :: Qhop_stride_m2v => null ()  
  procedure (opt_stride_mat2vec),pointer :: vdm_NC_stride_m2v => null ()
  procedure (opt_stride_mat2vec),pointer :: vdm_NCoff_stride_m2v => null ()
  procedure (opt_stride_mat2vec),pointer :: vdm_AC_stride_m2v => null ()
  !
END MODULE GZ_VARS_GLOBAL
