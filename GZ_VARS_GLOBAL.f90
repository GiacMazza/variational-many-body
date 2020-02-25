MODULE GZ_VARS_GLOBAL
  USE GZ_VARS_INPUT
  USE SF_CONSTANTS
  USE MATRIX_SPARSE
  implicit none
  !
  !+- LOCAL FOCK SPACE -+!
  integer                            :: nFock               ! dimension of the local Fock space = 2^{2*Norb}   
  integer                            :: Ns                  ! number of local energy levels = 2*Norb
  integer                              :: state_dim           ! dimension of a single Fock state  |(\up,\dn)_1,...,(\up,\dn)_Norb> ===> Ns = 2*Norb  
  integer,dimension(:,:),allocatable :: index               ! ispin,iorb  to  istate=1,Ns

  integer :: Nh_2p
  integer,dimension(:),allocatable :: ifk_to_i2p
  integer,dimension(:),allocatable :: i2p_to_ifk

  integer,dimension(:,:,:),allocatable :: i_ios               ! (ispin,iorb,isite) to istate
  

  !
  !# Operators in the space (nFock X nFock) #!
  real(8),allocatable                :: CC(:,:,:),CA(:,:,:)          ! creation annhilation 

  real(8),allocatable                :: CC_(:,:,:),CA_(:,:,:)          ! creation annhilation 

  real(8),allocatable                :: local_hamiltonian(:,:)       ! Hamiltonian of the atomic problem
  real(8),allocatable                :: local_hamiltonian_free(:,:)       ! free Hamiltonian of the atomic problem
  !



  !+---------------------------+!
  !+- GUTZWILLER WAVEFUNCTION -+!
  !+---------------------------+!

  !+- Numbers of dynamical equations -+!
  integer :: nDynamics

  !+- PARAMETERS TO BE OPTIMIZED -+!
  integer :: NRhop_opt,NQhop_opt,Nvdm_NC_opt,Nvdm_NCoff_opt,Nvdm_AC_opt,Nopt_reduced
  integer :: Nsl_normal_opt,Nsl_anomalous_opt
  logical :: optimization_flag
  
  !+- PHI BASIS -+!
  integer                            :: Nphi                
  complex(8),dimension(:,:,:),allocatable :: phi_basis,phi_basis_dag
  type(sparse_matrix_csr_z),allocatable :: phi_basis_sp(:),phi_basis_dag_sp(:)

  
  !# Gutzwiller Traces Matrices #!
  ! hopping renormalization
  complex(8),allocatable :: phi_traces_basis_Rhop(:,:,:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  complex(8),allocatable :: phi_traces_basis_Rhop_hc(:,:,:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  complex(8),allocatable :: phi_traces_basis_Qhop(:,:,:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  complex(8),allocatable :: phi_traces_basis_Qhop_hc(:,:,:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  ! density constraints
  complex(8),allocatable :: phi_traces_basis_dens(:,:,:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  complex(8),allocatable :: phi_traces_basis_dens_hc(:,:,:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  complex(8),allocatable :: phi_traces_basis_dens_anomalous(:,:,:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  complex(8),allocatable :: phi_traces_basis_dens_anomalous_hc(:,:,:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  ! local operators
  complex(8),allocatable :: phi_traces_basis_Hloc(:,:)     ! Single matrix of dimension (nPhi X nPhi) 
  complex(8),allocatable :: phi_traces_basis_free_Hloc(:,:)     ! Single matrix of dimension (nPhi X nPhi) 
  complex(8),allocatable :: phi_traces_basis_local_dens(:,:,:,:)  !Ns x Ns matrices (c+a cb)
  complex(8),allocatable :: phi_traces_basis_dens_dens(:,:,:,:)
  complex(8),allocatable :: phi_traces_basis_spin_flip(:,:,:,:)
  complex(8),allocatable :: phi_traces_basis_pair_hopping(:,:,:,:)
  complex(8),allocatable :: phi_traces_basis_sc_order(:,:,:,:)
  complex(8),allocatable :: phi_traces_basis_spin2(:,:)
  complex(8),allocatable :: phi_traces_basis_spinZ(:,:)
  complex(8),allocatable :: phi_traces_basis_isospin2(:,:)
  complex(8),allocatable :: phi_traces_basis_isospinZ(:,:)
  !
  !sparse_verision
  !
  type(sparse_matrix_csr_z),allocatable :: phi_spTraces_basis_Rhop(:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  type(sparse_matrix_csr_z),allocatable :: phi_spTraces_basis_Rhop_hc(:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  type(sparse_matrix_csr_z),allocatable :: phi_spTraces_basis_Qhop(:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  type(sparse_matrix_csr_z),allocatable :: phi_spTraces_basis_Qhop_hc(:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  ! density constraints
  type(sparse_matrix_csr_z),allocatable :: phi_spTraces_basis_dens(:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  type(sparse_matrix_csr_z),allocatable :: phi_spTraces_basis_dens_hc(:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  type(sparse_matrix_csr_z),allocatable :: phi_spTraces_basis_dens_anomalous(:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  type(sparse_matrix_csr_z),allocatable :: phi_spTraces_basis_dens_anomalous_hc(:,:) ! Ns X Ns matrices of dimension (nPhi X nPhi) 
  ! local operators
  type(sparse_matrix_csr_z),allocatable :: phi_spTraces_basis_Hloc     ! Single matrix of dimension (nPhi X nPhi) 
  type(sparse_matrix_csr_z),allocatable :: phi_spTraces_basis_free_Hloc     ! Single matrix of dimension (nPhi X nPhi) 
  type(sparse_matrix_csr_z),allocatable :: phi_spTraces_basis_local_dens(:,:)  !Ns x Ns matrices (c+a cb)
  type(sparse_matrix_csr_z),allocatable :: phi_spTraces_basis_dens_dens(:,:)
  type(sparse_matrix_csr_z),allocatable :: phi_spTraces_basis_spin_flip(:,:)
  type(sparse_matrix_csr_z),allocatable :: phi_spTraces_basis_pair_hopping(:,:)
  type(sparse_matrix_csr_z),allocatable :: phi_spTraces_basis_sc_order(:,:)
  type(sparse_matrix_csr_z) :: phi_spTraces_basis_spin2
  type(sparse_matrix_csr_z) :: phi_spTraces_basis_spinZ
  type(sparse_matrix_csr_z) :: phi_spTraces_basis_isospin2
  type(sparse_matrix_csr_z) :: phi_spTraces_basis_isospinZ
  
  !+- OPTIMIZED QUANTITIES -+!
  complex(8),allocatable                    :: GZ_vector(:),GZ_init_equ(:)  
  !# Slater Determinant #!
  complex(8),dimension(:,:,:),allocatable   :: GZ_opt_slater   ! (Ns,Ns,Lk) <c^+_{k\alpha} c_{k\beta}> \alpha,\beta=1,Ns
  complex(8),dimension(:,:,:,:),allocatable :: GZ_opt_slater_superc   ! (2,Ns,Ns,Lk) [<c^+_{k\alpha} c_{k\beta}>,<c^+_{k\alpha} c^+_{k\beta}>] \alpha,\beta=1,Ns
  !
  real(8)                                   :: GZ_opt_energy,GZ_opt_kinetic,GZ_opt_Eloc,GZ_opt_free_energy,GZ_opt_entropy
  !
  complex(8),allocatable                    :: GZ_opt_Rhop(:,:)         
  complex(8),allocatable                    :: GZ_opt_Qhop(:,:)         
  !
  complex(8),allocatable                       :: GZ_opt_vdm(:,:)           
  complex(8),dimension(:,:,:),allocatable   :: GZ_Opt_vdm_superc
  !
  complex(8),dimension(:,:),allocatable        :: GZ_opt_slater_lgr(:,:)
  complex(8),dimension(:,:,:),allocatable   :: GZ_opt_slater_lgr_superc(:,:,:)
  complex(8),dimension(:,:),allocatable   :: slater_ground_state_deriv   ! (Ns,Ns) aka d<H*>/R_{\alpha,\beta} \alpha,\beta=1,Ns
  
  !+- projectors lgr_multipliers
  complex(8),dimension(:,:),allocatable     :: GZ_opt_proj_lgr(:,:)
  complex(8),dimension(:,:,:),allocatable   :: GZ_opt_proj_lgr_superc(:,:,:)
  
  
  !# Observables #!
  real(8),allocatable                :: gz_dens_dens(:,:)  !  (1:Ns,1:Ns)
  real(8),allocatable                :: gz_dens(:)           ! local density    (1:Ns)
  real(8),allocatable                :: gz_dens_matrix(:,:)           ! local density    (1:Ns)
  complex(8),allocatable             :: gz_sc_order(:,:)
  real(8)                            :: gz_spin2,gz_spinZ,gz_isospin2,gz_isospinZ
  !
  real(8),allocatable                :: op_dens(:,:,:)              ! local density    (1:Ns,:,:)
  real(8),allocatable                :: op_local_dens(:,:,:,:)              ! local density matrix    (1:Ns,1:Ns,:,:)
  real(8),allocatable                :: op_local_dens_anomalous(:,:,:,:)              ! local density matrix    (1:Ns,1:Ns,:,:)
  !
  real(8),allocatable                :: op_dens_dens(:,:,:,:)     ! density-density    (Ns,Ns,:,:)
  real(8),allocatable                :: op_spin_flip(:,:,:,:)     ! spin_flip    (Norb,Norb,:,:)
  real(8),allocatable                :: op_pair_hopping(:,:,:,:)     ! pair_hopping    (Norb,Norb,:,:)
  real(8),allocatable                :: op_sc_order(:,:,:,:)     ! SC-order parameter    (Norb,Norb,:,:)

  complex(8),allocatable             :: op_spin2(:,:)
  complex(8),allocatable             :: op_spinZ(:,:)
  complex(8),allocatable             :: op_isospin2(:,:)
  complex(8),allocatable             :: op_isospinZ(:,:)


  integer :: opt_energy_unit,opt_rhop_unit,opt_qhop_unit
  integer :: GZmin_unit,GZmin_unit_,opt_GZ_unit

  !+------------+!
  !+- DYNAMICS -+!
  !+------------+!
  !
  !# time grids
  integer                                :: Nt_aux
  integer                                :: Ntgf,Nt0,Nttgf
  real(8),dimension(:),allocatable       :: t_grid
  real(8),dimension(:),allocatable       :: t_grid_aux
  !# time dependent hamiltonian parameters 
  real(8),dimension(:,:),allocatable     :: Uloc_t
  real(8),dimension(:),allocatable       :: Ust_t
  real(8),dimension(:),allocatable       :: Jh_t
  real(8),dimension(:),allocatable       :: Jsf_t,Jph_t
  real(8),dimension(:,:,:,:),allocatable :: Hk_tb_t 
  complex(8),dimension(:,:,:,:),allocatable :: read_neq_lgr
  real(8),dimension(:,:),allocatable     :: eLevels_t 
  
  !
  real(8),dimension(:),allocatable :: Ubcs_t,kdiss_t
  !
  integer,dimension(:,:),allocatable     :: print_grid_Rhop
  integer,dimension(:,:),allocatable     :: print_grid_Qhop
  integer,dimension(:,:),allocatable     :: print_grid_SC
  !

  logical,dimension(:,:),allocatable :: Rgrid,Qgrid,Ngrid


  
  !+-------------------------+!
  !+- MODEL DETAILS DETAILS -+!
  !+-------------------------+!
  real(8)                          :: Wband
  integer :: Lk
  real(8),dimension(:),allocatable :: wtk
  real(8),dimension(:),allocatable :: atomic_energy_levels
  real(8),dimension(:),allocatable :: eLevels
  real(8)                          :: Cfield !+- to be defined in the driver
  real(8),dimension(:,:,:),allocatable :: Hk_tb ! tight binding input matrix
  real(8),dimension(:),allocatable :: wr


  !+---> STRIDES <-----+!
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
  procedure (opt_stride_vec2mat),pointer :: sl_normal_stride_v2m => null()
  procedure (opt_stride_vec2mat),pointer :: sl_anomalous_stride_v2m => null ()
  procedure (opt_stride_vec2mat),pointer :: Rhop_stride_v2m => null ()
  procedure (opt_stride_vec2mat),pointer :: Qhop_stride_v2m => null ()  
  procedure (opt_stride_vec2mat),pointer :: vdm_NC_stride_v2m => null ()
  procedure (opt_stride_vec2mat),pointer :: vdm_NCoff_stride_v2m => null ()
  procedure (opt_stride_vec2mat),pointer :: vdm_AC_stride_v2m => null ()
  !
  !
  procedure (opt_stride_mat2vec),pointer :: sl_normal_stride_m2v => null ()
  procedure (opt_stride_mat2vec),pointer :: sl_anomalous_stride_m2v => null ()
  procedure (opt_stride_mat2vec),pointer :: Rhop_stride_m2v => null ()
  procedure (opt_stride_mat2vec),pointer :: Qhop_stride_m2v => null ()  
  procedure (opt_stride_mat2vec),pointer :: vdm_NC_stride_m2v => null ()
  procedure (opt_stride_mat2vec),pointer :: vdm_NCoff_stride_m2v => null ()
  procedure (opt_stride_mat2vec),pointer :: vdm_AC_stride_m2v => null ()
  !
  integer,dimension(:),allocatable :: IS_vdmAC,JS_vdmAC

  abstract interface
     subroutine opt_stride_iv2im(iv,imI,imJ)
       integer :: iv
       integer :: imI,imJ
     end subroutine opt_stride_iv2im
  end interface
  abstract interface
     subroutine opt_stride_im2iv(imI,imJ,iv)
       integer :: iv
       integer :: imI,imJ
     end subroutine opt_stride_im2iv
  end interface
  !
  procedure(opt_stride_iv2im),pointer :: slNi_v2m => null()
  procedure(opt_stride_iv2im),pointer :: slAi_v2m => null()
  procedure(opt_stride_im2iv),pointer :: slNi_m2v => null()
  procedure(opt_stride_im2iv),pointer :: slAi_m2v => null()
  !
  abstract interface
     subroutine zeros_stride(x1,x2)
       real(8),dimension(:) :: x1
       real(8),dimension(:) :: x2
     end subroutine zeros_stride
  end interface
  procedure (zeros_stride),pointer :: stride_zeros_red2orig => null ()
  procedure (zeros_stride),pointer :: stride_zeros_orig2red => null ()
  
  abstract interface
     subroutine get_Hk(Hkmat,ik,time)
       complex(8),dimension(:,:) :: Hkmat
       integer                   :: ik
       real(8)                   :: time
     end subroutine get_Hk
  end interface
  procedure (get_Hk),pointer  :: get_Hk_t => null ()
 
END MODULE GZ_VARS_GLOBAL
