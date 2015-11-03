MODULE GZ_VARS_GLOBAL
  implicit none
  !  private

  !+----------------+!
  !+- BAND DETAILS -+!
  !+----------------+!
  integer                          :: Nx,Lk
  real(8),dimension(:),allocatable :: epsik,hybik
  real(8),dimension(:),allocatable :: wtk
  real(8)                          :: Wband,Vhyb
  real(8)                          :: xmu
  real(8),dimension(:),allocatable :: Eloc
  real(8),parameter                :: beta=500.d0
  
  real(8) :: U,Cfield,Nread
  
  !+--------------------+!
  !+- LOCAL FOCK SPACE -+!
  !+--------------------+!
  integer,public                            :: nFock ! dimension of the local Fock space
  integer,public                            :: nFock_indep ! dimension of independent fock space
  integer,public                            :: state_dim ! dimension of a single Fock state    
  integer,public                            :: Norb  ! number of orbitals
  integer,dimension(:,:),allocatable,public :: index
  integer,dimension(:),allocatable,public   :: fock_indep
  integer,dimension(:),allocatable,public   :: full2indep_fock
  integer,dimension(:,:),allocatable,public :: indep2full_fock

  !# Operators #!
  real(8),public,allocatable                :: CC(:,:,:),CA(:,:,:) ! creation annhilation 
  real(8),public,allocatable                :: UHubbard(:,:)       ! Hubbard interaction
  real(8),public,allocatable                :: docc(:,:,:)         ! Double occupancy (Norb,:,:)
  real(8),public,allocatable                :: dens(:,:,:)         ! local density
  !# Gutzwiller stuff #!
  complex(8),public,allocatable             :: Rhop_c(:,:)         ! renormalization factors real
  real(8),public,allocatable                :: Rhop_r(:,:)         ! reormalization factors imag  
  real(8),public,allocatable                :: Rhop_diag(:)         ! reormalization factors diagonal imag  
  complex(8),public,allocatable             :: Phi_c(:,:)          ! gz_projectors -full- real
  real(8),public,allocatable                :: Phi_r(:,:)          ! gz_projectors -full- imag
  real(8),public,allocatable                :: phi_gz(:)              ! gz_projectors diagonal

  
  !# Observables #!
  real(8),public,allocatable                :: gz_docc(:)         ! Double occupancy (1:Norb)
  real(8),public,allocatable                :: gz_dens(:)         ! local density    (1:Norb*Nspin)

  
  real(8),dimension(:),allocatable,public :: mu_slater
  

  real(8),public :: Estar_iteration
  real(8),public,dimension(:,:,:),allocatable :: slater_matrix_elements


  logical :: fix_density_minimization

  integer :: lancelot_verbose
  logical :: amoeba_verbose

END MODULE GZ_VARS_GLOBAL
