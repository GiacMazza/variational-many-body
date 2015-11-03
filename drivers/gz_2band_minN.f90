program GUTZ_mb
  USE SCIFOR
  USE DMFT_MISC
  USE DMFT_PARSE_INPUT
!  USE SF_OPTIMIZE
  USE GZ_AUX_FUNX
  USE GZ_VARS_GLOBAL
  USE GZ_PROJECTORS
  USE GZ_ENERGY_FUNCTIONAL
  USE GZ_MINIMIZE
  USE MIN_AMOEBA
  implicit none
  !
  !+- hamiltonian details -+!
  integer :: ispin,iorb,i,istate,ifock
  integer,dimension(:),allocatable :: fock_vec
  real(8),dimension(3)                :: GZene  

  !+- amoeba variables -+!
  !  real(8),allocatable,dimension(:,:)  :: p
  !  real(8),allocatable,dimension(:)    :: y
  !  real(8)                             :: ftol
  !  integer                             :: np,mp,iter,j,ndim,i_vertex,i_dim,ipol
  !  integer                             :: jorb  

  real(8),allocatable,dimension(:)    :: local_density,local_dens_min
  real(8),allocatable,dimension(:)    :: kx
  !real(8),allocatable,dimension(:)    :: pol_test,pol_out
  !real(8)                             :: GZene,pol
  real(8)                             :: ene_min
  real(8)                             :: ts,test_k,kx_,ky_,kz_
  integer                             :: ix,iy,iz,ik,Nk
  !  integer                             :: Ndens
  integer                             :: out_unit,iter
  !integer                             :: minimization
  integer                             :: lattice ! 2=square;3=cubic



  call parse_input_variable(Norb,"Norb","inputGZ.conf",default=1)
  call parse_input_variable(U,"U","inputGZ.conf",default=0.5d0)
  call parse_input_variable(xmu,"XMU","inputGZ.conf",default=0.d0)
  !call parse_input_variable(n0,"N0","inputGZ.conf",default=1.d0)
  call parse_input_variable(Cfield,"Cfield","inputGZ.conf",default=0.d0)
  call parse_input_variable(Wband,"Wband","inputGZ.conf",default=1.d0)
  call parse_input_variable(Vhyb,"Vhyb","inputGZ.conf",default=0.5d0)
  call parse_input_variable(Nx,"Nx","inputGZ.conf",default=10)
  call parse_input_variable(lancelot_verbose,"LANCELOT_VERBOSE","inputGZ.conf",default=1)
  call parse_input_variable(amoeba_verbose,"AMOEBA_VERBOSE","inputGZ.conf",default=.false.)
  call parse_input_variable(fix_density_minimization,"MIN_FIX_DENSITY","inputGZ.conf",default=.false.)
  call parse_input_variable(lattice,"LAT_DIMENSION","inputGZ.conf",default=3)
  call save_input_file("inputGZ.conf")
  !+- build lattice dispersion -+!

  select case(lattice)
  case(2)
     ts = Wband/8.d0
     Vhyb = Vhyb/8.d0
     !
     nk=Nx+1
     Lk=Nk*Nk
     allocate(epsik(Lk),hybik(Lk),wtk(Lk))
     test_k=0.d0
     ik=0
     do ix=0,Nx
        do iy=0,Nx
           ik=ik+1
           kx_=dble(ix)/dble(Nx)*pi
           ky_=dble(iy)/dble(Nx)*pi
           epsik(ik) = -2.d0*ts*(cos(kx_)+cos(ky_))
           hybik(ik) = Vhyb*(cos(kx_)-cos(ky_))
           wtk(ik)   = 1.d0/dble(Lk)           
           test_k = test_k + fermi(epsik(ik),beta)*wtk(ik)*epsik(ik)
        end do
     end do
  case(3)
     ts = Wband/6.d0
     Vhyb = Vhyb/6.d0
     !
     nk=Nx+1
     !Lk=Nk*(Nk+1)/2
     Lk=Nk*Nk
     Lk=Lk*(Nk)
     allocate(epsik(Lk),hybik(Lk),wtk(Lk))
     test_k=0.d0
     ik=0
     do ix=0,Nx
        do iy=0,Nx
           do iz=0,Nx
              ik=ik+1
              kx_=dble(ix)/dble(Nx)*pi
              ky_=dble(iy)/dble(Nx)*pi
              kz_=dble(iz)/dble(Nx)*pi
              epsik(ik) = -2.d0*ts*(cos(kx_)+cos(ky_)+cos(kz_))
              hybik(ik) = Vhyb*(cos(kx_)-cos(ky_))*cos(kz_)
              !hybik(ik) = Vhyb*(sin(kx_)*sin(ky_))*cos(kz_)
              wtk(ik)   = 1.d0/dble(Lk)           
              test_k = test_k + fermi(epsik(ik),beta)*wtk(ik)*epsik(ik)
           end do
        end do
     end do
  end select
  !
  call get_free_dos(epsik,wtk,file='DOS_free.kgrid')
  !  
  !+- build local fock space -+!
  call build_local_fock
  allocate(slater_matrix_elements(state_dim,state_dim,Lk))  
  allocate(fock_vec(state_dim))
  call get_spin_indep_states
  do i=1,nFock_indep
     call bdecomp(fock_indep(i),fock_vec)
     write(*,'(I2,A,20I3)') fock_indep(i),'|>',fock_vec(:)
  end do

  allocate(Rhop_r(state_dim,state_dim),Rhop_diag(state_dim))
  allocate(Eloc(state_dim))
  select case(Norb)
  case(2)
     do iorb=1,Norb
        do ispin=1,2
           istate=index(ispin,iorb)
           Eloc(istate) = Cfield*0.5d0
           if(iorb.eq.2) Eloc(istate) = -Cfield*0.5d0
        end do
     end do
  case default
     Eloc=0.d0        
  end select

  allocate(local_density(state_dim),local_dens_min(Norb))  

  !+- look for the local density configuration file and read -+!
  call initialize_local_density(local_dens_min)  
  out_unit=free_unit()
  open(out_unit,file='used.density_seed.conf')
  do iorb=1,Norb
     write(out_unit,'(6(F18.10))') local_dens_min(iorb)
  end do
  close(out_unit)
  !+----------------------------------------------------------+!

  !+- MINIMIZE ENERGY -+!
  ene_min=gz_optimization_VS_density(local_dens_min)
  
  !call fmin_cg(p=local_dens_min,f=gz_energy_local_density,iter=iter,fret=ene_min)
  
  
  !<DEBUG
  !  if(.not.allocated(phi_gz)) allocate(phi_gz(nFock))
  !  phi_gz=1.d0/sqrt(dble(nFock))
  !DEBUG>
  !write(*,*) size(phi_gz),nFock
  GZene=gz_ground_state_energy_estimation(phi_gz)
  out_unit=free_unit()  
  open(out_unit,file='energy_minimization.out')
  write(out_unit,'(15(F18.10))') local_dens_min,GZene(:)
  close(out_unit)
  
  out_unit=free_unit()
  open(out_unit,file='restart.density_seed.conf')
  do iorb=1,Norb
     write(out_unit,'(6(F18.10))') local_dens_min(iorb)
  end do
  close(out_unit)

  !+-------------------+!
  
  do iorb=1,Norb
     do ispin=1,2
        istate=index(ispin,iorb)
        local_density(istate) = local_dens_min(iorb)*0.5d0
     end do
  end do
  
  Rhop_diag=gz_Rhop_dens(phi_gz,local_density,CC,CA)
  out_unit=free_unit()  
  open(out_unit,file='GZ_hoppings.out')
  write(out_unit,'(10(F18.10))') Rhop_diag
  close(out_unit)

  out_unit=free_unit()  
  open(out_unit,file='GZ_minimization.out')
  do ifock=1,NFock
     call bdecomp(ifock,fock_vec)
     write(out_unit,'(F18.10,I2,A,10I3)') phi_gz(ifock),ifock,'|>',fock_vec(:)
  end do
  close(out_unit)


CONTAINS

  function stupid_function(x) result(f)
    implicit none
    real(8) :: x(:)
    real(8) :: f
    integer :: N,i
    N=size(x);f=0
    do i=1,N
       f=f+x(i)*x(i)
    end do
  end function stupid_function
  
  
end program GUTZ_mb



!AMOEBA TEST


