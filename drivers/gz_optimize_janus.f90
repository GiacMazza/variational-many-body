program GUTZ_mb
  USE SCIFOR
  !
  !USE DMFT_MISC
  !USE DMFT_PARSE_INPUT
  USE DMFT_TOOLS
  !
  USE GZ_AUX_FUNX
  USE GZ_VARS_GLOBAL
  USE GZ_LOCAL_FOCK
  USE GZ_OPTIMIZED_ENERGY
  USE GZ_ENERGY_MINIMIZATION
  USE GZ_EFFECTIVE_HOPPINGS
  !
  USE GZ_MATRIX_BASIS
  !
  implicit none
  !
  !+- hamiltonian details -+!
  integer                            :: ispin,iorb,i,istate,jstate,ifock,jorb
  integer,dimension(:),allocatable   :: fock_vec
  complex(8),dimension(:),allocatable               :: init_vec
  real(8),dimension(:),allocatable   :: variational_density_natural
  real(8),dimension(:),allocatable   :: vdm_init,vdm_out
  complex(8),dimension(:),allocatable   :: R_init,R_out
  complex(8),dimension(:),allocatable   :: Rhop_init,Rhop_out
  real(8),dimension(:,:),allocatable :: variational_density_natural_simplex

  integer                            :: out_unit,iter
  integer                            :: lattice ! 2=square;3=cubic
  real(8),dimension(:),allocatable :: epsik,hybik
  integer :: Nx,is,xmu_unit
  !
  real(8) :: tmp_emin,Uiter,tmp_ene,orb_pol,Jh_ratio,N_target,emin,search_mu
  real(8) :: xmu1,xmu2

  real(8),dimension(:),allocatable :: target_vdm,init_mu

  complex(8),dimension(:,:),allocatable :: Rhop_init_matrix
  real(8),dimension(:,:),allocatable    :: variational_density_matrix



  character(len=5) :: dir_suffix
  character(len=6) :: dir_iter

  !
  !+- PARSE INPUT DRIVER -+!
  call parse_input_variable(Nx,"Nx","inputGZ.conf",default=10)
  call parse_input_variable(Cfield,"Cfield","inputGZ.conf",default=0.d0)
  call parse_input_variable(Wband,"Wband","inputGZ.conf",default=1.d0)
  call parse_input_variable(lattice,"LAT_DIMENSION","inputGZ.conf",default=3)  
  call parse_input_variable(N_target,"N_TARGET","inputGZ.conf",default=3.d0)  
  !
  call read_input("inputGZ.conf")
  call save_input_file("inputGZ.conf")

  if(Norb.eq.1.and.wf_symmetry.eq.1) then
     write(*,*) 'WARNING THE O(1) x SU(2)c x ORBITAL_ROTATION = O(1) x SU(2)c for the Norb=1 case!'
     wf_symmetry=0
  end if
  !
  !
  !NOTE: ON HUNDS COUPLINGS: (see arXiv 1207.3033)
  !
  !NORB=3 ROTATIONAL INVARIANT HAMILTONIAN       :: Jsf=Jh, Jph=U-Ust-J (NO relation between Ust and U)
  !       FULLY ROTATIONAL INVARIANT HAMILTONIAN :: Jsf=Jh, Jph=J, Ust = U - 2J   
  !
  !NORB=2 FULLY ROTATIONAL INVARIANCE            :: Jsf=Jh, Jph=0, Ust = U - J  
  !       PARTIAL O(2) INVARIANCE                :: Jsf=Jh, Jph=0, Ust = U
  !       ONLY SPIN SU(2) INVARIANCE             :: Jsf=Jh, Jph=J, Ust = U - 2J  
  !  
  Jh = Jh*Uloc(1)
  Jsf = Jh
  Jph = Jh
  Ust = Uloc(1)-2.d0*Jh
  !
  call initialize_local_fock_space
  !
  Nopt_diag=1
  Nopt_odiag=0  
  allocate(opt_map(Ns,Ns))  
  opt_map = 0
  do ispin=1,2
     do iorb=1,Norb
        is=index(ispin,iorb)
        opt_map(is,is) = 1
     end do
  end do
  allocate(lgr_init_slater(Nopt_diag+Nopt_odiag)); lgr_init_slater=0.d0
  allocate(lgr_init_gzproj(Nopt_diag+Nopt_odiag)); lgr_init_gzproj=0.d0
  !
  call init_variational_matrices
  !  
  call build_lattice_model
  !
  allocate(target_vdm(Ns))
  target_vdm = N_target/dble(Ns) 
  !+- minimize with respect to XMU the energy at fixed variational density -+!

  allocate(init_mu(1))

  optimization_flag=.true.
  if(.not.allocated(GZ_vector)) allocate(GZ_vector(Nphi))     
  if(.not.allocated(GZ_opt_slater_lgr)) allocate(GZ_opt_slater_lgr(Ns,Ns))     
  !
  Uiter = 0.d0
  Jh_ratio=Jh

  init_mu=-0.85d0


  allocate(Rhop_init_matrix(Ns,Ns)); Rhop_init_matrix = 0.d0;
  allocate(variational_density_matrix(Ns,Ns)); variational_density_matrix = 0.d0;
  do is=1,Ns
     variational_density_matrix(is,is) = 0.5d0
  end do
  do is=1,Ns
     Rhop_init_matrix(is,is) = 1.d0
  end do
  

  xmu1=-0.9d0
  xmu2=-0.8d0
  ! xmu_unit=free_unit()
  ! open(xmu_unit,file='bracket_XMU.out') 
  ! call bracket_density(xmu1,xmu2)
  ! write(xmu_unit,*)  xmu1,xmu2
  ! close(xmu_unit)
  ! stop

  do i=1,40
     !
     do iorb=1,Norb
        Uloc(iorb) = Uiter
     end do
     Jh = Jh_ratio*Uiter
     Jsf = Jh
     Jph = Jh
     Ust = Uiter-2.d0*Jh
     !
     call build_local_hamiltonian     
     phi_traces_basis_Hloc = get_traces_basis_phiOphi(local_hamiltonian)
     phi_traces_basis_free_Hloc = get_traces_basis_phiOphi(local_hamiltonian_free)
     !
     write(dir_suffix,'(F4.2)') Uiter
     dir_iter="U"//trim(dir_suffix)
     call system('mkdir -v '//dir_iter)     
     
     !+- SEARCH XMU TO FIX DENSITY TO THE TARGET VALUE -+!
     !
     xmu_unit=free_unit()
     open(xmu_unit,file='bracket_XMU.out') 
     call bracket_density(xmu1,xmu2)     
     write(xmu_unit,*)  xmu1,xmu2
     close(xmu_unit)
     !
     xmu_unit=free_unit()
     open(xmu_unit,file='search_XMU.out') 
     search_mu=fzero_brentq(gz_optimized_density_VS_xmu,xmu1,xmu2)       
     close(xmu_unit)     
     !
     xmu_unit=free_unit()
     open(xmu_unit,file='optimized_XMU.out') 
     write(xmu_unit,*) search_mu,N_target;
     close(xmu_unit)     
     !
     xmu=search_mu; call build_local_hamiltonian       
     phi_traces_basis_Hloc = get_traces_basis_phiOphi(local_hamiltonian)
     phi_traces_basis_free_Hloc = get_traces_basis_phiOphi(local_hamiltonian_free)  
     !+-------------------------------------------------+!

     !+- GZ OPTIMIZATION -+!
     Rhop_init_matrix = GZ_opt_Rhop
     variational_density_matrix = GZ_opt_VDM
     call gz_optimization_vdm_Rhop(Rhop_init_matrix,variational_density_matrix)
     call get_gz_ground_state(GZ_vector)
     call print_output
     call system('cp * '//dir_iter)
     call system('rm *.out *.data fort* ')
     Uiter = Uiter + 0.2
     !+-------------------+!
  end do
  !
CONTAINS
  !
  subroutine energy_VS_xmu
    real(8) :: energy
    real(8),dimension(NS) :: vdm
    integer :: i
    vdm = N_target/dble(Ns) 
    xmu=0.0
    do i=1,20
       call build_local_hamiltonian       
       phi_traces_basis_Hloc = get_traces_basis_phiOphi(local_hamiltonian)
       phi_traces_basis_free_Hloc = get_traces_basis_phiOphi(local_hamiltonian_free)  
       
       energy = gz_energy_vdm(vdm)
       write(321,*) xmu,energy
       xmu = xmu - 0.2
    end do
  end subroutine energy_VS_xmu


  subroutine bracket_density(xmu1,xmu2)
    real(8),intent(inout) :: xmu1,xmu2
    real(8)               :: dens1,dens2,xmu,dens,sign_dens,xmu_max,xmu_min
    real(8),dimension(2)  :: xmu_sort
    integer :: i
    
    if(xmu2.gt.xmu1) then
       xmu_sort(1) = xmu2
       xmu_sort(2) = xmu1
    else
       xmu_sort(1) = xmu1
       xmu_sort(2) = xmu2
    end if
       
    xmu_max = xmu_sort(1)
    xmu_min = xmu_sort(2)

    dens1 = gz_optimized_density_VS_xmu(xmu1)
    dens2 = gz_optimized_density_VS_xmu(xmu2)

    sign_dens = dens1*dens2

    write(*,*) 'BRACKET IN',xmu1,xmu2,dens1,dens2
    
    if(sign_dens.gt.0.d0) then
       !+- bracket xmu -+!
       if(dens1.gt.0.d0) then          
          xmu = xmu_min
          do i=1,20
             xmu = xmu - 0.1
             dens = gz_optimized_density_VS_xmu(xmu)
             sign_dens = dens*dens1
             write(*,*) 'loop bracket',xmu,xmu_sort,dens,dens1,dens2,sign_dens
             if(sign_dens.lt.0.d0) then
                ! take the largest between xmu1 and xmu2 and <--- xmu
                xmu_sort(1) = xmu
                exit
             end if
          end do
       else
          xmu = xmu_max
          do i=1,20
             xmu = xmu + 0.1
             dens = gz_optimized_density_VS_xmu(xmu)
             sign_dens = dens*dens1
             write(*,*) 'loop bracket',xmu,xmu_sort,dens,dens1,dens2
             if(sign_dens.lt.0.d0) then
                ! take the smallest between xmu1 and xmu2 and <--- xmu
                xmu_sort(2) = xmu
                exit             
             end if
          end do
       end if
    end if
    xmu1=xmu_sort(1)
    xmu2=xmu_sort(2)
    write(*,*) 'BRACKET OUT',xmu1,xmu2,dens,dens1,dens2
  end subroutine bracket_density

  function gz_optimized_density_VS_xmu(chem_pot) result(delta_local_density)
    real(8),intent(in) :: chem_pot
    real(8)              :: delta_local_density,local_density
    real(8)              :: Energy
    complex(8),dimension(:,:),allocatable   :: Rhop_init_matrix
    real(8),dimension(:,:),allocatable   :: variational_density_matrix
    real(8),dimension(1)   :: vdm_in,vdm_out
    integer :: is,js
    !
    xmu = chem_pot
    call build_local_hamiltonian    
    phi_traces_basis_Hloc = get_traces_basis_phiOphi(local_hamiltonian)
    phi_traces_basis_free_Hloc = get_traces_basis_phiOphi(local_hamiltonian_free)

    allocate(Rhop_init_matrix(Ns,Ns)); Rhop_init_matrix = 0.d0;
    do is=1,Ns
       Rhop_init_matrix(is,is) = 1.d0
    end do
    allocate(variational_density_matrix(Ns,Ns)); variational_density_matrix = 0.d0;
    do is=1,Ns
       variational_density_matrix(is,is) = dble(N_target)/dble(2*Ns)
    end do
    !
    lgr_init_slater=0.d0
    lgr_init_gzproj=0.d0
    !
    !vdm_init=0.5d0
    !call gz_optimization_vdm_nlsq(vdm_in,vdm_out)
    !
    !call gz_optimization_vdm_Rhop(Rhop_init_matrix,variational_density_matrix)
    call gz_optimization_vdm_Rhop(Rhop_init_matrix,variational_density_matrix)
    call get_gz_ground_state(GZ_vector)
    !
    variational_density_matrix = GZ_opt_VDM
    local_density=0.d0
    do is=1,Ns
       local_density = local_density + variational_density_matrix(is,is)
    end do
    !
    delta_local_density = local_density-N_target
    !write(*,*) xmu,delta_local_density,local_density
    write(xmu_unit,*) xmu,delta_local_density,local_density,N_target

  end function gz_optimized_density_VS_xmu


  function gz_optimized_density_VS_xmu_(chem_pot) result(delta_local_density)
    real(8),dimension(:)                  :: chem_pot
    real(8),dimension(size(chem_pot))     :: delta_local_density,local_density
    real(8)                               :: Energy
    ! complex(8),dimension(:,:),allocatable :: Rhop_init_matrix
    ! real(8),dimension(:,:),allocatable    :: variational_density_matrix
    integer :: is,js
    !
    xmu = chem_pot(1)
    call build_local_hamiltonian    
    phi_traces_basis_Hloc = get_traces_basis_phiOphi(local_hamiltonian)
    phi_traces_basis_free_Hloc = get_traces_basis_phiOphi(local_hamiltonian_free)

    !
    lgr_init_slater=0.d0
    lgr_init_gzproj=0.d0

    write(*,*) 'enter xmu',xmu
    call gz_optimization_vdm_Rhop(Rhop_init_matrix,variational_density_matrix)
    write(*,*) 'NaN di venerdi sera..'
    call get_gz_ground_state(GZ_vector)
    write(*,*) 'shshshshs'

    Rhop_init_matrix = GZ_opt_Rhop
    !variational_density_matrix = GZ_opt_VDM

    variational_density_matrix = GZ_opt_VDM
    local_density=0.d0
    do is=1,Ns
       local_density = local_density + variational_density_matrix(is,is)
    end do
    !

    delta_local_density = local_density-N_target
    write(*,*) xmu,delta_local_density,local_density
    write(xmu_unit,*) xmu,delta_local_density,local_density,N_target

  end function gz_optimized_density_VS_xmu_





  subroutine build_lattice_model  
    !
    integer :: ix,iy,iz,ik,Nk
    real(8),allocatable,dimension(:)   :: kx
    real(8)                            :: ts,test_k,kx_,ky_,kz_,wini,wfin,de,test_n1,test_n2
    !

    Lk=Nx
    allocate(epsik(Lk),wtk(Lk),hybik(Lk))

    wini=-Wband/2.d0
    wfin= Wband/2.d0
    epsik=linspace(wini,wfin,Lk,mesh=de)
    !
    test_k=0.d0
    test_n1=0.d0;test_n2=0.d0
    do ix=1,Lk
       wtk(ix)=4.d0/Wband/pi*sqrt(1.d0-(2.d0*epsik(ix)/Wband)**2.d0)*de
       !wtk(ix) = 1.d0/Wband*de
       if(ix==1.or.ix==Lk) wtk(ix)=0.d0
       test_n1=test_n1+wtk(ix)*fermi(epsik(ix)+Cfield*0.5d0,beta)
       test_n2=test_n2+wtk(ix)*fermi(epsik(ix)-Cfield*0.5d0,beta)
       write(77,*) epsik(ix),wtk(ix)
    end do
    hybik=0.d0
    !write(*,*) test_n1,test_n2,Cfield; stop


    ! allocate(kx(Nx))
    ! kx = linspace(0.d0,pi,Nx,.true.,.true.)
    ! Lk=Nx*Nx*Nx
    ! allocate(epsik(Lk),wtk(Lk),hybik(Lk))
    ! ik=0
    ! do ix=1,Nx
    !    do iy=1,Nx
    !       do iz=1,Nx
    !          ik=ik+1
    !          !kx_=dble(ix)/dble(Nx)*pi
    !          epsik(ik) = -2.d0/6.d0*(cos(kx(ix))+cos(kx(iy))+cos(kx(iz))) 
    !          hybik(ik) = 0.1d0/6.d0*(cos(kx(ix))-cos(kx(iy)))*cos(kx(iz)) 
    !          wtk(ik) = 1.d0/dble(Lk)
    !       end do
    !    end do
    ! end do

    call get_free_dos(epsik,wtk,file='DOS_free.kgrid')
    !stop

    allocate(Hk_tb(Ns,Ns,Lk))    
    Hk_tb=0.d0
    do ik=1,Lk
       do iorb=1,Norb
          do jorb=1,Norb
             do ispin=1,2
                istate=index(ispin,iorb)
                jstate=index(ispin,jorb)
                Hk_tb(istate,istate,ik) = epsik(ik)
                if(iorb.eq.jorb)  then
                   Hk_tb(istate,jstate,ik) = epsik(ik)
                else
                   Hk_tb(istate,jstate,ik) = hybik(ik)
                end if
             end do
          end do
       end do
    end do

    !<EXTREMA RATIO TEST
    e0test=0.d0
    do ik=1,Lk
       e0test = e0test + fermi_zero(epsik(ik),0.d0)*epsik(ik)*wtk(ik)
    end do
    !EXTREMA RATIO TEST>



  end subroutine build_lattice_model





  subroutine print_output(vdm)
    real(8),dimension(Ns,Ns),optional :: vdm
    integer :: out_unit,istate,iorb,iphi,ifock,jfock
    integer,dimension(Ns) :: fock_state
    real(8),dimension(Ns) :: tmp
    real(8) :: deltani,delta_tmp,vdm_tmp

    real(8),dimension(nFock,nFock) :: test_full_phi

    out_unit=free_unit()
    open(out_unit,file='optimized_projectors.data')

    !+- CHANGE THE NAME OF GZ_opt_projector_diag -+!
    test_full_phi=0.d0
    do iphi=1,Nphi
       ! write(*,*) iphi
       ! do ifock=1,nFock
       !    write(*,'(20F7.1)') phi_basis(iphi,ifock,:)
       ! end do
       ! write(*,*)
       !
       test_full_phi = test_full_phi + GZ_vector(iphi)*phi_basis(iphi,:,:)
    end do
    do ifock=1,nFock
       do jfock=1,nFock
          write(out_unit,*) test_full_phi(ifock,jfock),ifock,jfock
       end do
    end do

    do iphi=1,Nphi
       write(*,*) GZ_vector(iphi)
    end do
    close(out_unit)    




    !
    out_unit=free_unit()
    open(out_unit,file='optimized_internal_energy.data')
    write(out_unit,'(5F18.10)') GZ_opt_energy,GZ_opt_kinetic,GZ_opt_Eloc
    close(out_unit)
    !
    out_unit=free_unit()
    open(out_unit,file='optimized_variational_density_matrix.data')
    !write(out_unit,'(20F18.10)') variational_density_natural(1:Ns)
    do istate=1,Ns
       tmp(istate)=GZ_opt_VDM(istate,istate)
       write(out_unit,'(20F18.10)') GZ_opt_VDM(istate,:)
    end do
    write(out_unit,*) ! on the last line store the diagonal elements
    write(out_unit,'(20F18.10)') tmp(1:Ns)
    close(out_unit)
    !
    out_unit=free_unit()
    open(out_unit,file='optimized_Rhop_matrix.data')
    do istate=1,Ns
       tmp(istate)=GZ_opt_Rhop(istate,istate)
       write(out_unit,'(20F18.10)') GZ_opt_Rhop(istate,1:Ns)
    end do
    write(out_unit,*) ! on the last line store the diagonal elements
    write(out_unit,'(20F18.10)') tmp(1:Ns)
    close(out_unit)
    !
    out_unit=free_unit()
    open(out_unit,file='optimized_density.data')
    write(out_unit,'(20F18.10)') gz_dens(1:Ns)
    close(out_unit)
    !
    out_unit=free_unit()
    open(out_unit,file='orbital_double_occupancy.data')
    write(out_unit,'(20F18.10)') gz_docc(1:Norb)
    close(out_unit)
    !
    out_unit=free_unit()
    open(out_unit,file='orbital_density_density.data')
    do iorb=1,Norb
       write(out_unit,'(20F18.10)') gz_dens_dens_orb(iorb,:)
    end do
    close(out_unit)
    !

    if(present(vdm)) then
       out_unit=free_unit()
       open(out_unit,file='vdm_seed.restart')
       do istate=1,Ns
          write(out_unit,'(20F18.10)') vdm(istate,istate)
       end do
       close(out_unit)
    end if
    !
  end subroutine print_output


end program GUTZ_mb



!AMOEBA TEST


