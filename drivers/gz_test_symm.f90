program GUTZ_mb
  USE SCIFOR
  !
  USE DMFT_MISC
  USE DMFT_PARSE_INPUT
  !
  USE GZ_AUX_FUNX
  USE GZ_VARS_GLOBAL
  USE GZ_LOCAL_FOCK
  USE GZ_OPTIMIZED_ENERGY
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
  complex(8),dimension(:),allocatable   :: Rhop_init,Rhop_out
  complex(8),dimension(:,:),allocatable   :: Rhop_init_matrix
  real(8),dimension(:,:),allocatable :: variational_density_natural_simplex
  !  real(8),allocatable,dimension(:)   :: local_density,local_dens_min
  !  integer                            :: ix,iy,iz,ik,Nk
  integer                            :: out_unit,iter
  integer                            :: lattice ! 2=square;3=cubic

  real(8),dimension(:),allocatable :: epsik,hybik
  integer :: Nx,is
  real(8)                          :: Wband
  !
  real(8) :: tmp_emin





  !+- PARSE INPUT DRIVER -+!
  call parse_input_variable(Nx,"Nx","inputGZ.conf",default=10)
  call parse_input_variable(Cfield,"Cfield","inputGZ.conf",default=0.d0)
  call parse_input_variable(Wband,"Wband","inputGZ.conf",default=1.d0)
  call parse_input_variable(lattice,"LAT_DIMENSION","inputGZ.conf",default=3)  
  !
  call read_input("inputGZ.conf")
  call save_input_file("inputGZ.conf")
  
  if(Norb.eq.1.and.wf_symmetry.eq.1) then
     write(*,*) 'WARNING THE O(1) x SU(2)c x ORBITAL_ROTATION = O(1) x SU(2)c for the Norb=1 case!'
     wf_symmetry=0
  end if
  

  !rescale Jhund couplings
  Jh = Jh*Uloc(1)
  Jsf = Jh
  Jph = Jh
  Ust = Uloc(1)-2.d0*Jh
  !
  
  call initialize_local_fock_space
  

  
  !+- CONSTRAINTS THE VARIATIONAL DENSITY MATRIX -+!
  Nvdm=1
  allocate(vdm_map(Ns),vdm_c_map(Ns,Ns))
  do is=1,Ns
     vdm_map(is) = 1
  end do
  Nvdm_c = 1
  vdm_c_map=0
  do is=1,Ns
     vdm_c_map(is,is) = 1
  end do
  !+-----------------------------------------------+!

  call init_variational_matrices

  
  
  !call build_gz_local_traces_diag  
  !
  allocate(variational_density_natural_simplex(Ns+1,Ns))
  allocate(variational_density_natural(Ns))
  call initialize_variational_density_simplex(variational_density_natural_simplex)

  call build_lattice_model

  !call gz_optimization_simplex(variational_density_natural_simplex,variational_density_natural)  


  allocate(vdm_init(1),vdm_out(1))

  vdm_init = variational_density_natural_simplex(Ns,:)
  !stop
  !call gz_optimization_vdm(vdm_init,vdm_out)
  call gz_optimization_vdm_nlsq(vdm_init,vdm_out)
  
  write(*,*) "DAJE; quanto e' vera la madonna"


  stop
  

  !variational_density_natural_simplex(1,:)=0.5d0
  !
  !tmp_emin=gz_energy_broyden(variational_density_natural_simplex(1,:))  
  !
  !tmp_emin=gz_energy_recursive_nlep(variational_density_natural_simplex(1,:))
  stop
  !
  !<TEST MINIMIZATION
  ! allocate(vdm_init(Ns),vdm_out(Ns))
  ! allocate(Rhop_init(Ns),Rhop_out(Ns),init_vec(Nphi),Rhop_init_matrix(Ns,Ns))
  ! vdm_init=variational_density_natural_simplex(1,1:Ns)
  ! init_vec=1.d0/sqrt(dble(Nphi))
  ! Rhop_init_matrix=hopping_renormalization_normal(init_vec,vdm_init)  
  ! do is=1,NS
  !    Rhop_init(is) = Rhop_init_matrix(is,is)
  ! end do
  ! Rhop_init=one
  ! call gz_optimization_vdm_Rhop(vdm_init,Rhop_init,vdm_out,Rhop_out)
  ! stop
  !TEST MINIMIZATION
  !

  !
  !
  call gz_optimization_simplex(variational_density_natural_simplex,variational_density_natural)  
  !
  call get_gz_ground_state_estimation(variational_density_natural)
  !
  call print_output(variational_density_natural_simplex)
  !
CONTAINS
  !
  subroutine build_lattice_model  
    !
    integer :: ix,iy,iz,ik,Nk
    real(8),allocatable,dimension(:)   :: kx
    real(8)                            :: ts,test_k,kx_,ky_,kz_,wini,wfin,de
    !

    Lk=Nx
    allocate(epsik(Lk),wtk(Lk),hybik(Lk))
    
    wini=-Wband/2.d0
    wfin= Wband/2.d0
    epsik=linspace(wini,wfin,Lk,mesh=de)
    !
    test_k=0.d0
    do ix=1,Lk
       wtk(ix)=4.d0/Wband/pi*sqrt(1.d0-(2.d0*epsik(ix)/Wband)**2.d0)*de
       !wtk(ix) = 1.d0/Wband*de
       if(ix==1.or.ix==Lk) wtk(ix)=0.d0
       test_k=test_k+wtk(ix)
       write(77,*) epsik(ix),wtk(ix)
    end do
    hybik=0.d0
    ! write(*,*) test_k,de


    ! allocate(kx(Nx))
    ! kx = linspace(0.d0,pi,Nx,.false.,.false.)
    ! Lk=Nx*Nx*Nx
    ! allocate(epsik(Lk),wtk(Lk))
    ! ik=0
    ! do ix=1,Nx
    !    do iy=1,Nx
    !       do iz=1,Nx
    !          ik=ik+1
    !          epsik(ik) = -2.d0/6.d0*(cos(kx(ix))+cos(kx(iy))+cos(kx(iz))) !+- cfr lanata
    !          wtk(ik) = 1.d0/dble(Lk)
    !       end do
    !    end do
    ! end do

    allocate(Hk_tb(Ns,Ns,Lk))    
    Hk_tb=0.d0
    do ik=1,Lk
       do iorb=1,Norb
          do ispin=1,2
             istate=index(ispin,iorb)
             Hk_tb(istate,istate,ik) = epsik(ik)
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




  
  subroutine print_output(vdm_simplex)
    real(8),dimension(Ns+1,Ns) :: vdm_simplex
    integer :: out_unit,istate,iorb,iphi,ifock,jfock
    integer,dimension(Ns) :: fock_state
    real(8),dimension(Ns) :: tmp
    real(8) :: deltani,delta_tmp,vdm_tmp

    real(8),dimension(nFock,nFock) :: test_full_phi

    out_unit=free_unit()
    open(out_unit,file='optimized_projectors.data')

    
    !<BUP
    ! do ifock=1,nFock
    !    call bdecomp(ifock,fock_state)
    !    write(out_unit,'(F18.10,A,I4,A,20I3)') GZ_opt_projector_diag(ifock),'#|',ifock,' >',fock_state(:)
    ! end do
    ! do iphi=1,Nphi
    !    ifock=indep2full_fock(iphi,1)
    !    call bdecomp(ifock,fock_state)
    !    write(out_unit,'(F18.10,A,I4,A,20I3)') GZ_opt_projector_diag(iphi),'#|',ifock,' >',fock_state(:)
    ! end do   
    !BUP>
    
    !+- CHANGE THE NAME OF GZ_opt_projector_diag -+!
    test_full_phi=0.d0
    do iphi=1,Nphi
       write(*,*) iphi
       do ifock=1,nFock
          write(*,'(20F7.1)') phi_basis(iphi,ifock,:)
       end do
       write(*,*)
       !
       test_full_phi = test_full_phi + GZ_opt_projector_diag(iphi)*phi_basis(iphi,:,:)
    end do
    do ifock=1,nFock
       do jfock=1,nFock
          write(out_unit,*) test_full_phi(ifock,jfock),ifock,jfock
       end do
    end do
    
    do iphi=1,Nphi
       write(*,*) GZ_opt_projector_diag(iphi)
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
    write(out_unit,'(20F18.10)') variational_density_natural(1:Ns)
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
    out_unit=free_unit()
    open(out_unit,file='vdm_simplex.restart')
    do jstate=1,Ns+1
       if(jstate.le.Ns) then
          do istate=1,Ns
             write(out_unit,'(20F18.10)') vdm_simplex(jstate,istate)
          end do
          if(jstate.le.Ns) write(out_unit,*)  'x'
       else
          do istate=1,Ns
             deltani=vdm_simplex(jstate,istate)-0.5
             if(deltani.gt.0.d0) then
                delta_tmp=0.9999-vdm_simplex(jstate,istate)
                vdm_tmp=vdm_simplex(jstate,istate)+delta_tmp*0.1
                write(out_unit,'(20F18.10)') vdm_tmp
             else
                delta_tmp=vdm_simplex(jstate,istate)-0.0001
                vdm_tmp=vdm_simplex(jstate,istate)-delta_tmp*0.1
                write(out_unit,'(20F18.10)') vdm_tmp
             end if
          end do
       end if
    end do
    close(out_unit)
    !
  end subroutine print_output


end program GUTZ_mb



!AMOEBA TEST


