program GUTZ_mb
  USE SCIFOR
  !
  USE DMFT_MISC
  USE DMFT_PARSE_INPUT
  !
  USE GZ_AUX_FUNX
  USE GZ_VARS_GLOBAL
  USE GZ_PROJECTORS
  USE GZ_OPTIMIZED_ENERGY
  !
  implicit none
  !
  !+- hamiltonian details -+!
  integer                            :: ispin,iorb,i,istate,jstate,ifock,jorb
  integer,dimension(:),allocatable   :: fock_vec
!  real(8),dimension(3)               :: GZene  
  real(8),dimension(:),allocatable   :: variational_density_natural
  real(8),dimension(:,:),allocatable :: variational_density_natural_simplex
!  real(8),allocatable,dimension(:)   :: local_density,local_dens_min
!  integer                            :: ix,iy,iz,ik,Nk
  integer                            :: out_unit,iter
  integer                            :: lattice ! 2=square;3=cubic


  !+- PARSE INPUT -+!
  call parse_input_variable(Norb,"Norb","inputGZ.conf",default=1)
  call parse_input_variable(U,"U","inputGZ.conf",default=0.5d0)
  call parse_input_variable(xmu,"XMU","inputGZ.conf",default=0.d0)
  call parse_input_variable(Cfield,"Cfield","inputGZ.conf",default=0.d0)
  call parse_input_variable(Wband,"Wband","inputGZ.conf",default=1.d0)
  !call parse_input_variable(Vhyb,"Vhyb","inputGZ.conf",default=0.5d0)
  call parse_input_variable(Nx,"Nx","inputGZ.conf",default=10)
  call parse_input_variable(lancelot_verbose,"LANCELOT_VERBOSE","inputGZ.conf",default=1)
  call parse_input_variable(amoeba_verbose,"AMOEBA_VERBOSE","inputGZ.conf",default=.false.)
  call parse_input_variable(GZmin_verbose,"GZMIN_VERBOSE","inputGZ.conf",default=.false.)
  call parse_input_variable(min_method,"MIN_METHOD","inputGZ.conf",default='nlep')
  call parse_input_variable(Rseed,"RSEED","inputGZ.conf",default=1.d0)
  call parse_input_variable(Rmix,"RMIX","inputGZ.conf",default=1.d0)
  call parse_input_variable(Niter_self,"NITER_SELF","inputGZ.conf",default=100)
  call parse_input_variable(err_self,"ERR_SELF","inputGZ.conf",default=1.d-10)
  call parse_input_variable(amoeba_min_tol,"AMOEBA_MIN_TOL","inputGZ.conf",default=1.d-12)
  call parse_input_variable(fix_density_minimization,"MIN_FIX_DENSITY","inputGZ.conf",default=.false.)
  call parse_input_variable(lattice,"LAT_DIMENSION","inputGZ.conf",default=3)
  call save_input_file("inputGZ.conf")

  
  call build_lattice_model
  
  call initialize_local_fock_space
  
  call build_local_operators_fock_space

  call build_local_hamiltonian

  !  allocate(slater_matrix_elements(state_dim,state_dim,Lk),slater_ground_state_deriv(state_dim,state_dim))  
  
  !  allocate(fock_vec(state_dim))
  !call get_spin_indep_states
!  do i=1,nFock_indep
!     call bdecomp(fock_indep(i),fock_vec)
!     write(*,'(A,I5,A,20I3)') '|',fock_indep(i),'>',fock_vec(:)
!  end do

  call build_gz_local_traces_diag
  !
  allocate(variational_density_natural_simplex(state_dim+1,state_dim))
  allocate(variational_density_natural(state_dim))
  call initialize_variational_density_simplex(variational_density_natural_simplex)
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
    write(*,*) test_k,de
!    stop
    

    !
    !call get_free_dos(epsik,wtk,file='DOS_free.kgrid')
    !  


  end subroutine build_lattice_model




  
  subroutine print_output(vdm_simplex)
    real(8),dimension(state_dim+1,state_dim) :: vdm_simplex
    integer :: out_unit,istate,iorb
    integer,dimension(state_dim) :: fock_state
    real(8),dimension(state_dim) :: tmp
    real(8) :: deltani,delta_tmp,vdm_tmp

    out_unit=free_unit()
    open(out_unit,file='optimized_projectors.data')
    do ifock=1,nFock
       call bdecomp(ifock,fock_state)
       write(out_unit,'(F18.10,A,I4,A,20I3)') GZ_opt_projector_diag(ifock),'#|',ifock,' >',fock_state(:)
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
    write(out_unit,'(20F18.10)') variational_density_natural(1:state_dim)
    close(out_unit)
    !
    out_unit=free_unit()
    open(out_unit,file='optimized_Rhop_matrix.data')
    do istate=1,state_dim
       tmp(istate)=GZ_opt_Rhop(istate,istate)
       write(out_unit,'(20F18.10)') GZ_opt_Rhop(istate,1:state_dim)
    end do
    write(out_unit,*) ! on the last line store the diagonal elements
    write(out_unit,'(20F18.10)') tmp(1:state_dim)
    close(out_unit)
    !
    out_unit=free_unit()
    open(out_unit,file='optimized_density.data')
    write(out_unit,'(20F18.10)') gz_dens(1:state_dim)
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
    do jstate=1,state_dim+1
       if(jstate.le.state_dim) then
          do istate=1,state_dim
             write(out_unit,'(20F18.10)') vdm_simplex(jstate,istate)
          end do
          if(jstate.le.state_dim) write(out_unit,*)  'x'
       else
          do istate=1,state_dim
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


