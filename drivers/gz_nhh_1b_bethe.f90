program GUTZ_mb
  USE SCIFOR
  !
  USE DMFT_MISC
  USE SF_PARSE_INPUT
  USE RK_IDE
  !
  USE GZ_AUX_FUNX
  USE GZ_neqAUX_FUNX
  USE GZ_DYNAMICS
  USE GZ_VARS_GLOBAL
  USE GZ_LOCAL_FOCK
  USE GZ_LOCAL_HAMILTONIAN
  USE GZ_ENERGY_MINIMIZATION
  USE GZ_OPTIMIZED_ENERGY
  USE GZ_MATRIX_BASIS
  !
  implicit none
  real(8),dimension(:),allocatable :: epsik,hybik
  real(8) :: t,r,s,tmpU
  integer :: Nx,out_unit,is,js,ik,it,itt,i,iorb,ispin,igz,j
  integer :: nprint
  !
  character(len=200) :: store_dir,read_dir,read_optWF_dir,read_finSC_dir
  complex(8),dimension(:,:,:),allocatable :: slater_init
  complex(8),dimension(:),allocatable     :: gz_proj_init
  !
  complex(8),dimension(:,:),allocatable :: bcs_wf
  !
  complex(8),dimension(:),allocatable     :: psi_t,psi_t_tmp
  real(8),dimension(:,:),allocatable      :: Ut 
  !
  integer :: unit_neq_hloc
  integer :: unit_neq_local_dens
  integer :: unit_neq_local_dens_dens
  integer :: unit_neq_ene
  integer :: unit_neq_dens_constrSL
  integer :: unit_neq_dens_constrGZ
  integer :: unit_neq_dens_constrSLa
  integer :: unit_neq_dens_constrGZa
  integer :: unit_neq_constrU
  integer :: unit_neq_Rhop
  integer :: unit_neq_Qhop
  integer :: unit_neq_AngMom
  integer :: unit_neq_sc_order
  integer :: unit_neq_nqp
  integer :: unit_neq_lgr_diss
  !
  !+- observables -+!
  complex(8),dimension(:),allocatable   :: Rhop
  complex(8),dimension(:,:),allocatable   :: Rhop_matrix
  complex(8),dimension(:,:,:),allocatable :: slater
  complex(8),dimension(:,:),allocatable :: local_density_matrix
  real(8),dimension(:,:),allocatable    :: local_dens_dens
  complex(8),dimension(:,:),allocatable :: dens_constrSL
  complex(8),dimension(:,:),allocatable :: dens_constrGZ
  real(8)                               :: unitary_constr
  real(8),dimension(4)                  :: local_angular_momenta
  real(8),dimension(3)                  :: energies
  !
  real(8),dimension(:),allocatable      :: dump_vect
  real(8) :: fin_sc_dir
  real(8) :: Uneq,Uneq0,tStart_neqU,tRamp_neqU,tSin_neqU,dUneq
  real(8) :: tStart_kdiss,tRamp_kdiss
  
  real(8),dimension(2) :: diss_lgr,diss_lgr_init
  integer :: iter,fix_qp_lgr
  real(8) :: delta_lgr,diss_lgr_newton
  character(len=3)  :: fix_lgr_method
  !
  call parse_input_variable(Wband,"WBAND","inputGZ.conf",default=2.d0)
  call parse_input_variable(Nx,"Nx","inputGZ.conf",default=1000)
  !call parse_input_variable(read_dir,"READ_GZ_BASIS_DIR","inputGZ.conf",default='~/etc_local/GZ_basis/')
  call parse_input_variable(read_optWF_dir,"EQWF_DIR","inputGZ.conf",default='./')
  call parse_input_variable(nprint,"NPRINT","inputGZ.conf",default=10)  
  !
  call parse_input_variable(Uneq,"Uneq","inputGZ.conf",default=0.d0) 
  call parse_input_variable(Uneq0,"Uneq0","inputGZ.conf",default=0.d0) 
  call parse_input_variable(tStart_neqU,"TSTART_NEQU","inputGZ.conf",default=0.d0)
  call parse_input_variable(tRamp_neqU,"TRAMP_NEQU","inputGZ.conf",default=0.d0)  
  call parse_input_variable(tSin_neqU,"TSIN_NEQU","inputGZ.conf",default=0.5d0)
  call parse_input_variable(dUneq,"DUneq","inputGZ.conf",default=0.d0) 
  !
  call parse_input_variable(tStart_kdiss,"TSTART_KDISS","inputGZ.conf",default=0.d0)
  call parse_input_variable(tRamp_kdiss,"TRAMP_KDISS","inputGZ.conf",default=0.d0)
  call parse_input_variable(fix_lgr_method,"FIX_LGR_METHOD","inputGZ.conf",default='MIN')
  call parse_input_variable(diss_lgr_init,"DISS_LGR_INIT","inputGZ.conf",default=[0.d0,0.1d0])
  call parse_input_variable(fix_qp_lgr,"FIX_QP_LGR","inputGZ.conf",default=1)  
  !
  call read_input("inputGZ.conf")
  call save_input_file("inputGZ.conf")
  !

  
  Norb=1
  wf_symmetry=0
  call initialize_local_fock_space  
  call init_variational_matrices(wf_symmetry)
  !

  

  call build_lattice_model; get_Hk_t => getHk
  allocate(eLevels(Ns)); eLevels=0.d0
  !
  !+- INITIALIZE TIME GRIDS -+!
  Nt_aux=2*Nt+1
  allocate(t_grid(Nt),t_grid_aux(Nt_aux))
  !
  t_grid = linspace(tstart,tstep*real(Nt-1,8),Nt)
  t_grid_aux = linspace(tstart,0.5d0*tstep*real(Nt_aux-1,8),Nt_aux)
  !

  allocate(Ut(3,Nt_aux))  
  do itt=1,Nt_aux
     t = t_grid_aux(itt) 
     !
     if(t.lt.tStart_neqU) then
        r=0.d0
     else
        if(t.lt.tStart_neqU+tRamp_neqU) then
           r = (1.d0 - 1.5d0*cos(pi*(t-tStart_neqU)/tRamp_neqU) + 0.5d0*(cos(pi*(t-tStart_neqU)/tRamp_neqU))**3)*0.5d0
        else
           r = 1.d0 
        end if
     end if
     !
     if(t.lt.tStart_neqU+tRamp_neqU+tSin_neqU) then
        s = 1.d0
     else
        s = 1.d0 + dUneq*dsin(2.d0*pi*t/tSin_neqU)
     end if
     !
     tmpU = Uneq0 + r*(Uneq*s-Uneq0)
     !
     Ut(:,itt) = Uneq0 + r*(Uneq*s-Uneq0)
  end do
  !
  call setup_neq_hamiltonian(Uloc_t_=Ut)
  !
  unit_neq_hloc = free_unit()
  open(unit_neq_hloc,file="neq_local_interaction_parameters.out")
  !
  do itt=1,Nt_aux
     if(mod(itt-1,2*nprint).eq.0) then        
        write(unit_neq_hloc,*) t_grid_aux(itt),Uloc_t(:,itt),Jh_t(itt),Jsf_t(itt),Jph_t(itt),Ust_t(itt)
     end if
  end do
  close(unit_neq_hloc)
  !


  !+- t-ramp of the k-dissipation -+!
  allocate(k2p_loss_t(Nt_aux),k2p_pump_t(Nt_aux),kpump_t(Nt_aux),kloss_t(Nt_aux))
  unit_neq_hloc = free_unit()
  open(unit_neq_hloc,file="neq_kdiss.out")
  do itt=1,Nt_aux
     t = t_grid_aux(itt) 
     !
     if(t.lt.tStart_kdiss) then
        r=0.d0
     else
        if(t.lt.tStart_kdiss+tRamp_kdiss) then
           r = (1.d0 - 1.5d0*cos(pi*(t-tStart_kdiss)/tRamp_kdiss) + 0.5d0*(cos(pi*(t-tStart_kdiss)/tRamp_kdiss))**3)*0.5d0
        else
           r = 1.d0 
        end if
     end if
     !
     k2p_loss_t(itt) = r*k_2p_loss
     k2p_pump_t(itt) = r*k_2p_pump

     kpump_t(itt) = r*k_1p_pump
     kloss_t(itt) = r*k_1p_loss
     if(mod(itt-1,2*nprint).eq.0) then        
        write(unit_neq_hloc,'(5F18.10)') t,k2p_loss_t(itt),k2p_pump_t(itt),kpump_t(itt),kloss_t(itt)
     end if
  end do
  close(unit_neq_hloc)

  

  !
  !+- READ EQUILIBRIUM AND SETUP DYNAMICAL VECTOR -+!
  nDynamics = Ns*Ns*Lk + Nphi
  allocate(psi_t(nDynamics))
  allocate(slater_init(Ns,Ns,Lk),gz_proj_init(Nphi))  
  !
  call read_optimized_variational_wf(read_optWF_dir,slater_init,gz_proj_init)
  call wfMatrix_2_dynamicalVector(slater_init,gz_proj_init,psi_t)

  allocate(psi_t_tmp(nDynamics)) ; psi_t_tmp = psi_t

  !+- STRIDE FOR NEQ-LGR DISS.
  ! integer :: Nlp_diss
  ! integer,dimension(:,:),allocatable :: ilpvec2ij
  ! integer,dimension(:,:,:),allocatable :: ij2ilpvec
  ! Nlp_diss=Ns*(2*Ns-1)
  ! allocate(ilpvec2ij(Nlp_diss,3)) ;  ilpvec2ij=0
  ! allocate(ij2ilpvec(Ns,Ns,2)) ; ij2ilpvec=0
  ! !
  ! i=0
  ! do is=1,Ns
  !    do js=1,Ns
  !       if(is.ne.js) then
  !          i=i+1
  !          ij2ilpvec(is,js,1) = i
  !          ilpvec2ij(i,1) = is
  !          ilpvec2ij(i,2) = js
  !          ilpvec2ij(i,3) = 1
  !       end if
  !    end do
  ! end do
  ! do is=1,Ns
  !    do js=1,Ns
  !       i=i+1
  !       ij2ilpvec(is,js,2) = i
  !       ilpvec2ij(i,1) = is
  !       ilpvec2ij(i,2) = js
  !       ilpvec2ij(i,3) = 2       
  !    end do
  ! end do

  


  
  !
  it=1
  Uloc=Uloc_t(:,it)
  Ust =Ust_t(it)
  Jh=Jh_t(it)
  Jsf=Jsf_t(it)
  Jph=Jph_t(it)
  eLevels = eLevels_t(:,it)
  !
  call get_local_hamiltonian_trace(eLevels)      
  !
  call setup_neq_dynamics
  !    
  unit_neq_Rhop = free_unit()
  open(unit_neq_Rhop,file='neq_Rhop_matrix.data')
  !
  unit_neq_local_dens = free_unit()
  open(unit_neq_local_dens,file='neq_local_density_matrix.data')
  !
  unit_neq_local_dens_dens = free_unit()
  open(unit_neq_local_dens_dens,file='neq_local_dens_dens.data')
  !
  unit_neq_ene = free_unit()
  open(unit_neq_ene,file='neq_energy.data')
  !
  unit_neq_dens_constrSL = free_unit()
  open(unit_neq_dens_constrSL,file='neq_dens_constrSL.data')
  !
  unit_neq_dens_constrGZ = free_unit()
  open(unit_neq_dens_constrGZ,file='neq_dens_constrGZ.data')
  !
  unit_neq_constrU = free_unit()
  open(unit_neq_constrU,file='neq_constrU.data')
  !
  unit_neq_AngMom = free_unit()
  open(unit_neq_AngMom,file='neq_AngMom.data')
  !
  unit_neq_sc_order = free_unit()
  open(unit_neq_sc_order,file='neq_sc_order.data')
  !
  unit_neq_nqp = free_unit()
  open(unit_neq_nqp,file='neq_slater.data')
  !
  unit_neq_lgr_diss = free_unit()
  open(unit_neq_lgr_diss,file='neq_lgr_diss.data')
  !
  allocate(Rhop_matrix(Ns,Ns))
  allocate(local_density_matrix(Ns,Ns))
  allocate(local_dens_dens(Ns,Ns))
  allocate(dens_constrSL(Ns,Ns))
  allocate(dens_constrGZ(Ns,Ns))  
  allocate(dump_vect(Ns*Ns))
  allocate(slater(Ns,Ns,Lk))

  !*) ACTUAL DYNAMICS (simple do loop measuring each nprint times)
  do it=1,Nt
     write(*,*) it,Nt
     !
     t=t_grid(it)
     !
     if(mod(it-1,nprint).eq.0) then        
        !
        call gz_neq_measure(psi_t,t,read_slater=slater)
        !

        !
        if(mod(it-1,10*nprint).eq.0) then
           do ik=1,Lk
              write(unit_neq_nqp,'(10F18.10)') t,epsik(ik),slater(1,1,ik)
           end do
           write(unit_neq_nqp,'(10F18.10)')
           write(unit_neq_nqp,'(10F18.10)')
        end if
        !

        do is=1,Ns
           do js=1,Ns
              call get_neq_Rhop(is,js,Rhop_matrix(is,js))              
              call get_neq_local_dens(is,js,local_density_matrix(is,js))              
              call get_neq_local_dens_dens(is,js,local_dens_dens(is,js))              
              call get_neq_dens_constr_slater(is,js,dens_constrSL(is,js))
              call get_neq_dens_constr_gzproj(is,js,dens_constrGZ(is,js))              
           end do
        end do
        !
        call get_neq_energies(energies)
        call get_neq_local_angular_momenta(local_angular_momenta)
        call get_neq_unitary_constr(unitary_constr)
        !
        call write_complex_matrix_grid(Rhop_matrix,unit_neq_Rhop,print_grid_Rhop,t)
        !
        call write_hermitean_matrix(local_density_matrix,unit_neq_local_dens,t)
        call write_hermitean_matrix(dens_constrSL,unit_neq_dens_constrSL,t)
        call write_hermitean_matrix(dens_constrGZ,unit_neq_dens_constrGZ,t)
        call write_hermitean_matrix(local_density_matrix,unit_neq_local_dens,t)
        call write_symmetric_matrix(local_dens_dens,unit_neq_local_dens_dens,t)
        write(unit_neq_AngMom,'(10F18.10)') t,local_angular_momenta
        write(unit_neq_ene,'(10F18.10)') t,energies
        write(unit_neq_constrU,'(10F18.10)') t,unitary_constr
        write(unit_neq_lgr_diss,'(10F18.10)') t,diss_lgr 
        !        
     end if
     !

     !+- here instead of a direct step I first need to determine the lagrange parameters at each
     ! iteration
     select case(fix_qp_lgr)
     case(0)
        if(fix_lgr_method.eq.'MIN') then
           diss_lgr=diss_lgr_init
           call fmin_cg(diss_lgr,fix_neq_diss_lgr_min,iter,delta_lgr)
        else
           diss_lgr=diss_lgr_init
           call fsolve(fix_neq_diss_lgr,diss_lgr,tol=1d-10,check=.false.)
        end if
     case(1)
        diss_lgr=diss_lgr_init
        diss_lgr(2)=-1.10d0
        ! do j=1,20
        !    diss_lgr(2) = diss_lgr(2) + 0.1d0
        !    write(*,*) fix_neq_diss_lgr_imag(diss_lgr(2))
        ! end do
        diss_lgr(2) = brentq(fix_neq_diss_lgr_imag,-1d0,1d0)
        !call newton(fix_neq_diss_lgr_imag,diss_lgr(2))
        write(*,*) diss_lgr(2)
        !stop

     case(2)

        ! diss_lgr=diss_lgr_init
        ! diss_lgr(1)=-1.10d0
        ! do j=1,20
        !    diss_lgr(1) = diss_lgr(1) + 0.1d0
        !    write(*,*) fix_neq_diss_lgr_real(diss_lgr(1))
        ! end do
        diss_lgr(1) = brentq(fix_neq_diss_lgr_real,-1d0,1d0)
        !call newton(fix_neq_diss_lgr_imag,diss_lgr(2))
        !write(*,*) diss_lgr(2)
        !stop

        
        ! diss_lgr=diss_lgr_init
        ! call newton(fix_neq_diss_lgr_real,diss_lgr(1))
     end select
     !
     !diss_lgr=0d0
     !write(*,*) fix_neq_diss_lgr(diss_lgr)
     lgr_diss_1b = diss_lgr(1)+xi*diss_lgr(2)
     psi_t = RK_step(nDynamics,4,tstep,t,psi_t,gz_equations_of_motion)
     !stop
     !
  end do
  !
  !
  !
CONTAINS

  function fix_neq_diss_lgr(neq_diss_lgr) result(delta)
    real(8),dimension(:),intent(in)            :: neq_diss_lgr
    real(8),dimension(size(neq_diss_lgr))      :: delta
    complex(8),dimension(Ns)  :: cSL,cGZ

    if(size(neq_diss_lgr).ne.2) stop "(size(neq_diss_lgr).ne.2)"
    !+- I need to broadcast the neq_diss to the global variables entering the EOMs
    lgr_diss_1b = neq_diss_lgr(1)+xi*neq_diss_lgr(2)
    
    !+- do the step
    psi_t_tmp = RK_step(nDynamics,4,tstep,t,psi_t,gz_equations_of_motion)
    
    !+- measure the difference
    call gz_neq_measure(psi_t_tmp,t)
    do is=1,Ns
       call get_neq_dens_constr_slater(is,is,cSL(is))
       call get_neq_dens_constr_gzproj(is,is,cGZ(is))              
    end do
    is=1
    delta(1)=dreal(cSL(is)-cGZ(is))
    delta(2)=dimag(cSL(is)-cGZ(is))
    ! write(400,*) neq_diss_lgr,delta
    ! write(500,'(10F18.10)') cSL,cGZ
    !
  end function fix_neq_diss_lgr

  function fix_neq_diss_lgr_min(neq_diss_lgr) result(delta)
    real(8),dimension(:)            :: neq_diss_lgr
    real(8)      :: delta
    complex(8),dimension(Ns)  :: cSL,cGZ

    if(size(neq_diss_lgr).ne.2) stop "(size(neq_diss_lgr).ne.2)"
    !+- I need to broadcast the neq_diss to the global variables entering the EOMs
    lgr_diss_1b = neq_diss_lgr(1)+xi*neq_diss_lgr(2)
    
    !+- do the step
    psi_t_tmp = RK_step(nDynamics,4,tstep,t,psi_t,gz_equations_of_motion)
    
    !+- measure the difference
    call gz_neq_measure(psi_t_tmp,t)
    delta=0d0
    do is=1,Ns
       call get_neq_dens_constr_slater(is,is,cSL(is))
       call get_neq_dens_constr_gzproj(is,is,cGZ(is))
       delta = delta + (cSL(is)-cGZ(is))**2d0/dble(Ns)
    end do
    ! is=1
    ! delta(1)=dreal(cSL(is)-cGZ(is))
    ! delta(2)=dimag(cSL(is)-cGZ(is))
    ! write(700,*) neq_diss_lgr,delta
    ! write(800,'(10F18.10)') cSL,cGZ
    !stop
    !
  end function fix_neq_diss_lgr_min
  !
  !
  !
  function fix_neq_diss_lgr_imag(neq_diss_lgr) result(delta)
    real(8),intent(in)            :: neq_diss_lgr
    real(8)      :: delta
    complex(8),dimension(Ns)  :: cSL,cGZ

    !if(size(neq_diss_lgr).ne.2) stop "(size(neq_diss_lgr).ne.2)"
    !+- I need to broadcast the neq_diss to the global variables entering the EOMs
    lgr_diss_1b = xi*neq_diss_lgr
    
    !+- do the step -+!
    psi_t_tmp = RK_step(nDynamics,4,tstep,t,psi_t,gz_equations_of_motion)
    
    !+- measure the difference -+!
    !slater=0d0
    !call gz_neq_measure(psi_t_tmp,t,read_slater=slater)
    call gz_neq_measure(psi_t_tmp,t)

    ! do ik=1,Lk
    !    write(750,'(10F18.10)') epsik(ik),slater(1,1,ik)
    ! end do
    ! write(750,'(10F18.10)')
    ! write(750,'(10F18.10)')

    
    delta=0d0
    do is=1,Ns
       call get_neq_dens_constr_slater(is,is,cSL(is))
       call get_neq_dens_constr_gzproj(is,is,cGZ(is))
       delta = delta + (cSL(is)-cGZ(is))/dble(Ns)
    end do
    !
    ! write(700,*) neq_diss_lgr,delta 
    ! write(800,'(10F18.10)') cSL,cGZ 
    !
  end function fix_neq_diss_lgr_imag


  function fix_neq_diss_lgr_real(neq_diss_lgr) result(delta)
    real(8),intent(in)            :: neq_diss_lgr
    real(8)      :: delta
    complex(8),dimension(Ns)  :: cSL,cGZ

    !if(size(neq_diss_lgr).ne.2) stop "(size(neq_diss_lgr).ne.2)"
    !+- I need to broadcast the neq_diss to the global variables entering the EOMs
    lgr_diss_1b = neq_diss_lgr
    
    !+- do the step -+!
    psi_t_tmp = RK_step(nDynamics,4,tstep,t,psi_t,gz_equations_of_motion)
    
    !+- measure the difference -+!
    call gz_neq_measure(psi_t_tmp,t)
    delta=0d0
    do is=1,Ns
       call get_neq_dens_constr_slater(is,is,cSL(is))
       call get_neq_dens_constr_gzproj(is,is,cGZ(is))
       delta = delta + (cSL(is)-cGZ(is))/dble(Ns)
    end do
    !
    !+- write(700,*) neq_diss_lgr,delta -+!
    !+- write(800,'(10F18.10)') cSL,cGZ -+!
    !
  end function fix_neq_diss_lgr_real

  
  !
  subroutine build_lattice_model  
    implicit none
    !
    integer                          :: ix,iy,iz,ik,Nk,iorb,jorb,ispin,istate,jstate
    real(8),allocatable,dimension(:) :: kx
    real(8)                          :: ts,test_k,kx_,ky_,kz_,wini,wfin,de,n1,n2
    !
    !
    Lk=Nx
    allocate(epsik(Lk),wtk(Lk),hybik(Lk))    
    wini=-Wband/2.d0
    wfin= Wband/2.d0
    ! wini=-5.d0
    ! wfin= 5.d0
    epsik=linspace(wini,wfin,Lk,mesh=de)
    !
    test_k=0.d0
    do ix=1,Lk
       if(epsik(ix).gt.-Wband/2.d0.and.epsik(ix).lt.Wband/2) then
          wtk(ix)=4.d0/Wband/pi*sqrt(1.d0-(2.d0*epsik(ix)/Wband)**2.d0)*de
       else
          wtk(ix) = 0.d0
       end if
       !wtk(ix) = 1.d0/Wband*de
       !if(ix==1.or.ix==Lk) wtk(ix)=0.d0
       test_k=test_k+wtk(ix)
       write(77,*) epsik(ix),wtk(ix)
    end do
    hybik=0.d0
    ! write(*,*) test_k,de
    !
    ! allocate(kx(Nx))
    ! kx = linspace(0.d0,pi,Nx,.true.,.true.)
    ! Lk=Nx*Nx*Nx
    ! allocate(epsik(Lk),wtk(Lk),hybik(Lk))
    ! ik=0
    ! test_k=0.d0;n1=0.d0;n2=0.d0
    ! do ix=1,Nx
    !    do iy=1,Nx
    !       do iz=1,Nx
    !          ik=ik+1
    !          !kx_=dble(ix)/dble(Nx)*pi
    !          epsik(ik) = -2.d0/6.d0*(cos(kx(ix))+cos(kx(iy))+cos(kx(iz))) 
    !          hybik(ik) = 1.d0/6.d0*(cos(kx(ix))-cos(kx(iy)))*cos(kx(iz)) 
    !          wtk(ik) = 1.d0/dble(Lk)
    !          n1=n1+fermi(epsik(ik)+cfield*0.5,beta)*wtk(ik)
    !          n2=n2+fermi(epsik(ik)-cfield*0.5,beta)*wtk(ik)
    !       end do
    !    end do
    ! end do
    ! !
    ! write(*,*) 'n1/n2'
    ! write(*,*) n1,n2,n1+n2
    call get_free_dos(epsik,wtk,file='DOS_free.kgrid')
    !
    allocate(Hk_tb(Ns,Ns,Lk))
    !
    Hk_tb=0.d0
    do ik=1,Lk
       do iorb=1,Norb
          do jorb=1,Norb
             do ispin=1,2
                istate=index(ispin,iorb)
                jstate=index(ispin,jorb)
                if(iorb.eq.jorb)  then
                   Hk_tb(istate,jstate,ik) = epsik(ik)
                else
                   Hk_tb(istate,jstate,ik) = hybik(ik)
                end if
             end do
          end do
       end do
    end do
    !
  end subroutine build_lattice_model


  subroutine getHk(Hk,ik,time)
    complex(8),dimension(:,:) :: Hk
    integer                   :: ik
    real(8)                   :: time
    if(size(Hk,1).ne.size(Hk,2)) stop "wrong dimenions in getHk"
    if(size(Hk,1).ne.Ns) stop "wrong dimenions in getHk"
    Hk = Hk_tb(:,:,ik)
  end subroutine getHk







  subroutine print_output(vdm_simplex,vdm_opt)
    real(8),dimension(Ns+1,Ns),optional :: vdm_simplex
    real(8),dimension(Nvdm_NC_opt-Nvdm_NCoff_opt),optional :: vdm_opt
    integer :: out_unit,istate,jstate,iorb,iphi,ifock,jfock,is,js,ik
    integer,dimension(Ns) :: fock_state
    complex(8),dimension(Ns) :: tmp
    real(8) :: deltani,delta_tmp,vdm_tmp

    real(8),dimension(nFock,nFock) :: test_full_phi

    !+- STORE GZ PROJECTORS -+!
    out_unit=free_unit()
    open(out_unit,file='optimized_projectors.data')
    test_full_phi=0.d0
    do iphi=1,Nphi
       !
       test_full_phi = test_full_phi + GZ_vector(iphi)*phi_basis(iphi,:,:)
       write(out_unit,'(2F18.10)') GZ_vector(iphi)
    end do
    write(out_unit,*) '!+-----------------------------+!'
    write(out_unit,*) '!+-----------------------------+!'
    write(out_unit,*) '!+-----------------------------+!'
    do ifock=1,nFock
       do jfock=1,nFock
          write(out_unit,*) test_full_phi(ifock,jfock),ifock,jfock
       end do
    end do
    close(out_unit)    

    !+- STORE SLATER DETERMINANT -+!
    do is=1,Ns
       do js=1,Ns
          out_unit=free_unit()
          open(out_unit,file='optimized_slater_normal_IS'//reg(txtfy(is))//'_JS'//reg(txtfy(js))//'.data')
          do ik=1,Lk
             write(out_unit,'(2F18.10)') dreal(GZ_opt_slater(is,js,ik)),dimag(GZ_opt_slater(is,js,ik))
          end do
          close(out_unit)
       end do
    end do


    !
    out_unit=free_unit()
    open(out_unit,file='optimized_internal_energy.data')
    write(out_unit,'(5F18.10)') GZ_opt_energy,GZ_opt_kinetic,GZ_opt_Eloc
    close(out_unit)
    !
    out_unit=free_unit()
    open(out_unit,file='optimized_variational_density_matrix.data')
    write(out_unit,*) 'NORMAL VDM'
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
    do istate=1,Ns
       write(out_unit,'(20F18.10)') gz_dens_matrix(istate,1:Ns)
    end do
    write(out_unit,*) ! on the last line store the diagonal elements
    write(out_unit,'(20F18.10)') gz_dens(1:Ns)    
    close(out_unit)
    !
    out_unit=free_unit()
    open(out_unit,file='local_density_density.data')
    do is=1,Ns
       write(out_unit,'(20F18.10)') gz_dens_dens(is,:)
    end do
    close(out_unit)
    !
    out_unit=free_unit()
    open(out_unit,file='local_angular_momenta.data')
    write(out_unit,'(20F18.10)') gz_spin2,gz_spinZ,gz_isospin2,gz_isospinZ
    close(out_unit)
    !
    if(present(vdm_simplex)) then
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
    end if
    !
    if(present(vdm_opt)) then
       out_unit=free_unit()
       open(out_unit,file='vdm_seed.restart')
       do istate=1,Nvdm_NC_opt-Nvdm_NCoff_opt
          write(out_unit,'(10F18.10)')  vdm_opt(istate)
       end do
       close(out_unit)
    end if

  end subroutine print_output







  subroutine Rhop_vec2mat(Rhop_indep,Rhop_mat)
    complex(8),dimension(:)   :: Rhop_indep
    complex(8),dimension(:,:) :: Rhop_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(Rhop_mat,1).ne.size(Rhop_mat,2)) stop "wrong stride"
    if(size(Rhop_mat,1).ne.Ns) stop "wrong stride"
    if(size(Rhop_indep).ne.NRhop_opt) stop "wrong stride!"    
    Rhop_mat = zero
    do iorb=1,Norb
       do ispin=1,2
          is=index(ispin,iorb)
          Rhop_mat(is,is) = Rhop_indep(iorb)
       end do
    end do
    !
  end subroutine Rhop_vec2mat
  subroutine Rhop_mat2vec(Rhop_mat,Rhop_indep)
    complex(8),dimension(:,:) :: Rhop_mat
    complex(8),dimension(:)   :: Rhop_indep
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    complex(8) :: test_stride
    real(8) :: test
    if(size(Rhop_mat,1).ne.size(Rhop_mat,2)) stop "wrong stride"
    if(size(Rhop_mat,1).ne.Ns) stop "wrong stride"
    if(size(Rhop_indep).ne.NRhop_opt) stop "wrong stride!"    
    !
    ispin=1;iorb=1;is=index(ispin,iorb)
    Rhop_indep(1)=Rhop_mat(is,is)
    ispin=1;iorb=2;is=index(ispin,iorb)
    Rhop_indep(2)=Rhop_mat(is,is)
    !
  end subroutine Rhop_mat2vec



  subroutine vdm_NC_vec2mat(vdm_NC_indep,vdm_NC_mat)
    complex(8),dimension(:)   :: vdm_NC_indep
    complex(8),dimension(:,:) :: vdm_NC_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(vdm_NC_mat,1).ne.size(vdm_NC_mat,2)) stop "wrong stride"
    if(size(vdm_NC_mat,1).ne.Ns) stop "wrong stride"
    if(size(vdm_NC_indep).ne.Nvdm_NC_opt) stop "wrong stride!"    
    !
    vdm_NC_mat = zero
    do iorb=1,Norb
       do ispin=1,2
          is=index(ispin,iorb)
          vdm_NC_mat(is,is) = vdm_NC_indep(iorb)
       end do
    end do
    !
  end subroutine vdm_NC_vec2mat
  subroutine vdm_NC_mat2vec(vdm_NC_mat,vdm_NC_indep)
    complex(8),dimension(:)   :: vdm_NC_indep
    complex(8),dimension(:,:) :: vdm_NC_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(vdm_NC_mat,1).ne.size(vdm_NC_mat,2)) stop "wrong stride"
    if(size(vdm_NC_mat,1).ne.Ns) stop "wrong stride"
    if(size(vdm_NC_indep).ne.Nvdm_NC_opt) stop "wrong stride!"    
    !
    ispin=1;iorb=1;is=index(ispin,iorb)
    vdm_NC_indep(1)=vdm_NC_mat(is,is)
    ispin=1;iorb=2;is=index(ispin,iorb)
    vdm_NC_indep(2)=vdm_NC_mat(is,is)
    !
  end subroutine vdm_NC_mat2vec



  subroutine vdm_NCoff_vec2mat(vdm_NC_indep,vdm_NC_mat)
    complex(8),dimension(:)   :: vdm_NC_indep
    complex(8),dimension(:,:) :: vdm_NC_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(vdm_NC_mat,1).ne.size(vdm_NC_mat,2)) stop "wrong stride"
    if(size(vdm_NC_mat,1).ne.Ns) stop "wrong stride"
    if(size(vdm_NC_indep).ne.Nvdm_NCoff_opt) stop "wrong stride!"    
    !
    vdm_NC_mat = zero
    !
  end subroutine vdm_NCoff_vec2mat
  subroutine vdm_NCoff_mat2vec(vdm_NC_mat,vdm_NC_indep)
    complex(8),dimension(:)   :: vdm_NC_indep
    complex(8),dimension(:,:) :: vdm_NC_mat
    integer                   :: i,j,is,js,iorb,jorb,ispin,jspin
    if(size(vdm_NC_mat,1).ne.size(vdm_NC_mat,2)) stop "wrong stride"
    if(size(vdm_NC_mat,1).ne.Ns) stop "wrong stride"
    if(size(vdm_NC_indep).ne.Nvdm_NCoff_opt) stop "wrong stride!"    
    !
    vdm_NC_indep = zero
    !
  end subroutine vdm_NCoff_mat2vec




end program GUTZ_mb



!AMOEBA TEST


