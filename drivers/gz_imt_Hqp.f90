program GUTZ_mb
  USE SCIFOR
  !
  USE DMFT_MISC
  USE SF_PARSE_INPUT
  USE RK_IDE
  !
  USE GZ_AUX_FUNX
  USE GZ_neqAUX_FUNX
!  USE GZ_DYNAMICS
  USE GZ_imtDYNAMICS
  USE GZ_VARS_GLOBAL
  USE GZ_LOCAL_FOCK
  USE GZ_LOCAL_HAMILTONIAN
  USE GZ_ENERGY_MINIMIZATION
  USE GZ_OPTIMIZED_ENERGY
  USE GZ_MATRIX_BASIS
  !
  implicit none
  real(8),dimension(:),allocatable :: epsik,hybik,auxt
  real(8) :: t,r,s,tmpU,ndens
  integer :: Nx,out_unit,is,js,ik,im_it,itt,i,iorb,ispin,it
  integer :: nprint
  !
  character(len=200) :: store_dir,read_dir,read_optWF_dir,read_finSC_dir
  complex(8),dimension(:,:,:),allocatable :: slater_init,HK_save,Hqp_out,Hqp_in
  complex(8),dimension(:),allocatable     :: gz_proj_init
  !
  !
  complex(8),dimension(:),allocatable     :: psi_t,psi_save,psi_t_,psi_tmp
  complex(8),dimension(:,:),allocatable     :: lgr_NC
  real(8) :: lgrU
  !
  integer :: unit_imt_hloc
  integer :: unit_imt_local_dens
  integer :: unit_imt_local_dens_dens
  integer :: unit_imt_ene
  integer :: unit_imt_dens_constrSL
  integer :: unit_imt_dens_constrGZ
  integer :: unit_imt_dens_constrSLa
  integer :: unit_imt_dens_constrGZa
  integer :: unit_imt_constrU
  integer :: unit_imt_Rhop
  integer :: unit_imt_Qhop
  integer :: unit_imt_AngMom
  integer :: unit_imt_sc_order
  integer :: unit_imt_nqp
  !
  !+- observables -+!
  complex(8),dimension(:),allocatable   :: Rhop
  complex(8),dimension(:,:),allocatable   :: Rhop_matrix
  complex(8),dimension(:,:),allocatable :: local_density_matrix
  real(8),dimension(:,:),allocatable    :: local_dens_dens
  complex(8),dimension(:,:),allocatable :: dens_constrSL
  complex(8),dimension(:,:),allocatable :: dens_constrGZ
  real(8)                               :: unitary_constr
  real(8),dimension(4)                  :: local_angular_momenta
  real(8),dimension(3)                  :: energies
  !
  real(8),dimension(:),allocatable      :: dump_vect
  !
  real(8) :: itstart,itstop
  real(8) :: imt_dene,ene_save,imt_s,imt_f
  real(8),dimension(:),allocatable :: imt_entropy,imt_ene
  !
  
  !
  call parse_input_variable(Wband,"WBAND","inputGZ.conf",default=2.d0)  
  call parse_input_variable(Nx,"Nx","inputGZ.conf",default=1000)
  call parse_input_variable(read_optWF_dir,"EQWF_DIR","inputGZ.conf",default='./')
  call parse_input_variable(ndens,"NDENS","inputGZ.conf",default=1.d0)
  call parse_input_variable(nprint,"NPRINT","inputGZ.conf",default=10)

  !
  call read_input("inputGZ.conf")
  call save_input_file("inputGZ.conf")
  !
  Norb=1
  wf_symmetry=0
  call initialize_local_fock_space  
  call init_variational_matrices(wf_symmetry)
  !
  call build_lattice_model
  Nvdm_NC_opt=1; vdm_NC_stride_v2m => vdm_NC_vec2mat ; vdm_NC_stride_m2v => vdm_NC_mat2vec
  !
  allocate(eLevels(Ns)); eLevels=0.d0
  !
  !+- INITIALIZE TIME GRIDS -+!
  Nit_aux=2*Nit+1
  allocate(t_grid(Nt),t_grid_aux(Nt_aux))
  !
  itstart = beta_init
  itstop = beta_init + itstep*real(Nit-1)
  t_grid = linspace(itstart,itstop,Nit)
  !
  itstart = beta_init
  itstop = beta_init + itstep*0.5d0*real(Nit_aux-1)
  t_grid_aux = linspace(itstart,itstop,Nit_aux)
  
  !
  !
  call setup_imt_hamiltonian
  !
  !
  
  !
  !+- READ EQUILIBRIUM AND SETUP DYNAMICAL VECTOR -+!
  nDynamics = Nphi + Ns*Ns*Lk
  allocate(psi_t(nDynamics),psi_save(nDynamics),psi_t_(nDynamics),psi_tmp(nDynamics))
  
  allocate(slater_init(Ns,Ns,Lk),gz_proj_init(Nphi),Hqp_in(Ns,Ns,Lk))  
  !
  allocate(lgr_NC(Ns,Ns));lgr_NC=zero
  lgrU=0.d0
  !
  if(beta_init.eq.0.d0) then
     call read_optimized_variational_wf_normal_imt(read_optWF_dir,gz_proj_init)
     Hqp_in = zero
     ! gz_proj_init(1)=1./sqrt(2.d0)
     ! gz_proj_init(2)=1./2.d0
     ! gz_proj_init(3)=1./2.d0
     call beta0_init_imt_qpH(gz_proj_init,Hqp_in)
  else
     call read_optimized_variational_wf_normal_imt(read_optWF_dir,gz_proj_init,Hqp_in)
     !call beta_init_imt_qpH(gz_proj_init,Hqp_in)
  end if
  call wfMatrix_2_dynamicalVector(Hqp_in,gz_proj_init,psi_t)  
  !
  !psi_t = gz_proj_init  
  !
  !
  ! it=1
  ! Uloc=Uloc_t(:,it)
  ! Ust =Ust_t(it)
  ! Jh=Jh_t(it)
  ! Jsf=Jsf_t(it)
  ! Jph=Jph_t(it)
  eLevels = 0.d0
  ! !
  call get_local_hamiltonian_trace(eLevels)      
  !
  call setup_imt_dynamics
  !    
  unit_imt_Rhop = free_unit()
  open(unit_imt_Rhop,file='imt_Rhop_matrix.data')
  !
  unit_imt_local_dens = free_unit()
  open(unit_imt_local_dens,file='imt_local_density_matrix.data')
  !
  unit_imt_local_dens_dens = free_unit()
  open(unit_imt_local_dens_dens,file='imt_local_dens_dens.data')
  !
  unit_imt_ene = free_unit()
  open(unit_imt_ene,file='imt_energy.data')
  !
  unit_imt_dens_constrSL = free_unit()
  open(unit_imt_dens_constrSL,file='imt_dens_constrSL.data')
  !
  unit_imt_dens_constrGZ = free_unit()
  open(unit_imt_dens_constrGZ,file='imt_dens_constrGZ.data')
  !
  unit_imt_constrU = free_unit()
  open(unit_imt_constrU,file='imt_constrU.data')
  !
  unit_imt_AngMom = free_unit()
  open(unit_imt_AngMom,file='imt_AngMom.data')
  !
  unit_imt_sc_order = free_unit()
  open(unit_imt_sc_order,file='imt_sc_order.data')
  !
  allocate(Rhop_matrix(Ns,Ns))
  allocate(local_density_matrix(Ns,Ns))
  allocate(local_dens_dens(Ns,Ns))
  allocate(dens_constrSL(Ns,Ns))
  allocate(dens_constrGZ(Ns,Ns))  
  allocate(dump_vect(Ns*Ns))

  allocate(imt_ene(Nit),imt_entropy(Nit))

  !*) ACTUAL DYNAMICS (simple do loop measuring each nprint times)
  imt_s=0.d0
  do im_it=1,Nit
     !
     t=t_grid(im_it)
     !
     if(mod(im_it-1,nprint).eq.0) then        
        !
        call gz_imt_measure(psi_t,t)
        !
        do is=1,Ns
           do js=1,Ns
              call get_imt_Rhop(is,js,Rhop_matrix(is,js))              
              call get_imt_local_dens(is,js,local_density_matrix(is,js))              
              call get_imt_local_dens_dens(is,js,local_dens_dens(is,js))              
              call get_imt_dens_constr_slater(is,js,dens_constrSL(is,js))
              call get_imt_dens_constr_gzproj(is,js,dens_constrGZ(is,js))              
           end do
        end do
        !
        call get_imt_energies(energies)
        imt_ene(im_it) = energies(1)

        call get_imt_local_angular_momenta(local_angular_momenta)
        call get_imt_unitary_constr(unitary_constr)
        !
        call write_complex_matrix_grid(Rhop_matrix,unit_imt_Rhop,print_grid_Rhop,t)
        !
        call write_hermitean_matrix(local_density_matrix,unit_imt_local_dens,t)
        call write_hermitean_matrix(dens_constrSL,unit_imt_dens_constrSL,t)
        call write_hermitean_matrix(dens_constrGZ,unit_imt_dens_constrGZ,t)
        call write_hermitean_matrix(local_density_matrix,unit_imt_local_dens,t)
        call write_symmetric_matrix(local_dens_dens,unit_imt_local_dens_dens,t)
        write(unit_imt_AngMom,'(10F18.10)') t,local_angular_momenta
        if(im_it.lt.1) write(unit_imt_ene,'(10F18.10)') t,energies!,imt_dene,imt_dene*t_grid(im_it-1)
        write(unit_imt_constrU,'(10F18.10)') t,unitary_constr
        !        
     end if
     !
     psi_save=psi_t
     if(im_it.lt.Nit) then
        write(*,*) im_it,Nit
        !call  step_imt_dynamics(nDynamics,itstep,t,psi_t,lgr_NC,lgrU,gz_imt_eom)        
        call  step_imt_dynamics(nDynamics,itstep,t,psi_t,lgr_NC,lgrU,gz_imt_equations_of_motion,ndens,fix_constr_=.true.)        
        !psi_t = RK_step(nDynamics,4,itstep,t,psi_t,gz_imt_equations_of_motion)     !
        !
        !psi_t = trpz_implicit(nDynamics,4,itstep,t,psi_t,gz_imt_equations_of_motion_)     !
        !psi_t = mp_step(nDynamics,4,itstep,t,psi_t,gz_imt_equations_of_motion_)     !
        !psi_tmp = psi_t
        !psi_t = mp_symm_step(nDynamics,4,itstep,t,psi_t,psi_t_,gz_imt_equations_of_motion_)     !
        !psi_t_ = psi_tmp
        !
     end if
     !

  end do

  close(unit_imt_ene)
  open(unit_imt_ene,file='imt_free_energy.data')
  ene_save=imt_ene(Nit)
  imt_s=0.d0
  do im_it=1,Nit-1     
     !
     imt_dene = (imt_ene(Nit+1-im_it) - imt_ene(Nit-im_it))/itstep
     imt_s = imt_s - t_grid(Nit+1-im_it)*imt_dene*itstep
     !
     imt_f = imt_ene(im_it) - 1./t_grid(im_it)*imt_s
     write(unit_imt_ene,'(10F18.10)') t_grid(Nit+1-im_it),imt_f,imt_s
  end do


  allocate(Hqp_out(Ns,Ns,Lk))
  call gz_imt_measure(psi_t,t,slater_init,gz_proj_init,Hqp_out)

  
  
  call print_output(slater_init,gz_proj_init,Hqp_out)
  
  

  
  ! stop
  
  ! call system('cp * FORWARD_IMT_DYN')

  ! close(unit_imt_local_dens_dens)
  ! unit_imt_local_dens_dens = free_unit()
  ! open(unit_imt_local_dens_dens,file='imt_local_dens_dens.data')


  ! allocate(auxt(Nit))
  ! auxt=t_grid
  ! do it=1,Nit
  !    t_grid(it) = auxt(Nit-(it-1))+itstep     
  ! end do
  ! deallocate(auxt)
  ! allocate(auxt(Nit_aux))
  ! auxt=t_grid_aux
  ! do it=1,Nit
  !    t_grid_aux(it) = auxt(Nit_aux-(it-1))+itstep*0.5d0     
  ! end do
  ! deallocate(auxt)

  ! itstep=-1.d0*itstep
  ! !call init_imt_qpH(,psi_t)  
  ! !gz_imt_qpH=Hk_save  

  ! do im_it=1,Nit
  !    write(*,*) im_it,Nit
  !    !
  !    t=t_grid(im_it)
  !    !
  !    if(mod(im_it-1,nprint).eq.0) then        
  !       !
  !       call gz_imt_measure(psi_t,t)
  !       !
  !       do is=1,Ns
  !          do js=1,Ns
  !             call get_imt_Rhop(is,js,Rhop_matrix(is,js))              
  !             call get_imt_local_dens(is,js,local_density_matrix(is,js))              
  !             call get_imt_local_dens_dens(is,js,local_dens_dens(is,js))              
  !             call get_imt_dens_constr_slater(is,js,dens_constrSL(is,js))
  !             call get_imt_dens_constr_gzproj(is,js,dens_constrGZ(is,js))              
  !          end do
  !       end do
  !       !
  !       call get_imt_energies(energies)
  !       call get_imt_local_angular_momenta(local_angular_momenta)
  !       call get_imt_unitary_constr(unitary_constr)
  !       !
  !       call write_complex_matrix_grid(Rhop_matrix,unit_imt_Rhop,print_grid_Rhop,t)
  !       !
  !       call write_hermitean_matrix(local_density_matrix,unit_imt_local_dens,t)
  !       call write_hermitean_matrix(dens_constrSL,unit_imt_dens_constrSL,t)
  !       call write_hermitean_matrix(dens_constrGZ,unit_imt_dens_constrGZ,t)
  !       call write_hermitean_matrix(local_density_matrix,unit_imt_local_dens,t)
  !       call write_symmetric_matrix(local_dens_dens,unit_imt_local_dens_dens,t)
  !       write(unit_imt_AngMom,'(10F18.10)') t,local_angular_momenta
  !       write(unit_imt_ene,'(10F18.10)') t,energies
  !       write(unit_imt_constrU,'(10F18.10)') t,unitary_constr
  !       !        
  !    end if
  !    !

     
  !    call  step_imt_dynamics(nDynamics,itstep,t,psi_t,lgr_NC,lgrU,gz_imt_eom) 
  !    !
  !    !psi_t = RK_step(nDynamics,4,tstep,t,psi_t,gz_equations_of_motion)
  !    !
  ! end do

  
  

  ! write(800,*) psi_t

  !
  !
  !
  !
  !
  !
  !
CONTAINS
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
       test_k=test_k+wtk(ix)
       write(77,*) epsik(ix),wtk(ix)
    end do
    hybik=0.d0
    write(*,*) test_k,de
    !stop
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



  ! out_unit=free_unit()
  ! open(out_unit,file='optimized_projectors.data')
  ! test_full_phi=0.d0
    ! do iphi=1,Nphi
    !    !
    !    test_full_phi = test_full_phi + GZ_vector(iphi)*phi_basis(iphi,:,:)
    !    write(out_unit,'(2F18.10)') GZ_vector(iphi)
    ! end do
    ! close(out_unit)    
    ! open(out_unit,file='optimized_phi_matrix.data')    
    ! do ifock=1,nFock
    !    do jfock=1,nFock
    !       write(out_unit,'(2F18.10,2I2)') test_full_phi(ifock,jfock),ifock,jfock
    !    end do
    ! end do
    ! close(out_unit)    
    ! !
    ! do is=1,Ns
    !    do js=1,Ns
    !       out_unit=free_unit()
    !       open(out_unit,file='optimized_slater_normal_IS'//reg(txtfy(is))//'_JS'//reg(txtfy(js))//'.data')
    !       do ik=1,Lk
    !          write(out_unit,'(2F18.10)') dreal(GZ_opt_slater(is,js,ik)),dimag(GZ_opt_slater(is,js,ik))
    !       end do
    !       close(out_unit)
    !    end do
    ! end do






  subroutine print_output(slater,gzproj,Hqp_out)
    complex(8),dimension(Ns,Ns,Lk) :: slater
    complex(8),dimension(Ns,Ns,Lk) :: Hqp_out
    complex(8),dimension(Nphi) :: gzproj
    integer :: out_unit,istate,jstate,iorb,iphi,ifock,jfock,is,js,ik
    integer,dimension(Ns) :: fock_state
    complex(8),dimension(Ns) :: tmp
    real(8) :: deltani,delta_tmp,vdm_tmp

    complex(8),dimension(nFock,nFock) :: test_full_phi

    !+- STORE GZ PROJECTORS -+!
    out_unit=free_unit()
    open(out_unit,file='optimized_projectors.data')
    test_full_phi=0.d0
    do iphi=1,Nphi
       !
       test_full_phi = test_full_phi + gzproj(iphi)*phi_basis(iphi,:,:)
       write(out_unit,'(2F18.10)') gzproj(iphi)
    end do
    close(out_unit)

    open(out_unit,file='optimized_phi_matrix.data')
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
             write(out_unit,'(2F18.10)') dreal(slater(is,js,ik)),dimag(slater(is,js,ik))
          end do
          close(out_unit)
       end do
    end do
    !+- Hqp
    do is=1,Ns
       do js=1,Ns
          out_unit=free_unit()
          open(out_unit,file='HQP_IS'//reg(txtfy(is))//'_JS'//reg(txtfy(js))//'.data')
          do ik=1,Lk
             write(out_unit,'(2F18.10)') dreal(Hqp_out(is,js,ik)),dimag(Hqp_out(is,js,ik))
          end do
          close(out_unit)
       end do
    end do
    
  end subroutine print_output

  !   !
  !   out_unit=free_unit()
  !   open(out_unit,file='optimized_internal_energy.data')
  !   write(out_unit,'(5F18.10)') GZ_opt_energy,GZ_opt_kinetic,GZ_opt_Eloc
  !   close(out_unit)
  !   !
  !   out_unit=free_unit()
  !   open(out_unit,file='optimized_variational_density_matrix.data')
  !   write(out_unit,*) 'NORMAL VDM'
  !   do istate=1,Ns
  !      tmp(istate)=GZ_opt_VDM(istate,istate)
  !      write(out_unit,'(20F18.10)') GZ_opt_VDM(istate,:)
  !   end do
  !   write(out_unit,*) ! on the last line store the diagonal elements
  !   write(out_unit,'(20F18.10)') tmp(1:Ns)
  !   close(out_unit)
  !   !
  !   out_unit=free_unit()
  !   open(out_unit,file='optimized_Rhop_matrix.data')
  !   do istate=1,Ns
  !      tmp(istate)=GZ_opt_Rhop(istate,istate)
  !      write(out_unit,'(20F18.10)') GZ_opt_Rhop(istate,1:Ns)
  !   end do
  !   write(out_unit,*) ! on the last line store the diagonal elements
  !   write(out_unit,'(20F18.10)') tmp(1:Ns)
  !   close(out_unit)
  !   !
  !   out_unit=free_unit()
  !   open(out_unit,file='optimized_density.data')
  !   do istate=1,Ns
  !      write(out_unit,'(20F18.10)') gz_dens_matrix(istate,1:Ns)
  !   end do
  !   write(out_unit,*) ! on the last line store the diagonal elements
  !   write(out_unit,'(20F18.10)') gz_dens(1:Ns)    
  !   close(out_unit)
  !   !
  !   out_unit=free_unit()
  !   open(out_unit,file='local_density_density.data')
  !   do is=1,Ns
  !      write(out_unit,'(20F18.10)') gz_dens_dens(is,:)
  !   end do
  !   close(out_unit)
  !   !
  !   out_unit=free_unit()
  !   open(out_unit,file='local_angular_momenta.data')
  !   write(out_unit,'(20F18.10)') gz_spin2,gz_spinZ,gz_isospin2,gz_isospinZ
  !   close(out_unit)
  !   !
  !   if(present(vdm_simplex)) then
  !      out_unit=free_unit()
  !      open(out_unit,file='vdm_simplex.restart')
  !      do jstate=1,Ns+1
  !         if(jstate.le.Ns) then
  !            do istate=1,Ns
  !               write(out_unit,'(20F18.10)') vdm_simplex(jstate,istate)
  !            end do
  !            if(jstate.le.Ns) write(out_unit,*)  'x'
  !         else
  !            do istate=1,Ns
  !               deltani=vdm_simplex(jstate,istate)-0.5
  !               if(deltani.gt.0.d0) then
  !                  delta_tmp=0.9999-vdm_simplex(jstate,istate)
  !                  vdm_tmp=vdm_simplex(jstate,istate)+delta_tmp*0.1
  !                  write(out_unit,'(20F18.10)') vdm_tmp
  !               else
  !                  delta_tmp=vdm_simplex(jstate,istate)-0.0001
  !                  vdm_tmp=vdm_simplex(jstate,istate)-delta_tmp*0.1
  !                  write(out_unit,'(20F18.10)') vdm_tmp
  !               end if
  !            end do
  !         end if
  !      end do
  !      close(out_unit)
  !   end if
  !   !
  !   if(present(vdm_opt)) then
  !      out_unit=free_unit()
  !      open(out_unit,file='vdm_seed.restart')
  !      do istate=1,Nvdm_NC_opt-Nvdm_NCoff_opt
  !         write(out_unit,'(10F18.10)')  vdm_opt(istate)
  !      end do
  !      close(out_unit)
  !   end if

  ! end subroutine print_output







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
          Rhop_mat(is,is) = Rhop_indep(1)
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
          vdm_NC_mat(is,is) = vdm_NC_indep(1)
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


