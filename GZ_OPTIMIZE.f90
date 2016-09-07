MODULE GZ_OPTIMIZED_ENERGY
  USE SCIFOR
  !
  USE GZ_AUX_FUNX
  USE GZ_VARS_GLOBAL
  USE GZ_EFFECTIVE_HOPPINGS
  USE GZ_LOCAL_FOCK  
  !
  USE LANCELOT_simple_double
  !
  USE GZ_MATRIX_BASIS
  USE GZ_ENERGY_MINIMIZATION
  USE MIN_AMOEBA
  implicit none
  private
  !

  !+- OPTIMIZATION ROUTINES w/ respect the VariationalDensityMatrix (VDM) -+!
  public :: gz_optimization_vdm_simplex   ! Simplex method (AMOEBA.f90)  !
  public :: gz_optimization_vdm_nlsq      !+- constrained NonLinearLeastSquare method (GALHAD)

  !+- OPTIMIZATION considering VDM and Renormalization_matrices as free parameters -+!
  public :: gz_optimization_vdm_Rhop
  public :: gz_optimization_vdm_Rhop_reduced
  public :: gz_optimization_vdm_Rhop_superc
  public :: gz_optimization_vdm_Rhop_superc_reduced

  public ::  get_gz_optimized_vdm_Rhop_superc
  !
  public :: dump2mats_superc,dump2vec_superc
  public :: dump2mats,dump2vec
  integer :: root_unit,xmin_unit
  !
CONTAINS
  !+-----------------------------------------------------------------------------------------------------+!
  !+- PURPOSE: Minimize the GUTZWILLER ENERGY FUNCTIONAL WITH RESPECT TO THE VARIATIONAL DENSITY MATRIX -+!
  !+-----------------------------------------------------------------------------------------------------+!
  include 'gz_optimization.f90'
  !
  subroutine gz_optimization_vdm_Rhop(init_Rhop,init_lgr_slater,init_lgr_proj,init_vdm) 
    complex(8),dimension(Ns,Ns),intent(inout) :: init_Rhop
    complex(8),dimension(Ns,Ns),intent(inout) :: init_lgr_slater,init_lgr_proj
    !
    real(8),dimension(Ns,Ns),optional :: init_vdm
    !complex(8),dimension(Ns,Ns) :: tmp_vdm
    real(8),dimension(Ns) :: vdm
    complex(8),dimension(Ns,Ns) :: Rhop
    real(8) :: delta,tmp_ene,out_err
    integer :: iter    !
    !
    integer                           :: iunit,err_unit,ene_unit,Nsuccess  
    real(8),dimension(:),allocatable  :: xmin,xout
    integer :: Nopt,iopt,Nslater_lgr,NRhop,NQhop,Nproj_lgr
    integer :: is,js,imap
    !
    NRhop = 2*NRhop_opt
    !
    Nslater_lgr = 2*Nvdm_NC_opt 
    Nproj_lgr = 2*Nvdm_NCoff_opt 
    !
    Nopt = NRhop + Nslater_lgr + Nproj_lgr
    allocate(xmin(Nopt),xout(Nopt))    

    if(present(init_vdm)) then
       do is=1,Ns
          vdm(is) = init_vdm(is,is)
       end do
       call slater_minimization_lgr(init_Rhop,vdm,tmp_ene,init_lgr_slater,iverbose=.true.)    
    end if

    if(GZmin_verbose) then
       write(*,*) 'Rn0_minimization: INPUT PARAMETERS'
       write(*,*) 'Rhop'
       do is=1,Ns
          write(*,'(10F8.4)') init_Rhop(is,:)
       end do
       write(*,*) 'init_lgr_slater'
       do is=1,Ns
          write(*,'(10F8.4)') init_lgr_slater(is,:)
       end do
       write(*,*) 'init_lgr_proj'
       do is=1,Ns
          write(*,'(10F8.4)') init_lgr_proj(is,:)
       end do
    end if
    !
    !
    call dump2vec(xmin,init_Rhop,init_lgr_slater,init_lgr_proj) !
    !
    !  
    if(GZmin_verbose) then
       write(*,*) 'DUMPED OPTIMIZATION VECTOR'
       do is=1,Nopt
          write(*,*) is,xmin(is)
       end do
    end if
    !
    !
    if(GZmin_verbose) then
       root_unit=free_unit()
       open(unit=root_unit,file='Rn0_root_finding.out')
       xmin_unit=free_unit()
       open(unit=xmin_unit,file='Rn0_root_finding_xmin.out')
    end if
    !
    call fsolve(R_VDM_free_zeros_,xmin,tol=1.d-10,info=iter)
    !    
    !+- once optimization is achieved store the ground state results -+!
    optimization_flag=.true.
    if(allocated(GZ_vector)) deallocate(GZ_vector)
    allocate(GZ_vector(Nphi))
    if(allocated(GZ_opt_slater_lgr)) deallocate(GZ_opt_slater_lgr)
    allocate(GZ_opt_slater_lgr(Ns,Ns))
    if(allocated(GZ_opt_proj_lgr)) deallocate(GZ_opt_proj_lgr)
    allocate(GZ_opt_proj_lgr(Ns,Ns))

    !
    xout = R_VDM_free_zeros_(xmin)
    out_err=0.d0
    do is=1,Nopt
       out_err = out_err + xout(is)**2.d0
    end do
    !
    if(GZmin_verbose) then       
       write(root_unit,*)
       write(root_unit,*)
       write(root_unit,*) out_err
       do is=1,Nopt
          write(root_unit,*) is,xout(is),xmin(is)
       end do
       close(root_unit)
       close(xmin_unit)
    end if
    call dump2mats(xmin,init_Rhop,init_lgr_slater,init_lgr_proj)
    !
  end subroutine gz_optimization_vdm_Rhop
  !



  subroutine gz_optimization_vdm_Rhop_reduced(init_Rhop,init_lgr_slater,init_lgr_proj,init_vdm) 
    complex(8),dimension(Ns,Ns),intent(inout) :: init_Rhop
    complex(8),dimension(Ns,Ns),intent(inout) :: init_lgr_slater,init_lgr_proj
    !
    real(8),dimension(Ns,Ns),optional :: init_vdm
    !complex(8),dimension(Ns,Ns) :: tmp_vdm
    real(8),dimension(Ns) :: vdm
    complex(8),dimension(Ns,Ns) :: Rhop
    real(8) :: delta,tmp_ene,out_err
    integer :: iter    !
    !
    integer                           :: iunit,err_unit,ene_unit,Nsuccess  
    real(8),dimension(:),allocatable  :: xmin,xout
    real(8),dimension(:),allocatable  :: xmin_,xout_
    integer :: Nopt,iopt,Nslater_lgr,NRhop,NQhop,Nproj_lgr
    integer :: is,js,imap
    !
    NRhop = 2*NRhop_opt
    !
    Nslater_lgr = 2*Nvdm_NC_opt 
    Nproj_lgr = 2*Nvdm_NCoff_opt 
    !
    Nopt = NRhop + Nslater_lgr + Nproj_lgr
    allocate(xmin(Nopt),xout(Nopt))    

    if(present(init_vdm)) then
       do is=1,Ns
          vdm(is) = init_vdm(is,is)
       end do
       call slater_minimization_lgr(init_Rhop,vdm,tmp_ene,init_lgr_slater,iverbose=.true.)    
    end if

    if(GZmin_verbose) then
       write(*,*) 'Rn0_minimization: INPUT PARAMETERS'
       write(*,*) 'Rhop'
       do is=1,Ns
          write(*,'(10F8.4)') init_Rhop(is,:)
       end do
       write(*,*) 'init_lgr_slater'
       do is=1,Ns
          write(*,'(10F8.4)') init_lgr_slater(is,:)
       end do
       write(*,*) 'init_lgr_proj'
       do is=1,Ns
          write(*,'(10F8.4)') init_lgr_proj(is,:)
       end do
    end if
    !
    !
    call dump2vec(xmin,init_Rhop,init_lgr_slater,init_lgr_proj) !
    !
    !  
    if(GZmin_verbose) then
       write(*,*) 'DUMPED OPTIMIZATION VECTOR'
       do is=1,Nopt
          write(*,*) is,xmin(is)
       end do
    end if
    !
    !
    if(GZmin_verbose) then
       root_unit=free_unit()
       open(unit=root_unit,file='Rn0_root_finding.out')
       xmin_unit=free_unit()
       open(unit=xmin_unit,file='Rn0_root_finding_xmin.out')
    end if
    !
    
    !
    allocate(xmin_(Nopt_reduced))
    call stride_zeros_orig2red(xmin,xmin_)
    call fsolve(R_VDM_free_zeros__,xmin_,tol=1.d-10,info=iter)
    call stride_zeros_red2orig(xmin_,xmin)
    !
    !call fsolve(R_VDM_free_zeros_,xmin,tol=1.d-10,info=iter)
    !    
    !+- once optimization is achieved store the ground state results -+!
    optimization_flag=.true.
    if(allocated(GZ_vector)) deallocate(GZ_vector)
    allocate(GZ_vector(Nphi))
    if(allocated(GZ_opt_slater_lgr)) deallocate(GZ_opt_slater_lgr)
    allocate(GZ_opt_slater_lgr(Ns,Ns))
    if(allocated(GZ_opt_proj_lgr)) deallocate(GZ_opt_proj_lgr)
    allocate(GZ_opt_proj_lgr(Ns,Ns))
    !
    xout = R_VDM_free_zeros_(xmin)
    out_err=0.d0
    do is=1,Nopt
       out_err = out_err + xout(is)**2.d0
    end do
    !
    if(GZmin_verbose) then       
       write(root_unit,*)
       write(root_unit,*)
       write(root_unit,*) out_err
       do is=1,Nopt
          write(root_unit,*) is,xout(is),xmin(is)
       end do
       close(root_unit)
       close(xmin_unit)
    end if
    call dump2mats(xmin,init_Rhop,init_lgr_slater,init_lgr_proj)
    !
  end subroutine gz_optimization_vdm_Rhop_reduced







  !
  !
  subroutine gz_optimization_vdm_Rhop_superc(init_Rhop,init_Qhop,init_lgr_slater,init_lgr_proj,init_vdm) 
    complex(8),dimension(Ns,Ns),intent(inout) :: init_Rhop
    complex(8),dimension(Ns,Ns),intent(inout) :: init_Qhop
    complex(8),dimension(2,Ns,Ns),intent(inout) :: init_lgr_slater,init_lgr_proj
    !
    real(8),dimension(Ns,Ns),optional :: init_vdm
    complex(8),dimension(2,Ns,Ns) :: tmp_vdm
    real(8),dimension(Ns) :: vdm
    complex(8),dimension(Ns,Ns) :: Rhop
    real(8) :: delta,tmp_ene
    integer :: iter    !
    !
    integer                           :: n_min,neq,nin,maxit,print_level,exit_code
    real(8)                           :: gradtol,feastol,Ephi,nsite,err_iter,Ephi_,out_err
    real(8),allocatable               :: bL(:),bU(:),cx(:),y(:),phi_optimize(:),phi_optimize_(:)
    integer                           :: iunit,err_unit,ene_unit,Nsuccess  
    real(8),dimension(:),allocatable  :: xmin,xout
    integer :: Nopt,iopt,Nslater_lgr,NRhop,NQhop,Nproj_lgr
    integer :: is,js,imap
    !
    NRhop = 2*NRhop_opt
    NQhop = 2*NQhop_opt
    !
    Nslater_lgr = 2*Nvdm_NC_opt + 2*Nvdm_AC_opt
    Nproj_lgr = 2*Nvdm_NCoff_opt + 2*Nvdm_AC_opt
    !
    !Nopt = 2*NRhop_opt + 2*NQhop_opt + 2*Nvdm_NC_opt + 2*Nvdm_AC_opt + 2*Nvdm_NCoff_opt + 2*Nvdm_AC_opt
    Nopt = NRhop + NQhop + Nslater_lgr + Nproj_lgr
    allocate(xmin(Nopt),xout(Nopt))    

    if(present(init_vdm)) then
       tmp_vdm = zero
       tmp_vdm(1,:,:) = init_vdm
       do is=1,Ns
          vdm(is) = init_vdm(is,is)
       end do
       call slater_minimization_lgr_superc(init_Rhop,init_Qhop,vdm,tmp_ene,init_lgr_slater,iverbose=.true.)    
    end if

    if(GZmin_verbose) then
       write(*,*) 'RQn0_minimization: INPUT PARAMETERS'
       write(*,*) 'Rhop'
       do is=1,Ns
          write(*,'(10F8.4)') init_Rhop(is,:)
       end do
       write(*,*) 'Qhop'
       do is=1,Ns
          write(*,'(10F8.4)') init_Qhop(is,:)
       end do
       write(*,*) 'init_lgr_slater Normal'
       do is=1,Ns
          write(*,'(10F8.4)') init_lgr_slater(1,is,:)
       end do
       write(*,*) 'init_lgr_slater Anomalous'
       do is=1,Ns
          write(*,'(10F8.4)') init_lgr_slater(2,is,:)
       end do
       write(*,*) 'init_lgr_proj Normal'
       do is=1,Ns
          write(*,'(10F8.4)') init_lgr_proj(1,is,:)
       end do
       write(*,*) 'init_lgr_proj Anomalous'
       do is=1,Ns
          write(*,'(10F8.4)') init_lgr_proj(2,is,:)
       end do
    end if
    !
    call dump2vec_superc(xmin,init_Rhop,init_Qhop,init_lgr_slater,init_lgr_proj)
    !  
    if(GZmin_verbose) then
       write(*,*) 'DUMPED OPTIMIZATION VECTOR'
       do is=1,Nopt
          write(*,*) is,xmin(is)
       end do
    end if
    !
    !
    if(GZmin_verbose) then
       root_unit=free_unit()
       open(unit=root_unit,file='RQn0_root_finding.out')
       xmin_unit=free_unit()
       open(unit=xmin_unit,file='RQn0_root_finding_xmin.out')
    end if
    !
    call fsolve(R_Q_VDM_free_zeros_superc,xmin,tol=1.d-10,info=iter)
    !call fmin_cg(xmin,R_Q_VDM_free_opt_superc,iter,out_err,ftol=10d-15,istop=2)
    ! write(*,*) 'ITER',iter
    !    
    !+- once optimization is achieved store the ground state results -+!
    optimization_flag=.true.
    if(allocated(GZ_vector)) deallocate(GZ_vector)
    allocate(GZ_vector(Nphi))
    if(allocated(GZ_opt_slater_lgr_superc)) deallocate(GZ_opt_slater_lgr_superc)
    allocate(GZ_opt_slater_lgr_superc(2,Ns,Ns))
    if(allocated(GZ_opt_proj_lgr_superc)) deallocate(GZ_opt_proj_lgr_superc)
    allocate(GZ_opt_proj_lgr_superc(2,Ns,Ns))
    !
    xout = R_Q_VDM_free_zeros_superc(xmin)
    out_err=0.d0
    do is=1,Nopt
       out_err = out_err + xout(is)**2.d0
    end do
    !
    if(GZmin_verbose) then       
       write(root_unit,*)
       write(root_unit,*)
       write(root_unit,*) out_err
       do is=1,Nopt
          write(root_unit,*) is,xout(is),xmin(is)
       end do
       close(root_unit)
       close(xmin_unit)
    end if
    call dump2mats_superc(xmin,init_Rhop,init_Qhop,init_lgr_slater,init_lgr_proj)
    !
  end subroutine gz_optimization_vdm_Rhop_superc
  !
  !
  !
  !
  subroutine gz_optimization_vdm_Rhop_superc_reduced(init_Rhop,init_Qhop,init_lgr_slater,init_lgr_proj,init_vdm) 
    complex(8),dimension(Ns,Ns),intent(inout) :: init_Rhop
    complex(8),dimension(Ns,Ns),intent(inout) :: init_Qhop
    complex(8),dimension(2,Ns,Ns),intent(inout) :: init_lgr_slater,init_lgr_proj
    !
    real(8),dimension(Ns,Ns),optional :: init_vdm
    complex(8),dimension(2,Ns,Ns) :: tmp_vdm
    real(8),dimension(Ns) :: vdm
    complex(8),dimension(Ns,Ns) :: Rhop
    real(8) :: delta,tmp_ene
    integer :: iter    !
    !
    integer                           :: n_min,neq,nin,maxit,print_level,exit_code
    real(8)                           :: gradtol,feastol,Ephi,nsite,err_iter,Ephi_,out_err
    real(8),allocatable               :: bL(:),bU(:),cx(:),y(:),phi_optimize(:),phi_optimize_(:)
    integer                           :: iunit,err_unit,ene_unit,Nsuccess  
    real(8),dimension(:),allocatable  :: xmin,xout
    real(8),dimension(:),allocatable  :: xmin_,xout_
    integer :: Nopt,iopt,Nslater_lgr,NRhop,NQhop,Nproj_lgr
    integer :: is,js,imap
    !
    NRhop = 2*NRhop_opt
    NQhop = 2*NQhop_opt
    !
    Nslater_lgr = 2*Nvdm_NC_opt + 2*Nvdm_AC_opt
    Nproj_lgr = 2*Nvdm_NCoff_opt + 2*Nvdm_AC_opt
    !
    Nopt = NRhop + NQhop + Nslater_lgr + Nproj_lgr
    allocate(xmin(Nopt),xout(Nopt))    

    if(present(init_vdm)) then
       tmp_vdm = zero
       tmp_vdm(1,:,:) = init_vdm
       do is=1,Ns
          vdm(is) = init_vdm(is,is)
       end do
       call slater_minimization_lgr_superc(init_Rhop,init_Qhop,vdm,tmp_ene,init_lgr_slater,iverbose=.true.)    
    end if

    if(GZmin_verbose) then
       write(*,*) 'RQn0_minimization: INPUT PARAMETERS'
       write(*,*) 'Rhop'
       do is=1,Ns
          write(*,'(10F8.4)') init_Rhop(is,:)
       end do
       write(*,*) 'Qhop'
       do is=1,Ns
          write(*,'(10F8.4)') init_Qhop(is,:)
       end do
       write(*,*) 'init_lgr_slater Normal'
       do is=1,Ns
          write(*,'(10F8.4)') init_lgr_slater(1,is,:)
       end do
       write(*,*) 'init_lgr_slater Anomalous'
       do is=1,Ns
          write(*,'(10F8.4)') init_lgr_slater(2,is,:)
       end do
       write(*,*) 'init_lgr_proj Normal'
       do is=1,Ns
          write(*,'(10F8.4)') init_lgr_proj(1,is,:)
       end do
       write(*,*) 'init_lgr_proj Anomalous'
       do is=1,Ns
          write(*,'(10F8.4)') init_lgr_proj(2,is,:)
       end do
    end if
    !
    call dump2vec_superc(xmin,init_Rhop,init_Qhop,init_lgr_slater,init_lgr_proj)
    !  
    if(GZmin_verbose) then
       write(*,*) 'DUMPED OPTIMIZATION VECTOR @ REDUCED OPTIMIZATION'
       do is=1,Nopt
          write(*,*) is,xmin(is)
       end do
    end if
    !
    !
    if(GZmin_verbose) then
       root_unit=free_unit()
       open(unit=root_unit,file='RQn0_root_finding.out')
       xmin_unit=free_unit()
       open(unit=xmin_unit,file='RQn0_root_finding_xmin.out')
    end if
    !    
    allocate(xmin_(Nopt_reduced))
    call stride_zeros_orig2red(xmin,xmin_)
    call fsolve(R_Q_VDM_free_zeros_superc_,xmin_,tol=1.d-10,info=iter)
    call stride_zeros_red2orig(xmin_,xmin)
    !
    !+- once optimization is achieved store the ground state results -+!
    optimization_flag=.true.
    if(allocated(GZ_vector)) deallocate(GZ_vector)
    allocate(GZ_vector(Nphi))
    if(allocated(GZ_opt_slater_lgr_superc)) deallocate(GZ_opt_slater_lgr_superc)
    allocate(GZ_opt_slater_lgr_superc(2,Ns,Ns))
    if(allocated(GZ_opt_proj_lgr_superc)) deallocate(GZ_opt_proj_lgr_superc)
    allocate(GZ_opt_proj_lgr_superc(2,Ns,Ns))
    !
    xout = R_Q_VDM_free_zeros_superc(xmin)
    out_err=0.d0
    do is=1,Nopt
       out_err = out_err + xout(is)**2.d0
    end do
    !
    if(GZmin_verbose) then       
       write(root_unit,*)
       write(root_unit,*)
       write(root_unit,*) out_err
       do is=1,Nopt
          write(root_unit,*) is,xout(is),xmin(is)
       end do
       close(root_unit)
       close(xmin_unit)
    end if
    call dump2mats_superc(xmin,init_Rhop,init_Qhop,init_lgr_slater,init_lgr_proj)
    !
  end subroutine gz_optimization_vdm_Rhop_superc_reduced






  !+- reconstruct solution
  subroutine get_gz_optimized_vdm_Rhop_superc(init_Rhop,init_Qhop,init_lgr_slater,init_lgr_proj,init_vdm) 
    complex(8),dimension(Ns,Ns),intent(inout) :: init_Rhop
    complex(8),dimension(Ns,Ns),intent(inout) :: init_Qhop
    complex(8),dimension(2,Ns,Ns),intent(inout) :: init_lgr_slater,init_lgr_proj
    !
    real(8),dimension(Ns,Ns),optional :: init_vdm
    complex(8),dimension(2,Ns,Ns) :: tmp_vdm
    real(8),dimension(Ns) :: vdm
    complex(8),dimension(Ns,Ns) :: Rhop
    real(8) :: delta,tmp_ene
    integer :: iter    !
    !
    integer                           :: n_min,neq,nin,maxit,print_level,exit_code
    real(8)                           :: gradtol,feastol,Ephi,nsite,err_iter,Ephi_,out_err
    real(8),allocatable               :: bL(:),bU(:),cx(:),y(:),phi_optimize(:),phi_optimize_(:)
    integer                           :: iunit,err_unit,ene_unit,Nsuccess  
    real(8),dimension(:),allocatable  :: xmin,xout,xmin_,xout_
    integer :: Nopt,iopt,Nslater_lgr,NRhop,NQhop,Nproj_lgr
    integer :: is,js,imap
    !
    NRhop = 2*NRhop_opt
    NQhop = 2*NQhop_opt
    !
    Nslater_lgr = 2*Nvdm_NC_opt + 2*Nvdm_AC_opt
    Nproj_lgr = 2*Nvdm_NCoff_opt + 2*Nvdm_AC_opt
    !
    !Nopt = 2*NRhop_opt + 2*NQhop_opt + 2*Nvdm_NC_opt + 2*Nvdm_AC_opt + 2*Nvdm_NCoff_opt + 2*Nvdm_AC_opt
    Nopt = NRhop + NQhop + Nslater_lgr + Nproj_lgr
    allocate(xmin(Nopt),xout(Nopt))    


    if(GZmin_verbose) then
       write(*,*) 'RQn0 OPTIMIZED PARAMETERS'
       write(*,*) 'Rhop'
       do is=1,Ns
          write(*,'(10F8.4)') init_Rhop(is,:)
       end do
       write(*,*) 'Qhop'
       do is=1,Ns
          write(*,'(10F8.4)') init_Qhop(is,:)
       end do
       write(*,*) 'init_lgr_slater Normal'
       do is=1,Ns
          write(*,'(10F8.4)') init_lgr_slater(1,is,:)
       end do
       write(*,*) 'init_lgr_slater Anomalous'
       do is=1,Ns
          write(*,'(10F8.4)') init_lgr_slater(2,is,:)
       end do
       write(*,*) 'init_lgr_proj Normal'
       do is=1,Ns
          write(*,'(10F8.4)') init_lgr_proj(1,is,:)
       end do
       write(*,*) 'init_lgr_proj Anomalous'
       do is=1,Ns
          write(*,'(10F8.4)') init_lgr_proj(2,is,:)
       end do
    end if
    !
    call dump2vec_superc(xmin,init_Rhop,init_Qhop,init_lgr_slater,init_lgr_proj)
    !  
    if(GZmin_verbose) then
       write(*,*) 'DUMPED OPTIMIZED VECTOR'
       do is=1,Nopt
          write(*,*) is,xmin(is)
       end do
    end if
    !
    !
    ! allocate(xmin_(Nopt_reduced),xout_(Nopt_reduced))
    ! call stride_zeros_orig2red(xmin,xmin_)
    ! xout_ = R_Q_VDM_free_zeros_superc_(xmin_)
    ! !    call fsolve(R_Q_VDM_free_zeros_superc_,xmin_,tol=1.d-10,info=iter)
    ! call stride_zeros_red2orig(xout_,xout)

    !+- once optimization is achieved store the ground state results -+!
    optimization_flag=.true.
    if(allocated(GZ_vector)) deallocate(GZ_vector)
    allocate(GZ_vector(Nphi))
    if(allocated(GZ_opt_slater_lgr_superc)) deallocate(GZ_opt_slater_lgr_superc)
    allocate(GZ_opt_slater_lgr_superc(2,Ns,Ns))
    if(allocated(GZ_opt_proj_lgr_superc)) deallocate(GZ_opt_proj_lgr_superc)
    allocate(GZ_opt_proj_lgr_superc(2,Ns,Ns))
    !
    xout = R_Q_VDM_free_zeros_superc(xmin)
    out_err=0.d0
    do is=1,Nopt
       out_err = out_err + xout(is)**2.d0
    end do
    !
    if(GZmin_verbose) then       
       write(*,*)
       write(*,*)
       write(*,*) out_err
       do is=1,Nopt
          write(*,*) is,xout(is),xmin(is)
       end do
    end if
    call dump2mats_superc(xmin,init_Rhop,init_Qhop,init_lgr_slater,init_lgr_proj)
    !
  end subroutine get_gz_optimized_vdm_Rhop_superc
 


  



  !
  !+---------------------------------------------------------------------------+!
  !+- CONSTRAINED MINIMIZATION WITH RESPECT TO THE VARIATIONAL DENSITY MATRIX -+!
  !+---------------------------------------------------------------------------+!
  !
  subroutine gz_optimization_vdm_nlsq(init_vdm,optimized_vdm)
    real(8),dimension(:),intent(in)  :: init_vdm
    real(8),dimension(size(init_vdm)),intent(out) :: optimized_vdm
    integer                             :: iter,unit_vdm_opt,icall
    real(8)                             :: GZ_energy
    !
    !LANCELOT VARIABLES
    !
    integer                           :: n_min,neq,nin,maxit,print_level,exit_code
    real(8)                           :: gradtol,feastol,Ephi,nsite,err_iter,Ephi_
    real(8),allocatable               :: bL(:),bU(:),cx(:),y(:),phi_optimize(:),phi_optimize_(:)
    integer                           :: iunit,err_unit,ene_unit,Nsuccess    
    !
    if(size(init_vdm).ne.Nvdm_NC_opt-Nvdm_NCoff_opt) stop "wrong dimensions @ gz_optimizatoin_vdm_nlsq"
    optimized_vdm=init_vdm
    !
    unit_vdm_opt=free_unit()
    open(unit_vdm_opt,file='vdm_optimization.out'); icall=0
    !
    opt_energy_unit=free_unit()
    open(opt_energy_unit,file='GZ_OptEnergy_VS_vdm.out')
    opt_rhop_unit=free_unit()
    open(opt_rhop_unit,file='GZ_OptRhop_VS_vdm.out')
    if(gz_superc) then
       opt_qhop_unit=free_unit()
       open(opt_qhop_unit,file='GZ_OptQhop_VS_vdm.out')
    end if
    opt_GZ_unit=free_unit()
    open(opt_GZ_unit,file='GZ_OptProj_VS_vdm.out')
    if(GZmin_verbose) then
       GZmin_unit=free_unit()
       open(GZmin_unit,file='GZ_SelfCons_min_verbose.out')
       GZmin_unit_=free_unit()
       open(GZmin_unit_,file='GZ_proj_min.out')
    end if



    ! LANCELOT configuration parameters 
    n_min       = Nvdm_NC_opt-Nvdm_NCoff_opt  ! number of minimization parameters  
    neq         = 0   ! number of equality constraints                   
    nin         = 0           ! number of in-equality constraints                   
    maxit       = 1000        ! maximum iteration number 
    gradtol     = 1.d-7       ! maximum norm of the gradient at convergence 
    feastol     = 1.d-7       ! maximum violation of parameters at convergence  
    print_level = lancelot_verbose           ! verbosity
    allocate(bl(n_min),bu(n_min),cx(neq+nin),y(neq+nin))
    bL = 1.d-4               ! lower bounds for minimization parameters
    bU = 1.d0-bL                 ! upper bounds for minimization parameters          
    !    
    call lancelot_simple(n_min,optimized_vdm,GZ_energy,exit_code,my_fun=gz_get_energy_vdm, &
         bl = bl, bu = bu,                                                                      &
         neq = neq, nin = nin,                                                                  &
         cx = cx, y = y, iters  = iter, maxit = maxit,                                          &
         gradtol = gradtol, feastol = feastol,                                                  &
         print_level = print_level )
    !+--------------------------------------------------------------------------------------+!    

    optimization_flag=.true.
    if(.not.allocated(GZ_vector)) allocate(GZ_vector(Nphi))
    call gz_get_energy_vdm(optimized_vdm,GZ_energy)

  end subroutine gz_optimization_vdm_nlsq








  !+---------------------------------------------------------------------+!
  !+- SIMPLEX MINIMIZATION W/ RESPECT TO THE VARIATIONAL DENSITY MATRIX -+!
  !+---------------------------------------------------------------------+!
  subroutine gz_optimization_vdm_simplex(simplex_init,optimized_vdm) 
    real(8),dimension(Ns+1,Ns),intent(inout) :: simplex_init
    real(8),dimension(Ns),intent(out)               :: optimized_vdm
    !real(8),intent(out)                                    :: optimized_energy    
    !+- amoeba_variables-+!
    real(8),allocatable,dimension(:,:)                     :: p
    real(8),allocatable,dimension(:)                       :: y
    real(8)                                                :: ftol
    integer                                                :: np,mp,i_vertex,j_vertex,i_dim,iter
    !integer,allocatable,dimension(:)                       :: idum
    !integer                                                :: tmp_dum
    !real(8)                                                :: rnd,tot_dens
    real(8)                                                :: GZ_energy
    !
    !
    integer                                                :: amoeba_unit    
    !+-------------------+!
    optimization_flag=.false.
    amoeba_unit=free_unit()
    open(amoeba_unit,file='Amoeba_minimization.out')
    opt_energy_unit=free_unit()
    open(opt_energy_unit,file='GZ_OptEnergy_VS_vdm.out')
    opt_rhop_unit=free_unit()
    open(opt_rhop_unit,file='GZ_OptRhop_VS_vdm.out')
    if(gz_superc) then
       opt_qhop_unit=free_unit()
       open(opt_qhop_unit,file='GZ_OptQhop_VS_vdm.out')
    end if
    opt_GZ_unit=free_unit()
    open(opt_GZ_unit,file='GZ_OptProj_VS_vdm.out')
    if(GZmin_verbose) then
       GZmin_unit=free_unit()
       open(GZmin_unit,file='GZ_SelfCons_min_verbose.out')
       GZmin_unit_=free_unit()
       open(GZmin_unit_,file='GZ_proj_min.out')
    end if
    NP=Ns
    MP=NP+1  
    allocate(y(MP),p(MP,NP))
    !+- initialize simplex -+!
    p=simplex_init
    write(amoeba_unit,*) 'Initializing simplex verteces'
    !+----------------------+!  
    do i_vertex=1,MP
       y(i_vertex) = gz_energy_vdm(p(i_vertex,:))
       write(amoeba_unit,*) p(i_vertex,:),y(i_vertex)
    end do
    !<DEBUG
    !stop
    !DEBUG>
    !
    ftol=amoeba_min_tol
    call amoeba(p(1:MP,1:NP),y(1:MP),ftol,gz_energy_vdm,iter,amoeba_verbose)
    !+- Optimized Variational Density Matrix -+!
    optimized_vdm=p(1,:) 
    !+- Optimized Variational Energy -+!
    simplex_init=p
    write(amoeba_unit,*) 
    write(amoeba_unit,*) 'Last loop simplex verteces'
    write(amoeba_unit,*) 
    do i_vertex=1,MP
       write(amoeba_unit,*) p(i_vertex,:),y(i_vertex)
    end do
    write(amoeba_unit,*) 
    write(amoeba_unit,*) 
    deallocate(y,p)
    close(amoeba_unit)
    close(opt_energy_unit)    
    !
    optimization_flag=.true.
    if(.not.allocated(GZ_vector)) allocate(GZ_vector(Nphi))
    GZ_energy = gz_energy_vdm(optimized_vdm)
    !
  end subroutine gz_optimization_vdm_simplex












END MODULE GZ_OPTIMIZED_ENERGY
