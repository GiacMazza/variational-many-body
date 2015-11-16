MODULE GZ_ENERGY_FUNCTIONAL_SELF
  USE SF_LINALG
  USE SF_IOTOOLS
  USE SF_SPECIAL
  USE SF_OPTIMIZE
  USE GZ_VARS_GLOBAL
  USE GZ_AUX_FUNX
  USE LANCELOT_simple_double
  USE GZ_PROJECTORS
  implicit none
  private
  !
  public :: gz_energy_self
  !
  real(8),dimension(4) :: store_lgr_multip
contains
  !
  function gz_energy_self(n0)   result(GZ_energy)
    real(8),dimension(:),intent(in)           :: n0 !INPUT: Variational Density Matrix (VDM) (diagonal in istate)    
    real(8)                                           :: GZ_energy !INPUT: Optimized GZ energy at fixed 

    real(8) :: GZ_energy_old,energy_err     ! Value of the GZ energy functional

    real(8),dimension(state_dim,state_dim)            :: R_init        ! initial guess for the hopping renormalization matrix    
    real(8),dimension(state_dim,state_dim)            :: slater_derivatives    

    real(8),dimension(state_dim,state_dim)            :: R_iter,tmpR ! hopping matrix renormalization (during iterations)
    real(8),dimension(state_dim)            :: slater_lgr_multip
    real(8),dimension(state_dim,state_dim)            :: GZproj_lgr_multip  ! GZ vector (during iterations)

    real(8)                                           :: E_Hstar,E_Hloc
    real(8),dimension(nFock)                          :: GZvect_iter  ! GZ vector (during iterations)

    integer :: istate,iter,jstate,ifock,jfock
    integer :: unit
    logical :: bound
    !
    write(*,*) '*************************'
    write(*,*) 'INPUT DENSITY',n0(:)
    bound=.false.
    do istate=1,state_dim
       if(n0(istate).le.1.d-10.or.n0(istate).ge.1.d0-1.d-10) bound=.true.
    end do
    !
    ! GET phi_vect compatible with n0 ...mmm I don't think it's possible...
    
    !

    if(.not.bound) then
       R_init=0.d0 
       do istate=1,state_dim
          R_init(istate,istate)=Rseed
       end do
       !
       GZ_energy=0.d0    
       R_iter=R_init
       store_lgr_multip=0.d0
       GZvect_iter=0.d0
       do iter=1,Niter_self
          ! update hopping matrices
          GZ_energy_old=GZ_energy
          tmpR=R_iter

          GZvect_iter(11)=1.d0
          R_iter=Rhop_matrix(GZvect_iter,n0)
          !
          write(*,*) R_iter
          
          !+----------------------------+!
          !+- SLATER STEP MINIMIZATION -+!
          !+----------------------------+!    
          call slater_determinant_minimization_step(R_iter,n0,E_Hstar,slater_lgr_multip,slater_derivatives,GZmin_verbose)       

          !+----------------------------+!
          !+- GZproj STEP MINIMIZATION -+!
          !+----------------------------+!    
          call gz_projectors_minimization_step(slater_derivatives,n0,E_Hloc,GZvect_iter,GZproj_lgr_multip,GZmin_verbose)

          R_iter=0.1*Rhop_matrix(GZvect_iter,n0)+0.9*tmpR
          !R_iter=Rhop_matrix(GZvect_iter,n0)
          
          GZ_energy=E_Hstar+E_Hloc
          if(iter.lt.2) then 
             energy_err=1.d0
          else
             energy_err=abs(GZ_energy-GZ_energy_old)
          end if
          if(GZmin_verbose) then
             write(GZmin_unit,*) dble(iter),E_Hloc+E_Hstar,energy_err,E_Hloc,E_Hstar,R_iter(1,1),R_iter(2,2)
             write(100,*) dble(iter),E_Hloc+E_Hstar,energy_err,E_Hloc,E_Hstar,R_iter(1,1),R_iter(2,2)
          end if
          if(energy_err.lt.err_self) exit
       end do
       if(GZmin_verbose) write(GZmin_unit,*) 
       if(iter-1.eq.Niter_self) then
          write(*,*) 'Self consistent Gutzwiller minimization'
          write(*,*) 'Input VDM',n0
          write(*,*) 'final error',energy_err
          write(*,*) "Not converged after",Niter_self,'iterations: exiting'
          stop 
       end if
    else
       GZ_energy=100.d0
    end if
    write(opt_energy_unit,*) n0,GZ_energy,E_Hstar,E_Hloc
    !
  end function gz_energy_self
  !  
  subroutine gz_projectors_minimization_step(slater_derivatives,n0_target,E_Hloc,GZvect,lgr_multip,iverbose)
    real(8),dimension(state_dim,state_dim),intent(in) :: slater_derivatives !input:  Slater Deter GZ energy derivatives
    real(8),dimension(state_dim),intent(in)           :: n0_target          !input:  Variational density matrix
    ! real(8),dimension(state_dim),intent(out):: lgr_multip                 
    ! real(8),dimension(state_dim) :: lgr
    real(8),dimension(state_dim,state_dim),intent(out):: lgr_multip         !output: GZprojectors Lagrange Multipliers -diagonal-
    real(8),dimension(nFock)                          :: GZvect   !output: GZvector
    real(8)                                           :: E_Hloc   !output: optimized local energy
    real(8)                                           :: lgr_symm(1)
    logical,optional                                  :: iverbose
    
    real(8),dimension(state_dim) :: lgr
    real(8),dimension(state_dim)            :: err_dens
    real(8) :: delta_out
    logical                                 :: iverbose_
    ! real(8),dimension(state_dim,state_dim)            :: dens_test
    ! real(8),dimension(state_dim*state_dim)            :: dens_err    
    !    real(8),dimension(state_dim)            :: dens_err_diag
    !    real(8),dimension(state_dim) :: lgr_diag
    integer                                           :: info,istate,i,j
    
    iverbose_=.false.;if(present(iverbose)) iverbose_=iverbose    
    lgr=store_lgr_multip


    lgr_symm=-10
    do i=1,20
       lgr_symm(1)=lgr_symm(1)+1
       delta_out=get_delta_proj_variational_density_symm(lgr_symm)
       write(90,*) lgr_symm(1),delta_out
    end do
    !stop


    lgr_symm(1)=0.1d0
    !call fsolve(get_delta_proj_variational_density_symm,lgr,tol=1.d-15,info=info)    
    !call fmin_cgminimize(lgr_symm,get_delta_proj_variational_density_symm,info,delta_out)

    !call fzero_brentq(lgr_symm,get_delta_proj_variational_density_symm)

    ! !
    ! write(*,*)
    ! write(*,*) "DEBUGGING TEST"
    ! write(*,*)
    !
    ! lgr=store_lgr_multip

    !store_lgr_multip=lgr
    !write(*,*) "DELTA_OUT",delta_out
    !stop
    lgr_multip=0.d0
    do istate=1,state_dim
       !lgr_multip(istate,istate)=lgr(istate)
       lgr_multip(1,1)= lgr_symm(1)
       lgr_multip(2,2)=-lgr_symm(1)
       lgr_multip(3,3)= lgr_symm(1)
       lgr_multip(4,4)=-lgr_symm(1)
    end do
    call get_GZproj_ground_state(n0_target,slater_derivatives,lgr_multip,E_Hloc,GZvect)
    
    if(iverbose_) then
       write(*,*)
       write(*,*) "GZ projectors: Lagrange Parameters -diagonal case-",info
       !write(*,'(10F18.10)') lgr(1:state_dim)
       write(*,'(10F18.10)') lgr_symm(1)
       delta_out=get_delta_proj_variational_density_symm(lgr_symm)
       write(*,*) "GZ projectors: Variational density matrix error"
       write(*,'(10F18.10)') delta_out
       write(*,*) "GZ projectors: Optimized Local Energy"
       write(*,'(10F18.10)') E_Hloc
       write(*,*)
    end if
    !stop    
    !<DEBUG
    ! write(*,*)
    ! write(*,*) "GA PROJECTORS LAGRANGE PARAMETERS",info
    ! write(*,*)
    ! do istate=1,state_dim
    !    write(*,'(10F18.10)') lgr_multip(istate,1:state_dim)
    ! end do
    ! write(*,*)
    ! dens_test=get_proj_variational_density_full(lgr)
    ! write(*,*)
    ! write(*,*) "GA PROJECTORS VARIATIONAL DENSITY"
    ! write(*,*)
    ! do istate=1,state_dim
    !    write(*,'(10F18.10)') dens_test(istate,1:state_dim)
    ! end do
    ! write(*,*)
    ! dens_err=get_delta_proj_variational_density_full(lgr)
    ! write(*,*) "GA PROJECTORS VARIATIONAL DENSITY ERROR"
    ! write(*,*)
    ! call vec2mat_stride(dens_err,dens_test)
    ! do istate=1,state_dim
    !    write(*,'(10F18.10)') dens_test(istate,1:state_dim)
    ! end do
    ! write(*,*)
    ! !DEBUG>
    ! !<DEBUG
    ! lgr_diag=0.0d0
    ! call fsolve(get_delta_proj_variational_density_diag,lgr_diag,tol=1.d-15,info=info)
    ! write(*,*)
    ! write(*,*) "GA PROJECTORS LAGRANGE PARAMETERS DIAG",info
    ! write(*,*)
    ! write(*,'(10F18.10)') lgr_diag(1:state_dim)
    ! write(*,*)
    ! dens_test_diag=get_proj_variational_density_diag(lgr_diag)
    ! write(*,*)
    ! write(*,*) "GA PROJECTORS VARIATIONAL DENSITY DIAG"
    ! write(*,*)
    ! write(*,'(10F18.10)') dens_test_diag(1:state_dim)
    ! write(*,*)
    ! dens_test_diag=get_delta_proj_variational_density_diag(lgr_diag)
    ! write(*,*) "GA PROJECTORS VARIATIONAL DENSITY ERROR DIAG"
    ! write(*,*)
    ! ! call vec2mat_stride(dens_err,dens_test)
    ! ! do istate=1,state_dim
    ! write(*,'(10F18.10)') dens_test_diag(1:state_dim)
    ! !end do
    ! write(*,*)

    ! write(*,'(A,20F6.2)') "GA PROJECTORS VECTOR",GZvect
    ! write(*,*) "GA PROJECTORS GROUND STATE ENERGY",E_Hloc
  contains
    !
    include 'self_minimization_GZproj_routines.f90'
    !
  end subroutine gz_projectors_minimization_step


  subroutine get_GZproj_ground_state(n0,slater_derivatives,lgr_multip,E_Hloc,GZvect) 
    real(8),dimension(state_dim)           :: n0
    real(8),dimension(state_dim,state_dim) :: slater_derivatives,lgr_multip
    real(8)                                :: E_Hloc
    real(8),dimension(nFock)               :: GZvect
    !
    real(8),dimension(nFock,nFock)         :: H_projectors
    real(8),dimension(nFock)               :: H_eigens
    integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
    !
    !+- build up the local H_projectors -+!
    H_projectors=phi_traces_basis_Hloc
    do istate=1,state_dim
       do jstate=1,state_dim
          H_projectors = H_projectors + &
               slater_derivatives(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)/sqrt(n0(jstate)*(1.d0-n0(jstate)))
          H_projectors = H_projectors + lgr_multip(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
       end do
    end do
    !
    call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
    GZvect=H_projectors(1:nFock,1)
    !
    E_Hloc=0.d0
    do ifock=1,nFock
       do jfock=1,nFock
          E_Hloc=E_Hloc+GZvect(ifock)*phi_traces_basis_Hloc(ifock,jfock)*GZvect(jfock)
       end do
    end do
  end subroutine get_GZproj_ground_state
  
  
  !+-------------------------------+!
  !+- SLATER DETERMINANT ROUTINES -+!
  !+-------------------------------+!
  subroutine slater_determinant_minimization_step(Rhop,n0_target,Estar,lgr_multip,slater_derivatives,iverbose) 
    real(8),dimension(state_dim,state_dim),intent(in)  :: Rhop         !input:  renrmalization matrix
    real(8),dimension(state_dim),intent(in)            :: n0_target    !input:  variational density matrix
    real(8),intent(out)                                :: Estar        !output: Slater Deter GS energy
    real(8),dimension(state_dim),intent(out)           :: lgr_multip   !output: Slater Deter lagrange multipliers
    real(8),dimension(state_dim,state_dim),intent(out) :: slater_derivatives !output: Slater Deter GS energy derivatives
    logical,optional                                   :: iverbose     !input:  Verbosity level
    !
    real(8),dimension(state_dim)                       :: lgr
    real(8),dimension(state_dim,state_dim)             :: Hk,test_dens
    real(8),dimension(state_dim)                       :: ek,tmp_lgr,err_dens
    integer                                            :: iorb,jorb,ik,ispin,jspin,istate,jstate,info,korb
    real(8),dimension(state_dim)                       :: lgr_multip_vec
    logical                                            :: iverbose_

    iverbose_=.false.;if(present(iverbose)) iverbose_=iverbose

    lgr=0.d0
    call fsolve(get_delta_local_density_matrix_diag,lgr,tol=1.d-12,info=info)
    call store_slater_ground_state(Rhop,lgr,Estar,slater_derivatives)
    lgr_multip=lgr

    if(iverbose_) then
       write(*,*)
       write(*,*) "Slater Determinant: Lagrange Multipliers - diagonal form -"
       write(*,'(20F18.10)') lgr_multip
       write(*,*) "Slater Determinant: Variational density error"
       err_dens=get_delta_local_density_matrix_diag(lgr_multip)
       write(*,"(20F18.10)") err_dens
       write(*,*) "Slater Determinant: Ground State energy"
       write(*,"(20F18.10)") Estar
       write(*,*)
    end if
    ! Get slater ground state informations : Estar & slater_derivatives
    !HERE IS THE MISTERY!!!!! NOW I'M NOT ABLE TO FIGURE-OUT!!!!
    ! DAI CHE e' QUASI FINITA!!
    !<DEBUG
    ! write(*,*)
    ! write(*,*) "SLATER DETERMINANT LAGRANGE MULTIPLIERS"
    ! write(*,*)
    ! write(*,'(20F18.10)') lgr_multip
    ! write(*,*)
    ! write(*,*) "SLATER DETERMINANT VARIATIONAL DENSITY"
    ! write(*,*)
    !test_dens=get_local_density_matrix_diag(lgr_multip)
    ! do istate=1,state_dim
    !    write(*,"(20F18.10)") test_dens(istate,1:state_dim)
    ! end do
    ! write(*,*)
    ! write(*,*) "SLATER DETERMINANT VARIATIONAL DENSITY ERROR"
    ! write(*,*)
    !    do istate=1,state_dim
    !       write(*,"(20F18.10)") slater_derivatives(istate,:)
    !    end do
    !    write(*,*)
    !err_dens=get_delta_local_density_matrix_diag(lgr_multip)
    ! write(*,"(20F18.10)") err_dens
    ! write(*,*)
    ! write(*,*) 'SLATER DETERMINANT GROUND STATE ENERGY',Estar
    ! write(*,*) 
    !DEBUG>        
    !<DEBUG de le DEBUG
    do istate=1,state_dim
       write(*,"(20F18.10)") slater_derivatives(istate,:)
    end do
    write(*,*)
    ! do istate=1,state_dim
    !    write(*,"(20F18.10)") slater_derivatives(istate,:)
    ! end do
    !DEBUG de le DEBUG>
    !stop
    !
  contains    
    !
    include 'self_minimization_slater_routines.f90'

  end subroutine slater_determinant_minimization_step


  subroutine store_slater_ground_state(Rhop,lm,Estar,slater_derivatives)     
    real(8),dimension(state_dim,state_dim),intent(in) :: Rhop
    real(8),dimension(state_dim),intent(in)           :: lm
    real(8)                                           :: Estar
    real(8),dimension(state_dim,state_dim)            :: slater_derivatives
    real(8),dimension(state_dim,state_dim)            :: Hk,tmp,Hk_bare,Hstar
    real(8),dimension(state_dim)                      :: ek
    integer                                           :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik
    Estar=0.d0
    slater_derivatives=0.d0
    do ik=1,Lk
       Hk=0.d0
       ek=0.d0
       do iorb=1,Norb
          do ispin=1,2
             do jorb=1,Norb
                do jspin=1,2
                   istate=index(ispin,iorb)
                   jstate=index(jspin,jorb)               
                   ! build up the hopping hamiltonian !
                   if(ispin.eq.jspin) then
                      if(iorb.eq.jorb) then
                         Hk(istate,jstate)=epsik(ik)
                      else
                         Hk(istate,jstate)=hybik(ik)
                      end if
                   end if
                end do
             end do
          end do
       end do
       Hk_bare=Hk
       ! hopping renormalization !
       Hk=matmul(Hk,Rhop)
       Hk=matmul(Rhop,Hk)
       Hstar=Hk
       ! add Lagrange multipliers !
       do istate=1,state_dim
          Hk(istate,istate)=Hk(istate,istate)+lm(istate)             
       end do
       ! diagonalize hamiltonian !
       call  matrix_diagonalize(Hk,ek,'V','L')
       ! get ground state energy !       
       ! do kstate=1,state_dim
       !    Estar = Estar + fermi(ek(kstate),beta)*ek(kstate)*wtk(ik)
       ! end do
       ! store slater determinant matrix elements
       do iorb=1,Norb
          do ispin=1,2
             do jorb=1,Norb
                do jspin=1,2
                   istate=index(ispin,iorb)
                   jstate=index(jspin,jorb)               
                   tmp(istate,jstate)=0.d0
                   do kstate=1,state_dim
                      tmp(istate,jstate) = tmp(istate,jstate) + Hk(istate,kstate)*fermi(ek(kstate),beta)*Hk(jstate,kstate)
                      Estar = Estar + Hk(istate,kstate)*Hk(jstate,kstate)*Hstar(istate,jstate)*fermi(ek(kstate),beta)*wtk(ik)
                   end do
                end do
             end do
          end do
       end do
       ! store slater ground state derivatives
       tmp=matmul(Rhop,tmp)
       tmp=matmul(Hk_bare,tmp)             
       do istate=1,state_dim
          do jstate=1,state_dim             
             slater_derivatives(istate,jstate) = &
                  slater_derivatives(istate,jstate) + 2.d0*tmp(istate,jstate)*wtk(ik)
          end do
       end do
    end do
  end subroutine store_slater_ground_state


  !<TEST/DEBUG
  ! function delta_n0_vect(GZvect) result delta_n0
  !   real(8),dimension(:) :: GZvect
  !   real(8)
  ! end function delta_n0_vect
  !TEST/DEBUG>


END MODULE GZ_ENERGY_FUNCTIONAL_SELF

