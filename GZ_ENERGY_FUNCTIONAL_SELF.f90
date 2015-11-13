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

  public :: gz_energy_self!,fix_slater_matrix_elements
  !real(8),dimension(:),allocatable,public :: local_dens_target
  !real(8),dimension(:),allocatable :: phi_tmp


contains
  function gz_energy_self(n0)   result(GZ_energy)
    real(8),dimension(state_dim),intent(in)           :: n0            ! Variational density matrix (diagonal in istate index by definition)    
    real(8)                                           :: GZ_energy
    real(8) :: GZ_energy_old,energy_err     ! Value of the GZ energy functional

    real(8),dimension(state_dim,state_dim)            :: R_init        ! initial guess for the hopping renormalization matrix    
    real(8),dimension(state_dim,state_dim)            :: slater_derivatives    

    real(8),dimension(state_dim,state_dim)            :: R_iter,tmpR ! hopping matrix renormalization (during iterations)
    real(8),dimension(state_dim)            :: slater_lgr_multip
    real(8),dimension(state_dim,state_dim)            :: GZproj_lgr_multip  ! GZ vector (during iterations)

    real(8)                                           :: E_Hstar,E_Hloc
    real(8),dimension(nFock)                          :: GZvect_iter  ! GZ vector (during iterations)

    integer :: istate,iter,jstate,ifock,jfock
    !
    R_init=0.d0 ! it would be useful if this matrix is a global variable which can be updated as an input
    do istate=1,state_dim
       R_init(istate,istate)=1.d0
    end do
    
    GZ_energy=0.d0
    do iter=1,50
       GZ_energy_old=GZ_energy
       if(iter.eq.1) then
          R_iter=R_init
       else
          ! recompute new hopping parameters
          do istate=1,state_dim
             do jstate=1,state_dim
                R_iter(istate,jstate)=0.d0
                do ifock=1,nFock
                   do jfock=1,nFock
                      R_iter(istate,jstate)= &
                           R_iter(istate,jstate) + GZvect_iter(ifock)*GZvect_iter(jfock)*phi_traces_basis_Rhop(istate,jstate,ifock,jfock)
                   end do
                end do
                R_iter(istate,jstate)=R_iter(istate,jstate)/sqrt(n0(jstate)*(1.d0-n0(jstate)))
             end do
          end do
       end if
       !                              !
       !+----------------------------+!
       !+- SLATER STEP MINIMIZATION -+!
       !+----------------------------+!    
       !                              !
       call slater_determinant_minimization_step(R_iter,n0,E_Hstar,slater_lgr_multip,slater_derivatives)
       
       !<DEBUG
       ! do istate=1,state_dim
       !    write(*,*) slater_derivatives(istate,:)
       ! end do
       !DEBUG>

       !                              !
       !+----------------------------+!
       !+- GZproj STEP MINIMIZATION -+!
       !+----------------------------+!    
       !                              !
       call gz_projectors_minimization_step(slater_derivatives,n0,E_Hloc,GZvect_iter,GZproj_lgr_multip)
       !
       GZ_energy=E_Hstar+E_Hloc
       if(iter.lt.2) then 
          energy_err=1.d0
       else
          energy_err=abs(GZ_energy-GZ_energy_old)
       end if
       write(*,*) dble(iter),E_Hloc+E_Hstar,energy_err
       write(777,*) dble(iter),E_Hloc+E_Hstar,energy_err,E_Hloc,E_Hstar,R_iter(1,1)
       if(energy_err.lt.1.d-10) exit
    end do
    !
  end function gz_energy_self


  subroutine gz_projectors_minimization_step(slater_derivatives,n0_target,E_Hloc,GZvect,lgr_multip)
    real(8),dimension(state_dim,state_dim),intent(in) :: slater_derivatives
    real(8),dimension(state_dim),intent(in)           :: n0_target
    real(8),dimension(state_dim,state_dim) :: lgr_multip
    real(8),dimension(state_dim*state_dim) :: lgr_multip_vec
    real(8),dimension(nFock)               :: GZvect
    real(8)                                :: E_Hloc
    real(8),dimension(state_dim,state_dim) :: dens_test
    real(8),dimension(state_dim*state_dim) :: dens_err
    integer :: info,istate

    lgr_multip_vec=0.1d0
    call mat2vec_stride(lgr_multip,lgr_multip_vec)
    call fsolve(get_delta_proj_variational_density,lgr_multip_vec,tol=1.d-15,info=info)
    call vec2mat_stride(lgr_multip_vec,lgr_multip)

    !<DEBUG
    ! write(*,*)
    ! write(*,*) "GA PROJECTORS LAGRANGE PARAMETERS",info
    ! write(*,*)
    ! do istate=1,state_dim
    !    write(*,'(10F18.10)') lgr_multip(istate,1:state_dim)
    ! end do
    ! write(*,*)
    ! dens_test=get_proj_variational_density(lgr_multip_vec)
    ! write(*,*)
    ! write(*,*) "GA PROJECTORS VARIATIONAL DENSITY"
    ! write(*,*)
    ! do istate=1,state_dim
    !    write(*,'(10F18.10)') dens_test(istate,1:state_dim)
    ! end do
    ! write(*,*)
    ! dens_err=get_delta_proj_variational_density(lgr_multip_vec)
    ! write(*,*) "GA PROJECTORS VARIATIONAL DENSITY ERROR"
    ! write(*,*)
    ! call vec2mat_stride(dens_err,dens_test)
    ! do istate=1,state_dim
    !    write(*,'(10F18.10)') dens_test(istate,1:state_dim)
    ! end do
    ! write(*,*)
    !DEBUG>
    
    call get_GZproj_ground_state(n0_target,slater_derivatives,lgr_multip,E_Hloc,GZvect)
    ! write(*,'(A,20F6.2)') "GA PROJECTORS VECTOR",GZvect
    ! write(*,*) "GA PROJECTORS GROUND STATE ENERGY",E_Hloc

  contains


    ! function whose zeros determine the projectors lagrange multipliers !
    function get_delta_proj_variational_density(lm_) result(delta_proj_variational_density_vec)
      real(8),dimension(:)                   :: lm_
      real(8),dimension(state_dim*state_dim) :: delta_proj_variational_density_vec
      real(8),dimension(state_dim,state_dim) :: lm
      real(8),dimension(state_dim,state_dim) :: delta_proj_variational_density,proj_variational_density
      real(8),dimension(nFock,nFock)         :: H_projectors
      real(8),dimension(nFock)               :: H_eigens
      integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
      !
      call vec2mat_stride(lm_,lm)
      proj_variational_density=0.d0
      !+- build up the local H_projectors -+!
      H_projectors=phi_traces_basis_Hloc
      do istate=1,state_dim
         do jstate=1,state_dim
            H_projectors = H_projectors + slater_derivatives(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)/sqrt(n0_target(istate)*(1.d0-n0_target(istate)))
            H_projectors = H_projectors + lm(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
         end do
      end do
      call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
      do istate=1,state_dim
         do jstate=1,state_dim
            proj_variational_density(istate,jstate)=0.d0
            do ifock=1,nFock
               do jfock=1,nFock
                  proj_variational_density(istate,jstate)= &
                       proj_variational_density(istate,jstate) + &
                       H_projectors(ifock,1)*phi_traces_basis_dens(istate,jstate,ifock,jfock)*H_projectors(jfock,1)
               end do
            end do
         end do
      end do
      delta_proj_variational_density=proj_variational_density
      do istate=1,state_dim
         delta_proj_variational_density(istate,istate) = delta_proj_variational_density(istate,istate) - n0_target(istate)      
      end do
      call mat2vec_stride(delta_proj_variational_density,delta_proj_variational_density_vec)
    end function get_delta_proj_variational_density
    !
    function get_proj_variational_density(lm_) result(proj_variational_density)
      real(8),dimension(:)                   :: lm_
      real(8),dimension(state_dim,state_dim) :: proj_variational_density
      real(8),dimension(state_dim,state_dim) :: lm
      real(8),dimension(nFock,nFock)         :: H_projectors
      real(8),dimension(nFock)               :: H_eigens
      integer                                :: iorb,jorb,ispin,jspin,istate,jstate,ifock,jfock
      !
      call vec2mat_stride(lm_,lm)
      proj_variational_density=0.d0
      !+- build up the local H_projectors -+!
      H_projectors=phi_traces_basis_Hloc
      do istate=1,state_dim
         do jstate=1,state_dim
            H_projectors = H_projectors + slater_derivatives(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)/sqrt(n0_target(jstate)*(1.d0-n0_target(jstate)))
            H_projectors = H_projectors + lm(istate,jstate)*phi_traces_basis_dens(istate,jstate,:,:)
         end do
      end do
      !DEBUG
      ! do istate=1,nFock
      !    write(*,'(20F6.2)') H_projectors(istate,:)
      ! end do
      !DEBGU
      call matrix_diagonalize(H_projectors,H_eigens,'V','L')         
      !DEBUG
      ! write(*,'(20F6.2)') H_eigens
      ! write(*,*)
      ! do istate=1,nFock
      !    write(*,'(20F6.2)') H_projectors(istate,:)
      ! end do
      !DEBUG
      do istate=1,state_dim
         do jstate=1,state_dim
            proj_variational_density(istate,jstate)=0.d0
            do ifock=1,nFock
               do jfock=1,nFock
                  proj_variational_density(istate,jstate)= &
                       proj_variational_density(istate,jstate) + &
                       H_projectors(ifock,1)*phi_traces_basis_dens(istate,jstate,ifock,jfock)*H_projectors(jfock,1)
               end do
            end do
         end do
      end do
    end function get_proj_variational_density
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
          H_projectors = H_projectors + slater_derivatives(istate,jstate)*phi_traces_basis_Rhop(istate,jstate,:,:)/sqrt(n0(jstate)*(1.d0-n0(jstate)))
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
  subroutine slater_determinant_minimization_step(Rhop,n0_target,Estar,lgr_multip,slater_derivatives) 
    real(8),dimension(state_dim,state_dim),intent(in)    :: Rhop
    real(8),dimension(state_dim),intent(in)              :: n0_target
    real(8),intent(inout)                                :: Estar
    real(8),dimension(state_dim),intent(inout) :: lgr_multip
    real(8),dimension(state_dim,state_dim),intent(inout) :: slater_derivatives
    !
    real(8),dimension(state_dim,state_dim)               :: Hk,test_dens
    real(8),dimension(state_dim)                         :: ek,tmp_lgr,err_dens
    integer                                              :: iorb,jorb,ik,ispin,jspin,istate,jstate,info,korb
    real(8),dimension(state_dim)                         :: lgr_multip_vec

    ! initialize lagrange multipliers 
    lgr_multip=0.d0
    do istate=1,state_dim
       lgr_multip(istate)=0.01d0
    end do
    ! fix lagrange multipliers !
    lgr_multip=0.d0
    call fsolve(get_delta_local_density_matrix_diag,lgr_multip,tol=1.d-12,info=info)
    write(*,*)
    write(*,*) "SLATER DETERMINANT LAGRANGE MULTIPLIERS"
    write(*,*)
    write(*,'(20F18.10)') lgr_multip
    write(*,*)


    ! Get slater ground state informations : Estar & slater_derivatives

    !HERE IS THE MISTERY!!!!! NOW I'M NOT ABLE TO FIGURE-OUT!!!!
    ! DAI CHE e' QUASI FINITA!!

    call store_slater_ground_state(Rhop,lgr_multip,Estar,slater_derivatives)
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
    do istate=1,state_dim
!       write(*,"(20F18.10)") slater_derivatives(istate,:)
    end do
    write(*,*)
    !err_dens=get_delta_local_density_matrix_diag(lgr_multip)
    ! write(*,"(20F18.10)") err_dens
    ! write(*,*)
    ! write(*,*) 'SLATER DETERMINANT GROUND STATE ENERGY',Estar
    ! write(*,*) 
    !DEBUG>        
    !<DEBUG de le DEBUG
    do istate=1,state_dim
!       write(*,"(20F18.10)") slater_derivatives(istate,:)
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
    function get_delta_local_density_matrix_diag(lm_) result(delta_local_density_matrix_vec)
      real(8),dimension(:)  :: lm_
      real(8),dimension(state_dim)  :: delta_local_density_matrix_vec
      real(8),dimension(state_dim*state_dim)  :: delta_local_density_matrix_vec_
      real(8),dimension(state_dim,state_dim)  :: lm
      real(8),dimension(state_dim,state_dim)  :: delta_local_density_matrix,local_density_matrix
      real(8),dimension(state_dim,state_dim) :: Hk,tmp
      real(8),dimension(state_dim)          :: ek
      integer                              :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik
      !
      lm=0.d0
      do istate=1,state_dim
         lm(istate,istate)=lm_(istate)
      end do
      local_density_matrix=0.d0
      do ik=1,Lk
         Hk=0.d0
         ek=0.d0
         ! build-up the hopping hamiltonian !         
         do iorb=1,Norb
            do jorb=1,Norb
               do ispin=1,2
                  do jspin=1,2
                     istate=index(ispin,iorb)
                     jstate=index(jspin,jorb)                                            
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
         ! hopping renormalization !
         Hk=matmul(Hk,Rhop)
         Hk=matmul(Rhop,Hk)
         ! add Lagrange multipliers !
         Hk=Hk+lm                     
         ! diagonalize hamiltonian !
         call  matrix_diagonalize(Hk,ek,'V','L')
         !compute local density matrix
         do istate=1,state_dim
            do jstate=1,state_dim
               do kstate=1,state_dim
                  local_density_matrix(istate,jstate) = &
                       local_density_matrix(istate,jstate) + fermi(ek(kstate),beta)*Hk(istate,kstate)*Hk(jstate,kstate)*wtk(ik)                        
               end do
            end do
         end do
      end do
      ! return variation of local density matrix with respect to the target values
      delta_local_density_matrix = local_density_matrix
      do istate=1,state_dim
         delta_local_density_matrix(istate,istate) = delta_local_density_matrix(istate,istate) - n0_target(istate)      
      end do
      do istate=1,state_dim
         delta_local_density_matrix_vec(istate)=delta_local_density_matrix(istate,istate)
      end do
    end function get_delta_local_density_matrix_diag
    !
    function get_delta_local_density_matrix_full(lm_) result(delta_local_density_matrix_vec)
      real(8),dimension(:)  :: lm_
      real(8),dimension(state_dim,state_dim)  :: lm_full
      real(8),dimension(state_dim*state_dim)  :: delta_local_density_matrix_vec
      real(8),dimension(state_dim,state_dim)  :: lm
      real(8),dimension(state_dim,state_dim)  :: delta_local_density_matrix,local_density_matrix
      real(8),dimension(state_dim,state_dim) :: Hk,tmp
      real(8),dimension(state_dim)          :: ek
      integer                              :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik
      !
      call vec2mat_stride(lm_,lm)
      local_density_matrix=0.d0
      do ik=1,Lk
         Hk=0.d0
         ek=0.d0
         ! build-up the hopping hamiltonian !         
         do iorb=1,Norb
            do jorb=1,Norb
               do ispin=1,2
                  do jspin=1,2
                     istate=index(ispin,iorb)
                     jstate=index(jspin,jorb)                                            
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
         ! hopping renormalization !
         Hk=matmul(Hk,Rhop)
         Hk=matmul(Rhop,Hk)
         ! add Lagrange multipliers !
         Hk=Hk+lm                     
         ! diagonalize hamiltonian !
         call  matrix_diagonalize(Hk,ek,'V','L')         
         !compute local density matrix
         do istate=1,state_dim
            do jstate=1,state_dim
               do kstate=1,state_dim
                  local_density_matrix(istate,jstate) = &
                       local_density_matrix(istate,jstate) + fermi(ek(kstate),beta)*Hk(istate,kstate)*Hk(jstate,kstate)*wtk(ik)                        
               end do
            end do
         end do
      end do
      ! return variation of local density matrix with respect to the target values
      delta_local_density_matrix = local_density_matrix
      do istate=1,state_dim
         delta_local_density_matrix(istate,istate) = delta_local_density_matrix(istate,istate) - n0_target(istate)      
      end do
      call mat2vec_stride(delta_local_density_matrix,delta_local_density_matrix_vec)
    end function get_delta_local_density_matrix_full
    !



    function get_local_density_matrix_diag(lm_) result(local_density_matrix)
      real(8),dimension(state_dim) :: lm_
      real(8),dimension(state_dim,state_dim) :: lm
      real(8),dimension(state_dim,state_dim) :: local_density_matrix
      real(8),dimension(state_dim,state_dim) :: Hk
      real(8),dimension(state_dim)           :: ek
      integer                                :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik
      lm=0.d0
      do istate=1,state_dim
         lm(istate,istate)=lm_(istate)
      end do
      local_density_matrix=0.d0
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
         ! hopping renormalization !
         Hk=matmul(Hk,Rhop)
         Hk=matmul(Rhop,Hk)
         ! add Lagrange multipliers !
         Hk=Hk+lm                     
         ! diagonalize hamiltonian !
         call  matrix_diagonalize(Hk,ek,'V','L')         
         !compute local density matrix
         do iorb=1,Norb
            do ispin=1,2
               do jorb=1,Norb
                  do jspin=1,2
                     istate=index(ispin,iorb)
                     jstate=index(jspin,jorb)               
                     do kstate=1,state_dim
                        local_density_matrix(istate,jstate) = &
                             local_density_matrix(istate,jstate) + fermi(ek(kstate),beta)*Hk(istate,kstate)*Hk(jstate,kstate)*wtk(ik)                        
                     end do
                  end do
               end do
            end do
         end do
      end do
    end function get_local_density_matrix_diag
    !
    function get_local_density_matrix_full(lm) result(local_density_matrix)
      real(8),dimension(state_dim,state_dim)  :: lm
      real(8),dimension(state_dim,state_dim)  :: local_density_matrix
      real(8),dimension(state_dim,state_dim) :: Hk
      real(8),dimension(state_dim)          :: ek
      integer                              :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik
      local_density_matrix=0.d0
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
         ! hopping renormalization !
         Hk=matmul(Hk,Rhop)
         Hk=matmul(Rhop,Hk)
         ! add Lagrange multipliers !
         Hk=Hk+lm                     
         ! diagonalize hamiltonian !
         call  matrix_diagonalize(Hk,ek,'V','L')         
         !compute local density matrix
         do iorb=1,Norb
            do ispin=1,2
               do jorb=1,Norb
                  do jspin=1,2
                     istate=index(ispin,iorb)
                     jstate=index(jspin,jorb)               
                     do kstate=1,state_dim
                        local_density_matrix(istate,jstate) = &
                             local_density_matrix(istate,jstate) + fermi(ek(kstate),beta)*Hk(istate,kstate)*Hk(jstate,kstate)*wtk(ik)                        
                     end do
                  end do
               end do
            end do
         end do
      end do
    end function get_local_density_matrix_full
  end subroutine slater_determinant_minimization_step


  subroutine store_slater_ground_state(Rhop,lm,Estar,slater_derivatives)     
    real(8),dimension(state_dim,state_dim) :: Rhop
    real(8),dimension(state_dim,state_dim) :: lm
    real(8)                                :: Estar
    real(8),dimension(state_dim,state_dim),intent(inout) :: slater_derivatives
    real(8),dimension(state_dim,state_dim) :: Hk,tmp,Hk_bare
    real(8),dimension(state_dim)           :: ek
    integer                                :: iorb,jorb,ispin,jspin,istate,jstate,kstate,ik
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
       ! add Lagrange multipliers !
       Hk=Hk+lm                     
       ! diagonalize hamiltonian !
       call  matrix_diagonalize(Hk,ek,'V','L')
       ! get ground state energy !       
       do kstate=1,state_dim
          Estar= Estar + fermi(ek(kstate),beta)*ek(kstate)*wtk(ik)
       end do
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
END MODULE GZ_ENERGY_FUNCTIONAL_SELF

