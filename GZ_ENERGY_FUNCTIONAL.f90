MODULE GZ_ENERGY_FUNCTIONAL
  USE SF_LINALG
  USE SF_IOTOOLS
  USE SF_SPECIAL
  USE SF_OPTIMIZE
  USE GZ_VARS_GLOBAL
  USE LANCELOT_simple_double
  USE GZ_PROJECTORS
  implicit none
  private

  public :: gz_energy_local_density,gz_energy_local_polar,gz_ground_state_energy_estimation
  real(8),dimension(:),allocatable,public :: local_dens_target
  real(8),dimension(:),allocatable :: phi_tmp

contains




  function gz_energy_local_density(n0)  result(GZ_energy)
    real(8),intent(in)           :: n0(:)
    !real(8)           :: n0(:)  
    real(8)                      :: GZ_energy
    real(8),dimension(state_dim) :: ni
    real(8),dimension(nFock)     :: phi_
    logical                      :: converged
    integer                      :: iter,iter_max,istate,i,iter_
    integer                      :: iorb,ispin,i_ind,ifock
    integer                      :: n_min,neq,nin,maxit,print_level,exit_code
    real(8)                      :: gradtol,feastol,Ephi,nsite,err_iter,Ephi_
    real(8),allocatable          :: bL(:),bU(:),cx(:),y(:),phi_optimize(:),phi_optimize_(:)
    integer                      :: iunit,err_unit,ene_unit,Nsuccess
    logical                      :: bound

    !external  :: energy_gz_projectors
    external  :: energy_gz_projectors_iteration

    !+- initialize GZ_projectors -+!
    if(size(n0).ne.Norb) stop "size(n0) /= Norb"
    if(.not.allocated(phi_gz)) then
       allocate(phi_gz(nFock))
       phi_gz=1.d0/sqrt(1.d0*nFock)
       write(*,*) phi_gz
    end if
    if(size(phi_gz).ne.nFock) stop "size(Phi_gz) /= nFock"
    allocate(local_dens_target(Norb))    
    local_dens_target=n0    
    bound=.false.
    do iorb=1,Norb
       if(local_dens_target(iorb).le.1.d-10.or.local_dens_target(iorb).ge.2.d0-1.d-10) bound=.true.
    end do
    !+------------------------------------------------------------------------------------------------+!
    !+- optimize GZ_energy functional with respect to the GZ projectors - LANCELOT ROUTINE -         -+!    
    !+- inside the routine energy_gz_projectors the slater determinant is optimized for each values  -+!
    !+- of the GZ projectors                                                                         -+!
    !+------------------------------------------------------------------------------------------------+!
    if(.not.bound) then
       allocate(phi_optimize(nFock_indep),phi_optimize_(nFock_indep))
       do i_ind=1,nFock_indep
          phi_optimize(i_ind) = phi_gz(fock_indep(i_ind))
       end do
       iunit=free_unit()
       open(iunit,file='Lancelot_call.out',access='append')
       err_unit=free_unit()
       open(err_unit,file='error.err',access='append')       
       ene_unit=free_unit()
       open(ene_unit,file='Lancelot_energy.out',access='append')       
       Nsuccess=0
       do iter_=1,100
          if(iter_.gt.1) Ephi_=Ephi
          phi_optimize_=phi_optimize          
          do ifock=1,nFock
             phi_(ifock)=phi_optimize(full2indep_fock(ifock))     
          end do
          call fix_slater_matrix_elements(phi_,local_dens_target)          
          n_min       = nFock_indep ! number of minimization parameters
          neq         = Norb + 1    ! number of equality constraints                   
          nin         = 0           ! number of in-equality constraints                   
          maxit       = 1000        ! maximum iteration number 
          gradtol     = 1.d-7       ! maximum norm of the gradient at convergence 
          feastol     = 1.d-7       ! maximum violation of parameters at convergence  
          print_level = lancelot_verbose           ! verbosity
          allocate(bl(n_min),bu(n_min),cx(neq+nin),y(neq+nin))
          bL = 0.d0               ! lower bounds for minimization parameters
          bU = 1.d0-bL                 ! upper bounds for minimization parameters          
          write(iunit,*) '!+---------------------------------------------------+!'
          write(iunit,*) '!+---------------------------------------------------+!'
          write(iunit,*) 'Number of minimization parameters', n_min
          write(iunit,*) 'Number of equality constrains', neq
          write(iunit,*) 'Number of inquality constrains', nin
          write(iunit,*) 'Maximum number of iterations', maxit
          write(iunit,*) 'Gradient tolerance', gradtol
          write(iunit,*) 'Constraints tolerance', feastol
          write(iunit,*) 
          write(iunit,*) 
          write(iunit,*) 'INPUT PARAMETERS'
          do i=1,n_min
             write(iunit,*) phi_optimize(i)
          end do
          write(iunit,*) 
          write(iunit,*)           
          call lancelot_simple(n_min,phi_optimize,Ephi,exit_code,my_fun=energy_gz_projectors_iteration, &
               bl = bl, bu = bu,                                                              &
               neq = neq, nin = nin,                                                          &
               cx = cx, y = y, iters  = iter, maxit = maxit,                                 &
               gradtol = gradtol, feastol = feastol,                                          &
               print_level = print_level )
          !+--------------------------------------------------------------------------------------+!    
          write(iunit,*) 'LANCELOT EXIT STATUS'
          write(iunit,*)
          write(iunit,*) exit_code
          write(iunit,*) 'OPTIMIZED PARAMETERS'
          do i=1,n_min
             write(iunit,*) phi_optimize(i)
          end do
          write(iunit,*) 
          write(iunit,*) 
          write(iunit,*) 'FINAL VALUE'
          write(iunit,*) Ephi
          write(iunit,*) 
          write(iunit,*) 
          write(iunit,*)  'CONSTRAINTS'
          do i=1,neq+nin
             write(iunit,*) cx(i)
          end do
          write(iunit,*) 
          write(iunit,*) 
          write(iunit,*)  'LAGRANGE MULTIPLIERS'
          do i=1,neq+nin
             write(iunit,*) y(i)
          end do
          write(iunit,*) 
          write(iunit,*) 
          write(iunit,*)  'NUMBER OF ITERATION'
          write(iunit,*) iter,'/',maxit       
          deallocate(bl,bu,cx,y)
          write(iunit,*) 'Eiteration',Ephi
          write(ene_unit,'(10(F18.10))') Ephi,local_dens_target,dble(iter_)
          if(iter_.eq.1) then
             err_iter=1.d0
          else
             err_iter=abs(Ephi_-Ephi)
          end if
          if(err_iter.lt.1.d-6) Nsuccess = Nsuccess + 1
          if(err_iter.lt.1.d-6.and.Nsuccess.gt.10) exit
          write(err_unit,*) dble(iter_),err_iter,Nsuccess          
       end do
       close(iunit)
       close(err_unit)
       close(ene_unit)
       do ifock=1,nFock
          phi_gz(ifock) = phi_optimize(full2indep_fock(ifock))
       end do
       GZ_energy=Ephi
    else
       GZ_energy=100.d0
    end if
    deallocate(local_dens_target)    
  end function gz_energy_local_density

  
  function gz_energy_local_polar(pol)  result(GZ_energy)
    real(8),intent(in)           :: pol(:)  
    real(8)                      :: GZ_energy
    real(8),dimension(state_dim) :: ni
    real(8),dimension(nFock)     :: phi_
    logical                      :: converged
    integer                      :: iter,iter_max,istate,i,iter_
    integer                      :: iorb,ispin,i_ind,ifock
    integer                      :: n_min,neq,nin,maxit,print_level,exit_code
    real(8)                      :: gradtol,feastol,Ephi,nsite,err_iter,Ephi_
    real(8),allocatable          :: bL(:),bU(:),cx(:),y(:),phi_optimize(:),phi_optimize_(:)
    integer                      :: iunit,ene_unit,err_unit,Nsuccess,energy_unit
    logical                      :: bound
    external                     :: energy_gz_projectors_iteration
    
    !+- initialize GZ_projectors -+!
    if(size(pol).ne.1) stop "size(pol) /= 1"
    if(.not.allocated(phi_gz)) then
       allocate(phi_gz(nFock))
       phi_gz=1.d0/sqrt(1.d0*nFock)
       write(*,*) phi_gz
    end if
    if(size(phi_gz).ne.nFock) stop "size(Phi_gz) /= nFock"
    allocate(local_dens_target(Norb))    
    do iorb=1,Norb
       local_dens_target(iorb) = Nread + (-1)**dble(iorb)*pol(1)*0.5d0
    end do
    bound=.false.
    do iorb=1,Norb
       if(local_dens_target(iorb).le.1.d-10.or.local_dens_target(iorb).ge.2.d0-1.d-10) bound=.true.
    end do
    !+------------------------------------------------------------------------------------------------+!
    !+- optimize GZ_energy functional with respect to the GZ projectors - LANCELOT ROUTINE -         -+!    
    !+- inside the routine energy_gz_projectors the slater determinant is optimized for each values  -+!
    !+- of the GZ projectors                                                                         -+!
    !+------------------------------------------------------------------------------------------------+!
    
    open(energy_unit,file='energy_fixed_polarization.out',access='append')
    if(.not.bound) then
       allocate(phi_optimize(nFock_indep),phi_optimize_(nFock_indep))
       do i_ind=1,nFock_indep
          phi_optimize(i_ind) = phi_gz(fock_indep(i_ind))
       end do
       iunit=free_unit()
       open(iunit,file='Lancelot_call.out',access='append')
       err_unit=free_unit()
       open(err_unit,file='error.err',access='append')       
       ene_unit=free_unit()
       open(ene_unit,file='Lancelot_energy.out',access='append')       
       Nsuccess=0
       do iter_=1,100
          if(iter_.gt.1) Ephi_=Ephi
          phi_optimize_=phi_optimize          
          do ifock=1,nFock
             phi_(ifock)=phi_optimize(full2indep_fock(ifock))     
          end do
          call fix_slater_matrix_elements(phi_,local_dens_target)          
          n_min       = nFock_indep ! number of minimization parameters
          neq         = Norb + 1    ! number of equality constraints                   
          nin         = 0           ! number of in-equality constraints                   
          maxit       = 1000        ! maximum iteration number 
          gradtol     = 1.d-7       ! maximum norm of the gradient at convergence 
          feastol     = 1.d-7       ! maximum violation of parameters at convergence  
          print_level = lancelot_verbose           ! verbosity
          allocate(bl(n_min),bu(n_min),cx(neq+nin),y(neq+nin))
          bL = 0.d0               ! lower bounds for minimization parameters
          bU = 1.d0                 ! upper bounds for minimization parameters          
          write(iunit,*) '!+---------------------------------------------------+!'
          write(iunit,*) '!+---------------------------------------------------+!'
          write(iunit,*) 'Number of minimization parameters', n_min
          write(iunit,*) 'Number of equality constrains', neq
          write(iunit,*) 'Number of inquality constrains', nin
          write(iunit,*) 'Maximum number of iterations', maxit
          write(iunit,*) 'Gradient tolerance', gradtol
          write(iunit,*) 'Constraints tolerance', feastol
          write(iunit,*) 
          write(iunit,*) 
          write(iunit,*) 'INPUT PARAMETERS'
          do i=1,n_min
             write(iunit,*) phi_optimize(i)
          end do
          write(iunit,*) 
          write(iunit,*)           
          call lancelot_simple(n_min,phi_optimize,Ephi,exit_code,my_fun=energy_gz_projectors_iteration, &
               bl = bl, bu = bu,                                                              &
               neq = neq, nin = nin,                                                          &
               cx = cx, y = y, iters  = iter, maxit = maxit,                                 &
               gradtol = gradtol, feastol = feastol,                                          &
               print_level = print_level )
          !+--------------------------------------------------------------------------------------+!    
          write(iunit,*) 'LANCELOT EXIT STATUS'
          write(iunit,*)
          write(iunit,*) exit_code
          write(iunit,*) 'OPTIMIZED PARAMETERS'
          do i=1,n_min
             write(iunit,*) phi_optimize(i)
          end do
          write(iunit,*) 
          write(iunit,*) 
          write(iunit,*) 'FINAL VALUE'
          write(iunit,*) Ephi
          write(iunit,*) 
          write(iunit,*) 
          write(iunit,*)  'CONSTRAINTS'
          do i=1,neq+nin
             write(iunit,*) cx(i)
          end do
          write(iunit,*) 
          write(iunit,*) 
          write(iunit,*)  'LAGRANGE MULTIPLIERS'
          do i=1,neq+nin
             write(iunit,*) y(i)
          end do
          write(iunit,*) 
          write(iunit,*) 
          write(iunit,*)  'NUMBER OF ITERATION'
          write(iunit,*) iter,'/',maxit       
          deallocate(bl,bu,cx,y)
          write(iunit,*) 'Eiteration',Ephi
          write(ene_unit,*) Ephi,iter
          if(iter_.eq.1) then
             err_iter=1.d0
          else
             err_iter=abs(Ephi_-Ephi)
          end if
          write(err_unit,*) dble(iter_),err_iter          
          if(err_iter.lt.1.d-6) Nsuccess = Nsuccess + 1
          if(err_iter.lt.1.d-6.and.Nsuccess.gt.4) then
             write(err_unit,*)
             write(err_unit,*) err_iter,Nsuccess
             write(err_unit,*)
             exit
          end if
          phi_optimize=0.5d0*phi_optimize_+0.5d0*phi_optimize          
       end do
       close(iunit)
       close(err_unit)

       do ifock=1,nFock
          phi_gz(ifock) = phi_optimize(full2indep_fock(ifock))
       end do
       GZ_energy=Ephi
       write(energy_unit,*) GZ_energy,local_dens_target
    else
       GZ_energy=100.d0
       write(energy_unit,*) GZ_energy,local_dens_target
    end if
    deallocate(local_dens_target)    
    close(energy_unit)
  end function gz_energy_local_polar
  
  

  subroutine fix_slater_matrix_elements(phi,local_dens) 
    real(8),dimension(nFock)     :: phi
    real(8)                      :: energy_star,energy_star_,energy_star_test
    real(8),dimension(2)         :: mu,local_dens,dens_test,dens_test_
    real(8),dimension(state_dim) :: local_dens_
    real(8),dimension(size(mu),size(mu)) :: Hk,Hk_
    real(8),dimension(size(mu))          :: ek
    real(8),dimension(state_dim) :: Rhop    
    real(8)                      :: eps
    integer                      :: iorb,jorb,ik,ispin,jspin,istate,jstate,info,korb
    !mu=mu_slater   
    do iorb=1,Norb
       do ispin=1,2
          istate=index(ispin,iorb)
          local_dens_(istate)=local_dens(iorb)*0.5d0
       end do
    end do
    Rhop=gz_Rhop_dens(phi,local_dens_,CC,CA)            
    mu(1)=-1.d-6
    mu(2)=1.d-6
    !+- find mu_star -+!
    call fsolve(fix_density,mu,tol=1.d-10,info=info)
    ! mu(1)=-0.4d0
    ! mu(2)=0.4d0
    write(95,'(10F18.10)') Rhop
    dens_test_=get_density(mu)
    dens_test=fix_density(mu)
    write(96,'(10F18.10)') dens_test_,dens_test,mu,dble(info),local_dens_target
    !stop
    
    energy_star=0.d0
    energy_star_=0.d0
    energy_star_test=0.d0
    do ik=1,Lk
       ek=0.d0
       Hk=0.d0
       do iorb=1,2
          do jorb=1,2
             istate=index(1,iorb)
             jstate=index(1,jorb)
             if(iorb.eq.jorb) then
                Hk(iorb,jorb)=epsik(ik)*Rhop(istate)*Rhop(jstate) - mu(iorb)
             else
                Hk(iorb,jorb)=hybik(ik)*Rhop(istate)*Rhop(jstate)
             end if
          end do
       end do
       Hk_=Hk
       call  matrix_diagonalize(Hk,ek,'V','L')                

       !<DEBUG
       do iorb=1,Norb
          do ispin=1,2
             energy_star_ = energy_star_ + fermi(ek(iorb)-mu(iorb),beta)*ek(iorb)*wtk(ik)
          end do
       end do
       !DEBUG

!       if(ik.eq.15) write(*,'(10(F18.10))') Hk

       do iorb=1,Norb
          do jorb=1,Norb
             do ispin=1,2
                do jspin=1,2
                   istate=index(ispin,iorb)
                   jstate=index(jspin,jorb)                   
                   slater_matrix_elements(istate,jstate,ik)=0.d0
                   do korb=1,Norb
                      slater_matrix_elements(istate,jstate,ik) = slater_matrix_elements(istate,jstate,ik) + Hk(jorb,korb)*fermi(ek(korb),beta)*Hk(iorb,korb)*0.5d0
                      !NOTE: this 0.5d0 factor is needed because  fermi(ek(korb)- mu(korb),beta) is already summed over spin, a consequence
                      !      of the fact that the hamiltonian Hk has dimension (Norb,Norb) instaed of (state_dim,state_dim) 
                      !      where state_dim=Norb*Nspin. Pay attention! This factor should be removed when I will eventually make uniform the
                      !      code with respect to the standard structure (state_dim,state_dim) with nFock = 2^state_dim.
                   end do
                end do
             end do
          end do
       end do


       !<DEBUG
       ! if(ik.eq.15) then
       !    write(*,*)
       !    do istate=1,state_dim
       !       write(*,'(16(F18.10))') slater_matrix_elements(istate,:,ik)
       !    end do
       ! end if



       ! do iorb=1,Norb
       !    do jorb=1,Norb
       !       do ispin=1,2
       !          do jspin=1,2
       !             istate=index(ispin,iorb)
       !             jstate=index(jspin,jorb)
       !             energy_star=energy_star+slater_matrix_elements(istate,jstate,ik)*Hk_(jorb,iorb)*wtk(ik)
       !             !write(*,'(10(F18.10))') slater_matrix_elements(istate,jstate,ik),Hk_(iorb,jorb),dble(iorb),dble(jorb)
       !          end do
       !       end do
       !    end do
       ! end do
       !DEBUG>

    end do
    !<DEBUG
    ! write(*,*) 'energy_star'
    ! write(*,*) energy_star,energy_star_,energy_star_test
    ! stop
    !DEBUG>
  contains    
    function fix_density(mu) result(local_dens)
      real(8),dimension(:)                 :: mu
      real(8),dimension(size(mu))          :: local_dens
      real(8),dimension(size(mu),size(mu)) :: Hk
      real(8),dimension(size(mu))          :: ek
      integer                              :: iorb,jorb,ik,ispin,istate,jstate
      local_dens=0.d0
      do ik=1,Lk
         Hk=0.d0
         ek=0.d0
         do iorb=1,2
            do jorb=1,2
               istate=index(1,iorb)
               jstate=index(1,jorb)               
               if(iorb.eq.jorb) then
                  Hk(iorb,jorb)=epsik(ik)*Rhop(istate)*Rhop(jstate)  - mu(iorb)
               else
                  Hk(iorb,jorb)=hybik(ik)*Rhop(istate)*Rhop(jstate)
               end if
            end do
         end do
         call  matrix_diagonalize(Hk,ek,'V','L')         
         
         ! do ispin=1,2
         !    do iorb=1,Norb
         !       local_dens(iorb) = local_dens(iorb) + fermi(ek(iorb)   ,beta)*wtk(ik)
         !    end do
         ! end do

         do ispin=1,2
            do iorb=1,Norb
               do jorb=1,Norb
                  local_dens(iorb) = local_dens(iorb) + fermi(ek(jorb),beta)*Hk(iorb,jorb)*Hk(iorb,jorb)*wtk(ik)
               end do
            end do
         end do
      end do
      ! write(*,*) local_dens
      ! write(*,*) mu
      local_dens = local_dens - local_dens_target      
    end function fix_density
    function get_density(mu) result(local_dens)
      real(8),dimension(:)                 :: mu
      real(8),dimension(size(mu))          :: local_dens
      real(8),dimension(size(mu),size(mu)) :: Hk
      real(8),dimension(size(mu))          :: ek
      integer                              :: iorb,jorb,ik,ispin,istate,jstate
      local_dens=0.d0
      do ik=1,Lk
         Hk=0.d0
         ek=0.d0
         do iorb=1,2
            do jorb=1,2
               istate=index(1,iorb)
               jstate=index(1,jorb)               
               if(iorb.eq.jorb) then
                  Hk(iorb,jorb)=epsik(ik)*Rhop(istate)*Rhop(jstate) - mu(iorb)
               else
                  Hk(iorb,jorb)=hybik(ik)*Rhop(istate)*Rhop(jstate)
               end if
            end do
         end do
         call  matrix_diagonalize(Hk,ek,'V','L')  
         do ispin=1,2
            do iorb=1,Norb
               do jorb=1,Norb
                  local_dens(iorb) = local_dens(iorb) + fermi(ek(jorb),beta)*Hk(iorb,jorb)*Hk(iorb,jorb)*wtk(ik)
               end do
            end do
         end do
         ! do ispin=1,2
         !    do iorb=1,Norb
         !       local_dens(iorb) = local_dens(iorb) + fermi(ek(iorb) ,beta)*wtk(ik)
         !    end do
         ! end do
      end do
    end function get_density

  end subroutine fix_slater_matrix_elements




  function gz_ground_state_energy_estimation(x) result(Egz)
    implicit none
    !+- routine variables -+!
    real(8),dimension(:)           :: x
    real(8),dimension(3)          :: Egz
    real(8)                       :: f

    real(8),allocatable           :: phi_(:)
    real(8)                       :: nsite
    real(8),allocatable           :: niorb(:),niorb_(:)
    real(8)                       :: Estar
    real(8),dimension(state_dim) :: local_dens_
    real(8),dimension(state_dim) :: Rhop    
    real(8),dimension(Norb,Norb) :: Hk

    integer                       :: iorb,ispin,istate,jstate,ik,ifock,jorb,jspin
    !
    !write(*,*) x(:)    
    allocate(phi_(nFock))
    !

    do ifock=1,nFock
       phi_(ifock)=x(ifock)     
    end do
    ! write(*,*) phi_
    ! stop    
    !
    
    !    nsite=sum(local_dens_target)
    ! write(*,*) local_dens_target
    ! stop
    ! !
    ! do iorb=1,Norb
    !    niorb(iorb) =  local_dens_target(iorb)
    ! end do

    ! !
    ! do iorb=1,Norb
    !    do ispin=1,2
    !       istate=index(ispin,iorb)
    !       local_dens_(istate)=local_dens_target(iorb)*0.5d0
    !    end do
    ! end do

    allocate(niorb(Norb))
    do istate=1,state_dim
       local_dens_=gz_local_diag(phi_,dens(istate,:,:))
    end do
    nsite=sum(local_dens_)
    write(*,'(4(F15.8))') local_dens_
    Rhop=gz_Rhop_dens(phi_,local_dens_,CC,CA)            
    !

    !+- FREE ENERGY ESTIMATION WITHIN GZ APPROX -+!

    !- quasiparticles kinetic energy -
    Estar=0.d0
    do ik=1,Lk
       Hk=0.d0
       do iorb=1,Norb
          do jorb=1,Norb
             istate=index(1,iorb)
             jstate=index(1,jorb)
             if(iorb.eq.jorb) then
                Hk(iorb,jorb)=epsik(ik)*Rhop(istate)*Rhop(jstate) 
             else
                Hk(iorb,jorb)=hybik(ik)*Rhop(istate)*Rhop(jstate)
             end if
          end do
       end do

       do iorb=1,Norb
          do jorb=1,Norb
             do ispin=1,2
                do jspin=1,2
                   istate=index(ispin,iorb)
                   jstate=index(jspin,jorb)
                   Estar = Estar + slater_matrix_elements(istate,jstate,ik)*wtk(ik)*Hk(iorb,jorb)
                end do
             end do
          end do
       end do
    end do

    ! - local hamiltonian average -
    f=0.d0     
    f = f + gz_local_diag(phi_,UHubbard)*U*0.5d0
    do iorb=1,Norb
       do ispin=1,2
          istate=index(ispin,iorb)
          f = f + gz_local_diag(phi_,dens(istate,:,:))*Eloc(iorb)
          f = f - U*gz_local_diag(phi_,dens(istate,:,:))
          f = f - gz_local_diag(phi_,dens(istate,:,:))*xmu
       end do
    end do
    
    Egz(1)=f+Estar
    Egz(2)=Estar
    Egz(3)=f

  end function gz_ground_state_energy_estimation
  



END MODULE GZ_ENERGY_FUNCTIONAL



! subroutine energy_gz_projectors(x,f,i)
!   USE GZ_VARS_GLOBAL
!   USE GZ_ENERGY_FUNCTIONAL
!   USE GZ_PROJECTORS
!   implicit none
!   !+- routine variables -+!
!   real(8), intent(in)           :: x(:)
!   real(8), intent(out)          :: f
!   integer, intent(in), optional :: i
!   real(8),allocatable           :: phi_(:)
!   real(8)                       :: nsite
!   real(8),allocatable           :: niorb(:)
!   real(8)                          :: Estar
!   integer                       :: iorb,ispin,istate,ifock
!   !
!   allocate(phi_(nFock))
!   allocate(niorb(Norb))
!   !
!   do ifock=1,nFock
!      phi_(ifock)=x(full2indep_fock(ifock))     
!   end do
!   !
!   nsite=sum(local_dens_target)
!   !
!   do iorb=1,Norb
!      niorb(iorb) =  local_dens_target(iorb)
!   end do

!   !Estar=Hstar_energy(phi_,niorb)
!   !Estar=slater_energy(phi_,niorb)
!   Estar=Hstar(phi_,niorb)


!   ! stop

!   if (.not.present(i)) then
!      !+- FREE ENERGY ESTIMATION WITHIN GZ APPROX -+!
!      f=0.d0     
!      f = f + gz_local_diag(phi_,UHubbard)*U*0.5d0
!      do iorb=1,Norb
!         do ispin=1,2
!            istate=index(ispin,iorb)
!            f = f + gz_local_diag(phi_,dens(istate,:,:))*Eloc(iorb)
!            f = f - U*gz_local_diag(phi_,dens(istate,:,:))
!            f = f - gz_local_diag(phi_,dens(istate,:,:))*xmu
!         end do
!      end do
!      f = f + Estar     
!      write(94,*) f,Estar
!   else
!      !+- CONSTRAINTS ON GUTZWILLER PARAMETERS -+!

!      if(Norb.eq.2) then

!         select case(i)
!         case(1)
!            f=0.d0
!            iorb=1
!            do ispin=1,2
!               istate=index(ispin,iorb)
!               f=f+gz_local_diag(phi_,dens(istate,:,:))
!            end do
!            f=f-niorb(iorb)
!         case(2)
!            f=0.d0
!            iorb=2
!            do ispin=1,2
!               istate=index(ispin,iorb)
!               f=f+gz_local_diag(phi_,dens(istate,:,:))
!            end do
!            f=f-niorb(iorb)
!         case(3)
!            f=0.d0
!            do ifock=1,nFock
!               f=f+phi_(ifock)*phi_(ifock)
!            end do
!            f=f-1.d0
!         end select
!      else
!         select case(i)
!         case(1)
!            f=0.d0
!            iorb=1
!            do ispin=1,2
!               istate=index(ispin,iorb)
!               f=f+gz_local_diag(phi_,dens(istate,:,:))
!            end do
!            f=f-niorb(iorb)
!         case(2)
!            f=0.d0
!            do ifock=1,nFock
!               f=f+phi_(ifock)*phi_(ifock)
!            end do
!            f=f-1.d0
!         end select
!      end if



!   end if


! end subroutine energy_gz_projectors




subroutine energy_gz_projectors_iteration(x,f,i)
  USE GZ_VARS_GLOBAL
  USE GZ_ENERGY_FUNCTIONAL
  USE GZ_PROJECTORS
  implicit none
  !+- routine variables -+!
  real(8), intent(in)           :: x(:)
  real(8), intent(out)          :: f
  integer, intent(in), optional :: i
  real(8),allocatable           :: phi_(:)
  real(8)                       :: nsite
  real(8),allocatable           :: niorb(:)
  real(8)                          :: Estar
  real(8),dimension(state_dim) :: local_dens_
  real(8),dimension(state_dim) :: Rhop    
  real(8),dimension(Norb,Norb) :: Hk

  integer                       :: iorb,ispin,istate,jstate,ik,ifock,jorb,jspin
  !
  allocate(phi_(nFock))
  allocate(niorb(Norb))
  !
  do ifock=1,nFock
     phi_(ifock)=x(full2indep_fock(ifock))     
  end do
  !
  nsite=sum(local_dens_target)
  !
  do iorb=1,Norb
     niorb(iorb) =  local_dens_target(iorb)
  end do

  !
  do iorb=1,Norb
     do ispin=1,2
        istate=index(ispin,iorb)
        local_dens_(istate)=local_dens_target(iorb)*0.5d0
     end do
  end do
  Rhop=gz_Rhop_dens(phi_,local_dens_,CC,CA)            
  !



  if (.not.present(i)) then
     !+- FREE ENERGY ESTIMATION WITHIN GZ APPROX -+!

     Estar=0.d0

     do ik=1,Lk
        Hk=0.d0
        do iorb=1,Norb
           do jorb=1,Norb
              istate=index(1,iorb)
              jstate=index(1,jorb)
              if(iorb.eq.jorb) then
                 Hk(iorb,jorb)=epsik(ik)*Rhop(istate)*Rhop(jstate) 
              else
                 Hk(iorb,jorb)=hybik(ik)*Rhop(istate)*Rhop(jstate)
              end if
           end do
        end do
        
        do iorb=1,Norb
           do jorb=1,Norb
              do ispin=1,2
                 do jspin=1,2
                    istate=index(ispin,iorb)
                    jstate=index(jspin,jorb)
                    Estar = Estar + slater_matrix_elements(istate,jstate,ik)*wtk(ik)*Hk(iorb,jorb)
                 end do
              end do
           end do
        end do
     end do
     
     f=0.d0     
     f = f + gz_local_diag(phi_,UHubbard)*U*0.5d0
     do iorb=1,Norb
        do ispin=1,2
           istate=index(ispin,iorb)
           f = f + gz_local_diag(phi_,dens(istate,:,:))*Eloc(iorb)
           f = f - U*gz_local_diag(phi_,dens(istate,:,:))
           f = f - gz_local_diag(phi_,dens(istate,:,:))*xmu
        end do
     end do
     f = f + Estar     
  else
     !+- CONSTRAINTS ON GUTZWILLER PARAMETERS -+!
     if(Norb.eq.2) then
        select case(i)
        case(1)
           f=0.d0
           iorb=1
           do ispin=1,2
              istate=index(ispin,iorb)
              f=f+gz_local_diag(phi_,dens(istate,:,:))
           end do
           f=f-niorb(iorb)
        case(2)
           f=0.d0
           iorb=2
           do ispin=1,2
              istate=index(ispin,iorb)
              f=f+gz_local_diag(phi_,dens(istate,:,:))
           end do
           f=f-niorb(iorb)
        case(3)
           f=0.d0
           do ifock=1,nFock
              f=f+phi_(ifock)*phi_(ifock)
           end do
           f=f-1.d0
        end select
     else
        select case(i)
        case(1)
           f=0.d0
           iorb=1
           do ispin=1,2
              istate=index(ispin,iorb)
              f=f+gz_local_diag(phi_,dens(istate,:,:))
           end do
           f=f-niorb(iorb)
        case(2)
           f=0.d0
           do ifock=1,nFock
              f=f+phi_(ifock)*phi_(ifock)
           end do
           f=f-1.d0
        end select
     end if
  end if

end subroutine energy_gz_projectors_iteration







subroutine fun ( x, f, i )
  !.............................................................................                                                                               
  real(8), intent( in )   :: x( : )
  real(8), intent( out )  :: f
  integer, intent( in ), optional   :: i

  if ( .not. present( i ) ) then
     !       the objective function value (user defined)                                                                                                       
     !==============================================================================                                                                           
     f = 100.0d0*(x(2)-x(1)**2)**2 +(1.0d0-x(1))**2                      !                                                                                     
     !==============================================================================                                                                           
  else
     select case ( i )
     case ( 1 )
        !               the equality constraint value (user defined)                                                                                           
	!==============================================================================                                                                        
        f = x(1)+3.0d0*x(2)-3.0d0                                   !                                                                                          
        !==============================================================================                                                                        
     case ( 2 )
	!               the inequality constraint value (user defined)                                                                                         
        !==============================================================================                                                                        
        f = x(1)**2+x(2)**2-4.0d0                                    !                                                                                         
        !==============================================================================                                                                        
     end select
  end if
  return
end subroutine fun


! function Hstar_energy(phi,local_dens) result(Estar)
!     real(8),dimension(nFock)     :: phi
!     real(8),dimension(Norb)      :: mstar,ntest,local_dens
!     real(8)                      :: Estar
!     real(8),dimension(state_dim) :: Rhop,Z,local_dens_
!     integer                      :: iorb,ispin,ie,istate,info
!     real(8)                      :: eps
!     do iorb=1,Norb
!        do ispin=1,2
!           istate=index(ispin,iorb)
!           local_dens_(istate)=local_dens(iorb)*0.5d0
!        end do
!     end do
!     Rhop=gz_Rhop_dens(phi,local_dens_,CC,CA)        
!     Z=Rhop**2.d0
!     do iorb=1,Norb
!        do ispin=1,2
!           istate=index(ispin,iorb)
!           mstar(iorb)=(local_dens(iorb)-1.d0)*Wband
!        end do
!     end do
!     Estar=0.d0
!     do iorb=1,Norb
!        do ispin=1,2
!           istate=index(ispin,iorb)
!           Estar=Estar+0.25d0*Z(istate)*Wband*(mstar(iorb)**2.d0-Wband**2.d0)
!        end do
!     end do
!   end function Hstar_energy


!+- output =====> GZ_energy as a function of the local density matrix <===== 
!+- SIMPLE TESTING FUNCTIONS -+!
! if(Norb.eq.2) then
!    GZ_energy=0.d0
!    do iorb=1,Norb
!       GZ_energy = 3.24d0 + (n0(1)-0.432d0)**2.d0 + (n0(2)-0.2d0)**2.d0
!    end do
! else
!    tmp=n0(1)
!    if(n0(1).lt.0.d0) tmp=0.d0
!    if(n0(1).gt.1.d0) tmp=1.d0
!    GZ_energy=tmp*(tmp-1.d0)
! end if
! write(77,*) GZ_energy
