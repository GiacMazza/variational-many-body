MODULE GZ_MINIMIZE
  USE SCIFOR
  !USE SF_OPTIMIZE
  USE GZ_AUX_FUNX
  USE GZ_VARS_GLOBAL
  USE GZ_PROJECTORS
  USe GZ_ENERGY_FUNCTIONAL
  USE MIN_AMOEBA
  implicit none
  private

  public ::  gz_optimization_VS_density,gz_optimization_VS_polarization


CONTAINS




  function gz_optimization_VS_density(density_seed) result(optimized_energy)
    real(8),dimension(2),intent(inout) :: density_seed
    real(8)                            :: optimized_energy    
    !+- amoeba_variables-+!
    real(8),allocatable,dimension(:,:) :: p
    real(8),allocatable,dimension(:)   :: y
    real(8)                            :: ftol
    integer                            :: np,mp,i_vertex,j_vertex,i_dim,iter
    integer,allocatable,dimension(:)   :: idum
    integer                            :: tmp_dum
    real(8)                            :: rnd,tot_dens,tmp_pol
    integer                            :: amoeba_unit    
    !+-------------------+!

    iter=100

    call random_seed(size=tmp_dum)
    allocate(idum(tmp_dum))
    idum=1234567
    call random_seed(put=[idum])
    
    amoeba_unit=free_unit()
    open(amoeba_unit,file='Amoeba_call.out',access='append')
    !+- AMOEBA MINIMIZATION WITH RESPECT TO LOCAL DENSITY n0=(n1+n2) -+!
    NP=2
    MP=NP+1  
    allocate(y(MP),p(MP,NP))
    !+- initialize simplex -+!
    do i_dim = 1,NP
       p(:,i_dim) = density_seed(i_dim)
    end do
    do i_vertex=2,MP
       if(.not.fix_density_minimization) then
          !+- with this choise each vertex of the simplex has different number of particles -+!
          tmp_pol = density_seed(2)-density_seed(1)
          do j_vertex=1,NP
             call random_number(rnd)
             rnd=rnd*0.2
             if(tmp_pol.gt.1.d0) then
                p(i_vertex,j_vertex) = density_seed(j_vertex)-rnd*(-1.d0)**dble(j_vertex)
             else
                p(i_vertex,j_vertex) = density_seed(j_vertex)+rnd*(-1.d0)**dble(j_vertex)
             end if
          end do
       else
          !+- this choise conserve the total number of particles -+!
          call random_number(rnd)
          rnd=rnd*0.2d0
          tmp_pol = density_seed(2)-density_seed(1)
          do j_vertex=1,NP
             if(tmp_pol.gt.1.d0) then
                !+- decrease a bit polarization -+!
                p(i_vertex,j_vertex) = density_seed(j_vertex)-rnd*(-1.d0)**dble(j_vertex)
             else
                !+- increase a bit polarization -+!
                p(i_vertex,j_vertex) = density_seed(j_vertex)+rnd*(-1.d0)**dble(j_vertex)
             end if

          end do
       end if
    end do
    write(amoeba_unit,*) 'Initializing simplex'
    !+----------------------+!  
    do i_vertex=1,MP
       y(i_vertex)=gz_energy_local_density(p(i_vertex,:))
       write(amoeba_unit,*) 'vertex',p(i_vertex,:),'gz_energy_functional',y(i_vertex)
    end do
    ftol=1.d-6
    call amoeba(p(1:MP,1:NP),y(1:MP),FTOL,gz_energy_local_density,iter,amoeba_verbose)       
    density_seed=p(1,:)
    optimized_energy=y(1)     
    write(amoeba_unit,*) 
    write(amoeba_unit,*) 'Result'
    write(amoeba_unit,*) 
    do i_vertex=1,MP
       write(amoeba_unit,*) p(i_vertex,:),y(i_vertex)
    end do
    write(amoeba_unit,*) 
    write(amoeba_unit,*) 
    deallocate(y,p)
    close(amoeba_unit)       
  end function gz_optimization_VS_density





  
  function gz_optimization_VS_density_self(density_seed) result(optimized_energy)
    real(8),dimension(2),intent(inout) :: density_seed
    real(8)                            :: optimized_energy    
    !+- amoeba_variables-+!
    real(8),allocatable,dimension(:,:) :: p
    real(8),allocatable,dimension(:)   :: y
    real(8)                            :: ftol
    integer                            :: np,mp,i_vertex,j_vertex,i_dim,iter
    integer,allocatable,dimension(:)   :: idum
    integer                            :: tmp_dum
    real(8)                            :: rnd,tot_dens,tmp_pol
    integer                            :: amoeba_unit    
    !+-------------------+!

    iter=100

    call random_seed(size=tmp_dum)
    allocate(idum(tmp_dum))
    idum=1234567
    call random_seed(put=[idum])
    
    amoeba_unit=free_unit()
    open(amoeba_unit,file='Amoeba_call.out',access='append')
    !+- AMOEBA MINIMIZATION WITH RESPECT TO LOCAL DENSITY n0=(n1+n2) -+!
    NP=2
    MP=NP+1  
    allocate(y(MP),p(MP,NP))
    !+- initialize simplex -+!
    do i_dim = 1,NP
       p(:,i_dim) = density_seed(i_dim)
    end do
    do i_vertex=2,MP
       if(.not.fix_density_minimization) then
          !+- with this choise each vertex of the simplex has different number of particles -+!
          tmp_pol = density_seed(2)-density_seed(1)
          do j_vertex=1,NP
             call random_number(rnd)
             rnd=rnd*0.2
             if(tmp_pol.gt.1.d0) then
                p(i_vertex,j_vertex) = density_seed(j_vertex)-rnd*(-1.d0)**dble(j_vertex)
             else
                p(i_vertex,j_vertex) = density_seed(j_vertex)+rnd*(-1.d0)**dble(j_vertex)
             end if
          end do
       else
          !+- this choise conserve the total number of particles -+!
          call random_number(rnd)
          rnd=rnd*0.2d0
          tmp_pol = density_seed(2)-density_seed(1)
          do j_vertex=1,NP
             if(tmp_pol.gt.1.d0) then
                !+- decrease a bit polarization -+!
                p(i_vertex,j_vertex) = density_seed(j_vertex)-rnd*(-1.d0)**dble(j_vertex)
             else
                !+- increase a bit polarization -+!
                p(i_vertex,j_vertex) = density_seed(j_vertex)+rnd*(-1.d0)**dble(j_vertex)
             end if

          end do
       end if
    end do
    write(amoeba_unit,*) 'Initializing simplex'
    !+----------------------+!  
    do i_vertex=1,MP
       y(i_vertex)=gz_energy_local_density(p(i_vertex,:))
       write(amoeba_unit,*) 'vertex',p(i_vertex,:),'gz_energy_functional',y(i_vertex)
    end do
    ftol=1.d-6
    call amoeba(p(1:MP,1:NP),y(1:MP),FTOL,gz_energy_local_density,iter,amoeba_verbose)       
    density_seed=p(1,:)
    optimized_energy=y(1)     
    write(amoeba_unit,*) 
    write(amoeba_unit,*) 'Result'
    write(amoeba_unit,*) 
    do i_vertex=1,MP
       write(amoeba_unit,*) p(i_vertex,:),y(i_vertex)
    end do
    write(amoeba_unit,*) 
    write(amoeba_unit,*) 
    deallocate(y,p)
    close(amoeba_unit)       
  end function gz_optimization_VS_density_self
  
  function gz_optimization_VS_polarization(polarization_seed) result(optimized_energy)
    real(8),intent(inout)              :: polarization_seed
    real(8)                            :: optimized_energy    
    !+- amoeba_variables-+!
    real(8),allocatable,dimension(:,:) :: p
    real(8),allocatable,dimension(:)   :: y
    real(8)                            :: ftol
    integer                            :: np,mp,i_vertex,j_vertex,i_dim,iter
    integer,allocatable,dimension(:)   :: idum
    integer                            :: tmp_dum

    real(8)                            :: rnd
    integer                            :: amoeba_unit    
    !+-------------------+!
    iter=100
    call random_seed(size=tmp_dum)
    allocate(idum(tmp_dum))
    idum=1234567
    call random_seed(put=[idum])
    
    amoeba_unit=free_unit()
    open(amoeba_unit,file='Amoeba_call.out')
    !+- AMOEBA MINIMIZATION WITH RESPECT TO LOCAL DENSITY n0=(n1+n2) -+!
    NP=1
    MP=NP+1  
    allocate(y(MP),p(MP,NP))
    !+- initialize simplex -+!
    do i_dim = 1,NP
       p(:,i_dim) = polarization_seed
    end do
    do i_vertex=2,MP
       do j_vertex=1,NP
          call random_number(rnd)
          p(i_vertex,j_vertex) = dble(rnd)
       end do
    end do
    write(amoeba_unit,*) 'Initializing simplex'
    !+----------------------+!  
    do i_vertex=1,MP
       y(i_vertex)=gz_energy_local_polar(p(i_vertex,:))
       write(amoeba_unit,*) 'vertex',p(i_vertex,:),'gz_energy_functional',y(i_vertex)
    end do
    ftol=1.d-6
    call amoeba(p(1:MP,1:NP),y(1:MP),FTOL,gz_energy_local_polar,iter)       
    polarization_seed=p(1,1)
    optimized_energy=y(1)     
    write(amoeba_unit,*) 
    write(amoeba_unit,*) 'Result'
    write(amoeba_unit,*) 
    do i_vertex=1,MP
       write(amoeba_unit,*) p(i_vertex,:),y(i_vertex)
    end do
    write(amoeba_unit,*) 
    write(amoeba_unit,*) 
    deallocate(y,p)
    close(amoeba_unit)       
  end function gz_optimization_VS_polarization

  
  
END MODULE GZ_MINIMIZE
  
