subroutine enforce_su2_rotational_symmetry  
  complex(8),dimension(nFock,nFock,3) :: Svec
  complex(8),dimension(2,2,3) :: sigma_pauli

  complex(8),dimension(nFock,nFock) :: Splus,Sminus
  complex(8),dimension(nFock,nFock) :: S2,tmp_matrix,test,test_
  real(8),dimension(nFock,nFock) :: S2diag,Ntest
  real(8),dimension(nFock) :: tmp
  real(8),dimension(nFock) :: S2eigen,SZeigen
  real(8),dimension(nFock) :: Svalue,MSvalue
  real(8),dimension(nFock,2) :: spin_state
  real(8),dimension(2) :: tmp_spin_state
  integer,dimension(nFock) :: MS,search_index
  real(8),dimension(:),allocatable :: get_Svalue
  integer,dimension(:),allocatable :: count_NSstates,tmp_vec
  integer,dimension(:,:),allocatable :: irreducible_states,tmp_irreducible_states,SZ_states
  integer :: i,j,k,iorb,jorb,ispin,jspin,istate,jstate,is,iss,jfock,ifock
  real(8) :: storeS,tmp_sz,deltaS

  integer :: map,NS,NMS

  sigma_pauli=0.d0
  !
  sigma_pauli(1,2,1) = 1.d0!one
  sigma_pauli(2,1,1) = 1.d0!one
  !
  sigma_pauli(1,2,2) = -xi
  sigma_pauli(2,1,2) =  xi
  !
  sigma_pauli(1,1,3) = 1.d0!one
  sigma_pauli(2,2,3) = -1.d0!*one
  !
  S2=0.d0
  do i=1,3
     Svec(:,:,i)=0.d0
     do iorb=1,Norb
        do ispin=1,2
           do jspin=1,2
              istate=index(ispin,iorb)
              jstate=index(jspin,iorb)
              Svec(:,:,i) = Svec(:,:,i) + &
                   0.5d0*sigma_pauli(ispin,jspin,i)*matmul(CC(istate,:,:),CA(jstate,:,:))
           end do
        end do
     end do
     S2 = S2 + matmul(Svec(:,:,i),Svec(:,:,i))
  end do

  Ntest=0.d0
  do iorb=1,Norb
     do ispin=1,2
        istate=index(ispin,iorb)
        Ntest = Ntest + matmul(CC(istate,:,:),CA(istate,:,:))
     end do
  end do
  


  Splus  = Svec(:,:,1)+xi*Svec(:,:,2)
  Sminus = Svec(:,:,1)-xi*Svec(:,:,2)

  tmp_matrix = S2
  !call matrix_diagonalize(tmp_matrix,S2eigen,'V','U')
  call matrix_diagonalize(tmp_matrix,S2eigen)
  S2diag=0.d0
  do i=1,nFock
     S2diag(i,i)=S2eigen(i)
  end do
     
  !

  !<TEST
  do ifock=1,nFock
     do jfock=1,nFock
        test(ifock,jfock)=0.d0        
        test_(ifock,jfock)=0.d0        
        !
        do i=1,nFock
           do j=1,nFock
              test(ifock,jfock) = test(ifock,jfock) +&
                                !     tmp_matrix(ifock,i)*conjg(tmp_matrix(jfock,j))*Svec(i,j,3)
                   tmp_matrix(j,jfock)*conjg(tmp_matrix(i,ifock))*Ntest(i,j)
           end do
        end do
        ! do i=1,nFock
        !    test(ifock,jfock) = test(ifock,jfock) + &
        !         S2(ifock,i)*Svec(i,jfock,3)

        !    test_(ifock,jfock) = test_(ifock,jfock) + &
        !         Svec(ifock,i,3)*S2(i,jfock)
        ! end do        
        if(abs(test(ifock,jfock)-test_(ifock,jfock)).gt.1.d-10) then
           write(*,'(2I4,4F18.10)')  ifock,jfock,test(ifock,jfock),test_(ifock,jfock)
        end if

     end do
  end do
  !TEST>

  stop

  call matrix_diagonalize(test,SZeigen)

  do i=1,nFock
     write(*,*) SZeigen(i)
  end do


  do ifock=1,nFock
     do jfock=1,nFock
        test_(ifock,jfock)=0.d0        
        !
        do i=1,nFock
           do j=1,nFock
              test_(ifock,jfock) = test_(ifock,jfock) +&
                   test(i,ifock)*conjg(test(j,jfock))*S2diag(i,j)
           end do
        end do
        if(abs(test_(ifock,jfock)).gt.1.d-10) then
           write(*,'(2I4,4F18.10)')  ifock,jfock,test_(ifock,jfock)
        end if

     end do
  end do




  stop


  do i=1,nFock     
     Svalue(i)=-0.5d0+0.5d0*sqrt(1+4.d0*S2eigen(i))
     MS(i)=int(2.d0*Svalue(i)+1.d0)
     write(*,*) Svalue(i)
  end do
  !
  do k=1,nfock
     !
     do i=1,nFock
        tmp(i)=0.d0
        do j=1,nfock
           tmp(i) = tmp(i) + dreal(Svec(i,j,3))*tmp_matrix(j,k)
        end do
     end do

     write(*,'(60F5.2)') tmp_matrix(:,k)     
     write(*,'(60F5.2)') tmp(:)

     !
     MSvalue(k) = 0.d0
     do j=1,nFock
        MSvalue(k) = MSvalue(k) + tmp(j)*tmp_matrix(j,k)
     end do
     !
     spin_state(k,1) = Svalue(k)
     spin_state(k,2) = MSvalue(k)
     !
     write(*,*) Svalue(k),MSvalue(k)
     !
  end do
  stop
  !  

  ! group states according to S
  NS = get_mask_equal_values(Svalue,irreducible_states)

  ! group states according to MS
  NMS = get_mask_equal_values(MSvalue,SZ_states)

  write(*,*)
  write(*,*) 'S masks'
  write(*,*)

  do is=1,NS
     write(*,'(20I3)') irreducible_states(is,:)
     do ifock=1,nFock
        if(irreducible_states(is,ifock).eq.1) then           
           write(*,*) Svalue(ifock)           
        end if
     end do
  end do

  write(*,*)
  write(*,*) 'Sz masks'
  write(*,*)

  do is=1,NMS
     write(*,'(20I3)') Sz_states(is,:)
     do ifock=1,nFock
        if(Sz_states(is,ifock).eq.1) then           
           write(*,*) MSvalue(ifock)           
        end if
     end do
  end do




  stop




  ! !
  ! NS=1
  ! do jfock=1,Nfock
  !    tmp_spin_state=spin_state(jfock,:)
  !    do ifock=1,nFock     
  !       !
  !       deltaS = vec_dist(spin_state(ifock,:),tmp_spin_state)
  !       if(deltaS.lt.1.d-8) then

  !       end if
  !       !
  !    end do
  ! end do





  storeS=Svalue(1)
  NS=1
  do i=1,Nfock
     if(abs(Svalue(i)-storeS).gt.1.d-6) then
        storeS = Svalue(i)
        NS = NS + 1
     end if
  end do
  write(*,*) NS
  !
  allocate(count_NSstates(NS),get_Svalue(NS))
  count_NSstates=0
  is=1;  storeS=Svalue(1)
  get_Svalue(1) = storeS
  do i=1,Nfock
     if(abs(Svalue(i)-storeS).gt.1.d-6) then
        storeS = Svalue(i)
        is = is + 1
        count_NSstates(is) = 1        
        get_Svalue(is) = storeS
     else
        count_NSstates(is)=count_NSstates(is) + 1
     end if
  end do
  write(*,*) count_NSstates(:)

  do is=1,NS     
     !+- find kernel of subspace S -+!
     allocate(tmp_vec(count_NSstates(is)))
     !
     tmp_vec=0
     iss = 0
     do i=1,nFock
        if(abs(Svalue(i) - get_Svalue(is)).lt.1.d-6) then
           iss = iss + 1
           tmp_vec(iss) = i
        end if
     end do
     !
     !write(*,*) tmp_vec(:)
     !the vector tmp_vec contains 

     deallocate(tmp_vec)
  end do


  !stop
  !
  !<CHECK
  ! do i=1,nFock
  !    write(*,*) S2eigen(i)
  !    write(*,*)
  !    write(*,'(20F7.3)') S2(:,i)
  !    write(*,*)
  ! end do
  !CHECK>

  !+- FIND THE | Gamma; (S,Sz) > eigenstates

  !+- group the states with the same S and Sz

  !+- build the new transformation matrix

  !+- check the diagonalization of the spin operator  


end subroutine enforce_su2_rotational_symmetry


function vec_dist(a,b) result(r)
  real(8),dimension(:) :: a,b
  real(8) :: r
  integer :: Na,Nb,i
  Na=size(a);Nb=size(b)
  if(Na.ne.Nb) then 
     write(*,*) 'vec_dist_error'; stop
  end if
  r=0.d0
  do i=1,Na
     r = r+ (a(i)-b(i))**2.d0
  end do
  r=sqrt(r)
end function vec_dist




function get_mask_equal_values(xValues,mask,eps) result(Ns)
  !
  real(8),dimension(:)               :: xValues
  integer,dimension(:,:),allocatable :: mask
  real(8),optional                   :: eps
  integer                            :: Ns
  real(8)                            :: eps_
  real(8),dimension(:,:),allocatable :: tmp_mask
  integer,dimension(:),allocatable   :: search_index
  integer                            :: Nx,ix,jx,is
  real(8) :: deltaX
  !
  eps_=1.d-8
  if(present(eps)) eps_=eps
  if(allocated(mask)) deallocate(mask)
  Nx=size(xValues)
  allocate(tmp_mask(Nx,Nx),search_index(Nx))  
  tmp_mask=0
  search_index=0
  NS=0
  do ix=1,Nx
     if(search_index(ix).ge.0) then
        NS = NS + 1
        do jx=1,Nx
           deltaX = abs(Xvalues(jx)-Xvalues(ix))
           if(deltaX.lt.eps_) then
              search_index(jx) = -1
              tmp_mask(NS,jx) = 1
           end if
        end do
     end if
  end do
  allocate(mask(NS,Nx))
  do is=1,NS
     mask(is,:)=tmp_mask(is,:)
  end do
  deallocate(tmp_mask)
  !
end function get_mask_equal_values






