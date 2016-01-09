subroutine enforce_su2_rotational_symmetry  
  !+- SPIN AND ORBITAL ANGULAR MOMENTUM OPERATOR -+!
  complex(8),dimension(nFock,nFock,3)     :: Svec
  complex(8),dimension(nFock,nFock)       :: S2
  !
  complex(8),dimension(nFock,nFock,3)     :: isoSvec
  complex(8),dimension(nFock,nFock)       :: isoS2
  !
  complex(8),dimension(2,2,3)             :: sigma_pauli
  complex(8),dimension(3,3,3)             :: levi_civita
  complex(8),dimension(nFock,nFock)       :: Splus,Sminus
  complex(8),dimension(nFock,nFock)       :: tmp_matrix,test,test_
  !
  complex(8),dimension(:,:,:),allocatable :: test_joint_diag
  complex(8),dimension(:,:),allocatable   :: test_jointV
  real(8),dimension(:,:),allocatable      :: test_jointD
  !
  real(8),dimension(nFock,nFock)          :: S2diag,Ntest
  real(8),dimension(nFock)                :: tmp
  real(8),dimension(nFock)                :: S2eigen,SZeigen
  real(8),dimension(nFock)                :: Svalue,MSvalue
  real(8),dimension(nFock,2)              :: spin_state
  real(8),dimension(2)                    :: tmp_spin_state
  integer,dimension(nFock)                :: MS,search_index
  real(8),dimension(:),allocatable        :: get_Svalue
  integer,dimension(:),allocatable        :: count_NSstates,tmp_vec
  integer,dimension(:,:),allocatable      :: irreducible_states,tmp_irreducible_states,SZ_states
  integer                                 :: i,j,k,iorb,jorb,ispin,jspin,istate,jstate,is,iss,jfock,ifock
  real(8)                                 :: storeS,tmp_sz,deltaS,tmp_test

  integer                                 :: map,NS,NMS


  !+- build sigma pauli and levi-civita tensor -+!
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
  !
  !
  levi_civita=0.d0
  levi_civita(1,2,3) =  1.d0
  levi_civita(1,3,2) = -1.d0
  !
  levi_civita(2,3,1) =  1.d0
  levi_civita(3,2,1) = -1.d0
  !
  levi_civita(3,1,2) =  1.d0
  levi_civita(3,2,1) = -1.d0
  !

  S2=zero
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

  isoS2=0.d0
  do i=1,3
     isoSvec(:,:,i)=0.d0
     select case(Norb)
     case(1)        
        forall(ifock=1:nFock) isoSvec(ifock,ifock,i) = 1.d0
     case(2)
        do iorb=1,Norb
           do jorb=1,Norb
              do ispin=1,2
                 istate=index(ispin,iorb)
                 jstate=index(ispin,jorb)
                 isoSvec(:,:,i) = isoSvec(:,:,i) + &
                      0.5d0*sigma_pauli(iorb,jorb,i)*matmul(CC(istate,:,:),CA(jstate,:,:))
              end do
           end do
        end do
     case(3) 
        do iorb=1,Norb
           do jorb=1,Norb
              do ispin=1,2
                 istate=index(ispin,iorb)
                 jstate=index(ispin,jorb)
                 isoSvec(:,:,i) = isoSvec(:,:,i) + &
                      xi*levi_civita(i,iorb,jorb)*matmul(CC(istate,:,:),CA(jstate,:,:))
              end do
           end do
        end do
     end select
     isoS2 = isoS2 + matmul(isoSvec(:,:,i),isoSvec(:,:,i))
  end do

  !check isoS2 commutes with S2
  test=matmul(S2,isoS2)
  test_=matmul(isoS2,S2)
  do ifock=1,nFock
     do jfock=1,nFock
        if(abs(test(ifock,jfock)-test_(ifock,jfock)).gt.1.d-8) write(*,*) ifock,jfock
     end do
  end do
  !stop

  !< TEST JACOBI JOINT DIAGONALIZATION 
  tmp_matrix = Svec(:,:,2)
  !call matrix_diagonalize(tmp_matrix,S2eigen,'V','U')
  call matrix_diagonalize(tmp_matrix,S2eigen)
  S2diag=0.d0
  do i=1,nFock
     S2diag(i,i)=S2eigen(i)
     write(*,*) i,S2eigen(i),Ntest(i,i)
  end do



  !  allocate(test_joint_diag(nFock,nFock,1),test_jointV(nFock,nFock),test_jointD(nFock,nFock))

  allocate(test_joint_diag(nFock,nFock,5),test_jointV(nFock,nFock),test_jointD(nFock,5))
  test_joint_diag(:,:,1)=S2
  test_joint_diag(:,:,2)=Svec(:,:,3)
  test_joint_diag(:,:,3)=isoS2
  test_joint_diag(:,:,4)=isoSvec(:,:,3)
  test_joint_diag(:,:,5)=Ntest

  write(*,*)
  call simultaneous_diag(test_joint_diag,test_jointV,test_jointD,eps=1.d-10) 

  do ifock=1,nFock
     do jfock=ifock+1,nFock
        !if(abs(test_jointD(ifock,jfock)).gt.1.d-10) stop 'diagonalization failed'
     end do
     write(778,*) ifock,test_jointD(ifock,:)
  end do

  ! stop
  !  TEST JACOBI JOINT DIAGONALIZATION>




  Splus  = Svec(:,:,1)+xi*Svec(:,:,2)
  Sminus = Svec(:,:,1)-xi*Svec(:,:,2)

  stop
  !


  !+- FIND THE | Gamma; (S,Sz) > eigenstates

  !+- group the states with the same S and Sz

  !+- build the new transformation matrix

  !+- check the diagonalization of the spin operator  


end subroutine enforce_su2_rotational_symmetry



subroutine simultaneous_diag(A,V,diagA,eps)
  implicit none
  complex(8),dimension(:,:,:) :: A
  complex(8),dimension(:,:)   :: V
  complex(8),dimension(:,:),allocatable   :: Rgivens
  complex(8),dimension(:,:,:),allocatable   :: tmp
  real(8),dimension(:,:)      :: diagA
  real(8),dimension(:,:),allocatable      :: tmp_diag
  real(8),dimension(3,3)      :: G
  complex(8),dimension(3)      :: Gvec
  complex(8),dimension(2,2)      :: tmp_rotate,tmpA,tmp_rotate_dag
  real(8),dimension(3) :: Geigen,Jacobi_vec
  real(8)                     :: eps,C,off
  complex(8)                  :: S
  integer                     :: N,M,NM
  integer                     :: iloop,imax,i,j,k
  !
  real(8) :: min,max
  integer,dimension(:),allocatable :: sort,tmp_sort,ivec
  integer :: imin
  !
  integer                     :: ii,jj,kk,kkk

  integer,dimension(:,:),allocatable :: maskN,maskS,maskisoS
  integer,dimension(:,:),allocatable :: tmp_group_states,group_states
  !integer,dimension(:,:,:),allocatable ::  cross_mask
  integer,dimension(:,:,:,:),allocatable ::  cross_mask
  integer :: count_Sstates,count_Nstates,count_isoSstates,count,test_count

  imax=200
  !
  if(size(A,1).ne.size(A,2)) stop "simultaneous diag: wrong dimensions" 
  N=size(A,1);M=size(A,3);NM=N*M
  if(size(V,1).ne.size(V,2)) stop "simultaneous diag: wrong dimensions" 
  if(size(V,1).ne.N) stop "simultaneous diag: wrong dimensions" 

  if(size(diagA,1).ne.N) stop "simultaneous diag: wrong dimensions"
  if(size(diagA,2).ne.M) stop "simultaneous diag: wrong dimensions"
  !

  write(*,*) N,M,NM
  !  stop

  V=0.d0
  do i=1,N
     V(i,i) = 1.d0
  end do

  allocate(tmp(N,N,M),Rgivens(N,N));

  imax=20
  do iloop=1,imax

     off = 0.d0
     do k=1,M
        off = off + off_diag(A(:,:,k))
     end do
     write(777,*) iloop,off


     do i=1,N
        do j=1,N           
           if(i.ne.j) then
              !
              G=0.d0     
              do k=1,M
                 !
                 Gvec(1) = A(i,i,k) - A(j,j,k)
                 Gvec(2) = A(i,j,k) + A(j,i,k)
                 Gvec(3) = xi*(A(j,i,k) - A(i,j,k))
                 G = G + get_gmatrix(Gvec)
                 !
              end do
              !              
              call matrix_diagonalize(G,Geigen,'V','U')
              Jacobi_vec = G(:,3)
              !
              !Givens rotations
              if(Jacobi_vec(1).lt.0.d0) Jacobi_vec=-Jacobi_vec
              C = sqrt(Jacobi_vec(1)*0.5d0+0.5d0)
              S = Jacobi_vec(2)-xi*Jacobi_vec(3)
              S = 0.5d0*S/C
              !
              Rgivens=0.d0
              do k=1,N
                 Rgivens(k,k) = 1.d0
              end do
              !+- R A R^+ rotation -+!
              Rgivens(i,i) = C
              Rgivens(i,j) = conjg(S)
              Rgivens(j,i) = -S
              Rgivens(j,j) = C
              !
              do k=1,M              
                 !+-  compute rotations R A R^+  -+!
                 tmp(:,:,k)=A(:,:,k)              
                 do kk=1,N
                    tmp(kk,i,k) = A(kk,i,k)*conjg(Rgivens(i,i)) + A(kk,j,k)*conjg(Rgivens(i,j))
                    tmp(KK,j,k) = A(kk,i,k)*conjg(Rgivens(j,i)) + A(kk,j,k)*conjg(Rgivens(j,j))
                 end do
                 A(:,:,k)=tmp(:,:,k)
                 do kk=1,N
                    A(i,kk,k) =  Rgivens(i,i)*tmp(i,kk,k) + Rgivens(i,j)*tmp(j,kk,k)
                    A(j,kk,k) =  Rgivens(j,i)*tmp(i,kk,k) + Rgivens(j,j)*tmp(j,kk,k)
                 end do
                 !
              end do
              !
              V(i,i) = C*V(i,i)+conjg(S)*V(j,i)
              V(i,j) = C*V(i,j)+conjg(S)*V(j,j)
              V(j,i) = -S*V(i,i)+C*V(j,i)
              V(j,j) = -S*V(i,j)+C*V(j,j)
              !
           end if
        end do
     end do
  end do
  !
  do k=1,M
     do i=1,N
        diagA(i,k)=A(i,i,k)
     end do
  end do
  !
  do i=1,N
     write(*,'(20F7.3)') diagA(i,:)
  end do
  !+- sort eigenvalues for k=1 -+!
  max=diagA(1,1)
  min=diagA(N,1)
  imax=1
  imin=N  
  allocate(sort(N),tmp_sort(N),tmp_diag(N,2)); forall(i=1:N) sort(i)=0; tmp_sort=sort  
  do i=1,N
     max = -10     
     do j=1,N
        if(tmp_sort(j).eq.0.and.diagA(j,1).gt.max) then
           max = diagA(j,1)
           imax = j
        end if
     end do
     tmp_sort(imax) = 1
     sort(i) = imax
     write(*,'(I3,A,20I3)') i-1,'   ',sort(:)
  end do
  !
  tmp_diag=diagA
  do i=1,N
     diagA(i,:) = tmp_diag(sort(i),:)
     write(*,'(I4,20F8.3)') sort(i),diagA(i,1),diagA(i,3),diagA(i,5)
  end do

  !+- match and group (N,S) states -+!
  count_Nstates=get_mask_equal_values(diagA(:,5),maskN)  
  write(*,*) count_Nstates
  do i=1,N
     write(*,*) maskN(:,i)
  end do
  !
  count_Sstates=get_mask_equal_values(diagA(:,1),maskS)  
  write(*,*) count_Sstates
  do i=1,N
     write(*,*) maskS(:,i)
  end do
  !
  count_isoSstates=get_mask_equal_values(diagA(:,3),maskisoS)  
  write(*,*) count_isoSstates
  do i=1,N
     write(*,*) maskisoS(:,i)
  end do  
  !
  allocate(cross_mask(count_Nstates,count_Sstates,count_isoSstates,N))
  allocate(tmp_group_states(3,N))
  count=0
  write(*,*)
  do i=1,count_Nstates
     do j=1,count_Sstates
        do ii=1,count_isoSstates
           test_count=0
           do k=1,N
              cross_mask(i,j,ii,k) = maskN(i,k)*maskS(j,k)*maskisoS(ii,k)
              test_count = test_count + cross_mask(i,j,ii,k)
           end do
           if(test_count.gt.0) then
              count = count+1  
              tmp_group_states(1,count) = i; tmp_group_states(2,count) = j           
              tmp_group_states(3,count) = ii
              write(*,*) i,j,ii,test_count
           end if
        end do
     end do
  end do
  write(*,*) count
  !
  allocate(group_states(3,count))
  do i=1,count
     group_states(:,i) = tmp_group_states(:,i)
  end do
  !
  do i=1,count
     do k=1,N
        if(cross_mask(group_states(1,i),group_states(2,i),group_states(3,i),k).eq.1) then
!           write(*,'(5F7.3,10I3)') diagA(k,1),diagA(k,2),diagA(k,5)
           write(*,'(5F7.3,10I3)') diagA(k,:)
        end if
     end do
  end do



!   allocate(cross_mask(count_Nstates,count_Sstates,N))
!   allocate(tmp_group_states(2,N))
!   count=0
!   write(*,*)
!   do i=1,count_Nstates
!      do j=1,count_Sstates
!         test_count=0
!         do k=1,N
!            cross_mask(i,j,k) = maskN(i,k)*maskS(j,k)
!            test_count = test_count + cross_mask(i,j,k)
!         end do
!         if(test_count.gt.0) then
!            count = count+1  
!            tmp_group_states(1,count) = i; tmp_group_states(2,count) = j           
!            write(*,*) i,j,test_count
!         end if
!      end do
!   end do
!   write(*,*) count
!   !
!   allocate(group_states(2,count))
!   do i=1,count
!      group_states(:,i) = tmp_group_states(:,i)
!   end do
!   !
!   do i=1,count
!      do k=1,N
!         if(cross_mask(group_states(1,i),group_states(2,i),k).eq.1) then
! !           write(*,'(5F7.3,10I3)') diagA(k,1),diagA(k,2),diagA(k,5)
!            write(*,'(5F7.3,10I3)') diagA(k,:)
!         end if
!      end do
!   end do


  !
  !+- TO DO LIST -+!
  !+- match and group the equal entries of diagA
  !+- obtain map
  !+- sort the corresponding diagA according to one of the column values
  !  !
contains
  !
  function get_Gmatrix(h) result(G)
    implicit none
    complex(8),dimension(:) :: h
    real(8),dimension(size(h),size(h)) :: G
    integer :: N,i,j
    N=size(h)
    ! if(size(G,1).ne.size(G,2)) stop 'wrong dimensions'
    ! if(size(G,1).ne.N) stop 'wrong dimensions'
    !
    do i=1,N
       do j=1,N
          G(i,j) = dreal(conjg(h(i))*h(j))
       end do
    end do
    !

  end function get_Gmatrix


  !
  function off_diag(A) result(off)
    complex(8),dimension(:,:) :: A
    real(8) :: off
    integer :: i,j,N
    N=size(A,1)
    !
    off = 0.d0
    do i=1,N
       do j=1,N
          if(i.ne.j) off = off + abs(A(i,j))**2.d0
       end do
    end do
    !
  end function off_diag
  !
end subroutine simultaneous_diag

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


! function get_equal_states(xValues,mask,eps) result(Ns)
!   !
!   real(8),dimension(:)               :: xValues
!   integer,dimension(:,:),allocatable :: mask
!   real(8),optional                   :: eps
!   integer                            :: Ns
!   real(8)                            :: eps_
!   real(8),dimension(:,:),allocatable :: tmp_mask
!   integer,dimension(:),allocatable   :: search_index
!   integer                            :: Nx,ix,jx,is
!   real(8) :: deltaX
!   !
!   eps_=1.d-8
!   if(present(eps)) eps_=eps
!   if(allocated(mask)) deallocate(mask)
!   Nx=size(xValues)
!   allocate(tmp_mask(Nx,Nx),search_index(Nx))  
!   tmp_mask=0
!   search_index=0
!   NS=0
!   do ix=1,Nx
!      if(search_index(ix).ge.0) then
!         NS = NS + 1
!         do jx=1,Nx
!            deltaX = abs(Xvalues(jx)-Xvalues(ix))
!            if(deltaX.lt.eps_) then
!               search_index(jx) = -1
!               tmp_mask(NS,jx) = 1
!            end if
!         end do
!      end if
!   end do
!   allocate(mask(NS,Nx))
!   do is=1,NS
!      mask(is,:)=tmp_mask(is,:)
!   end do
!   deallocate(tmp_mask)  
!   !
! end function get_equal_states













  ! !<TEST
  ! tmp=0.d0
  ! do ifock=1,nFock
  !    tmp(ifock)=0.d0
  !    do jfock=1,nFock
  !       test(ifock,jfock)=0.d0        
  !       test_(ifock,jfock)=0.d0        
  !       !
  !       do i=1,nFock
  !          do j=1,nFock
  !             test(ifock,jfock) = test(ifock,jfock) + &
  !                  tmp_matrix(j,jfock)*conjg(tmp_matrix(i,ifock))*Ntest(i,j)
  !             !tmp_matrix(j,jfock)*conjg(tmp_matrix(i,ifock))*S2(i,j)
  !          end do
  !       end do
  !       !
  !       !
  !       ! test(ifock,jfock)=0.d0        
  !       ! do i=1,nFock
  !       !    do j=1,nFock
  !       !       test(ifock,jfock) = test(ifock,jfock) + &
  !       !            tmp_matrix(j,jfock)*conjg(tmp_matrix(i,ifock))*S2(i,j)
  !       !    end do
  !       ! end do        
  !       !
  !       !
  !       if(abs(test(ifock,jfock)).gt.1.d-10) then
  !          write(*,'(2I4,4F18.10)')  ifock,jfock,test(ifock,jfock)
  !          if(ifock.ne.jfock) EXIT
  !       end if

  !       ! if(abs(test(ifock,jfock)).gt.1.d-10) then
  !       !    write(*,'(2I4,4F18.10)')  ifock,jfock,test(ifock,jfock),Ntest(ifock,jfock)
  !       !    if(ifock.ne.jfock) stop 'test ifock /= jfock'
  !       ! end if
  !       !tmp(ifock) = tmp(ifock) + Ntest(ifock,jfock)*tmp_matrix(jfock,ifock)
  !    end do
  !    !write(*,*) Ntest(:,ifock)
  !    !     write(*,*) dreal(tmp_matrix(:,ifock))
  ! end do

  ! tmp_test=0.d0
  ! do jfock=1,nFock
  !    tmp(jfock) = 0.d0
  !    tmp_test = 0.d0
  !    do ifock=1,nFock
  !       tmp(jfock) = tmp(jfock) + dreal(tmp_matrix(ifock,jfock))**2.d0
  !       write(77,*) ifock,dreal(tmp_matrix(ifock,jfock)),dreal(tmp_matrix(ifock,jfock))**2.d0
  !       tmp_test = tmp_test + tmp_matrix(ifock,jfock)*tmp_matrix(ifock,jfock)*Ntest(ifock,ifock)
  !    end do
  !    !
  !    write(77,*)
  !    write(77,*)
  !    !
  !    write(78,*) jfock,Ntest(jfock,jfock)
  !    !
  !    write(79,*) tmp(jfock)
  !    write(*,*) jfock,tmp_test
  ! end do


  ! !TEST>

  ! stop

  ! call matrix_diagonalize(test,SZeigen)

  ! do i=1,nFock
  !    write(*,*) SZeigen(i)
  ! end do


  ! do ifock=1,nFock
  !    do jfock=1,nFock
  !       test_(ifock,jfock)=0.d0        
  !       !
  !       do i=1,nFock
  !          do j=1,nFock
  !             test_(ifock,jfock) = test_(ifock,jfock) +&
  !                  test(i,ifock)*conjg(test(j,jfock))*S2diag(i,j)
  !          end do
  !       end do
  !       if(abs(test_(ifock,jfock)).gt.1.d-10) then
  !          write(*,'(2I4,4F18.10)')  ifock,jfock,test_(ifock,jfock)
  !       end if

  !    end do
  ! end do




  ! stop


  ! do i=1,nFock     
  !    Svalue(i)=-0.5d0+0.5d0*sqrt(1+4.d0*S2eigen(i))
  !    MS(i)=int(2.d0*Svalue(i)+1.d0)
  !    write(*,*) Svalue(i)
  ! end do
  ! !
  ! do k=1,nfock
  !    !
  !    do i=1,nFock
  !       tmp(i)=0.d0
  !       do j=1,nfock
  !          tmp(i) = tmp(i) + dreal(Svec(i,j,3))*tmp_matrix(j,k)
  !       end do
  !    end do

  !    write(*,'(60F5.2)') tmp_matrix(:,k)     
  !    write(*,'(60F5.2)') tmp(:)

  !    !
  !    MSvalue(k) = 0.d0
  !    do j=1,nFock
  !       MSvalue(k) = MSvalue(k) + tmp(j)*tmp_matrix(j,k)
  !    end do
  !    !
  !    spin_state(k,1) = Svalue(k)
  !    spin_state(k,2) = MSvalue(k)
  !    !
  !    write(*,*) Svalue(k),MSvalue(k)
  !    !
  ! end do
  ! stop
  ! !  

  ! ! group states according to S
  ! NS = get_mask_equal_values(Svalue,irreducible_states)

  ! ! group states according to MS
  ! NMS = get_mask_equal_values(MSvalue,SZ_states)

  ! write(*,*)
  ! write(*,*) 'S masks'
  ! write(*,*)

  ! do is=1,NS
  !    write(*,'(20I3)') irreducible_states(is,:)
  !    do ifock=1,nFock
  !       if(irreducible_states(is,ifock).eq.1) then           
  !          write(*,*) Svalue(ifock)           
  !       end if
  !    end do
  ! end do

  ! write(*,*)
  ! write(*,*) 'Sz masks'
  ! write(*,*)

  ! do is=1,NMS
  !    write(*,'(20I3)') Sz_states(is,:)
  !    do ifock=1,nFock
  !       if(Sz_states(is,ifock).eq.1) then           
  !          write(*,*) MSvalue(ifock)           
  !       end if
  !    end do
  ! end do




  ! stop




  ! ! !
  ! ! NS=1
  ! ! do jfock=1,Nfock
  ! !    tmp_spin_state=spin_state(jfock,:)
  ! !    do ifock=1,nFock     
  ! !       !
  ! !       deltaS = vec_dist(spin_state(ifock,:),tmp_spin_state)
  ! !       if(deltaS.lt.1.d-8) then

  ! !       end if
  ! !       !
  ! !    end do
  ! ! end do





  ! storeS=Svalue(1)
  ! NS=1
  ! do i=1,Nfock
  !    if(abs(Svalue(i)-storeS).gt.1.d-6) then
  !       storeS = Svalue(i)
  !       NS = NS + 1
  !    end if
  ! end do
  ! write(*,*) NS
  ! !
  ! allocate(count_NSstates(NS),get_Svalue(NS))
  ! count_NSstates=0
  ! is=1;  storeS=Svalue(1)
  ! get_Svalue(1) = storeS
  ! do i=1,Nfock
  !    if(abs(Svalue(i)-storeS).gt.1.d-6) then
  !       storeS = Svalue(i)
  !       is = is + 1
  !       count_NSstates(is) = 1        
  !       get_Svalue(is) = storeS
  !    else
  !       count_NSstates(is)=count_NSstates(is) + 1
  !    end if
  ! end do
  ! write(*,*) count_NSstates(:)

  ! do is=1,NS     
  !    !+- find kernel of subspace S -+!
  !    allocate(tmp_vec(count_NSstates(is)))
  !    !
  !    tmp_vec=0
  !    iss = 0
  !    do i=1,nFock
  !       if(abs(Svalue(i) - get_Svalue(is)).lt.1.d-6) then
  !          iss = iss + 1
  !          tmp_vec(iss) = i
  !       end if
  !    end do
  !    !
  !    !write(*,*) tmp_vec(:)
  !    !the vector tmp_vec contains 

  !    deallocate(tmp_vec)
  ! end do


  ! !stop
  ! !
  ! !<CHECK
  ! ! do i=1,nFock
  ! !    write(*,*) S2eigen(i)
  ! !    write(*,*)
  ! !    write(*,'(20F7.3)') S2(:,i)
  ! !    write(*,*)
  ! ! end do
  ! !CHECK>

