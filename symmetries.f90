subroutine basis_O1xSU2_irr_reps(irr_reps,equ_reps,Virr_reps)  
  !+-BASIS STRUCTURE FOR THE IRREDUCIBLE REPS OF THE GROUP O(1)c x SU(2)s on the local Fock space-+!
  complex(8),dimension(nFock,nFock)               :: Virr_reps ! trasnformation to the irreducible reps
  integer,dimension(:,:),allocatable              :: irr_reps !irreducible reps info: block-structure and equivalent reps
  integer,dimension(:,:),allocatable              :: equ_reps
  complex(8),dimension(nFock,nFock,3)             :: Svec
  complex(8),dimension(nFock,nFock)               :: S2
  real(8),dimension(nFock,nFock)                  :: Ncharge
  !
  complex(8),dimension(2,2,3)                     :: sigma_pauli
  complex(8),dimension(nFock,nFock)               :: Splus,Sminus
  !
  complex(8),dimension(:,:,:),allocatable         :: joint_diag  
  complex(8),dimension(:,:),allocatable           :: jointV 
  real(8),dimension(:,:),allocatable              :: joint_eigen 
  !
  complex(8),dimension(:),allocatable             :: tmp_vec
  real(8)                                         :: modV,modV_  
  integer,dimension(:,:),allocatable              :: eigen_labels
  !
  integer                                         :: i,j,k,iorb,jorb,ispin,jspin,istate,jstate
  integer                                         :: ifock,jfock
  integer                                         :: imin,imax,ii,jj,dim_irr
  integer                                         :: map
  integer,dimension(nFock)                        :: ker_map
  !
  type(local_multiplets),dimension(:),allocatable :: mult_list
  integer,dimension(:,:),allocatable              :: irr_reps_
  !
  integer                                         :: Nirr_reps,jtmp,Nineq_reps
  integer,dimension(:,:),allocatable              :: equ_reps_
  

  !+- build sigma pauli and levi-civita tensor -+!
  sigma_pauli=0.d0
  !
  sigma_pauli(1,2,1) = one
  sigma_pauli(2,1,1) = one
  !
  sigma_pauli(1,2,2) = -xi
  sigma_pauli(2,1,2) =  xi
  !

  sigma_pauli(1,1,3) =  one
  sigma_pauli(2,2,3) = -one
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
  Splus  = Svec(:,:,1) + xi*Svec(:,:,2)
  Sminus = Svec(:,:,1) - xi*Svec(:,:,2)
  !
  Ncharge=0.d0
  do iorb=1,Norb
     do ispin=1,2
        istate=index(ispin,iorb)
        Ncharge = Ncharge + matmul(CC(istate,:,:),CA(istate,:,:))
     end do
  end do
  !


  
  !< TEST JACOBI JOINT DIAGONALIZATION 
  allocate(joint_diag(nFock,nFock,3),jointV(nFock,nFock),joint_eigen(nFock,3))
  joint_diag(:,:,1)=S2
  joint_diag(:,:,2)=Svec(:,:,3)
  joint_diag(:,:,3)=Ncharge
  call simultaneous_diag(joint_diag,jointV,joint_eigen,eps=1.d-10) 
  !
  allocate(eigen_labels(3,1)); eigen_labels = 0
  eigen_labels(1,1) = 1
  eigen_labels(3,1) = 1
  !
  call get_multiplets_list(joint_eigen,eigen_labels,mult_list)

  !+- I obtained the basis for irreducible representation of total-spin rotations -+!
  allocate(tmp_vec(nFock))
  ifock=0
  ker_map = 0
  do i=1,mult_list(1)%N_mult
     do j=1,mult_list(1)%Nequiv_mult(i)
        !
        map = mult_list(1)%Maps(i)%index(j)
        !apply S+
        tmp_vec = jointV(:,map)
        tmp_vec = matmul(Splus,tmp_vec)
        !
        modV = sqrt(dot_product(tmp_vec,tmp_vec))
        !
        if(abs(modV).lt.1.d-10) then
           ker_map(map) = 1
        end if
        !
     end do
  end do

  ifock=0  
  Nirr_reps=0;Nineq_reps=0

  write(*,*) ker_map

  allocate(irr_reps_(nFock,4))
  allocate(equ_reps_(nFock,nFock)); equ_reps_=0
  jtmp=0
  do jj=1,mult_list(1)%N_mult
     do ii=1,mult_list(1)%Nequiv_mult(jj)
        !
        i=mult_list(1)%Maps(jj)%index(ii)        
        if(ker_map(i).eq.1) then
           !
           tmp_vec = jointV(:,i)                   
           modV = sqrt(dot_product(tmp_vec,tmp_vec))
           !
           imin = ifock+1
           !
           dim_irr=0
           do while(modV.gt.1.d-10) 
              !
              dim_irr = dim_irr+1
              ifock = ifock + 1
              Virr_reps(:,ifock) = tmp_vec/modV
              tmp_vec = matmul(Sminus,tmp_vec)
              modV = sqrt(dot_product(tmp_vec,tmp_vec))
              !
           end do
           imax = ifock
           j=mult_list(1)%inv_map(i)
           !
           Nirr_reps = Nirr_reps+1
           !
           if(j.ne.jtmp) Nineq_reps = Nineq_reps+1
           equ_reps_(Nirr_reps,Nineq_reps) = 1
           !
           jtmp=j
           !
           irr_reps_(Nirr_reps,1) = imin 
           irr_reps_(Nirr_reps,2) = imax
           irr_reps_(Nirr_reps,3) = dim_irr
           irr_reps_(Nirr_reps,4) = mult_list(1)%inv_map(i)
           !
        end if
     end do
  end do
  !
  if(allocated(irr_reps)) deallocate(irr_reps)
  if(allocated(equ_reps)) deallocate(equ_reps)
  allocate(irr_reps(Nirr_reps,4))
  allocate(equ_reps(Nirr_reps,Nineq_reps))
  !
  do i=1,Nirr_reps
     irr_reps(i,:) = irr_reps_(i,:)
     do j=1,Nineq_reps
        equ_reps(i,j) = equ_reps_(i,j)
     end do
  end do
  !
end subroutine basis_O1xSU2_irr_reps









subroutine get_matrix_basis_irr_reps(irr_reps,equ_reps,phi_irr)
  integer,dimension(:,:),allocatable :: irr_reps
  integer,dimension(:,:),allocatable :: equ_reps
  integer,dimension(:),allocatable :: map_equ_reps
  integer :: Nirr_reps,Nineq_reps,i,j,Neq,ireps,jreps,ieq,dim_phi,jeq
  integer :: ifock,jfock,idim,jdim,reps_dim,imin,jmin,iphi
  complex(8),dimension(:,:,:),allocatable :: phi_irr


  if(allocated(phi_irr)) deallocate(phi_irr)
  !
  Nirr_reps=size(irr_reps,1)
  Nineq_reps=size(equ_reps,2)
  !
  dim_phi=0
  do i=1,Nineq_reps
     Neq=0
     do j=1,Nirr_reps
        Neq=Neq+equ_reps(j,i)
     end do
     dim_phi = dim_phi + Neq*Neq
  end do

  allocate(phi_irr(dim_phi,nFock,nFock)) ; phi_irr=0.d0

  iphi=0
  do i=1,Nineq_reps
     ieq=0
     do j=1,Nirr_reps
        ieq=ieq+equ_reps(j,i)
     end do
     Neq=ieq

     ieq=0
     allocate(map_equ_reps(Neq))
     do j=1,Nirr_reps
        ieq=ieq+equ_reps(j,i)
        if(equ_reps(j,i).eq.1) map_equ_reps(ieq)=j
     end do
     !write(*,*) 'MAP',map_equ_reps

     do ieq=1,Neq
        do jeq=1,Neq
           iphi=iphi+1
           !write(*,*) iphi

           ireps=map_equ_reps(ieq)
           jreps=map_equ_reps(jeq)           
           imin=irr_reps(ireps,1);jmin=irr_reps(jreps,1)
           idim=irr_reps(ireps,3);jdim=irr_reps(jreps,3)
           if(idim.ne.jdim) stop "Error: equivalent representations with different dimensions"
           reps_dim=idim
           !
           do j=1,reps_dim
              ifock=imin + (j-1)
              jfock=jmin + (j-1)
              !write(*,*) 'blocks identities',ifock,jfock              
              phi_irr(iphi,ifock,jfock) = 1.d0
           end do
           !
        end do
     end do
     ! write(*,*)
     ! write(*,*)
     deallocate(map_equ_reps)
  end do
  
end subroutine get_matrix_basis_irr_reps




subroutine get_matrix_basis_original_fock(phi_irr,phi_fock,Virr_reps)
  complex(8),dimension(:,:,:) :: phi_irr
  complex(8),dimension(nFock,nFock) :: Virr_reps
  complex(8),dimension(:,:,:),allocatable :: phi_fock,phi_fock_
  integer :: Nphi,out_phi
  integer,dimension(:),allocatable :: mask_out
  logical :: flag
  integer :: ifock,jfock,ii,jj,iphi,Nfin,i
  real(8) :: tmp

  Nphi = size(phi_irr,1)
  if(allocated(phi_fock)) deallocate(phi_fock)
  allocate(phi_fock_(Nphi,nFock,nFock))
  allocate(mask_out(Nphi)) ; mask_out=0
  !
  out_phi=0
  do iphi=1,nphi
     tmp=0.d0
     Nfin=0
     flag=.false.
     do ifock=1,nFock
        do jfock=1,nFock
           phi_fock_(iphi,ifock,jfock)=0.d0
           do ii=1,nFock
              do jj=1,nFock
                 phi_fock_(iphi,ifock,jfock)= &
                      phi_fock_(iphi,ifock,jfock) + &
                                !Virr_reps(ii,ifock)*conjg(Virr_reps(jj,jfock))*phi_irr(iphi,ii,jj)
                      Virr_reps(ifock,ii)*conjg(Virr_reps(jfock,jj))*phi_irr(iphi,ii,jj)
              end do
           end do
           tmp = tmp + abs(dreal(phi_fock_(iphi,ifock,jfock)))**2.d0
           if(abs(dreal(phi_fock_(iphi,ifock,jfock))).gt.1.d-10.and.ifock.ne.jfock) flag=.true.
           if(abs(dreal(phi_fock_(iphi,ifock,jfock))).gt.1.d-10) then
              Nfin=Nfin+1
           end if
        end do
     end do

     if(tmp.gt.1.d-8)  then
        out_phi = out_phi+1
        mask_out(iphi) = 1
     end if


     ! do ifock=1,nFock
     !    write(*,'(20F7.2)') dreal(phi_fock(iphi,ifock,:))
     ! end do
     ! write(*,*)
     ! write(*,*)

  end do

  !write(*,*) mask_out
  !stop


  allocate(phi_fock(out_phi,nFock,nFock))
  out_phi=0
  do iphi=1,Nphi
     !
     if(mask_out(iphi).eq.1) then
        out_phi=out_phi+1
        write(*,*) out_phi
        phi_fock(out_phi,:,:)  = phi_fock_(iphi,:,:)
        !
     end if
  end do



end subroutine get_matrix_basis_original_fock







subroutine basisirr_reps  
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
  complex(8),dimension(:,:),allocatable   :: test_jointV,jointV
  complex(8),dimension(:),allocatable   :: tmp_vec
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
  integer,dimension(:),allocatable        :: count_NSstates
  integer,dimension(:,:),allocatable      :: irreducible_states,tmp_irreducible_states,SZ_states
  integer                                 :: i,j,k,iorb,jorb,ispin,jspin,istate,jstate,is,iss,jfock,ifock
  real(8)                                 :: storeS,tmp_sz,deltaS,tmp_test,modV,modV_

  integer                                 :: map,NS,NMS
  integer,dimension(nFock) :: ker_map


  integer :: imap,jmap,imin,imax,jmin,jmax,ii,jj,dim_irr

  type(local_multiplets),dimension(:),allocatable :: mult_list
  type(intarray),dimension(:),allocatable :: irr_reps,irr_reps_
  integer :: Nirr_reps,jtmp,Nineq_reps
  integer,dimension(:,:),allocatable :: equ_reps,equ_reps_


  integer,dimension(:,:),allocatable :: eigen_labels
  logical :: ker_flag


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


  Splus  = Svec(:,:,1) + xi*Svec(:,:,2)
  Sminus = Svec(:,:,1) - xi*Svec(:,:,2)

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
  ! tmp_matrix = Svec(:,:,2)
  ! !call matrix_diagonalize(tmp_matrix,S2eigen,'V','U')
  ! call matrix_diagonalize(tmp_matrix,S2eigen)
  ! S2diag=0.d0
  ! do i=1,nFock
  !    S2diag(i,i)=S2eigen(i)
  !    write(*,*) i,S2eigen(i),Ntest(i,i)
  ! end do

  allocate(test_joint_diag(nFock,nFock,3),test_jointV(nFock,nFock),test_jointD(nFock,3))
  test_joint_diag(:,:,1)=S2
  test_joint_diag(:,:,2)=Svec(:,:,3)
  test_joint_diag(:,:,3)=Ntest

  ! allocate(test_joint_diag(nFock,nFock,5),test_jointV(nFock,nFock),test_jointD(nFock,5))
  ! test_joint_diag(:,:,1)=S2
  ! test_joint_diag(:,:,2)=Svec(:,:,3)
  ! test_joint_diag(:,:,3)=isoS2
  ! test_joint_diag(:,:,4)=isoSvec(:,:,3)
  ! test_joint_diag(:,:,5)=Ntest


  write(*,*)
  call simultaneous_diag(test_joint_diag,test_jointV,test_jointD,eps=1.d-10) 



  allocate(eigen_labels(3,1)); eigen_labels = 0
  eigen_labels(1,1) = 1
  eigen_labels(3,1) = 1
  !
  ! eigen_labels(1,2) = 1
  ! eigen_labels(2,2) = 1
  ! eigen_labels(3,2) = 1
  !
  call get_multiplets_list(test_jointD,eigen_labels,mult_list)

  !+- I obtained the basis for irreducible representation of total-spin rotations -+!
  allocate(jointV(nFock,nFock))
  allocate(tmp_vec(nFock))
  ifock=0
  ker_map = 0
  do i=1,mult_list(1)%N_mult



     do j=1,mult_list(1)%Nequiv_mult(i)

        map = mult_list(1)%Maps(i)%index(j)
        !apply S+
        tmp_vec = test_jointV(:,map)

        ! write(*,*) map
        ! write(*,'(20F7.2)') dreal(tmp_vec)
        ! write(*,'(20F7.2)') dimag(tmp_vec)

        tmp_vec = matmul(Splus,tmp_vec)
        modV=0.d0
        do k=1,nFock
           modV=modV+tmp_vec(k)*conjg(tmp_vec(k))
        end do
        modV=sqrt(modV)

        if(abs(modV).lt.1.d-10) then
           ker_map(map) = 1
        end if

     end do
  end do

  ifock=0  
  Nirr_reps=0
  allocate(irr_reps_(nFock))
  do i=1,nFock
     allocate(irr_reps_(i)%index(4))
  end do
  allocate(equ_reps_(nFock,nFock)); equ_reps_=0
  jtmp=0
  do jj=1,mult_list(1)%N_mult
     do ii=1,mult_list(1)%Nequiv_mult(jj)
        !
        i=mult_list(1)%Maps(jj)%index(ii)        
        if(ker_map(i).eq.1) then
           !
           tmp_vec = test_jointV(:,i)        
           modV=0.d0
           do k=1,nFock
              modV=modV+tmp_vec(k)*conjg(tmp_vec(k))
           end do
           modV=sqrt(modV)                   
           !
           imin = ifock+1
           !
           dim_irr=0
           do while(modV.gt.1.d-10) 
              !
              dim_irr=dim_irr+1
              ifock = ifock + 1
              jointV(:,ifock) = tmp_vec/modV
              tmp_vec = matmul(Sminus,tmp_vec)
              !
              modV=0.d0
              do k=1,nFock
                 modV=modV+tmp_vec(k)*conjg(tmp_vec(k))
              end do
              modV=sqrt(modV)           
           end do
           imax = ifock
           j=mult_list(1)%inv_map(i)
           !
           Nirr_reps = Nirr_reps+1
           !
           if(j.eq.jtmp) then
              equ_reps_(Nirr_reps,Nineq_reps) = 1
           else
              Nineq_reps = Nineq_reps+1
              equ_reps_(Nirr_reps,Nineq_reps) = 1
           end if
           jtmp=j
           !
           irr_reps_(Nirr_reps)%index(1) = imin 
           irr_reps_(Nirr_reps)%index(2) = imax
           irr_reps_(Nirr_reps)%index(3) = dim_irr
           irr_reps_(Nirr_reps)%index(4) = mult_list(1)%inv_map(i)
           !
        end if
     end do
  end do
  write(*,*) Nirr_reps,Nineq_reps

  allocate(irr_reps(Nirr_reps))
  do i=1,Nirr_reps
     allocate(irr_reps(i)%index(4))
  end do
  allocate(equ_reps(Nirr_reps,Nineq_reps))
  !
  do i=1,Nirr_reps
     irr_reps(i)%index(:) = irr_reps_(i)%index(:)
     do j=1,Nineq_reps
        equ_reps(i,j) = equ_reps_(i,j)
     end do
  end do

  write(*,*) 'IRREDUCIBLE REPS'
  do i=1,Nirr_reps
     write(*,*) irr_reps(i)%index(:)
  end do

  write(*,*) 'EQUIVALENT IRREDUCIBLE REPS'
  do i=1,Nirr_reps
     write(*,*) equ_reps(i,:)
  end do


end subroutine basisirr_reps








subroutine get_multiplets_list(eigenvalues,eigen_label,multiplets)
  real(8),dimension(:,:)               :: eigenvalues
  integer,dimension(:,:)   :: eigen_label
  type(local_multiplets),dimension(:),allocatable :: multiplets
  !
  integer,dimension(:),allocatable :: Nlabels
  real(8),dimension(:,:),allocatable :: group_eigenvalues
  real(8),dimension(:,:),allocatable :: indep_eigen
  integer,dimension(:,:),allocatable :: group_mask
  integer :: N,M,Nl
  integer :: i,j,k,count_equal_states
  integer :: igroup,tmp,imap

  !
  N=size(eigenvalues,1)
  M=size(eigenvalues,2)
  if(size(eigen_label,1).ne.M) stop "get multiplets: wrong dimension"
  Nl=size(eigen_label,2)
  !
  if(allocated(multiplets)) deallocate(multiplets)
  allocate(multiplets(Nl))
  !
  allocate(Nlabels(Nl))  
  Nlabels=0
  do i=1,M
     Nlabels(:) = Nlabels(:) + eigen_label(i,:)
  end do

  do j=1,Nl
     !
     allocate(group_eigenvalues(N,Nlabels(j)))
     igroup=0
     do i=1,M     
        if(eigen_label(i,j).eq.1) then
           igroup=igroup+1
           group_eigenvalues(:,igroup) = eigenvalues(:,i)
        end if
     end do
     !
     count_equal_states = get_mask_equal_eigenvalues(group_eigenvalues,group_mask,indep_eigen)  
     ! !
     multiplets(j)%N_mult = count_equal_states
     multiplets(j)%size_mult = Nlabels(j)
     allocate(multiplets(j)%mult(count_equal_states,Nlabels(j)))
     allocate(multiplets(j)%Nequiv_mult(count_equal_states))
     allocate(multiplets(j)%Maps(count_equal_states))
     allocate(multiplets(j)%inv_map(N))
     ! !
     do i=1,count_equal_states
        multiplets(j)%mult(i,1:Nlabels(j)) = indep_eigen(i,1:Nlabels(j))        
        !
        tmp=0
        do k=1,N
           tmp = tmp + group_mask(i,k)           
        end do
        multiplets(j)%Nequiv_mult(i) = tmp
        ! !
        allocate(multiplets(j)%Maps(i)%index(tmp))
        !multiplets(j)%Maps(i)%index
        imap=0
        do k=1,N
           if(group_mask(i,k).eq.1) then
              imap = imap+1
              multiplets(j)%Maps(i)%index(imap) = k
              multiplets(j)%inv_map(k) = i
           end if
        end do
     end do
     deallocate(group_eigenvalues)

     write(*,*)
     write(*,*) 'Label',eigen_label(:,j)
     write(*,*)
     write(*,*) 'Size of multiplets',multiplets(j)%size_mult
     write(*,*)
     write(*,*) 'Number of different multiplets',multiplets(j)%N_mult
     write(*,*)
     write(*,*) 'Multiplets list - multeplicity - maps to fock index'
     do i=1,multiplets(j)%N_mult
        write(*,*) multiplets(j)%mult(i,:),multiplets(j)%Nequiv_mult(i),&
             multiplets(j)%Maps(i)%index(:)
     end do
  end do

end subroutine get_multiplets_list










subroutine group_statesNS(eigenvalues)
  real(8),dimension(:,:)               :: eigenvalues
  real(8),dimension(:,:),allocatable               :: indep_eigen_
  real(8),dimension(:,:,:),allocatable               :: indep_eigen
  integer,dimension(:,:),allocatable   :: eigen_label_mask
  integer,dimension(:,:),allocatable   :: map_group_
  integer,dimension(:),allocatable     :: map_group,mult,mult_
  integer,dimension(:,:),allocatable   :: group_mask_
  integer,dimension(:,:),allocatable   :: final_mask
  integer,dimension(:,:,:),allocatable :: group_mask
  real(8),dimension(:,:),allocatable   :: group_eigenvalues
  real(8),dimension(:),allocatable   :: NS_multiplets,NSSz_multiplets
  integer,dimension(2)                 :: count_equal_states,Ngroup
  integer                              :: igroup,i,j,k
  integer                              :: M,N
  integer                              :: imap,imap_,Ntmp,tmp,imult
  logical                              :: tmp_flag
  !
  N=size(eigenvalues,1)
  M=size(eigenvalues,2)
  !
  allocate(eigen_label_mask(M,2))
  eigen_label_mask=0
  !
  !+- S2 eigenvalues -+!
  eigen_label_mask(1,1)=1
  eigen_label_mask(5,1)=1
  !+- Sz eigenvalues -+!
  eigen_label_mask(1,2)=1
  eigen_label_mask(2,2)=1
  eigen_label_mask(5,2)=1
  !
  !
  Ngroup=0
  do i=1,M
     Ngroup(1:2) = Ngroup(1:2) + eigen_label_mask(i,1:2)
  end do
  !


  write(*,*) 'CREATE MASKS'
  
  allocate(group_mask(N,N,2)); group_mask=0
  allocate(indep_eigen(N,M,2)); indep_eigen=0.d0
  do j=1,2
     allocate(group_eigenvalues(N,Ngroup(j)))
     igroup=0
     do i=1,M     
        if(eigen_label_mask(i,j).eq.1) then
           igroup=igroup+1
           group_eigenvalues(:,igroup) = eigenvalues(:,i)
        end if
     end do
     write(*,*) j
     if(.not.allocated(group_mask_))allocate(group_mask_(N,N))
     count_equal_states(j) = get_mask_equal_eigenvalues(group_eigenvalues,group_mask_,indep_eigen_)  
     do i=1,count_equal_states(j)
        group_mask(i,:,j)  = group_mask_(i,:)        
        indep_eigen(i,1:Ngroup(j),j) = indep_eigen_(i,1:Ngroup(j))
     end do
     !<TEST
     do i=1,N
        write(*,*) i,group_mask(1:count_equal_states(j),i,j)
     end do
     write(*,*)
     !TEST>
     deallocate(group_eigenvalues)
  end do  
  !
  !allocate(map_group_(N,2),map_group(N))

  write(*,*) 'CREATE MAPS'

  allocate(map_group_(N,2))
  do j=1,2     
     imap=0
     do i=1,count_equal_states(j)
        do k=1,N
           if(group_mask(i,k,j).eq.1) then
              imap=imap+1
              map_group_(imap,j) = k
           end if
        end do
     end do
     !<TEST
     do i=1,N
        write(*,'(5F7.2,I4)') eigenvalues(map_group_(i,j),:),map_group_(i,j)
     end do
     write(*,*)     
     !TEST>
  end do


  !
  ! allocate(NS_mulitplets(Ngroup(1)),NSSz_mulitplets(Ngroup(2)))  
  ! tmp_flag=.true.
  ! j=1
  ! do i=1,Ngroup(j)
  !    tmp=0
  !    do k=1,N
  !       tmp = tmp + group_mask(i,k,j)
  !       if(group_mask(i,k,j).eq.1.and.tmp_flag) then
  !          NS_multiplets
  !       end if
  !    end do
     
  ! end do

  
  write(*,*) 'STORE MULTIPLETS'


  !
  allocate(map_group(N),mult_(N))
  imap=0
  mult_=0
  imult=0
  do i=1,count_equal_states(1)
     do j=1,count_equal_states(2)
        tmp=0
        tmp_flag=.true.
        do k=1,N
           tmp = tmp + group_mask(i,k,1)*group_mask(j,k,2)
           if(group_mask(i,k,1)*group_mask(j,k,2).eq.1) then
              if(tmp_flag) then
                 imap_=k; tmp_flag=.false.
              end if
              imap=imap+1
              map_group(imap) = k                            
           end if
        end do
        !
        if(tmp.ne.0) then
           imult=imult+1
           !write(*,*) i,j,indep_eigen(imult,1:Ngroup(2),2),tmp,imap_
           !write(*,*) i,j,eigenvalues(map_group(imap_),1),eigenvalues(map_group(imap_),2),eigenvalues(map_group(imap_),5),tmp,imap_
           write(*,*) eigenvalues(imap_,1),eigenvalues(imap_,2),eigenvalues(imap_,5),tmp
           mult_(imult)=tmp
        end if
        !
     end do
  end do

  allocate(mult(imult))
  do i=1,imult
     mult(imult) = mult_(imult)
  end do

  !<TEST
  do i=1,N
     write(*,'(5F7.2,2I4)') eigenvalues(map_group(i),:),map_group(i)
  end do
  write(*,*)     
  !TEST>




end subroutine group_statesNS




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


function get_mask_equal_eigenvalues(xValues,mask,indep_values,eps) result(Ns)
  !
  real(8),dimension(:,:)               :: xValues
  integer,dimension(:,:),allocatable :: mask
  real(8),dimension(:,:),allocatable :: indep_values
  real(8),optional                   :: eps
  integer                            :: Ns
  real(8)                            :: eps_
  real(8),dimension(:,:),allocatable :: tmp_mask
  integer,dimension(:),allocatable   :: search_index
  integer                            :: Nx,Nm,ix,jx,is,k,flag
  real(8) :: deltaX
  !
  eps_=1.d-8
  if(present(eps)) eps_=eps
  if(allocated(mask)) deallocate(mask)
  if(allocated(indep_values)) deallocate(indep_values)

  Nx=size(xValues,1); Nm=size(xValues,2)
  allocate(tmp_mask(Nx,Nx),search_index(Nx))  
  tmp_mask=0
  search_index=0
  NS=0
  do ix=1,Nx
     if(search_index(ix).ge.0) then
        NS = NS + 1
        do jx=1,Nx
           !
           deltaX=0.d0
           do k=1,Nm
              deltaX = deltaX + (Xvalues(jx,k)-Xvalues(ix,k))**2.d0
           end do
           !
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

  allocate(indep_values(NS,Nm))
  do is=1,NS
     flag=0
     do ix=1,Nx
        if(mask(is,ix).eq.1.and.flag.eq.0) then
           indep_values(is,:) = Xvalues(ix,:)
           !write(*,*) Xvalues(ix,:)
           flag=1
        end if
     end do
  end do

  !
  deallocate(tmp_mask)  
  !
end function get_mask_equal_eigenvalues




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

