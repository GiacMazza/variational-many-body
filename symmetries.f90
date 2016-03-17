subroutine basis_O1cXSU2sXSU2c_irr_reps(irr_reps,equ_reps,Virr_reps)  
  !+-BASIS STRUCTURE FOR THE IRREDUCIBLE REPS OF THE GROUP O(1)c x SU(2)s on the local Fock space-+!
  complex(8),dimension(nFock,nFock)               :: Virr_reps ! trasnformation to the irreducible reps
  integer,dimension(:,:),allocatable              :: irr_reps !irreducible reps info: block-structure and equivalent reps
  integer,dimension(:,:),allocatable              :: equ_reps
  complex(8),dimension(nFock,nFock,3)             :: Svec,isoSvec
  complex(8),dimension(nFock,nFock)               :: S2,isoS2
  real(8),dimension(nFock,nFock)                  :: Ncharge
  real(8),dimension(nFock)                  :: tmp_eigen
  complex(8),dimension(nFock,nFock)                  :: tmp_matrix
  !
  complex(8),dimension(2,2,3)                     :: sigma_pauli
  complex(8),dimension(3,3,3)                     :: levi_civita
  complex(8),dimension(nFock,nFock)               :: Splus,Sminus
  complex(8),dimension(nFock,nFock)               :: isoSplus,isoSminus
  !
  complex(8),dimension(:,:,:),allocatable         :: joint_diag  
  complex(8),dimension(:,:),allocatable           :: jointV 
  real(8),dimension(:,:),allocatable              :: joint_eigen 
  !
  complex(8),dimension(:),allocatable             :: tmp_vec,tmp_vec_
  complex(8),dimension(:),allocatable             :: tmp_vec_S,tmp_vec_isoS
  real(8)                                         :: modV,modV_  
  real(8)                                         :: modV_S,modV_isoS
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

  !
  levi_civita=zero
  levi_civita(1,2,3) =  one
  levi_civita(1,3,2) = -one
  !
  levi_civita(2,3,1) =  one
  levi_civita(2,1,3) = -one
  !
  levi_civita(3,1,2) =  one
  levi_civita(3,2,1) = -one
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

  isoS2=zero
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
                 isoSvec(:,:,i) = isoSvec(:,:,i) - &
                      xi*levi_civita(i,iorb,jorb)*matmul(CC(istate,:,:),CA(jstate,:,:))
              end do
           end do
        end do
     end select
     isoS2 = isoS2 + matmul(isoSvec(:,:,i),isoSvec(:,:,i))
  end do
  isoSplus  = isoSvec(:,:,1) + xi*isoSvec(:,:,2)
  isoSminus = isoSvec(:,:,1) - xi*isoSvec(:,:,2)  
  !
  Ncharge=0.d0
  do iorb=1,Norb
     do ispin=1,2
        istate=index(ispin,iorb)
        Ncharge = Ncharge + matmul(CC(istate,:,:),CA(istate,:,:))
     end do
  end do
  !

  !
  allocate(joint_diag(nFock,nFock,5),jointV(nFock,nFock),joint_eigen(nFock,5))
  joint_diag(:,:,1)=S2(:,:)
  joint_diag(:,:,2)=Svec(:,:,3)
  joint_diag(:,:,3)=isoS2
  joint_diag(:,:,4)=isoSvec(:,:,3)
  joint_diag(:,:,5)=Ncharge
  !
  call simultaneous_diag(joint_diag,jointV,joint_eigen,eps=1.d-10) 
  !
  allocate(eigen_labels(5,1)); eigen_labels = 0
  eigen_labels(1,1) = 1
  eigen_labels(3,1) = 1
  eigen_labels(5,1) = 1
  !
  call get_multiplets_list(joint_eigen,eigen_labels,mult_list)
  !


  !+- here I should obtain the simultanoues kernels of S+ and L+
  write(*,*)
  !+- I obtained the basis for irreducible representation of total-spin rotations -+!
  allocate(tmp_vec(nFock),tmp_vec_(nFock))
  ifock=0
  ker_map = 0
  do i=1,mult_list(1)%N_mult
     do j=1,mult_list(1)%Nequiv_mult(i)
        !
        map = mult_list(1)%Maps(i)%index(j)

        tmp_vec = jointV(:,map)
        ! apply S+
        !        write(*,*) map
        tmp_vec = matmul(Splus,tmp_vec)
        modV = sqrt(dot_product(tmp_vec,tmp_vec))
        !write(*,*) modV
        ! apply L+
        tmp_vec = jointV(:,map)
        tmp_vec = matmul(isoSplus,tmp_vec)
        modV_ = sqrt(dot_product(tmp_vec,tmp_vec))
        ! write(*,*) modV_
        ! write(*,*)
        !
        if(abs(modV).lt.1.d-10.and.abs(modV_).lt.1.d-10) then
           ker_map(map) = 1
        end if
        !
     end do
  end do

  ifock=0  
  Nirr_reps=0;Nineq_reps=0

  write(*,*) ker_map  


  allocate(tmp_vec_S(nFock),tmp_vec_isoS(nFock))
  allocate(irr_reps_(nFock,4))
  allocate(equ_reps_(nFock,nFock)); equ_reps_=0
  jtmp=0
  do jj=1,mult_list(1)%N_mult
     do ii=1,mult_list(1)%Nequiv_mult(jj)
        !
        i=mult_list(1)%Maps(jj)%index(ii)        
        if(ker_map(i).eq.1) then
           !
           imin = ifock+1
           dim_irr=0
           !
           write(*,*) jj,ii
           !

           !
           tmp_vec_S = jointV(:,i) 
           modV_S = sqrt(dot_product(tmp_vec_S,tmp_vec_S)) 
           !
           do while(modV_S.gt.1.d-10)
              !
              tmp_vec_isoS = tmp_vec_S
              modV_isoS = sqrt(dot_product(tmp_vec_isoS,tmp_vec_isoS))

              tmp_vec_S = matmul(Sminus,tmp_vec_S)
              modV_S = sqrt(dot_product(tmp_vec_S,tmp_vec_S))
              !
              do while(modV_isoS.gt.1.d-10)
                 !
                 dim_irr = dim_irr+1
                 ifock = ifock + 1
                 write(*,*) 'ifock',ifock
                 Virr_reps(:,ifock) = tmp_vec_isoS/modV_isoS
                 !
                 tmp_vec_isoS = matmul(isoSminus,tmp_vec_isoS)
                 modV_isoS = sqrt(dot_product(tmp_vec_isoS,tmp_vec_isoS))
                 !
              end do
           end do
           !
           imax = ifock
           ! write(*,*) 'ifock',ifock,dim_irr
           ! write(*,*) 'imin,imax',imin,imax
           !
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
end subroutine basis_O1cXSU2sXSU2c_irr_reps


subroutine basis_O1cXSU2sXisoZ_irr_reps(irr_reps,equ_reps,Virr_reps)  
  !+-BASIS STRUCTURE FOR THE IRREDUCIBLE REPS OF THE GROUP O(1)c x SU(2)s on the local Fock space-+!
  complex(8),dimension(nFock,nFock)               :: Virr_reps ! trasnformation to the irreducible reps
  integer,dimension(:,:),allocatable              :: irr_reps !irreducible reps info: block-structure and equivalent reps
  integer,dimension(:,:),allocatable              :: equ_reps
  complex(8),dimension(nFock,nFock,3)             :: Svec,isoSvec
  complex(8),dimension(nFock,nFock)               :: S2,isoS2
  real(8),dimension(nFock,nFock)                  :: Ncharge
  real(8),dimension(nFock)                  :: tmp_eigen
  complex(8),dimension(nFock,nFock)                  :: tmp_matrix
  !
  complex(8),dimension(2,2,3)                     :: sigma_pauli
  complex(8),dimension(3,3,3)                     :: levi_civita
  complex(8),dimension(nFock,nFock)               :: Splus,Sminus
  complex(8),dimension(nFock,nFock)               :: isoSplus,isoSminus
  !
  complex(8),dimension(:,:,:),allocatable         :: joint_diag  
  complex(8),dimension(:,:),allocatable           :: jointV 
  real(8),dimension(:,:),allocatable              :: joint_eigen 
  !
  complex(8),dimension(:),allocatable             :: tmp_vec,tmp_vec_
  complex(8),dimension(:),allocatable             :: tmp_vec_S,tmp_vec_isoS
  real(8)                                         :: modV,modV_  
  real(8)                                         :: modV_S,modV_isoS
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

  !
  levi_civita=zero
  levi_civita(1,2,3) =  one
  levi_civita(1,3,2) = -one
  !
  levi_civita(2,3,1) =  one
  levi_civita(2,1,3) = -one
  !
  levi_civita(3,1,2) =  one
  levi_civita(3,2,1) = -one
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

  isoS2=zero
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
                 isoSvec(:,:,i) = isoSvec(:,:,i) - &
                      xi*levi_civita(i,iorb,jorb)*matmul(CC(istate,:,:),CA(jstate,:,:))
              end do
           end do
        end do
     end select
     isoS2 = isoS2 + matmul(isoSvec(:,:,i),isoSvec(:,:,i))
  end do
  isoSplus  = isoSvec(:,:,1) + xi*isoSvec(:,:,2)
  isoSminus = isoSvec(:,:,1) - xi*isoSvec(:,:,2)  
  !
  Ncharge=0.d0
  do iorb=1,Norb
     do ispin=1,2
        istate=index(ispin,iorb)
        Ncharge = Ncharge + matmul(CC(istate,:,:),CA(istate,:,:))
     end do
  end do
  !
  allocate(joint_diag(nFock,nFock,4),jointV(nFock,nFock),joint_eigen(nFock,4))
  joint_diag(:,:,1)=S2(:,:)
  joint_diag(:,:,2)=Svec(:,:,3)
  joint_diag(:,:,3)=isoSvec(:,:,3)
  joint_diag(:,:,4)=Ncharge
  !
  call simultaneous_diag(joint_diag,jointV,joint_eigen,eps=1.d-10) 
  !
  allocate(eigen_labels(4,1)); eigen_labels = 0
  eigen_labels(1,1) = 1
  eigen_labels(3,1) = 1
  eigen_labels(4,1) = 1
  !
  call get_multiplets_list(joint_eigen,eigen_labels,mult_list)
  !
  !+- here I should obtain the simultanoues kernels of S+ -+!
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
end subroutine basis_O1cXSU2sXisoZ_irr_reps



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
  !+- JACOBI JOINT DIAGONALIZATION -+!
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
!
subroutine basis_SU2sXSU2c_irr_reps(irr_reps,equ_reps,Virr_reps)  
  !+-BASIS STRUCTURE FOR THE IRREDUCIBLE REPS OF THE GROUP O(1)c x SU(2)s on the local Fock space-+!
  complex(8),dimension(nFock,nFock)               :: Virr_reps ! trasnformation to the irreducible reps
  integer,dimension(:,:),allocatable              :: irr_reps !irreducible reps info: block-structure and equivalent reps
  integer,dimension(:,:),allocatable              :: equ_reps
  complex(8),dimension(nFock,nFock,3)             :: Svec,isoSvec
  complex(8),dimension(nFock,nFock)               :: S2,isoS2
  real(8),dimension(nFock,nFock)                  :: Ncharge
  real(8),dimension(nFock)                  :: tmp_eigen
  complex(8),dimension(nFock,nFock)                  :: tmp_matrix
  !
  complex(8),dimension(2,2,3)                     :: sigma_pauli
  complex(8),dimension(3,3,3)                     :: levi_civita
  complex(8),dimension(nFock,nFock)               :: Splus,Sminus
  complex(8),dimension(nFock,nFock)               :: isoSplus,isoSminus
  !
  complex(8),dimension(:,:,:),allocatable         :: joint_diag  
  complex(8),dimension(:,:),allocatable           :: jointV 
  real(8),dimension(:,:),allocatable              :: joint_eigen 
  !
  complex(8),dimension(:),allocatable             :: tmp_vec,tmp_vec_
  complex(8),dimension(:),allocatable             :: tmp_vec_S,tmp_vec_isoS
  real(8)                                         :: modV,modV_  
  real(8)                                         :: modV_S,modV_isoS
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

  !
  levi_civita=zero
  levi_civita(1,2,3) =  one
  levi_civita(1,3,2) = -one
  !
  levi_civita(2,3,1) =  one
  levi_civita(2,1,3) = -one
  !
  levi_civita(3,1,2) =  one
  levi_civita(3,2,1) = -one
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

  isoS2=zero
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
                 isoSvec(:,:,i) = isoSvec(:,:,i) - &
                      xi*levi_civita(i,iorb,jorb)*matmul(CC(istate,:,:),CA(jstate,:,:))
              end do
           end do
        end do
     end select
     isoS2 = isoS2 + matmul(isoSvec(:,:,i),isoSvec(:,:,i))
  end do
  isoSplus  = isoSvec(:,:,1) + xi*isoSvec(:,:,2)
  isoSminus = isoSvec(:,:,1) - xi*isoSvec(:,:,2)  
  !
  Ncharge=0.d0
  do iorb=1,Norb
     do ispin=1,2
        istate=index(ispin,iorb)
        Ncharge = Ncharge + matmul(CC(istate,:,:),CA(istate,:,:))
     end do
  end do
  !

  !
  allocate(joint_diag(nFock,nFock,4),jointV(nFock,nFock),joint_eigen(nFock,4))
  joint_diag(:,:,1)=S2(:,:)
  joint_diag(:,:,2)=Svec(:,:,3)
  joint_diag(:,:,3)=isoS2
  joint_diag(:,:,4)=isoSvec(:,:,3)
  !
  call simultaneous_diag(joint_diag,jointV,joint_eigen,eps=1.d-10) 
  !
  allocate(eigen_labels(4,1)); eigen_labels = 0
  eigen_labels(1,1) = 1
  eigen_labels(3,1) = 1
  !
  call get_multiplets_list(joint_eigen,eigen_labels,mult_list)
  !
  !+- here I should obtain the simultanoues kernels of S+ and L+
  !+- I obtained the basis for irreducible representation of total-spin rotations -+!
  allocate(tmp_vec(nFock),tmp_vec_(nFock))
  ifock=0
  ker_map = 0
  do i=1,mult_list(1)%N_mult
     do j=1,mult_list(1)%Nequiv_mult(i)
        !
        map = mult_list(1)%Maps(i)%index(j)

        tmp_vec = jointV(:,map)
        ! apply S+
        tmp_vec = matmul(Splus,tmp_vec)
        modV = sqrt(dot_product(tmp_vec,tmp_vec))
        ! apply L+
        tmp_vec = jointV(:,map)
        tmp_vec = matmul(isoSplus,tmp_vec)
        modV_ = sqrt(dot_product(tmp_vec,tmp_vec))
        !
        if(abs(modV).lt.1.d-10.and.abs(modV_).lt.1.d-10) then
           ker_map(map) = 1
        end if
        !
     end do
  end do
  !
  ifock=0  
  Nirr_reps=0;Nineq_reps=0
  allocate(tmp_vec_S(nFock),tmp_vec_isoS(nFock))
  allocate(irr_reps_(nFock,4))
  allocate(equ_reps_(nFock,nFock)); equ_reps_=0
  jtmp=0
  do jj=1,mult_list(1)%N_mult
     do ii=1,mult_list(1)%Nequiv_mult(jj)
        !
        i=mult_list(1)%Maps(jj)%index(ii)        
        if(ker_map(i).eq.1) then
           !
           imin = ifock+1
           dim_irr=0
           !
           tmp_vec_S = jointV(:,i) 
           modV_S = sqrt(dot_product(tmp_vec_S,tmp_vec_S)) 
           !
           do while(modV_S.gt.1.d-10)
              !
              tmp_vec_isoS = tmp_vec_S
              modV_isoS = sqrt(dot_product(tmp_vec_isoS,tmp_vec_isoS))

              tmp_vec_S = matmul(Sminus,tmp_vec_S)
              modV_S = sqrt(dot_product(tmp_vec_S,tmp_vec_S))
              !
              do while(modV_isoS.gt.1.d-10)
                 !
                 dim_irr = dim_irr+1
                 ifock = ifock + 1
                 write(*,*) 'ifock',ifock
                 Virr_reps(:,ifock) = tmp_vec_isoS/modV_isoS
                 !
                 tmp_vec_isoS = matmul(isoSminus,tmp_vec_isoS)
                 modV_isoS = sqrt(dot_product(tmp_vec_isoS,tmp_vec_isoS))
                 !
              end do
           end do
           !
           imax = ifock
           !
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
end subroutine basis_SU2sXSU2c_irr_reps
!
subroutine basis_SU2_irr_reps(irr_reps,equ_reps,Virr_reps)  
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
  complex(8),dimension(:,:,:),allocatable         :: test_joint_diag  
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
  !+- TEST JACOBI JOINT DIAGONALIZATION  -+!
  allocate(joint_diag(nFock,nFock,2),jointV(nFock,nFock),joint_eigen(nFock,2))
  joint_diag(:,:,1)=S2
  joint_diag(:,:,2)=Svec(:,:,3)
  !
  call simultaneous_diag(joint_diag,jointV,joint_eigen,eps=1.d-10) 
  !
  allocate(eigen_labels(2,1)); eigen_labels = 0
  eigen_labels(1,1) = 1
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
  !< TMP_DEBUG
  ! joint_diag(:,:,1)=S2
  ! joint_diag(:,:,2)=Svec(:,:,3)
  ! allocate(test_joint_diag(nFock,nFock,2)); test_joint_diag=zero
  ! do j=1,2     
  !    do ifock=1,nFock
  !       do jfock=1,nFock
  !          test_joint_diag(ifock,jfock,j) = zero
  !          do ii=1,nFock
  !             do jj=1,nFock
  !                test_joint_diag(ifock,jfock,j) = test_joint_diag(ifock,jfock,j) + &
  !                     conjg(Virr_reps(ii,ifock))*Virr_reps(jj,jfock)*joint_diag(ii,jj,j)
  !             end do
  !          end do
  !       end do
  !    end do

  !    do ifock=1,nFock
  !       write(*,'(20F6.2)') dreal(test_joint_diag(ifock,:,j))
  !    end do

  !    write(*,*)
  !    do ifock=1,nFock
  !       write(*,'(20F6.2)') dimag(test_joint_diag(ifock,:,j))
  !    end do

  !    write(*,*)
  !    write(*,*)
  ! end do
  !TMP_DEBUG>
end subroutine basis_SU2_irr_reps




subroutine basis_SZ_irr_reps(irr_reps,equ_reps,Virr_reps)  
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
  allocate(joint_diag(nFock,nFock,2),jointV(nFock,nFock),joint_eigen(nFock,2))
  joint_diag(:,:,1)=S2
  joint_diag(:,:,2)=Svec(:,:,3)
  call simultaneous_diag(joint_diag,jointV,joint_eigen,eps=1.d-10) 
  !
  allocate(eigen_labels(2,1)); eigen_labels = 0
  eigen_labels(1,1) = 1
  eigen_labels(2,1) = 1
  !
  call get_multiplets_list(joint_eigen,eigen_labels,mult_list)

  !+- I obtained the basis for irreducible representation of total-spin rotations -+!
  allocate(tmp_vec(nFock))
  allocate(irr_reps_(nFock,4))
  allocate(equ_reps_(nFock,nFock)); equ_reps_=0
  jtmp=0
  ifock=0  
  Nirr_reps=0;Nineq_reps=0  
  do jj=1,mult_list(1)%N_mult !+- sweep over all the multiplets
     Nineq_reps = Nineq_reps+1
     do ii=1,mult_list(1)%Nequiv_mult(jj) !+- sweep over the equivalents multiplets
        !
        i=mult_list(1)%Maps(jj)%index(ii) !+- get fock space index               
        imin = ifock + 1; dim_irr = 1; imax=imin
        ifock = ifock + 1
        Virr_reps(:,ifock) = jointV(:,i)
        !
        Nirr_reps = Nirr_reps+1
        equ_reps_(Nirr_reps,Nineq_reps) = 1
        !
        jtmp=j
        !
        write(*,*) Nirr_reps,'Nirr_reps'
        irr_reps_(Nirr_reps,1) = imin 
        irr_reps_(Nirr_reps,2) = imax
        irr_reps_(Nirr_reps,3) = dim_irr
        irr_reps_(Nirr_reps,4) = mult_list(1)%inv_map(i)
        !
     end do
     !
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
end subroutine basis_SZ_irr_reps
















subroutine get_matrix_basis_irr_reps(irr_reps,equ_reps,phi_irr)
  integer,dimension(:,:),allocatable :: irr_reps
  integer,dimension(:,:),allocatable :: equ_reps
  integer,dimension(:),allocatable :: map_equ_reps
  integer :: Nirr_reps,Nineq_reps,i,j,Neq,ireps,jreps,ieq,dim_phi,jeq
  integer :: ifock,jfock,idim,jdim,reps_dim,imin,jmin,iphi
  complex(8),dimension(:,:,:),allocatable :: phi_irr
  complex(8),dimension(:,:),allocatable :: test_phi

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
  !
  allocate(phi_irr(dim_phi,nFock,nFock)) ; phi_irr=zero

  iphi=0
  do i=1,Nineq_reps
     ieq=0
     do j=1,Nirr_reps
        ieq=ieq+equ_reps(j,i)
     end do
     Neq=ieq
     !
     ieq=0
     allocate(map_equ_reps(Neq))
     do j=1,Nirr_reps
        ieq=ieq+equ_reps(j,i)
        if(equ_reps(j,i).eq.1) map_equ_reps(ieq)=j
     end do
     !
     do ieq=1,Neq
        do jeq=1,Neq
           !
           iphi=iphi+1
           !
           ireps=map_equ_reps(ieq)
           jreps=map_equ_reps(jeq)           
           imin=irr_reps(ireps,1);jmin=irr_reps(jreps,1)
           idim=irr_reps(ireps,3);jdim=irr_reps(jreps,3)
           if(idim.ne.jdim) stop "Error: equivalent representations with different dimensions"
           reps_dim=idim
           !
           do j=1,reps_dim
              !
              ifock=imin + (j-1)
              jfock=jmin + (j-1)
              !
              phi_irr(iphi,ifock,jfock) = 1.d0
              !
           end do
        end do
     end do
     deallocate(map_equ_reps)
  end do
  !<TMP_DEBUG
  ! allocate(test_phi(nFock,nFock)); test_phi=zero
  ! do iphi=1,dim_phi
  !    test_phi = test_phi + phi_irr(iphi,:,:)
  ! end do
  
  ! do ifock=1,nFock
  !    write(*,'(20F7.3)') dreal(test_phi(ifock,:))
  ! end do
  ! write(*,*)
  ! do ifock=1,nFock
  !    write(*,'(20F7.3)') dimag(test_phi(ifock,:))
  ! end do
  !TMP_DEBUG>
end subroutine get_matrix_basis_irr_reps






subroutine get_matrix_basis_original_fock(phi_irr,phi_fock,Virr_reps)
  complex(8),dimension(:,:,:) :: phi_irr
  complex(8),dimension(nFock,nFock) :: Virr_reps
  complex(8),dimension(:,:,:),allocatable :: phi_fock,phi_fock_
  complex(8),dimension(:,:),allocatable :: test_phi
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
                      Virr_reps(ifock,ii)*conjg(Virr_reps(jfock,jj))*phi_irr(iphi,ii,jj)
              end do
           end do
           tmp = tmp + phi_fock_(iphi,ifock,jfock)*conjg(phi_fock_(iphi,ifock,jfock))
           if(abs(dreal(phi_fock_(iphi,ifock,jfock))).gt.1.d-10.and.ifock.ne.jfock) flag=.true.
           if(abs(dreal(phi_fock_(iphi,ifock,jfock))).gt.1.d-10) then
              Nfin=Nfin+1
           end if
        end do
     end do
     !
     if(tmp.gt.1.d-8)  then
        out_phi = out_phi+1
        mask_out(iphi) = 1
     end if
     !
  end do

  allocate(phi_fock(out_phi,nFock,nFock))
  out_phi=0
  do iphi=1,Nphi
     !
     if(mask_out(iphi).eq.1) then
        out_phi=out_phi+1
        phi_fock(out_phi,:,:)  = phi_fock_(iphi,:,:)
     end if
  end do
  !<TMP_DEBUG
  ! allocate(test_phi(nFock,nFock)); test_phi=zero
  ! do iphi=1,nphi
  !    test_phi = test_phi + phi_fock(iphi,:,:)
  ! end do
  ! do ifock=1,nFock
  !    write(*,'(20F7.3)') dreal(test_phi(ifock,:))
  ! end do
  ! write(*,*)
  ! do ifock=1,nFock
  !    write(*,'(20F7.3)') dimag(test_phi(ifock,:))
  ! end do
  !TMP_DEBUG>
end subroutine get_matrix_basis_original_fock

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
  Nl=size(eigen_label,2)   !+- number of different labels
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
!
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



