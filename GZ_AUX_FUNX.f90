MODULE GZ_AUX_FUNX
  USE GZ_VARS_GLOBAL
  USE SF_IOTOOLS
  implicit none
  private
  !
  public :: build_local_fock
  public :: bdecomp
  public :: get_spin_indep_states
  public :: initialize_local_density
  !
CONTAINS

  subroutine build_local_fock
    !    integer                            :: State_dim
    integer                            :: fock_state,ifock
    integer                            :: i,iup,jup,idw,iup_,idw_,np,nup_
    integer                            :: ispin,iorb,istate
    integer                            :: nup,ndw,irnd
    integer,dimension(:),allocatable   :: Fock,ivec
    integer,dimension(:,:),allocatable :: Fock_vec
    real(8)                            :: rnd    
    !
    State_dim = 2*Norb        
    NFock = 2**State_dim
    !
    allocate(ivec(State_dim))
    allocate(Fock(Nfock))
    do ifock=1,NFock
       Fock(ifock)=ifock
       call bdecomp(Fock(ifock),ivec)
       write(*,'(10I3)') ivec(:),ifock
    end do
    !

    !+- allocate and initialize stride  -+! 
    allocate(index(2,Norb))
    do ispin=1,2
       do iorb=1,Norb
          index(ispin,iorb)=iorb+(ispin-1)*Norb
       enddo
    end do

    !+- allocate and build creation/anhhilation operators -+!
    allocate(CC(state_dim,nFock,nFock),CA(state_dim,nFock,nFock))
    do ispin=1,2
       do iorb=1,Norb
          istate=index(ispin,iorb)
          call build_c_operator(ispin,iorb,CA(istate,:,:))
          call build_cdg_operator(ispin,iorb,CC(istate,:,:))          
       end do
    end do

    !+- allocate and build local operators -+!
    allocate(UHubbard(nFock,nFock),docc(Norb,nFock,nFock),dens(state_dim,nFock,nFock))
    Uhubbard=HubbInt_RI(CC,CA)
    docc = local_doubly(CC,CA)
    dens = local_density(CC,CA)    
  end subroutine build_local_fock
  !
  subroutine bdecomp(i,ivec)
    integer :: ivec(:)         
    integer :: l,i,Ntot
    logical :: busy
    Ntot=size(ivec)
    !this is the configuration vector |1,..,Ns,Ns+1,...,Ntot>
    !obtained from binary decomposition of the state/number i\in 2^Ntot
    do l=0,Ntot-1
       busy=btest(i-1,l)
       ivec(l+1)=0
       if(busy)ivec(l+1)=1
    enddo
  end subroutine bdecomp

  subroutine build_c_operator(ispin,iorb,Cmat)
    integer                        :: ispin,iorb
    real(8),dimension(nFock,nFock) :: Cmat
    integer                        :: imp,i,j
    integer                        :: ib(state_dim)
    real(8)                        :: c_
    !build <j|C|i>
    imp = index(ispin,iorb)
    Cmat=0d0
    do i=1,nFock
       call bdecomp(i,ib)
       if(ib(imp)==0)cycle
       call c(imp,i,j,c_)
       Cmat(j,i)=c_
    enddo
    return
  end subroutine build_c_operator
  !+-------------------------------------------------------------------+
  !PURPOSE  : input state |i> of the basis and calculates |j>=Cm|i>
  !the sign of j has the phase convention
  !m labels the sites
  !+-------------------------------------------------------------------+
  subroutine c(m,i,j,sgn)
    integer :: ib(state_dim)
    integer :: i,j,m,km
    integer :: isg
    real(8) :: sgn
    call bdecomp(i,ib)
    if (ib(m)==0)then
       j=0
    else
       if(m==1)then
          j=i-1
       else
          km=sum(ib(1:m-1))
          isg=(-1)**km
          j=(i-2**(m-1))*isg
       endif
    endif
    sgn=dfloat(j)/dfloat(abs(j))
    j=abs(j)
  end subroutine c


  !
  subroutine build_cdg_operator(ispin,iorb,Cmat)
    integer                        :: ispin,iorb
    real(8),dimension(nFock,nFock) :: Cmat
    integer                        :: imp,i,j
    integer                        :: ib(state_dim)
    real(8)                        :: c_
    !build <j|C^+|i>
    imp = index(ispin,iorb)
    Cmat=0d0
    do i=1,nFock
       call bdecomp(i,ib)
       if(ib(imp)==1)cycle
       call cdg(imp,i,j,c_)
       Cmat(j,i)=c_
    enddo
    return
  end subroutine build_cdg_operator
  !+-------------------------------------------------------------------+
  !PURPOSE  : input state |i> of the basis and calculates |j>=Cm+|i>
  !the sign of j has the phase convention
  !m labels the sites
  !+-------------------------------------------------------------------+
  subroutine cdg(m,i,j,sgn)
    integer :: ib(state_dim)
    integer :: i,j,m,km
    integer :: isg
    real(8) :: sgn
    call bdecomp(i,ib)
    if (ib(m)==1)then
       j=0
    else
       if(m==1)then
          j=i+1
       else
          km=sum(ib(1:m-1))
          isg=(-1)**km
          j=(i+2**(m-1))*isg
       endif
    endif
    sgn=dfloat(j)/dfloat(abs(j))
    j=abs(j)
  end subroutine cdg


  function sz_rotate(fock_in) result(fock_out)
    integer                      :: fock_in,fock_out
    integer,dimension(state_dim) :: state_in,state_out    
    integer                      :: iorb,istate
    call bdecomp(fock_in,state_in)    
    do iorb=1,Norb
       state_out(iorb) = state_in(iorb+Norb)
       state_out(iorb+Norb) = state_in(iorb)
    end do
    fock_out=1
    do istate=0,state_dim-1
       fock_out = fock_out + state_out(istate+1)*2**istate
    end do
  end function sz_rotate


  subroutine get_spin_indep_states
    integer :: i_ind,i,iorb,istate
    integer :: tmp_search(nFock),tmp_target(nFock)
    integer :: ifock,isymm
    integer :: check_maps
    integer :: test_vec(state_dim)
    !+- get independent states under sz rotation symmetry -+!        
    tmp_search=0
    i_ind=0
    do ifock=1,nFock       
       tmp_target(ifock)=sz_rotate(ifock)              
       if(tmp_search(ifock).ge.0) then
          i_ind=i_ind+1
          tmp_search(ifock)=ifock
          if(tmp_target(ifock).ne.ifock) tmp_search(tmp_target(ifock)) = -1
       end if
    end do
    nFock_indep=i_ind
    allocate(fock_indep(nFock_indep),full2indep_fock(nFock),indep2full_fock(nFock_indep,2))
    i_ind=0
    do i=1,nFock
       if(tmp_search(i).ge.0) then
          i_ind=i_ind+1
          fock_indep(i_ind) = tmp_search(i)
       end if
    end do    
    do i_ind=1,nFock_indep
       full2indep_fock(fock_indep(i_ind))=i_ind       
       full2indep_fock(tmp_target(fock_indep(i_ind)))=i_ind
    end do
    do i_ind=1,nFock_indep       
       indep2full_fock(i_ind,1) = fock_indep(i_ind)
       indep2full_fock(i_ind,2) = tmp_target(fock_indep(i_ind))
    end do
    !+- check maps +-!
    do i_ind=1,nFock_indep
       do isymm=1,2
          check_maps=indep2full_fock(i_ind,isymm)
          if(i_ind /= full2indep_fock(check_maps)) stop "WRONG MAPS"
       end do
    end do
  end subroutine get_spin_indep_states




  !+- build local operators +-!

  ! Rotationally invariant Hubbard interaction !
  function HubbInt_RI(cc,ca) result(Oi)
    real(8),dimension(state_dim,nFock,nFock) :: cc,ca
    real(8),dimension(nFock,nFock) :: Oi,Id
    real(8),dimension(nFock,nFock) :: ni
    integer                        :: i,ispin,iorb,istate
    !
    Id=0.d0                 
    do i=1,nFock
       Id(i,i)=1.d0
    end do
    !
    ni=0.d0    
    do ispin=1,2
       do iorb=1,Norb
          istate=index(ispin,iorb)
          ni=ni+matmul(cc(istate,:,:),ca(istate,:,:))
       end do
    end do
    !
    Oi=matmul(ni-Id,ni-Id)
  end function HubbInt_RI
  
  ! Local density !
  function local_density(cc,ca) result(ni)
    real(8),dimension(state_dim,nFock,nFock) :: cc,ca
    real(8),dimension(nFock,nFock) :: Id
    real(8),dimension(state_dim,nFock,nFock) :: ni
    integer                        :: i,ispin,iorb,istate
    !
    do ispin=1,2
       do iorb=1,Norb
          istate=index(ispin,iorb)
          ni(istate,:,:)=matmul(cc(istate,:,:),ca(istate,:,:))
       end do
    end do
    !
  end function local_density

  
  function local_doubly(cc,ca) result(di)
    real(8),dimension(state_dim,nFock,nFock) :: cc,ca
    real(8),dimension(nFock,nFock) :: Id,tmp_up,tmp_dw
    real(8),dimension(Norb,nFock,nFock) :: di
    integer                        :: i,ispin,iorb,istate
    !
    do iorb=1,Norb
       !
       ispin=1
       istate=index(ispin,iorb)
       tmp_up = matmul(cc(istate,:,:),ca(istate,:,:))
       ispin=2
       istate=index(ispin,iorb)
       tmp_dw = matmul(cc(istate,:,:),ca(istate,:,:))
       !
       di(iorb,:,:) = matmul(tmp_up,tmp_dw)
       !
    end do
    !    
  end function local_doubly



  subroutine initialize_local_density(local_density)
    real(8),dimension(Norb),intent(inout) :: local_density
    logical                 :: IOfile
    integer                 :: unit,flen,i
    inquire(file="restart.density_seed.conf",exist=IOfile)
    if(IOfile) then
       flen=file_length("restart.density_seed.conf")
       unit=free_unit()
       open(unit,file="restart.density_seed.conf")
       write(*,*) 'reading denisty seed from file density_seed.conf'
       if(flen.eq.Norb) then
          !+- read from file -+!
          do i=1,flen
             read(unit,*) local_density(i)
          end do
       else
          !+- initialize in the usual way -+!
          do i=1,Norb
             local_density(i)=1.d0+(-1.d0)**dble(i)*0.5d0
          end do
       end if
       close(unit)
    else
       do i=1,Norb
          local_density(i)=1.d0+(-1.d0)**dble(i)*0.5d0
       end do       
    end if
  end subroutine initialize_local_density


  
END MODULE GZ_AUX_FUNX
