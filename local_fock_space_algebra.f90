subroutine build_local_fock_algebra
  integer                            :: fock_state,ifock
  integer                            :: i,iup,jup,idw,iup_,idw_,np,nup_
  integer                            :: ispin,iorb,istate
  integer                            :: nup,ndw,irnd
  real(8)                            :: rnd    
  !+- allocate and build creation/anhhilation operators -+!
  allocate(CC(Ns,nFock,nFock),CA(Ns,nFock,nFock))
  do ispin=1,2
     do iorb=1,Norb
        istate=index(ispin,iorb)
        call build_c_operator(ispin,iorb,CA(istate,:,:))
        call build_cdg_operator(ispin,iorb,CC(istate,:,:))          
     end do
  end do
  !
contains
  !
  subroutine build_c_operator(ispin,iorb,Cmat)
    integer                        :: ispin,iorb
    real(8),dimension(nFock,nFock) :: Cmat
    integer                        :: imp,i,j
    integer                        :: ib(Ns)
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
    integer :: ib(Ns)
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
    integer                        :: ib(Ns)
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
    integer :: ib(Ns)
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
end subroutine build_local_fock_algebra
