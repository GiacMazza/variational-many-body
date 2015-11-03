MODULE GZ_PROJECTORS
  USE GZ_VARS_GLOBAL
  USE GZ_AUX_FUNX
  implicit none


  interface gz_trace
     module procedure gz_local_diag,gz_local_diag_
  end interface gz_trace

  public :: gz_trace
  public :: gz_hop_diag
  public :: gz_Rhop
  public :: gz_Rhop_dens

CONTAINS






  function gz_local_diag(phi,Oi) result(out)
    real(8),dimension(nFock)       :: phi  ! gz_projector     --> input
    real(8),dimension(nFock,nFock) :: Oi   ! local observable --> operator
    real(8)                        :: out
    integer                        :: i
    out=0.d0
    do i=1,Nfock
       out=out+phi(i)*phi(i)*Oi(i,i)
    end do

  end function gz_local_diag

  function gz_local_diag_(phi,Oi) result(out)
    complex(8),dimension(nFock)    :: phi  ! gz_projector     --> input
    real(8),dimension(nFock,nFock) :: Oi   ! local observable --> operator
    real(8)                        :: out
    integer                        :: i
    out=0.d0
    do i=1,Nfock
       out=out+abs(phi(i))*abs(phi(i))*Oi(i,i)
    end do
  end function gz_local_diag_



  function gz_local(phi,Oi) result(out)
    real(8),dimension(nFock,nFock) :: phi  ! gz_projector     --> input
    real(8),dimension(nFock,nFock) :: Oi   ! local observable --> operator
    real(8),dimension(nFock,nFock) :: tmp   ! local observable --> operator
    real(8)                        :: out    
    integer                        :: i,j,ispin,iorb,istate
    tmp=matmul(phi,Oi)
    tmp=matmul(tmp,phi)
    out=0.d0
    do i=1,nFock
       out=out+tmp(i,i)
    end do
  end function gz_local

  function gz_local_(phi,Oi) result(out)
    complex(8),dimension(nFock,nFock) :: phi  ! gz_projector     --> input
    complex(8),dimension(nFock,nFock) :: Oi   ! local observable --> operator
    complex(8),dimension(nFock,nFock) :: tmp,phi_   ! local observable --> operator
    real(8)                           :: out    
    integer                           :: i,j,ispin,iorb,istate
    do i=1,nFock
       do j=1,nFock
          phi_(i,j)=conjg(phi(j,i))
       end do
    end do
    tmp=matmul(phi_,Oi)
    tmp=matmul(tmp,phi)
    out=0.d0
    do i=1,nFock
       out=out+tmp(i,i)
    end do
  end function gz_local_


  function gz_hop_diag(phi,cc,ca) result(R)
    real(8),dimension(nFock) :: phi
    real(8),dimension(state_dim,nFock,nFock) :: ca,cc
    real(8),dimension(state_dim) :: ni
    real(8),dimension(state_dim,state_dim) :: R
    real(8),dimension(nFock,nFock) :: phi_,tmp
    integer                        :: iorb,ispin,istate,i
    tmp=0.d0
    phi_=0.d0
    do i=1,nFock
       phi_(i,i)=phi(i)
    end do
    do i=1,state_dim
       ni(i)=gz_local_diag(phi,dens(i,:,:))
    end do
    R=0.d0
    do ispin=1,2
       do iorb=1,Norb
          istate=index(ispin,iorb)
          tmp = matmul(phi_,cc(istate,:,:))
          tmp = matmul(ca(istate,:,:),tmp)
          tmp = matmul(phi_,tmp)
          do i=1,nFock
             R(istate,istate) = R(istate,istate) + tmp(i,i)
          end do
          R(istate,istate) = R(istate,istate)/sqrt(ni(istate)*(1.d0-ni(istate)))
       end do
    end do
  end function gz_hop_diag



  function gz_Rhop(phi,cc,ca) result(R)
    real(8),dimension(nFock)                 :: phi
    real(8),dimension(state_dim,nFock,nFock) :: ca,cc
    real(8),dimension(state_dim)             :: ni
    real(8),dimension(state_dim)             :: R
    real(8),dimension(nFock,nFock)           :: phi_,tmp
    integer                                  :: iorb,ispin,istate,i
    tmp=0.d0
    phi_=0.d0
    do i=1,nFock
       phi_(i,i)=phi(i)
    end do
    do i=1,state_dim
       ni(i)=gz_local_diag(phi,dens(i,:,:))
    end do
    R=0.d0
    do ispin=1,2
       do iorb=1,Norb
          istate=index(ispin,iorb)
          tmp = matmul(phi_,cc(istate,:,:))
          tmp = matmul(ca(istate,:,:),tmp)
          tmp = matmul(phi_,tmp)
          do i=1,nFock
             R(istate) = R(istate) + tmp(i,i)
          end do
          R(istate) = R(istate)/sqrt(ni(istate)*(1.d0-ni(istate)))
       end do
    end do
  end function gz_Rhop



  function gz_Rhop_dens(phi,ni,cc,ca) result(R)
    real(8),dimension(nFock)                 :: phi
    real(8),dimension(state_dim,nFock,nFock) :: ca,cc
    real(8),dimension(state_dim)             :: ni
    real(8),dimension(state_dim)             :: R
    real(8),dimension(nFock,nFock)           :: phi_,tmp
    real(8)                                  :: ntmp,eps
    integer                                  :: iorb,ispin,istate,i
    eps=1.d-8
    tmp=0.d0
    phi_=0.d0
    do i=1,nFock
       phi_(i,i)=phi(i)
    end do
    R=0.d0
    do ispin=1,2
       do iorb=1,Norb
          istate=index(ispin,iorb)
          if(ni(istate).gt.eps.and.ni(istate).lt.1.d0-eps) then
             ntmp=ni(istate)
          else
             if(ni(istate).lt.0.5d0) then
                ntmp=eps
             else
                ntmp=1.d0-eps
             end if
          end if
          istate=index(ispin,iorb)
          tmp = matmul(phi_,cc(istate,:,:))
          tmp = matmul(ca(istate,:,:),tmp)
          tmp = matmul(phi_,tmp)
          do i=1,nFock
             R(istate) = R(istate) + tmp(i,i)
          end do
          R(istate) = R(istate)/sqrt(ni(istate)*(1.d0-ni(istate)))
          !R(istate) = R(istate)/sqrt(ntmp*(1.d0-ntmp))
       end do
    end do
  end function gz_Rhop_dens  

  
END MODULE GZ_PROJECTORS
