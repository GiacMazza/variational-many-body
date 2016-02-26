MODULE GZ_EFFECTIVE_HOPPINGS
  USE GZ_VARS_GLOBAL
  USE GZ_MATRIX_BASIS
  implicit none

  public :: hopping_renormalization_normal

CONTAINS

  function hopping_renormalization_normal(phi,vdm_natural) result(Rhop)
    real(8),dimension(Ns)   :: vdm_natural
    complex(8),dimension(Nphi) :: phi
    complex(8),dimension(Ns,Ns) :: Rhop
    real(8) :: a
    complex(8) :: b
    integer :: is,js,iphi,jphi
    !
    do is=1,Ns
       do js=1,Ns
          a = sqrt(vdm_natural(js)*(1.d0-vdm_natural(js)))
          b = trace_phi_basis(phi,phi_traces_basis_Rhop(is,js,:,:))          
          Rhop(is,js) = b/a
          ! if(abs(a).gt.1.d-4) then
          !    Rhop(is,js) = b/a
          ! else
          !    Rhop = 0.d0
          !    if(is.eq.js) Rhop(is,js) = 1.d0
          ! end if
       end do
    end do

    
    !<DEBUG
    !write(*,*) b,a,abs(a),b/a
    ! do is=1,Ns
    !    write(*,'(10F7.4)') dreal(Rhop(is,:))
    ! end do
    ! write(*,*)
    !DEBUG>



  end function hopping_renormalization_normal

  !
END MODULE GZ_EFFECTIVE_HOPPINGS
