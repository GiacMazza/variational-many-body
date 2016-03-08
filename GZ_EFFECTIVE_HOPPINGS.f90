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
       end do
    end do
    !
  end function hopping_renormalization_normal
  !
END MODULE GZ_EFFECTIVE_HOPPINGS
