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
    integer :: is,js,iphi,jphi
    !
    do is=1,Ns
       do js=1,Ns
          Rhop(is,js) = trace_phi_basis(phi,phi_traces_basis_Rhop(is,js,iphi,jphi))
          Rhop(is,js) = Rhop(is,js)/sqrt(vdm_natural(js)*(1.d0-vdm_natural(js)))
       end do
    end do

    ! do is=1,Ns
    !    do js=1,Ns
    !       Rhop(is,js) = 0.d0
    !       do iphi=1,Nphi
    !          do jphi=1,Nphi
    !             Rhop(is,js) = &
    !                  Rhop(is,js) + phi(iphi)*phi(jphi)*phi_traces_basis_Rhop(is,js,iphi,jphi)/sqrt(vdm_natural(js)*(1.d0-vdm_natural(js)))                           
    !          end do
    !       end do
    !    end do
    ! end do

  end function hopping_renormalization_normal

  !
END MODULE GZ_EFFECTIVE_HOPPINGS
