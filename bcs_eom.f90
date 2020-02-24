function BCS_equations_of_motion(time,y,Nsys) result(f)
  implicit none
  !inputs                                                                                                                                                                                                  
  integer                                     :: Nsys ! nr of equations                                                                                                                                   
  real(8)                                     :: time ! time variable                                                                                                                                      
  complex(8),dimension(Nsys)                  :: y    ! argument array                                                                                                                                     
  complex(8),dimension(Nsys)                  :: f    ! result 
  !
  complex(8),dimension(3,Lk)               :: bcsWF,bcsWF_dot
  complex(8) :: ekt
  complex(8) :: delta_t
  complex(8),dimension(Lk) :: delta_tk
  real(8),dimension(Lk) :: n_tk
  complex(8),dimension(Ns,Ns)                 :: Hk
  complex(8) ::  nhh_dot 
  integer                                     :: is,js,ik,it,ks,kks,iis,jjs,iphi,jphi
  !
  real(8) :: n_t

  !HERE write the GZ EQUATIONS OF MOTION
  if(Nsys.ne.3*LK) stop "wrong dimensions in the BCS equations of motion"
  !
  call dynamicalVector_2_BCSwf(y,BCSwf)
  !
  BCSwf_dot=zero
  !
  it = t2it(time,tstep*0.5d0)
  !
  delta_t = zero
  n_t=0.d0
  do ik=1,Lk
     delta_t = delta_t + 0.5d0*(bcsWF(1,ik)+xi*bcsWF(2,ik))*wtk(ik)
     delta_tk(ik) = 0.5d0*(bcsWF(1,ik)+xi*bcsWF(2,ik))
     n_t = n_t +0.5d0*(bcsWF(3,ik)+1.d0)*wtk(ik)
     n_tk(ik) = 0.5d0*(bcsWF(3,ik)+1.d0)
  end do
  delta_t = (Ubcs_t(it)-xi*kdiss_t(it))*delta_t 
  !
  do ik=1,Lk
     call get_Hk_t(Hk,ik,time)
     ekt = Hk(1,1)
     !
     bcsWF_dot(1,ik) = -2.d0*ekt*bcsWF(2,ik) + 2.d0*dimag(delta_t)*bcsWF(3,ik)
     bcsWF_dot(1,ik) = bcsWF_dot(1,ik) - kdiss_t(it)*n_t*bcsWF(1,ik) 
     !
     bcsWF_dot(2,ik) =  2.d0*ekt*bcsWF(1,ik) - 2.d0*dreal(delta_t)*bcsWF(3,ik)
     bcsWF_dot(2,ik) = bcsWF_dot(2,ik) - kdiss_t(it)*n_t*bcsWF(2,ik) 
     !
     bcsWF_dot(3,ik) =  2.d0*dreal(delta_t)*bcsWF(2,ik) - 2.d0*dimag(delta_t)*bcsWF(1,ik)
     bcsWF_dot(3,ik) = bcsWF_dot(3,ik) - kdiss_t(it)*n_t*(bcsWF(3,ik)+1.d0)
     !
     !+ add here the non-hermitean part -+!
     nhh_dot = 2.d0*n_t*delta_tk(ik)*n_tk(ik)
     nhh_dot = nhh_dot +2.d0*conjg(delta_t)*delta_tk(ik)**2.d0+delta_t*n_tk(ik)**2.d0
     !nhh_dot = nhh_dot + 2.d0*(conjg(delta_t)*delta_tk(ik)**2.d0 - delta_t*n_tk(ik)**2.d0)
     !
     bcsWF_dot(1,ik) = bcsWF_dot(1,ik) + (1.d0-a_nhh)*kdiss_t(it)*2.d0*dreal(nhh_dot)
     bcsWF_dot(2,ik) = bcsWF_dot(2,ik) + (1.d0-a_nhh)*kdiss_t(it)*2.d0*dimag(nhh_dot)
     !
     nhh_dot = n_t*(abs(delta_tk(ik))**2.d0-n_tk(ik)**2.d0)
     nhh_dot = nhh_dot-4.d0*dreal(delta_t*conjg(delta_tk(ik))*n_tk(ik))
     bcsWF_dot(3,ik) = bcsWF_dot(3,ik) - 2.d0*(1.d0-a_nhh)*kdiss_t(it)*nhh_dot
     !
     ! nhh_dot = 2.d0*(abs(delta_tk(ik))**2.d0-n_tk(ik)**2.d0)
     ! nhh_dot = nhh_dot + 4.d0*2.d0*dreal(delta_t*n_tk(ik)*conjg(delta_tk(ik)))
     ! bcsWF_dot(3,ik) = bcsWF_dot(3,ik) - (1.d0-a_nhh)*kdiss_t(it)*nhh_dot
     !
  end do
  !
  call BCSwf_2_dynamicalVector(bcsWF_dot,f)
  !
end function BCS_equations_of_motion
