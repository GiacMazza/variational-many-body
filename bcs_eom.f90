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
  complex(8),dimension(Ns,Ns)                 :: Hk
  integer                                     :: is,js,ik,it,ks,kks,iis,jjs,iphi,jphi
  !

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
  do ik=1,Lk
     delta_t = delta_t + 0.5d0*(bcsWF(1,ik)-xi*bcsWF(2,ik))*wtk(ik)
  end do
  delta_t = Ubcs_t(it)*delta_t
  !
  do ik=1,Lk
     call get_Hk_t(Hk,ik,time)
     ekt = Hk(1,1)
     !
     bcsWF_dot(1,ik) = -2.d0*ekt*bcsWF(2,ik) - 2.d0*dimag(delta_t)*bcsWF(3,ik)
     !
     bcsWF_dot(2,ik) =  2.d0*ekt*bcsWF(1,ik) - 2.d0*dreal(delta_t)*bcsWF(3,ik)
     !
     bcsWF_dot(3,ik) =  2.d0*dreal(delta_t)*bcsWF(2,ik) + 2.d0*dimag(delta_t)*bcsWF(1,ik)
     !
  end do
  !
  call BCSwf_2_dynamicalVector(bcsWF_dot,f)
  !
end function BCS_equations_of_motion
