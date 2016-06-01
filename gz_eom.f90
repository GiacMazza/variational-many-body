function gz_equations_of_motion(time,y,Nsys) result(f)
  implicit none
  !inputs                                                                                                                                                                                                  
  integer                                     :: Nsys ! nr of equations                                                                                                                                    
  real(8)                                     :: time ! time variable                                                                                                                                      
  complex(8),dimension(Nsys)                  :: y    ! argument array                                                                                                                                     
  complex(8),dimension(Nsys)                  :: f    ! result 
  !
  complex(8),dimension(Ns,Ns,Lk)              :: slater,slater_dot
  complex(8),dimension(Nphi)                  :: gzproj,gzproj_dot
  complex(8),dimension(Nphi,Nphi)             :: Hproj
  integer                                     :: is,js,ik,it,ks,kks,iis,jjs,iphi,jphi
  complex(8),dimension(Ns,Ns)                 :: tmpHk,Rhop,slater_derivatives
  complex(8),dimension(Ns,Ns)                 :: vdm_natural
  real(8),dimension(Ns)                       :: vdm_diag
  !

  !HERE write the GZ EQUATIONS OF MOTION
  if(Nsys.ne.nDynamics) stop "wrong dimensions in the GZ_equations_of_motion"
  !
  call dynamicalVector_2_wfMatrix(y,slater,gzproj)
  do is=1,Ns
     do js=1,Ns
        vdm_natural(is,js) = trace_phi_basis(gzproj,phi_traces_basis_dens(is,js,:,:))
     end do
     vdm_diag(is) = dreal(vdm_natural(is,is))
  end do
  Rhop = hopping_renormalization_normal(gzproj,vdm_diag)
  !
  slater_dot=zero
  gzproj_dot=zero
  !
  it = t2it(time,tstep*0.5d0)
  !
  slater_derivatives=zero
  do ik=1,Lk
     !
     do is=1,Ns
        do js=1,Ns
           tmpHk(is,js)=zero           
           do iis=1,Ns
              do jjs=1,Ns                 
                 tmpHk(is,js) = tmpHk(is,js) + conjg(Rhop(iis,is))*Hk_tb_t(iis,jjs,ik,it)*Rhop(jjs,js)
              end do
           end do
        end do
     end do
     !
     do is=1,Ns
        do js=1,Ns
           slater_dot(is,js,ik) = zero
           do ks=1,Ns
              slater_dot(is,js,ik) = slater_dot(is,js,ik) + tmpHk(js,ks)*slater(is,ks,ik)
              slater_dot(is,js,ik) = slater_dot(is,js,ik) - tmpHk(ks,is)*slater(ks,js,ik)
           end do
        end do
     end do
     !     
     do is=1,Ns
        do js=1,Ns
           do ks=1,Ns
              do kks=1,Ns
                 slater_derivatives(is,js) = slater_derivatives(is,js) + conjg(Rhop(kks,ks))*Hk_tb_t(kks,is,ik,it)*slater(ks,js,ik)*wtk(ik)
              end do
           end do
        end do
     end do
     !
  end do
  slater_dot = -xi*slater_dot

  !+- create HLOC

  Uloc=Uloc_t(:,it)
  Ust =Ust_t(it)
  Jh=Jh_t(it)
  Jsf=Jsf_t(it)
  Jph=Jph_t(it)
  eLevels = eLevels_t(:,it)
  write(777,*) time,Uloc
  call get_local_hamiltonian_trace(eLevels)  
  call build_neqH_GZproj(Hproj,slater_derivatives,Rhop,vdm_diag)
  !
  do iphi=1,Nphi
     gzproj_dot(iphi) = zero
     do jphi=1,Nphi
        gzproj_dot(iphi) = gzproj_dot(iphi) + Hproj(iphi,jphi)*gzproj(jphi)
     end do
  end do
  gzproj_dot = -xi*gzproj_dot
  !
  call wfMatrix_2_dynamicalVector(slater_dot,gzproj_dot,f)
  !
end function gz_equations_of_motion





subroutine build_neqH_GZproj(Hproj,slater_derivatives,Rhop,n0)
  complex(8),dimension(Nphi,Nphi) :: Hproj
  complex(8),dimension(Ns,Ns)     :: slater_derivatives,Rhop
  real(8),dimension(Ns)           :: n0 
  integer                         :: is,js
  Hproj=zero
  Hproj=phi_traces_basis_Hloc
  do is=1,Ns
     do js=1,Ns
        Hproj = Hproj + slater_derivatives(is,js)*phi_traces_basis_Rhop(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))
        Hproj = Hproj + conjg(slater_derivatives(is,js))*phi_traces_basis_Rhop_hc(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))
        Hproj = Hproj + slater_derivatives(is,js)*Rhop(is,js)*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens(is,js,:,:)
        Hproj = Hproj + conjg(slater_derivatives(is,js)*Rhop(is,js))*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens_hc(is,js,:,:)
     end do
  end do
end subroutine build_neqH_GZproj
