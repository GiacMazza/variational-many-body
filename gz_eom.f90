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




function gz_equations_of_motion_superc(time,y,Nsys) result(f)
  implicit none
  !inputs
  integer                          :: Nsys ! nr of equations
  real(8)                          :: time ! time variable                                                                                                                                      
  complex(8),dimension(Nsys)       :: y    ! argument array                                                                                                                                     
  complex(8),dimension(Nsys)       :: f    ! result 
  !
  complex(8),dimension(2,Ns,Ns,Lk) :: slater_,slater_dot
  complex(8),dimension(3,Ns,Ns,Lk) :: slater
  complex(8),dimension(Nphi)       :: gzproj,gzproj_dot
  complex(8),dimension(Nphi,Nphi)  :: Hproj
  integer                          :: is,js,ik,it,ks,kks,iis,jjs,iphi,jphi
  complex(8),dimension(Ns,Ns)      :: tmpHk
  complex(8),dimension(2*Ns,2*Ns)  :: tmp
  complex(8),dimension(2,Ns,Ns)    :: slater_derivatives
  complex(8),dimension(Ns,Ns)      :: Rhop,Qhop
  complex(8),dimension(Ns,Ns)      :: vdm_natural
  real(8),dimension(Ns)            :: vdm_diag
  !

  !HERE write the GZ EQUATIONS OF MOTION
  if(Nsys.ne.nDynamics) stop "wrong dimensions in the GZ_equations_of_motion"
  !
  call dynamicalVector_2_wfMatrix_superc(y,slater_,gzproj)
  slater(1:2,:,:,:) = slater_(1:2,:,:,:)
  slater(3,:,:,:) = zero
  do is=1,Ns
     do js=1,Ns
        if(is.eq.js) slater(3,is,js,:) = 1.d0
        slater(3,is,js,:) = slater_(1,js,is,:)
     end do
  end do
  !
  do is=1,Ns
     do js=1,Ns
        vdm_natural(is,js) = trace_phi_basis(gzproj,phi_traces_basis_dens(is,js,:,:))
     end do
     vdm_diag(is) = dreal(vdm_natural(is,is))
  end do
  !
  Rhop = hopping_renormalization_normal(gzproj,vdm_diag)
  Qhop = hopping_renormalization_anomalous(gzproj,vdm_diag)
  !
  slater_dot=zero
  gzproj_dot=zero
  !
  it = t2it(time,tstep*0.5d0)
  !
  slater_derivatives = zero
  slater_dot = zero
  do ik=1,Lk
     !
     do is=1,Ns
        do js=1,Ns           
           slater_dot(1,is,js,ik) = zero
           do ks=1,Ns
              !
              do iis=1,Ns
                 do jjs=1,Ns

                    !
                    !
                    !

                    !tRR_{js,ks}SL1(is,ks) - tRR_{ks,is}SL1(ks,js) 
                    slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + &
                         conjg(Rhop(iis,js))*Hk_tb_t(iis,jjs,ik,it)*Rhop(jjs,ks)*slater(1,is,ks,ik)
                    slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) - &
                         conjg(Rhop(iis,ks))*Hk_tb_t(iis,jjs,ik,it)*Rhop(jjs,is)*slater(1,ks,js,ik)                    
                    !tRQ_{js,ks}SL2(is,ks) + tRQ_{ks,js}SL2(ks,is)
                    slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + &
                         conjg(Rhop(iis,js))*Hk_tb_t(iis,jjs,ik,it)*Qhop(jjs,ks)*slater(2,is,ks,ik)
                    slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + &
                         conjg(Rhop(iis,ks))*Hk_tb_t(iis,jjs,ik,it)*Qhop(jjs,js)*slater(2,ks,is,ik)                    
                    !tQR_{is,ks}CONJG_SL2(js,ks) + tQR_{ks,is}CONJG_SL2(ks,js)
                    slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + &
                         conjg(Qhop(iis,is))*Hk_tb_t(iis,jjs,ik,it)*Rhop(jjs,ks)*conjg(slater(2,js,ks,ik))
                    slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + &
                         conjg(Qhop(iis,ks))*Hk_tb_t(iis,jjs,ik,it)*Rhop(jjs,is)*conjg(slater(2,ks,js,ik))
                    !tQQ_{is,ks}SL1(ks,js) - tQQ_{ks,js}SL1(is,ks)
                    slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) + &
                         conjg(Qhop(iis,is))*Hk_tb_t(iis,jjs,ik,it)*Qhop(jjs,ks)*slater(1,ks,js,ik)
                    slater_dot(1,is,js,ik) = slater_dot(1,is,js,ik) - &
                         conjg(Qhop(iis,ks))*Hk_tb_t(iis,jjs,ik,it)*Qhop(jjs,js)*slater(1,is,ks,ik)                    

                    !
                    !
                    !

                    !tRR_{ks,js}SL2(ks,is) - tRR_{ks,is}SL2(ks,js)
                    slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + &
                         conjg(Rhop(iis,ks))*Hk_tb_t(iis,jjs,ik,it)*Rhop(jjs,js)*slater(2,ks,is,ik)
                    slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - &
                         conjg(Rhop(iis,ks))*Hk_tb_t(iis,jjs,ik,it)*Rhop(jjs,is)*slater(2,ks,js,ik)                    
                    !tQR_{js,ks}SL1(is,ks) - tQR_{is,ks}SL1(js,ks)
                    slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + &
                         conjg(Qhop(iis,js))*Hk_tb_t(iis,jjs,ik,it)*Rhop(jjs,ks)*slater(1,is,ks,ik)
                    slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - &
                         conjg(Qhop(iis,is))*Hk_tb_t(iis,jjs,ik,it)*Rhop(jjs,ks)*slater(1,js,ks,ik)                    
                    !tQR_{ks,js}SL3(ks,is) - tQR_{ks,is}SL3(ks,js)
                    slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + &
                         conjg(Qhop(iis,ks))*Hk_tb_t(iis,jjs,ik,it)*Rhop(jjs,js)*slater(3,ks,is,ik)
                    slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - &
                         conjg(Qhop(iis,ks))*Hk_tb_t(iis,jjs,ik,it)*Rhop(jjs,is)*slater(3,ks,js,ik)                    
                    !tQQ_{is,ks}SL2(ks,js) - tQQ_{js,ks}SL2(ks,is)
                    slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) + &
                         conjg(Qhop(iis,is))*Hk_tb_t(iis,jjs,ik,it)*Qhop(jjs,ks)*slater(2,ks,js,ik)
                    slater_dot(2,is,js,ik) = slater_dot(2,is,js,ik) - &
                         conjg(Qhop(iis,js))*Hk_tb_t(iis,jjs,ik,it)*Qhop(jjs,ks)*slater(2,ks,is,ik)                    

                    !
                    !
                    !

                 end do
              end do
              !
           end do
           !tmp(1:Ns,1:Ns) = <c^+ c>
           tmp(is,js) = slater(1,is,js,ik)     
           !tmp(1:Ns,Ns+1:2*Ns) = <c^+ c^+>
           tmp(is,js+Ns) = slater(2,is,js,ik)     
           !tmp(Ns+1:2*Ns,1:Ns) = <c c>
           tmp(is+Ns,js) = conjg(slater(2,js,is,ik))     
           !tmp(Ns+1:2*Ns,Ns+1:2*Ns) = <c c^+>
           tmp(is+Ns,js+Ns) = slater(3,is,js,ik)     
        end do
     end do
     !   

     do is=1,Ns
        do js=1,Ns
           do ks=1,Ns
              do kks=1,Ns
                 !
                 slater_derivatives(1,is,js) = slater_derivatives(1,is,js) + &
                      conjg(Rhop(kks,ks))*Hk_tb_t(kks,is,ik,it)*tmp(ks,js)*wtk(ik)    
                 !
                 slater_derivatives(1,is,js) = slater_derivatives(1,is,js) + &
                      conjg(Qhop(kks,ks))*Hk_tb_t(kks,is,ik,it)*tmp(ks+Ns,js)*wtk(ik)
                 !
                 slater_derivatives(2,is,js) = slater_derivatives(2,is,js) + &
                      conjg(Rhop(kks,ks))*Hk_tb_t(kks,is,ik,it)*tmp(ks,js+Ns)*wtk(ik) 
                 !
                 slater_derivatives(2,is,js) = slater_derivatives(2,is,js) + &
                      conjg(Qhop(kks,ks))*Hk_tb_t(kks,is,ik,it)*tmp(ks+Ns,js+Ns)*wtk(ik)
                 !
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
  call get_local_hamiltonian_trace(eLevels)  
  call build_neqH_GZproj_superc(Hproj,slater_derivatives,Rhop,Qhop,vdm_diag)
  !
  do iphi=1,Nphi
     gzproj_dot(iphi) = zero
     do jphi=1,Nphi
        gzproj_dot(iphi) = gzproj_dot(iphi) + Hproj(iphi,jphi)*gzproj(jphi)
     end do
  end do
  gzproj_dot = -xi*gzproj_dot
  ! !
  call wfMatrix_superc_2_dynamicalVector(slater_dot,gzproj_dot,f)
  !
end function gz_equations_of_motion_superc






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




subroutine build_neqH_GZproj_superc(Hproj,slater_derivatives,Rhop,Qhop,n0)
  complex(8),dimension(Nphi,Nphi) :: Hproj
  complex(8),dimension(2,Ns,Ns)   :: slater_derivatives
  complex(8),dimension(Ns,Ns)     :: Rhop,Qhop
  real(8),dimension(Ns)           :: n0 
  integer                         :: is,js
  Hproj=zero
  Hproj=phi_traces_basis_Hloc
  do is=1,Ns
     do js=1,Ns
        !
        Hproj = Hproj + slater_derivatives(1,is,js)*phi_traces_basis_Rhop(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))
        Hproj = Hproj + conjg(slater_derivatives(1,is,js))*phi_traces_basis_Rhop_hc(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))
        Hproj = Hproj + slater_derivatives(2,is,js)*phi_traces_basis_Qhop(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))
        Hproj = Hproj + conjg(slater_derivatives(2,is,js))*phi_traces_basis_Qhop_hc(is,js,:,:)/sqrt(n0(js)*(1.d0-n0(js)))
        !
        Hproj = Hproj + slater_derivatives(1,is,js)*Rhop(is,js)*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens(is,js,:,:)
        Hproj = Hproj + conjg(slater_derivatives(1,is,js)*Rhop(is,js))*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens_hc(is,js,:,:)
        Hproj = Hproj + slater_derivatives(2,is,js)*Qhop(is,js)*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens(is,js,:,:)
        Hproj = Hproj + conjg(slater_derivatives(2,is,js)*Qhop(is,js))*(2.d0*n0(js)-1.d0)/2.d0/(n0(js)*(1.d0-n0(js)))*phi_traces_basis_dens_hc(is,js,:,:)

     end do
  end do
end subroutine build_neqH_GZproj_superc
