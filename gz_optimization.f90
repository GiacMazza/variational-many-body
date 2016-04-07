subroutine R_VDM_free_opt_function(x,delta_opt,i)     
  real(8),dimension(:),intent(in) :: x      
  real(8),intent(out)             :: delta_opt
  integer,optional,intent(in)     :: i

  real(8)                         :: E_Hstar,E_Hloc
  complex(8),dimension(Nphi)      :: GZvect

  complex(8),dimension(Ns,Ns) :: Rhop,Rhop_GZproj
  complex(8),dimension(Ns,Ns) :: slater_derivatives
  real(8),dimension(Ns,Ns)    :: slater_lgr_multip
  real(8),dimension(Ns,Ns)    :: proj_lgr_multip
  complex(8),dimension(Ns,Ns) :: Rhop_lgr_multip
  real(8),dimension(Ns,Ns)    :: n0_slater,n0_GZproj
  real(8),dimension(Ns)       :: n0_diag
  real(8)                     :: tmp
  integer                     :: imap,is,js
  integer :: Nslater_lgr,Nrhop,Nproj_lgr
  !
  Nslater_lgr = Nopt_diag + Nopt_odiag
  NRhop = Nopt_diag + Nopt_odiag
  Nproj_lgr = Nopt_odiag
  !
  slater_lgr_multip=0.d0
  proj_lgr_multip=0.d0
  Rhop = zero
  do is=1,Ns
     do js=1,Ns
        imap = opt_map(is,js)
        if(imap.gt.0) then
           slater_lgr_multip(is,js)=x(imap)
           Rhop(is,js) = x(imap + Nslater_lgr) 
           Rhop(is,js) = Rhop(is,js) + xi*x(imap + Nslater_lgr + NRhop)
           if(is.ne.js) proj_lgr_multip(is,js) = x(imap + Nslater_lgr + 2*NRhop)
        end if
     end do
  end do
  !
  call slater_minimization_fixed_lgr(Rhop,slater_lgr_multip,E_Hstar,n0_slater,slater_derivatives)
  !
  do is=1,Ns
     do js=1,Ns
        if(abs(n0_slater(js,js)).gt.1.d-10.and.abs(n0_slater(js,js)).lt.1.d0-1.d-10) then
           Rhop_lgr_multip(is,js) = -1.d0*slater_derivatives(is,js)/sqrt(n0_slater(js,js)*(1.d0-n0_slater(js,js)))
        else
           Rhop_lgr_multip(is,js) = 0.d0
        end if
     end do
     !
  end do
  !
  do is=1,Ns
     tmp = 0.d0
     do js=1,Ns
        tmp = tmp + dreal(Rhop_lgr_multip(is,js)*Rhop(is,js))
     end do
     if(abs(n0_slater(is,is)).gt.1.d-10.and.abs(n0_slater(is,is)).lt.1.d0-1.d-10) then
        proj_lgr_multip(is,is) = -1.d0*slater_lgr_multip(is,is) + & 
             0.5d0*(1.d0-2.d0*n0_slater(is,is))/sqrt(n0_slater(is,is)*(1.d0-n0_slater(is,is)))*tmp 
     else
        proj_lgr_multip(is,is) = -1.d0*slater_lgr_multip(is,is)
     end if
  end do
  !
  do is=1,Ns
     n0_diag(is) = n0_slater(is,is)
  end do
  !
  call gz_proj_minimization_fixed_lgr(n0_diag,proj_lgr_multip,Rhop_lgr_multip,E_Hloc,GZvect)
  !
  n0_GZproj = 0.d0
  do is=1,Ns
     do js=1,Ns
        n0_GZproj(is,js)   = trace_phi_basis(GZvect,phi_traces_basis_dens(is,js,:,:))
        Rhop_GZProj(is,js) = trace_phi_basis(GZvect,phi_traces_basis_Rhop(is,js,:,:))
     end do
  end do
  !
  delta_opt = 0.d0
  do is=1,Ns
     do js=1,Ns
        if(is.ne.js) then
           delta_opt = delta_opt +  abs(n0_slater(is,js))**2.d0
           delta_opt = delta_opt +  abs(n0_GZproj(is,js))**2.d0
        else
           delta_opt = delta_opt + abs(n0_slater(is,js)-n0_GZproj(is,js))**2.d0
        end if
        delta_opt = delta_opt + (Rhop_GZproj(is,js)-Rhop(is,js)*sqrt(n0_diag(js)*(1.d0-n0_diag(js))))*conjg(Rhop_GZproj(is,js)-Rhop(is,js)*sqrt(n0_diag(js)*(1.d0-n0_diag(js)))) 
     end do
  end do
  if(optimization_flag) then
     !+- store final informations to global variables -+!                   
     GZ_vector      = GZvect
     GZ_opt_energy  = E_Hstar+E_Hloc
     GZ_opt_kinetic = E_Hstar
     GZ_opt_Eloc    = E_Hloc
     GZ_opt_slater_lgr = slater_lgr_multip
  end if
  !<DEBUG
  !write(*,*) delta_opt,x
  ! write(124,*) n0_slater
  ! write(125,*) n0_GZproj
  ! write(126,*) dreal(Rhop_GZproj(1,1)),dreal(Rhop_GZproj(2,2)),dimag(Rhop_GZproj(1,1)),dimag(Rhop_GZproj(2,2))
  ! write(127,*) dreal(Rhop_GZproj(1,1)),dreal(Rhop_GZproj(2,2)),dreal(Rhop(1,1)*sqrt(n0_diag(1)*(1.d0-n0_diag(1)))),dreal(Rhop(2,2)*sqrt(n0_diag(2)*(1.d0-n0_diag(2))))    
  !DEBUG>
end subroutine R_VDM_free_opt_function
!
function R_VDM_free_zeros(x)  result(Fout)
  real(8),dimension(:) :: x      
  real(8),dimension(size(x))      :: Fout
  real(8)                         :: E_Hstar,E_Hloc
  complex(8),dimension(Nphi)      :: GZvect

  complex(8),dimension(Ns,Ns) :: Rhop,Rhop_GZproj
  complex(8),dimension(Ns,Ns) :: slater_derivatives
  real(8),dimension(Ns,Ns)    :: slater_lgr_multip
  real(8),dimension(Ns,Ns)    :: proj_lgr_multip
  complex(8),dimension(Ns,Ns) :: Rhop_lgr_multip
  real(8),dimension(Ns,Ns)    :: n0_slater,n0_GZproj
  real(8),dimension(Ns)       :: n0_diag
  real(8)                     :: tmp
  complex(8)                     :: tmpR
  integer                     :: imap,is,js
  integer :: Nslater_lgr,Nrhop,Nproj_lgr

  !
  Nslater_lgr = Nopt_diag + Nopt_odiag
  NRhop = Nopt_diag + Nopt_odiag
  Nproj_lgr = Nopt_odiag
  !
  slater_lgr_multip=0.d0
  proj_lgr_multip=0.d0
  Rhop = zero
  do is=1,Ns
     do js=1,Ns
        imap = opt_map(is,js)
        if(imap.gt.0) then
           slater_lgr_multip(is,js)=x(imap)
           Rhop(is,js) = x(imap + Nslater_lgr) 
           Rhop(is,js) = Rhop(is,js) + xi*x(imap + Nslater_lgr + NRhop)
           if(is.ne.js) proj_lgr_multip(is,js) = x(imap + Nslater_lgr + 2*NRhop)
        end if
     end do
  end do
  !

  !<TMP DEBUG
  do is=1,NS
     write(545,*) slater_lgr_multip(is,:)
  end do
  write(545,*)
  !TMP DEBUG>


  call slater_minimization_fixed_lgr(Rhop,slater_lgr_multip,E_Hstar,n0_slater,slater_derivatives)
  !
  do is=1,Ns
     do js=1,Ns
        if(abs(n0_slater(js,js)).gt.1.d-10.and.abs(n0_slater(js,js)).lt.1.d0-1.d-10) then
           Rhop_lgr_multip(is,js) = -1.d0*slater_derivatives(is,js)/sqrt(n0_slater(js,js)*(1.d0-n0_slater(js,js)))
        else
           Rhop_lgr_multip(is,js) = 0.d0
        end if
     end do
     !
  end do
  !
  do is=1,Ns
     tmp = 0.d0
     do js=1,Ns
        tmp = tmp + dreal(Rhop_lgr_multip(is,js)*Rhop(is,js))
     end do
     if(abs(n0_slater(is,is)).gt.1.d-10.and.abs(n0_slater(is,is)).lt.1.d0-1.d-10) then
        proj_lgr_multip(is,is) = -0.5d0*slater_lgr_multip(is,is) + &  !+- YO, THIS 0.5 is VERY IMPORTANT BITCH!
             0.5d0*(1.d0-2.d0*n0_slater(is,is))/sqrt(n0_slater(is,is)*(1.d0-n0_slater(is,is)))*tmp 
     else
        proj_lgr_multip(is,is) = -1.d0*slater_lgr_multip(is,is)
     end if
  end do
  !
  do is=1,Ns
     n0_diag(is) = n0_slater(is,is)
  end do
  !
  call gz_proj_minimization_fixed_lgr(n0_diag,proj_lgr_multip,Rhop_lgr_multip,E_Hloc,GZvect)
  !
  n0_GZproj = 0.d0
  do is=1,Ns
     do js=1,Ns
        n0_GZproj(is,js)   = trace_phi_basis(GZvect,phi_traces_basis_dens(is,js,:,:))
        Rhop_GZProj(is,js) = trace_phi_basis(GZvect,phi_traces_basis_Rhop(is,js,:,:))
     end do
  end do
  !
  ! delta_opt = 0.d0
  ! do is=1,Ns
  !    do js=1,Ns
  !       if(is.ne.js) then
  !          delta_opt = delta_opt +  abs(n0_slater(is,js))**2.d0
  !          delta_opt = delta_opt +  abs(n0_GZproj(is,js))**2.d0
  !       else
  !          delta_opt = delta_opt + abs(n0_slater(is,js)-n0_GZproj(is,js))**2.d0
  !       end if
  !       delta_opt = delta_opt + (Rhop_GZproj(is,js)-Rhop(is,js)*sqrt(n0_diag(js)*(1.d0-n0_diag(js))))*conjg(Rhop_GZproj(is,js)-Rhop(is,js)*sqrt(n0_diag(js)*(1.d0-n0_diag(js)))) 
  !    end do
  ! end do



  do is=1,Ns
     do js=1,Ns
        imap = opt_map(is,js)
        if(imap.gt.0) then
           if(is.ne.js) then
              Fout(imap) = n0_slater(is,js)
           else
              Fout(imap) = n0_slater(is,js) - n0_GZproj(is,js)
           end if
           !
           tmpR = Rhop_GZproj(is,js) - Rhop(is,js)*sqrt(n0_diag(js)*(1.d0-n0_diag(js)))
           Fout(imap+Nslater_lgr) = dreal(tmpR)
           Fout(imap+Nslater_lgr+NRhop) = dimag(tmpR)
           !
           if(is.ne.js) then
              Fout(imap+Nslater_lgr+2*NRhop) =   n0_GZproj(is,js)
           end if
           !
        end if
     end do
  end do




  if(optimization_flag) then
     !+- store final informations to global variables -+!                   
     GZ_vector      = GZvect
     GZ_opt_energy  = E_Hstar+E_Hloc
     GZ_opt_kinetic = E_Hstar
     GZ_opt_Eloc    = E_Hloc
     GZ_opt_slater_lgr = slater_lgr_multip
  end if
  !<DEBUG
  !  write(*,*) delta_opt,x
  ! write(124,*) n0_slater
  ! write(125,*) n0_GZproj
  ! write(126,*) dreal(Rhop_GZproj(1,1)),dreal(Rhop_GZproj(2,2)),dimag(Rhop_GZproj(1,1)),dimag(Rhop_GZproj(2,2))
  ! write(127,*) dreal(Rhop_GZproj(1,1)),dreal(Rhop_GZproj(2,2)),dreal(Rhop(1,1)*sqrt(n0_diag(1)*(1.d0-n0_diag(1)))),dreal(Rhop(2,2)*sqrt(n0_diag(2)*(1.d0-n0_diag(2))))    
  !DEBUG>
end function R_VDM_free_zeros
!
!
!






!
!
function R_Q_VDM_free_zeros_superc(x)  result(Fout)
  real(8),dimension(:) :: x      
  real(8),dimension(size(x))      :: Fout
  real(8),dimension(size(x))      :: tmp_x
  real(8)                         :: E_Hstar,E_Hloc
  complex(8),dimension(Nphi)      :: GZvect

  complex(8),dimension(Ns,Ns)   :: Rhop,Rhop_GZproj
  complex(8),dimension(Ns,Ns)   :: Qhop,Qhop_GZproj
  complex(8),dimension(2,Ns,Ns) :: slater_derivatives
  complex(8),dimension(2,Ns,Ns)    :: slater_lgr_multip
  complex(8),dimension(2,Ns,Ns)    :: proj_lgr_multip
  complex(8),dimension(Ns,Ns) :: Rhop_lgr_multip
  complex(8),dimension(Ns,Ns) :: Qhop_lgr_multip
  !
  complex(8),dimension(2,Ns,Ns)    :: n0_slater,n0_GZproj
  real(8),dimension(Ns)       :: n0_diag
  real(8)                     :: tmp
  complex(8)                     :: tmpR
  integer                     :: imap,is,js,i,ix,i0
  integer :: Nslater_lgr,NQhop,NRhop,Nproj_lgr,Nopt

  complex(8),dimension(:),allocatable :: dump_input
  complex(8),dimension(:),allocatable :: dump_output
  !
  complex(8),dimension(Ns,Ns) :: Rhop_out,Qhop_out
  complex(8),dimension(2,Ns,Ns) :: slater_out,GZproj_out


  !
  NRhop = 2*NRhop_opt
  NQhop = 2*NQhop_opt
  !
  Nslater_lgr = 2*(Nvdm_NC_opt+Nvdm_AC_opt)
  Nproj_lgr = 2*(Nvdm_NCoff_opt+Nvdm_AC_opt) 
  !
  Nopt = NRhop + NQhop + Nslater_lgr + Nproj_lgr
  if(size(x).ne.Nopt) stop "what the fuck are you doing!?!?!?"
  call dump2mats_superc(x,Rhop,Qhop,slater_lgr_multip,proj_lgr_multip)
  !
  write(*,*) "DENTRO R_Q_VDM"
  do is=1,Nopt
     write(*,*) is,x(is)
  end do

  do is=1,Ns
     write(*,*) Rhop(is,:)
  end do
  write(*,*) 
  do is=1,Ns
     write(*,*) Qhop(is,:)
  end do
  write(*,*) 
  do is=1,Ns
     write(*,*) slater_lgr_multip(1,is,:)
  end do
  write(*,*) 
  do is=1,Ns
     write(*,*) slater_lgr_multip(2,is,:)
  end do
  write(*,*) 
  do is=1,Ns
     write(*,*) proj_lgr_multip(1,is,:)
  end do
  write(*,*) 
  do is=1,Ns
     write(*,*) proj_lgr_multip(2,is,:)
  end do
  write(*,*)   
  !
  call slater_minimization_fixed_lgr_superc(Rhop,Qhop,slater_lgr_multip,E_Hstar,n0_slater,slater_derivatives)
  !
  do is=1,Ns
     do js=1,Ns
        Rhop_lgr_multip(is,js) = -1.d0*slater_derivatives(1,is,js)/sqrt(n0_slater(1,js,js)*(1.d0-n0_slater(1,js,js)))
        Qhop_lgr_multip(is,js) = -1.d0*slater_derivatives(2,is,js)/sqrt(n0_slater(1,js,js)*(1.d0-n0_slater(1,js,js)))
     end do
  end do
  !
  !
  do is=1,Ns
     tmp = 0.d0
     do js=1,Ns
        tmp = tmp + dreal(Rhop_lgr_multip(is,js)*Rhop(is,js)) + dreal(Qhop_lgr_multip(is,js)*Qhop(is,js))
     end do
     proj_lgr_multip(1,is,is) = -0.5*slater_lgr_multip(1,is,is) + & 
          0.5d0*(1.d0-2.d0*n0_slater(1,is,is))/sqrt(n0_slater(1,is,is)*(1.d0-n0_slater(1,is,is)))*tmp 
  end do
  !
  do is=1,Ns
     n0_diag(is) = n0_slater(1,is,is)
  end do
  !  
  call gz_proj_minimization_fixed_lgr_superc_(n0_diag,proj_lgr_multip,Rhop_lgr_multip,Qhop_lgr_multip,E_Hloc,GZvect)

  n0_GZproj = 0.d0
  do is=1,Ns
     do js=1,Ns
        n0_GZproj(1,is,js) = trace_phi_basis(GZvect,phi_traces_basis_dens(is,js,:,:))
        n0_GZproj(2,is,js) = trace_phi_basis(GZvect,phi_traces_basis_dens_anomalous(is,js,:,:))
        Rhop_GZProj(is,js) = trace_phi_basis(GZvect,phi_traces_basis_Rhop(is,js,:,:))
        Qhop_GZProj(is,js) = trace_phi_basis(GZvect,phi_traces_basis_Qhop(is,js,:,:))
     end do
  end do
  !  
  do is=1,Ns
     do js=1,Ns
        Rhop_out(is,js) = Rhop_GZproj(is,js)-Rhop(is,js)*sqrt(n0_diag(js)*(1.d0-n0_diag(js)))
        Qhop_out(is,js) = Qhop_GZproj(is,js)-Qhop(is,js)*sqrt(n0_diag(js)*(1.d0-n0_diag(js)))
        if(is.eq.js) then
           slater_out(1,is,js) = n0_slater(1,is,js) - n0_GZproj(1,is,js)
        else
           slater_out(1,is,js) = n0_slater(1,is,js)
        end if
        slater_out(2,is,js) = n0_slater(2,is,js)
        GZproj_out(1,is,js) = n0_GZproj(1,is,js)
        GZproj_out(2,is,js) = n0_GZproj(2,is,js)
     end do
  end do
  !<DEBUG
  write(*,*) "DEBUG_N0"
  do is=1,Ns
     write(*,'(10F8.4)') GZproj_out(1,is,:)
  end do
  write(*,*)
  do is=1,Ns
     write(*,'(10F8.4)') GZproj_out(2,is,:)
  end do
  !DEBUG>


  call dump2vec_superc(tmp_x,Rhop_out,Qhop_out,slater_out,GZproj_out)
  Fout=tmp_x
  !<DEBUG
  do is=1,Nopt
     write(*,*) Fout(is)
  end do
  write(*,*) '!+------------------+!'
  !DEBUG>
  !
  if(optimization_flag) then
     !+- store final informations to global variables -+!                   
     GZ_vector      = GZvect
     GZ_opt_energy  = E_Hstar+E_Hloc
     GZ_opt_kinetic = E_Hstar
     GZ_opt_Eloc    = E_Hloc
     !+- here I should modify this
     GZ_opt_slater_lgr_superc = slater_lgr_multip
     !<TMP DEBUG
     ! if(allocated(GZ_opt_slater_lgr_superc)) write(*,*) 'allocated prima'
     ! write(*,*) '?????',size(GZ_opt_slater_lgr_superc,1),size(GZ_opt_slater_lgr_superc,2),size(GZ_opt_slater_lgr_superc,3)
     ! write(*,*) GZ_opt_slater_lgr_superc(1,1,1)
     ! if(allocated(GZ_opt_slater_lgr_superc)) write(*,*) 'allocated dopo'
     !TMP DEBUG>
  end if


end function R_Q_VDM_free_zeros_superc



! subroutine R_VDM_free_opt_function(x,delta_opt,i)     
!   real(8),dimension(:),intent(in) :: x      
!   real(8),intent(out)             :: delta_opt
!   integer,optional,intent(in)     :: i




function R_Q_VDM_free_opt_superc(x)  result(Fout)
  real(8),dimension(:) :: x      
  real(8)      :: Fout
  real(8),dimension(size(x))      :: tmp_x
  real(8)                         :: E_Hstar,E_Hloc
  complex(8),dimension(Nphi)      :: GZvect

  complex(8),dimension(Ns,Ns)   :: Rhop,Rhop_GZproj
  complex(8),dimension(Ns,Ns)   :: Qhop,Qhop_GZproj
  complex(8),dimension(2,Ns,Ns) :: slater_derivatives
  complex(8),dimension(2,Ns,Ns)    :: slater_lgr_multip
  complex(8),dimension(2,Ns,Ns)    :: proj_lgr_multip
  complex(8),dimension(Ns,Ns) :: Rhop_lgr_multip
  complex(8),dimension(Ns,Ns) :: Qhop_lgr_multip
  !
  complex(8),dimension(2,Ns,Ns)    :: n0_slater,n0_GZproj
  real(8),dimension(Ns)       :: n0_diag
  real(8)                     :: tmp
  complex(8)                     :: tmpR
  integer                     :: imap,is,js,i,ix,i0
  integer :: Nslater_lgr,NQhop,NRhop,Nproj_lgr,Nopt

  complex(8),dimension(:),allocatable :: dump_input
  complex(8),dimension(:),allocatable :: dump_output
  !
  complex(8),dimension(Ns,Ns) :: Rhop_out,Qhop_out
  complex(8),dimension(2,Ns,Ns) :: slater_out,GZproj_out


  !
  NRhop = 2*NRhop_opt
  NQhop = 2*NQhop_opt
  !
  Nslater_lgr = 2*(Nvdm_NC_opt+Nvdm_AC_opt)
  Nproj_lgr = 2*(Nvdm_NCoff_opt+Nvdm_AC_opt) 
  !
  Nopt = NRhop + NQhop + Nslater_lgr + Nproj_lgr
  if(size(x).ne.Nopt) stop "what the fuck are you doing!?!?!?"
  call dump2mats_superc(x,Rhop,Qhop,slater_lgr_multip,proj_lgr_multip)
  !
  write(*,*) "DENTRO R_Q_VDM"
  do is=1,Nopt
     write(*,*) is,x(is)
  end do

  do is=1,Ns
     write(*,*) Rhop(is,:)
  end do
  write(*,*) 
  do is=1,Ns
     write(*,*) Qhop(is,:)
  end do
  write(*,*) 
  do is=1,Ns
     write(*,*) slater_lgr_multip(1,is,:)
  end do
  write(*,*) 
  do is=1,Ns
     write(*,*) slater_lgr_multip(2,is,:)
  end do
  write(*,*) 
  do is=1,Ns
     write(*,*) proj_lgr_multip(1,is,:)
  end do
  write(*,*) 
  do is=1,Ns
     write(*,*) proj_lgr_multip(2,is,:)
  end do
  write(*,*)   
  !
  call slater_minimization_fixed_lgr_superc(Rhop,Qhop,slater_lgr_multip,E_Hstar,n0_slater,slater_derivatives)
  !
  do is=1,Ns
     do js=1,Ns
        Rhop_lgr_multip(is,js) = -1.d0*slater_derivatives(1,is,js)/sqrt(n0_slater(1,js,js)*(1.d0-n0_slater(1,js,js)))
        Qhop_lgr_multip(is,js) = -1.d0*slater_derivatives(2,is,js)/sqrt(n0_slater(1,js,js)*(1.d0-n0_slater(1,js,js)))
     end do
  end do
  !
  !
  do is=1,Ns
     tmp = 0.d0
     do js=1,Ns
        tmp = tmp + dreal(Rhop_lgr_multip(is,js)*Rhop(is,js)) + dreal(Qhop_lgr_multip(is,js)*Qhop(is,js))
     end do
     proj_lgr_multip(1,is,is) = -0.5*slater_lgr_multip(1,is,is) + & 
          0.5d0*(1.d0-2.d0*n0_slater(1,is,is))/sqrt(n0_slater(1,is,is)*(1.d0-n0_slater(1,is,is)))*tmp 
  end do
  !
  do is=1,Ns
     n0_diag(is) = n0_slater(1,is,is)
  end do
  !  
  call gz_proj_minimization_fixed_lgr_superc_(n0_diag,proj_lgr_multip,Rhop_lgr_multip,Qhop_lgr_multip,E_Hloc,GZvect)

  n0_GZproj = 0.d0
  do is=1,Ns
     do js=1,Ns
        n0_GZproj(1,is,js) = trace_phi_basis(GZvect,phi_traces_basis_dens(is,js,:,:))
        n0_GZproj(2,is,js) = trace_phi_basis(GZvect,phi_traces_basis_dens_anomalous(is,js,:,:))
        Rhop_GZProj(is,js) = trace_phi_basis(GZvect,phi_traces_basis_Rhop(is,js,:,:))
        Qhop_GZProj(is,js) = trace_phi_basis(GZvect,phi_traces_basis_Qhop(is,js,:,:))
     end do
  end do
  !  
  do is=1,Ns
     do js=1,Ns
        Rhop_out(is,js) = Rhop_GZproj(is,js)-Rhop(is,js)*sqrt(n0_diag(js)*(1.d0-n0_diag(js)))
        Qhop_out(is,js) = Qhop_GZproj(is,js)-Qhop(is,js)*sqrt(n0_diag(js)*(1.d0-n0_diag(js)))
        if(is.eq.js) then
           slater_out(1,is,js) = n0_slater(1,is,js) - n0_GZproj(1,is,js)
        else
           slater_out(1,is,js) = n0_slater(1,is,js)
        end if
        slater_out(2,is,js) = n0_slater(2,is,js)
        GZproj_out(1,is,js) = n0_GZproj(1,is,js)
        GZproj_out(2,is,js) = n0_GZproj(2,is,js)
     end do
  end do
  !<DEBUG
  ! write(*,*) "DEBUG_N0"
  ! do is=1,Ns
  !    write(*,'(10F8.4)') GZproj_out(1,is,:)
  ! end do
  ! write(*,*)
  ! do is=1,Ns
  !    write(*,'(10F8.4)') GZproj_out(2,is,:)
  ! end do
  !DEBUG>
  call dump2vec_superc(tmp_x,Rhop_out,Qhop_out,slater_out,GZproj_out)
  Fout=0.d0
  
  write(*,*) 'OUT_RQN_opt'
  do is=1,Nopt
     Fout=Fout+tmp_x(is)**2.d0
     !<DEBUG     
     write(*,*) tmp_x(is)
     !DEBUG>
  end do
  

  ! !<DEBUG
  ! do is=1,Nopt
  !    write(*,*) Fout(is)
  ! end do
  ! write(*,*) '!+------------------+!'
  !DEBUG>
  !
  if(optimization_flag) then
     !+- store final informations to global variables -+!                   
     GZ_vector      = GZvect
     GZ_opt_energy  = E_Hstar+E_Hloc
     GZ_opt_kinetic = E_Hstar
     GZ_opt_Eloc    = E_Hloc
     !+- here I should modify this
     GZ_opt_slater_lgr_superc = slater_lgr_multip
     !<TMP DEBUG
     ! if(allocated(GZ_opt_slater_lgr_superc)) write(*,*) 'allocated prima'
     ! write(*,*) '?????',size(GZ_opt_slater_lgr_superc,1),size(GZ_opt_slater_lgr_superc,2),size(GZ_opt_slater_lgr_superc,3)
     ! write(*,*) GZ_opt_slater_lgr_superc(1,1,1)
     ! if(allocated(GZ_opt_slater_lgr_superc)) write(*,*) 'allocated dopo'
     !TMP DEBUG>
  end if


end function R_Q_VDM_free_opt_superc




subroutine dump2mats_superc(x,Rhop,Qhop,slater_lgr_multip,proj_lgr_multip)
  real(8),dimension(:)                :: x
  complex(8),dimension(Ns,Ns)         :: Rhop
  complex(8),dimension(Ns,Ns)         :: Qhop
  complex(8),dimension(2,Ns,Ns)       :: slater_lgr_multip
  complex(8),dimension(2,Ns,Ns)       :: proj_lgr_multip
  integer                             :: NRhop,NQhop,Nslater_lgr,Nproj_lgr,Nopt  
  integer                             :: i,ix,i0,imap,is,js
  complex(8),dimension(:),allocatable :: dump_input
  !
  ! NRhop = 2*Nopt_normal
  ! NQhop = 2*Nopt_anomalous
  ! !
  ! Nslater_lgr = 2*Nopt_lgr  
  ! Nproj_lgr = 2*(Nopt_lgr - Nopt_diag)
  ! !
  ! Nopt = NRhop + NQhop + Nslater_lgr + Nproj_lgr
  ! write(*,*) 'NOPT',Nopt
  ! write(*,*) 'NRhop',NRhop
  ! write(*,*) 'NQhop',NQhop
  ! write(*,*) 'Nslater_lgr',Nslater_lgr
  ! write(*,*) 'Nproj_lgr',Nproj_lgr

  Nopt = 2*NRhop_opt + 2*NQhop_opt + 2*Nvdm_NC_opt + 2*Nvdm_AC_opt + 2*Nvdm_NCoff_opt + 2*Nvdm_AC_opt
  if(size(x).ne.Nopt) stop "wrong dimensions @ dump2mats_superc"
  !
  Rhop = zero
  Qhop = zero
  slater_lgr_multip = zero
  proj_lgr_multip = zero
  !
  i=0
  i0=i
  if(NRhop_opt.gt.0) then
     allocate(dump_input(NRhop_opt))
     do ix=1,NRhop_opt
        dump_input(ix) = x(i0+ix)+xi*x(i0+ix+NRhop_opt); i=i+2
     end do
     !+- HERE WE NEED A FUNCTION (user defined in the driver ->so probably a pointer to a function)
     !   so that ===> the number of Nopt_normal is automatically spread into the (Ns,Ns) matrix (or automatically defined by symmetry...)
     !   Rhop = symmetry_stride_vec2mat_Rhop(dump_input)
     call Rhop_stride_v2m(dump_input,Rhop)
     ! do is=1,Ns
     !    do js=1,Ns
     !       imap=opt_map(is,js)
     !       if(imap.gt.0) Rhop(is,js) = dump_input(imap)
     !    end do
     ! end do
     deallocate(dump_input)
  end if
  !
  i0=i
  if(NQhop_opt.gt.0) then
     allocate(dump_input(NQhop_opt))
     do ix=1,NQhop_opt
        dump_input(ix) = x(i0+ix)+xi*x(i0+ix+NQhop_opt); i=i+2
     end do
     call Qhop_stride_v2m(dump_input,Qhop)
     ! do is=1,Ns
     !    do js=1,Ns
     !       imap=opt_map_anomalous(is,js)
     !       if(imap.gt.0) Qhop(is,js) = dump_input(imap)  !AAAAAA here the SU(2) symmetry requires Q_{up,dn} = -Q_{dn,up} 
     !    end do
     ! end do
     deallocate(dump_input)
  end if
  !
  i0=i
  if(Nvdm_NC_opt.gt.0) then
     allocate(dump_input(Nvdm_NC_opt))
     do ix=1,Nvdm_NC_opt
        dump_input(ix) = x(i0+ix)+xi*x(i0+ix+Nvdm_NC_opt); i=i+2
     end do
     call vdm_NC_stride_v2m(dump_input,slater_lgr_multip(1,:,:))
     ! do is=1,Ns
     !    do js=1,Ns
     !       imap=opt_map(is,js)           
     !       if(imap.gt.0) slater_lgr_multip(1,is,js) = dump_input(imap)
     !    end do
     ! end do
     deallocate(dump_input)
  end if
  i0=i
  if(Nvdm_AC_opt.gt.0) then
     allocate(dump_input(Nvdm_AC_opt))
     do ix=1,Nvdm_AC_opt
        dump_input(ix) = x(i0+ix)+xi*x(i0+ix+Nvdm_AC_opt); i=i+2
     end do
     call vdm_AC_stride_v2m(dump_input,slater_lgr_multip(2,:,:))
     ! do is=1,Ns
     !    do js=1,Ns
     !       imap=opt_map_anomalous(is,js)           
     !       if(imap.gt.0) slater_lgr_multip(2,is,js) = dump_input(imap)
     !    end do
     ! end do
     deallocate(dump_input)
  end if
  i0=i
  if(Nvdm_NCoff_opt.gt.0) then
     allocate(dump_input(Nvdm_NCoff_opt))
     do ix=1,Nvdm_NCoff_opt
        dump_input(ix) = x(i0+ix)+xi*x(i0+ix+Nvdm_NCoff_opt); i=i+2
     end do
     call vdm_NCoff_stride_v2m(dump_input,proj_lgr_multip(1,:,:))
     ! do is=1,Ns
     !    do js=1,Ns
     !       imap=opt_map(is,js)           
     !       if(imap.gt.0.and.is.ne.js) proj_lgr_multip(1,is,js) = dump_input(imap)
     !    end do
     ! end do
     deallocate(dump_input)
  end if
  i0=i
  if(Nvdm_AC_opt.gt.0) then
     allocate(dump_input(Nvdm_AC_opt))
     do ix=1,Nvdm_AC_opt
        dump_input(ix) = x(i0+ix)+xi*x(i0+ix+Nvdm_AC_opt); i=i+2
     end do
     call vdm_AC_stride_v2m(dump_input,proj_lgr_multip(2,:,:))
     ! do is=1,Ns
     !    do js=1,Ns
     !       imap=opt_map_anomalous(is,js)           
     !       if(imap.gt.0) proj_lgr_multip(2,is,js) = dump_input(imap)
     !    end do
     ! end do
     deallocate(dump_input)
  end if
  !
end subroutine dump2mats_superc





subroutine dump2vec_superc(x,Rhop,Qhop,slater_lgr,GZproj_lgr)
  real(8),dimension(:)                :: x
  complex(8),dimension(Ns,Ns)         :: Rhop
  complex(8),dimension(Ns,Ns)         :: Qhop
  complex(8),dimension(2,Ns,Ns)       :: slater_lgr
  complex(8),dimension(2,Ns,Ns)       :: GZproj_lgr
  integer                             :: NRhop,NQhop,Nslater_lgr,Nproj_lgr,Nopt  
  integer                             :: i,ix,i0,imap,is,js
  complex(8),dimension(:),allocatable :: dump_output
  !
  Nopt = 2*NRhop_opt + 2*NQhop_opt + 2*Nvdm_NC_opt + 2*Nvdm_AC_opt + 2*Nvdm_NCoff_opt + 2*Nvdm_AC_opt
  if(size(x).ne.Nopt) stop "wrong dimensions @ dump2vec_superc"
  !
  ! write(*,*) 'NOPT',Nopt
  ! write(*,*) 'NRhop',NRhop
  ! write(*,*) 'NQhop',NQhop
  ! write(*,*) 'Nslater_lgr',Nslater_lgr
  ! write(*,*) 'Nproj_lgr',Nproj_lgr
  !
  i=0
  i0=i
  if(Nrhop_opt.gt.0) then
     allocate(dump_output(NRhop_opt))
     call Rhop_stride_m2v(Rhop,dump_output)
     ! do is=1,Ns
     !    do js=1,Ns
     !       imap=opt_map(is,js)
     !       if(imap.gt.0) dump_output(imap) = Rhop(is,js)
     !    end do
     ! end do
     do ix=1,NRhop_opt
        x(i0+ix) = dreal(dump_output(ix))
        x(i0+ix+NRhop_opt) = dimag(dump_output(ix))
        i=i+2
     end do
     deallocate(dump_output)
  end if
  !
  i0=i
  if(NQhop_opt.gt.0) then
     allocate(dump_output(NQhop_opt))
     call Qhop_stride_m2v(Qhop,dump_output)
     ! do is=1,Ns
     !    do js=1,Ns
     !       imap=opt_map_anomalous(is,js)
     !       if(imap.gt.0) dump_output(imap) = Qhop(is,js)
     !    end do
     ! end do
     do ix=1,NQhop_opt
        x(i0+ix) = dreal(dump_output(ix))
        x(i0+ix+NQhop_opt) = dimag(dump_output(ix))
        i=i+2
     end do
     deallocate(dump_output)     
  end if
  !
  i0=i
  if(Nvdm_NC_opt.gt.0) then
     allocate(dump_output(Nvdm_NC_opt))
     call vdm_NC_stride_m2v(slater_lgr(1,:,:),dump_output)
     ! do is=1,Ns
     !    do js=1,Ns
     !       imap=opt_map(is,js)
     !       if(imap.gt.0) dump_output(imap) = slater_lgr(1,is,js)
     !    end do
     ! end do
     do ix=1,Nvdm_NC_opt
        x(i0+ix) = dreal(dump_output(ix))
        x(i0+ix+Nvdm_NC_opt) = dimag(dump_output(ix))
        i=i+2
     end do
     deallocate(dump_output)
  end if
  i0=i
  if(Nvdm_AC_opt.gt.0) then
     allocate(dump_output(Nvdm_AC_opt))
     call vdm_AC_stride_m2v(slater_lgr(2,:,:),dump_output)
     ! do is=1,Ns
     !    do js=1,Ns
     !       imap=opt_map_anomalous(is,js)
     !       if(imap.gt.0) dump_output(imap) = slater_lgr(2,is,js)
     !    end do
     ! end do
     do ix=1,Nvdm_AC_opt
        x(i0+ix) = dreal(dump_output(ix))
        x(i0+ix+Nvdm_AC_opt) = dimag(dump_output(ix))
        i=i+2
     end do
     deallocate(dump_output)     
  end if
  i0=i
  if(Nvdm_NCoff_opt.gt.0) then
     allocate(dump_output(Nvdm_NCoff_opt))
     call vdm_NCoff_stride_m2v(GZproj_lgr(1,:,:),dump_output)
     ! do is=1,Ns
     !    do js=1,Ns
     !       imap=opt_map(is,js)
     !       if(imap.gt.0.and.is.ne.js) dump_output(imap) = GZproj_lgr(1,is,js)
     !    end do
     ! end do
     do ix=1,Nvdm_NCoff_opt
        x(i0+ix) = dreal(dump_output(ix))
        x(i0+ix+Nvdm_NCoff_opt) = dimag(dump_output(ix))
        i=i+2
     end do
     deallocate(dump_output)
  end if
  i0=i
  if(Nvdm_AC_opt.gt.0) then
     allocate(dump_output(Nvdm_AC_opt))
     call vdm_AC_stride_m2v(GZproj_lgr(2,:,:),dump_output)
     ! do is=1,Ns
     !    do js=1,Ns
     !       imap=opt_map_anomalous(is,js)
     !       if(imap.gt.0) dump_output(imap) = GZproj_lgr(2,is,js)
     !    end do
     ! end do
     do ix=1,Nvdm_AC_opt
        x(i0+ix) = dreal(dump_output(ix))
        x(i0+ix+Nvdm_AC_opt) = dimag(dump_output(ix))
        i=i+2
     end do
     deallocate(dump_output)     
  end if
  !
end subroutine dump2vec_superc

















function R_VDM_free_root(x) result(delta_opt)

  real(8),dimension(:),intent(in) :: x      
  real(8)             :: delta_opt
  !  integer,optional,intent(in)     :: i

  real(8)                         :: E_Hstar,E_Hloc
  complex(8),dimension(Nphi)      :: GZvect

  complex(8),dimension(Ns,Ns) :: Rhop,Rhop_GZproj
  complex(8),dimension(Ns,Ns) :: slater_derivatives
  real(8),dimension(Ns,Ns)    :: slater_lgr_multip
  real(8),dimension(Ns,Ns)    :: proj_lgr_multip
  complex(8),dimension(Ns,Ns) :: Rhop_lgr_multip
  real(8),dimension(Ns,Ns)    :: n0_slater,n0_GZproj
  real(8),dimension(Ns)       :: n0_diag
  real(8)                     :: tmp
  integer                     :: imap,is,js
  integer :: Nslater_lgr,Nrhop,Nproj_lgr
  logical :: bound_flag


  Nslater_lgr = Nopt_diag + Nopt_odiag
  NRhop = Nopt_diag + Nopt_odiag
  Nproj_lgr = Nopt_odiag

  bound_flag=.false.
  slater_lgr_multip=0.d0
  proj_lgr_multip=0.d0
  Rhop = zero
  do is=1,Ns
     do js=1,Ns
        imap = opt_map(is,js)
        if(imap.gt.0) then
           slater_lgr_multip(is,js)=x(imap)
           tmp=x(imap)
           if(tmp.gt.10.d0.or.tmp.lt.-10.d0) bound_flag=.true.
           write(*,*) 'density_con',bound_flag
           Rhop(is,js) = x(imap + Nslater_lgr) 
           tmp = x(imap + Nslater_lgr) 
           if(tmp.gt.1.d0.or.tmp.lt.0.d0) bound_flag=.true.
           write(*,*) 'Rhop',bound_flag
           Rhop(is,js) = Rhop(is,js) + xi*x(imap + Nslater_lgr + NRhop)
           tmp = x(imap + Nslater_lgr + NRhop) 
           if(tmp.gt.1.d0.or.tmp.lt.0.d0) bound_flag=.true.
           write(*,*) 'Rhop_imag',bound_flag
           if(is.ne.js) then
              proj_lgr_multip(is,js) = x(imap + Nslater_lgr + 2*NRhop)
              tmp = x(imap + Nslater_lgr + 2*NRhop) 
              if(tmp.gt.1.d0.or.tmp.lt.0.d0) bound_flag=.true.
              write(*,*) 'gz_lgr',bound_flag
           end if
        end if
     end do
  end do

  if(.not.bound_flag) then

     !
     call slater_minimization_fixed_lgr(Rhop,slater_lgr_multip,E_Hstar,n0_slater,slater_derivatives)
     !
     !
     !+-
     write(*,*) n0_slater
     !-+
     !
     !
     do is=1,Ns
        do js=1,Ns
           if(abs(n0_slater(js,js)).gt.1.d-10.and.abs(n0_slater(js,js)).lt.1.d0-1.d-10) then
              Rhop_lgr_multip(is,js) = slater_derivatives(is,js)/sqrt(n0_slater(js,js)*(1.d0-n0_slater(js,js)))
           else
              Rhop_lgr_multip(is,js) = 0.d0
           end if
        end do
        !
     end do
     !
     do is=1,Ns
        tmp = 0.d0
        do js=1,Ns
           tmp = tmp + dreal(Rhop_lgr_multip(js,is)*Rhop(js,is))
        end do
        if(abs(n0_slater(is,is)).gt.1.d-10.and.abs(n0_slater(is,is)).lt.1.d0-1.d-10) then
           proj_lgr_multip(is,is) = -1.d0*slater_lgr_multip(is,is) + &
                0.5d0*(1.d0-2.d0*n0_slater(is,is))/sqrt(n0_slater(is,is)*(1.d0-n0_slater(is,is)))*tmp 
        else
           proj_lgr_multip(is,is) = -1.d0*slater_lgr_multip(is,is)
        end if
     end do
     !
     do is=1,Ns
        n0_diag(is) = n0_slater(is,is)
     end do
     !

      
     call gz_proj_minimization_fixed_lgr(n0_diag,proj_lgr_multip,Rhop_lgr_multip,E_Hloc,GZvect)
     !
     n0_GZproj = 0.d0
     do is=1,Ns
        do js=1,Ns
           n0_GZproj(is,js) = trace_phi_basis(GZvect,phi_traces_basis_dens(is,js,:,:))
           Rhop_GZProj(is,js) = trace_phi_basis(GZvect,phi_traces_basis_Rhop(is,js,:,:))
           !Rhop_GZProj(is,js) = Rhop_GZProj(is,js)
           !
           write(*,*) is,js,n0_GZproj(is,js),n0_slater(is,js)
           !
        end do
     end do
     !
     delta_opt = 0.d0
     do is=1,Ns
        do js=1,Ns
           if(is.ne.js) then
              delta_opt = delta_opt +  abs(n0_slater(is,js))**2.d0
              delta_opt = delta_opt +  abs(n0_GZproj(is,js))**2.d0
           else
              delta_opt = delta_opt + abs(n0_slater(is,js)-n0_GZproj(is,js))**2.d0               
           end if
           delta_opt = delta_opt + abs(Rhop_GZproj(is,js)-Rhop(is,js)*sqrt(n0_diag(js)*(1.d0-n0_diag(js))))**2.d0
        end do
     end do
  else
     delta_opt = 30.d0
  end if

  if(optimization_flag) then
     !+- store final informations to global variables -+!              
     GZ_vector      = GZvect
     GZ_opt_energy  = E_Hstar+E_Hloc
     GZ_opt_kinetic = E_Hstar
     GZ_opt_Eloc    = E_Hloc
  end if

  write(123,*) delta_opt,x
  write(124,*) n0_slater
  write(125,*) n0_GZproj

  
end function R_VDM_free_root
