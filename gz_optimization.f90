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


  Nslater_lgr = Nopt_diag + Nopt_odiag
  NRhop = Nopt_diag + Nopt_odiag
  Nproj_lgr = Nopt_odiag

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
        if(abs(n0_slater(js,js)).gt.1.d-10) then
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
     if(abs(n0_slater(is,is)).gt.1.d-10) then
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
        Rhop_GZProj(is,js) = Rhop_GZProj(is,js)
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
  if(optimization_flag) then
     !+- store final informations to global variables -+!              
     GZ_vector      = GZvect
     GZ_opt_energy  = E_Hstar+E_Hloc
     GZ_opt_kinetic = E_Hstar
     GZ_opt_Eloc    = E_Hloc
  end if
end subroutine R_VDM_free_opt_function
