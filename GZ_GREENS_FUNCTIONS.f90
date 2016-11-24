MODULE GZ_neqGREENS_FUNCTIONS
  USE GZ_VARS_GLOBAL
  USE SF_LINALG
  implicit none  
  private


  
  !   parameters given by the driver:
  !   Hk_bare (use the same pointer_functions of the usual drivers)
  !   Rhop,Qhop
  !   initial_slater
  
  !   the main idea:
  !   computing the G_local
  ! 1 ) build up the time evolution operators Vt(Ntgf,Ns,Ns)
  !   !+----
  ! 2 ) loop on the momenta (actual calculation of Gloc)
  !   init gloc
  ! 2 a)  loop on the two times [for the Gloc is ok to use the parallelogramma config for the t,t plane (t+s,t) t,s=1,..,Ntgf] -> NO CHANGE, beacause trivially the Gret is known only the lower triangular
  ! 2 b)  create Gk_loc
  ! 2 c)  sum_up to Gloc
  !   !+----

  public :: get_relative_time_FT
  public :: gz_get_Gloc_ret
  public :: gz_get_Gloc_ret_superc
  public :: gz_get_Gloc_ret_superc_diag_hk
  

contains

  
  subroutine gz_get_Gloc_ret(Rhop,Gloc_ret)
    complex(8),dimension(2*Ntgf,Ns,Ns) :: Rhop
    complex(8),dimension(Ntgf,Ntgf,Ns,Ns) :: Gloc_ret
    !
    complex(8),dimension(Ns,Ns) :: Gk_ret,Gk_ret_00
    !
    real(8) :: t,tt    
    complex(8),dimension(2*Ntgf-1,Ns,Ns) :: Vt,Vt_dag
    complex(8),dimension(Ns,Ns) :: Vtk,Vttk,Rt,Rtt
    !
    integer :: ik,is,js,it,jt,iti,jtj
    !
    write(*,*) "building Gloc_ret "
    Gloc_ret=zero
    do ik=1,Lk
       !
       Gk_ret_00=zero
       do is=1,Ns
          Gk_ret_00(is,is) = -xi
       end do       
       !
       call get_time_evolution_operators(ik,Vt,Vt_dag,Rhop)
       !
       do it=1,Ntgf
          do jt=1,Ntgf
             !
             iti = it + Ntgf - 1
             jtj = iti + 1 - jt
             !
             Rt=Rhop(iti,:,:)
!             if(ik.eq.200) write(*,'(10F7.3)') Rt
             do is=1,Ns
                do js=1,Ns
                   Rtt(is,js)=conjg(Rhop(jtj,js,is))
                end do
             end do
             !             
             Vtk  = Vt(iti,:,:)
             Vttk = Vt_dag(jtj,:,:)
             !
             Gk_ret= matmul(Vtk,Gk_ret_00)
             Gk_ret=matmul(Gk_ret,Vttk)
             !
             !slave bosons like post production
             Gk_ret=matmul(Rt,Gk_ret)
             Gk_ret=matmul(Gk_ret,Rtt)
             !
             Gloc_ret(it,jt,:,:) = Gloc_ret(it,jt,:,:) + Gk_ret*wtk(ik)
             !             
          end do
       end do
       if(mod(ik,Lk/10).eq.0) write(*,'(F5.1)') dble(ik)/dble(Lk)*100       
    end do
  end subroutine gz_get_Gloc_ret
  



  subroutine gz_get_Gloc_ret_superc(Rhop,Qhop,Gloc_ret,sl_lgr_)
    implicit none
    complex(8),dimension(Nttgf,Ns,Ns) :: Rhop,Qhop
    complex(8),dimension(Nttgf,Ns,Ns),optional :: sl_lgr_
    complex(8),dimension(Nttgf,Ns,Ns) :: sl_lgr
    complex(8),dimension(Ntgf,Ntgf,2*Ns,2*Ns) :: Gloc_ret
    !
    complex(8),dimension(2*Ns,2*Ns) :: Gk_ret_00,tmpGk
    complex(8),dimension(2*Ns,2*Ns) :: Gk_ret,tmpGk_
    !
    real(8) :: t,tt    
    complex(8),dimension(Nttgf,2*Ns,2*Ns) :: Vt,Vt_dag,intHt
    complex(8),dimension(2*Ns,2*Ns) :: Vtk,Vttk,tmpVt,tmpVtt
    real(8),dimension(2*Ns) :: tmp_iHt,tmp_iHtt
    complex(8),dimension(Ns,Ns) :: Rt,Rtt,Qt,Qtt
    complex(8),dimension(2*Ns,2*Ns) :: sqZt,sqZtt
    !
    integer :: ik,is,js,it,jt,iti,jtj,iis
    !
    Gloc_ret=zero
    sqZt=zero
    sqZtt=zero
    Rt=zero;Rtt=zero;Qt=zero;Qtt=zero
    
    sl_lgr=zero
    if(present(sl_lgr_)) sl_lgr = sl_lgr_

    do ik=1,Lk
       !
       Gk_ret_00=zero
       do is=1,2*Ns
          Gk_ret_00(is,is) = -xi
       end do
       !
       call get_hamiltonian_time_int_superc(ik,intHt,Rhop,Qhop,sl_lgr)
       !
       do it=1,Ntgf
          do jt=1,Ntgf
             !
             iti = it + Nt0 - 1
             jtj = iti + 1 - jt
             !

             Rt=Rhop(iti,:,:)
             do is=1,Ns
                do js=1,Ns
                   Rtt(is,js)=conjg(Rhop(jtj,js,is))
                end do
             end do
             Qt=Qhop(iti,:,:)
             do is=1,Ns
                do js=1,Ns
                   Qtt(is,js)=conjg(Qhop(jtj,js,is))
                end do
             end do
             !
             do is=1,Ns
                do js=1,Ns
                   sqZt(is,js) = Rt(is,js)
                   sqZt(is,js+Ns) = Qt(is,js)
                   sqZt(is+Ns,js) = conjg(Qt(js,is))
                   sqZt(is+Ns,js+Ns) = conjg(Rt(js,is))
                   !
                   sqZtt(is,js) = Rtt(is,js)
                   sqZtt(is,js+Ns) = conjg(Qtt(js,is))
                   sqZtt(is+Ns,js) = Qtt(is,js)
                   sqZtt(is+Ns,js+Ns) = conjg(Rtt(js,is))
                end do
             end do
             !
             tmpVt = intHt(iti,:,:) - intHt(jtj,:,:)
             call matrix_diagonalize(tmpVt,tmp_iHt)
             !
             do is=1,2*Ns
                do js=1,2*Ns
                   tmpVtt(is,js)=-xi*exp(xi*tmp_iHt(is))*conjg(tmpVt(js,is))
                end do
             end do
             !
             !
             tmpGk = matmul(tmpVt,tmpVtt)
             !
             !
             ! do is=1,2*Ns
             !    do js=1,2*Ns
             !       Vttk(is,js) = zero
             !       do iis=1,2*Ns
             !          Vttk(is,js) = Vttk(is,js) + tmpVt(is,iis)*conjg(tmpVt(js,iis))*exp(xi*tmp_iHt(iis))
             !       end do
             !    end do
             ! end do
             ! !
             ! tmpGk = matmul(Vttk,Gk_ret_00)
             ! !
             Gk_ret = matmul(tmpGk,sqZtt)
             Gk_ret = matmul(sqZt,Gk_ret) 
             !
             Gloc_ret(it,jt,:,:) = Gloc_ret(it,jt,:,:) + Gk_ret*wtk(ik)
             !
          end do          
       end do
       write(*,*) ik,Lk
    end do
  end subroutine gz_get_Gloc_ret_superc


  subroutine gz_get_Gloc_ret_superc_diag_hk(Rhop,Qhop,Gloc_ret,sl_lgr_)
    implicit none
    complex(8),dimension(Nttgf,Ns,Ns) :: Rhop,Qhop
    complex(8),dimension(Nttgf,Ns,Ns),optional :: sl_lgr_
    complex(8),dimension(Nttgf,Ns,Ns) :: sl_lgr
    complex(8),dimension(Ntgf,Ntgf,2*Ns,2*Ns) :: Gloc_ret
    !
    complex(8),dimension(2*Ns,2*Ns) :: Gk_ret_00,tmpGk
    complex(8),dimension(2*Ns,2*Ns) :: Gk_ret,tmpGk_
    !
    real(8) :: t,tt,ek
    complex(8),dimension(Nttgf,2*Ns,2*Ns) :: Vt,Vt_dag,intHt
    complex(8),dimension(2*Ns,2*Ns) :: Vtk,Vttk,tmpVt,tmpVt_,tmpVtt
    real(8),dimension(2*Ns) :: tmp_iHt,tmp_iHtt
    complex(8),dimension(Ns,Ns) :: Rt,Rtt,Qt,Qtt,Hk
    complex(8),dimension(2*Ns,2*Ns) :: sqZt,sqZtt
    !
    complex(8),dimension(Ntgf,Ntgf,2*Ns,2*Ns) :: Vt_exp
    real(8),dimension(Ntgf,Ntgf,2*Ns) :: Ht_exp
    !
    integer :: ik,is,js,it,jt,iti,jtj,iis
    !
    Gloc_ret=zero
    sqZt=zero
    sqZtt=zero
    Rt=zero;Rtt=zero;Qt=zero;Qtt=zero
    sl_lgr = zero
    if(present(sl_lgr_)) sl_lgr = sl_lgr_
    !
    call get_hamiltonian_time_int_superc_(intHt,Rhop,Qhop,sl_lgr)
    !
    write(*,*) "diagonalizing t,t' evolution operators"
    !
    do it=1,Ntgf
       do jt=1,Ntgf
          !
          iti = it + Nt0 - 1
          jtj = iti + 1 - jt
          !
          tmpVt = intHt(iti,:,:) - intHt(jtj,:,:)
          call matrix_diagonalize(tmpVt,tmp_iHt)
          !
          Vt_exp(it,jt,:,:) = tmpVt
          Ht_exp(it,jt,:) = tmp_iHt
          !
       end do
       if(mod(it,Ntgf/10).eq.0) write(*,*) it,'/',Ntgf,'accomplished'
    end do
    !
    do ik=1,Lk
       !
       Gk_ret_00=zero
       do is=1,2*Ns
          Gk_ret_00(is,is) = -xi
       end do
       !
       do it=1,Ntgf
          do jt=1,Ntgf
             !
             iti = it + Nt0 - 1
             jtj = iti + 1 - jt
             !
             Rt=Rhop(iti,:,:)
             do is=1,Ns
                do js=1,Ns
                   Rtt(is,js)=conjg(Rhop(jtj,js,is))
                end do
             end do
             Qt=Qhop(iti,:,:)
             do is=1,Ns
                do js=1,Ns
                   Qtt(is,js)=conjg(Qhop(jtj,js,is))
                end do
             end do
             !
             do is=1,Ns
                do js=1,Ns
                   sqZt(is,js) = Rt(is,js)
                   sqZt(is,js+Ns) = Qt(is,js)
                   sqZt(is+Ns,js) = conjg(Qt(js,is))
                   sqZt(is+Ns,js+Ns) = conjg(Rt(js,is))
                   !
                   sqZtt(is,js) = Rtt(is,js)
                   sqZtt(is,js+Ns) = conjg(Qtt(js,is))
                   sqZtt(is+Ns,js) = Qtt(is,js)
                   sqZtt(is+Ns,js+Ns) = conjg(Rtt(js,is))
                end do
             end do
             !
             !
             call get_Hk_t(Hk,ik,0.d0);ek=dreal(Hk(1,1))
             !
             tmpVt=Vt_exp(it,jt,:,:)
             tmp_iHt = Ht_exp(it,jt,:)
             !
             do is=1,2*Ns
                do js=1,2*Ns
                   tmpVt_(is,js) = -xi*exp(xi*tmp_iHt(is)*ek)*conjg(Vt_exp(it,jt,js,is))
                end do
             end do
             !
             tmpGk = matmul(tmpVt,tmpVt_)
             !
             Gk_ret = matmul(tmpGk,sqZtt)
             Gk_ret = matmul(sqZt,Gk_ret) 
             !
             Gloc_ret(it,jt,:,:) = Gloc_ret(it,jt,:,:) + Gk_ret*wtk(ik)
             !
          end do
       end do
       write(*,*) ik,Lk
       !if(mod(ik,Lk/10).eq.0) write(*,'(F5.1)') dble(ik)/dble(Lk)*100       
    end do
  end subroutine gz_get_Gloc_ret_superc_diag_hk


  subroutine get_time_evolution_operators(ik,Vt,Vt_dag,Rhop)
    complex(8),dimension(2*Ntgf-1,Ns,Ns) :: Vt,Vt_dag,Rhop,Rhop_dag
    integer :: ik,is,js,iis,it
    complex(8),dimension(Ns,Ns) :: expHt,Hk
    complex(8),dimension(Ns,Ns) :: expHt_v
    real(8),dimension(Ns) :: expHt_e
    real(8) :: time
    !
    expHt=zero
    do it=1,2*Ntgf-1
       !
       do is=1,Ns
          do js=1,Ns
             Rhop_dag(it,is,js) = conjg(Rhop(it,js,is))
          end do
       end do
       !
       if(it.gt.1) then
          time = t_grid(it)          
          call get_Hk_t(Hk,ik,time)
          Hk = matmul(Rhop_dag(it,:,:),Hk)
          Hk = matmul(Hk,Rhop(it,:,:))
          expHt = expHt + Hk*tstep*0.5d0
          !
          time = t_grid(it-1)
          call get_Hk_t(Hk,ik,time)
          Hk = matmul(Rhop_dag(it-1,:,:),Hk)
          Hk = matmul(Hk,Rhop(it-1,:,:))
          expHt = expHt + Hk*tstep*0.5d0
       end if
       !
       expHt_v = expHt
       call matrix_diagonalize(expHt_v,expHt_e)
       !
       do is=1,2
          do js=1,2
             Vt(it,is,js)=zero
             do iis=1,Ns
                Vt(it,is,js)=Vt(it,is,js)+expHt_v(is,iis)*conjg(expHt_v(js,iis))*exp(xi*expHt_e(iis))
              end do
           end do
        end do
        do is=1,Ns
           do js=1,Ns
              Vt_dag(it,is,js) = conjg(Vt(it,js,is))
           end do
        end do
     end do
   end subroutine get_time_evolution_operators
   !
   !
   !
   subroutine get_hamiltonian_time_int_superc(ik,intHt,Rhop,Qhop,sl_lgr)
     implicit none
     complex(8),dimension(Nttgf,2*Ns,2*Ns) :: intHt
     complex(8),dimension(Nttgf,Ns,Ns) :: Rhop,Qhop,sl_lgr
     complex(8),dimension(Nttgf,Ns,Ns) :: Rhop_dag,Qhop_dag     
     integer :: ik,is,js,iis,it,iit
     complex(8),dimension(Ns,Ns) :: Hk,expHt_tmp
     complex(8),dimension(2*Ns,2*Ns) :: expHt,expHt_v,Ht,Vtmp
     complex(8),dimension(2*Ns,2*Ns) :: expHt_,expHt_v_
     real(8),dimension(2*Ns) :: expHt_e,expHt_e_
     real(8) :: time
     !
     !
     intHt=zero     
     Ht=zero
     expHt=zero
     !
     !
     do it=1,Nttgf
        !
        !
        do is=1,Ns
           do js=1,Ns
              Rhop_dag(it,is,js) = conjg(Rhop(it,js,is))
              Qhop_dag(it,is,js) = conjg(Qhop(it,js,is))             
           end do
        end do
        !
        intHt(it,:,:) = expHt
        
        if(it.gt.1) then
           time = t_grid(it)          
           call get_Hk_t(Hk,ik,time)
           !
           expHt_tmp = matmul(Hk,Rhop(it,:,:))
           expHt_tmp = matmul(Rhop_dag(it,:,:),expHt_tmp)
           Ht(1:Ns,1:Ns) = expHt_tmp 
           !
           expHt_tmp = matmul(Hk,Qhop(it,:,:))
           expHt_tmp = matmul(Rhop_dag(it,:,:),expHt_tmp)
           Ht(1:Ns,Ns+1:2*Ns) = expHt_tmp + sl_lgr(it,:,:)
           !
           expHt_tmp = matmul(Hk,Rhop(it,:,:))
           expHt_tmp = matmul(Qhop_dag(it,:,:),expHt_tmp)  
           Ht(Ns+1:2*Ns,1:Ns) = expHt_tmp + conjg(transpose(sl_lgr(it,:,:)))
           !
           expHt_tmp = matmul(Hk,Qhop(it,:,:))
           expHt_tmp = matmul(Qhop_dag(it,:,:),expHt_tmp)
           Ht(Ns+1:2*Ns,Ns+1:2*Ns) = expHt_tmp
           !
           expHt = expHt + Ht*tstep*0.5d0
           !
           time = t_grid(it-1)
           call get_Hk_t(Hk,ik,time)
           !
           expHt_tmp = matmul(Hk,Rhop(it-1,:,:))
           expHt_tmp = matmul(Rhop_dag(it-1,:,:),expHt_tmp)
           Ht(1:Ns,1:Ns) = expHt_tmp
           !
           expHt_tmp = matmul(Hk,Qhop(it-1,:,:))
           expHt_tmp = matmul(Rhop_dag(it-1,:,:),expHt_tmp)
           Ht(1:Ns,Ns+1:2*Ns) = expHt_tmp + sl_lgr(it-1,:,:)
          !
           expHt_tmp = matmul(Hk,Rhop(it-1,:,:))
           expHt_tmp = matmul(Qhop_dag(it-1,:,:),expHt_tmp)
           Ht(Ns+1:2*Ns,1:Ns) = expHt_tmp + conjg(transpose(sl_lgr(it-1,:,:)))
           !
           expHt_tmp = matmul(Hk,Qhop(it-1,:,:))
           expHt_tmp = matmul(Qhop_dag(it-1,:,:),expHt_tmp)
           Ht(Ns+1:2*Ns,Ns+1:2*Ns) = expHt_tmp
           !
           expHt = expHt + Ht*tstep*0.5d0
           !
        end if
     end do
     !
   end subroutine get_hamiltonian_time_int_superc
   !
   subroutine get_hamiltonian_time_int_superc_(intHt,Rhop,Qhop,sl_lgr)
     implicit none
     complex(8),dimension(Nttgf,2*Ns,2*Ns) :: intHt
     complex(8),dimension(Nttgf,Ns,Ns) :: Rhop,Qhop,sl_lgr
     complex(8),dimension(Nttgf,Ns,Ns) :: Rhop_dag,Qhop_dag
     integer :: ik,is,js,iis,it,iit
     complex(8),dimension(Ns,Ns) :: Hk,expHt_tmp
     complex(8),dimension(2*Ns,2*Ns) :: expHt,expHt_v,Ht,Vtmp
     complex(8),dimension(2*Ns,2*Ns) :: expHt_,expHt_v_
     real(8),dimension(2*Ns) :: expHt_e,expHt_e_
     real(8) :: time
     !
     !
     intHt=zero     
     Ht=zero
     expHt=zero
     !
     write(*,*) "computing time integrals"     
     !
     do it=1,Nttgf
        !
        !
        do is=1,Ns
           do js=1,Ns
              Rhop_dag(it,is,js) = conjg(Rhop(it,js,is))
              Qhop_dag(it,is,js) = conjg(Qhop(it,js,is))             
           end do
        end do
        !
        !
        intHt(it,:,:) = expHt
        !
        !
        if(it.gt.1) then
           !
           expHt_tmp = matmul(Rhop_dag(it,:,:),Rhop(it,:,:))
           Ht(1:Ns,1:Ns) = expHt_tmp
           !
           expHt_tmp = matmul(Rhop_dag(it,:,:),Qhop(it,:,:))
           Ht(1:Ns,Ns+1:2*Ns) = expHt_tmp + sl_lgr(it,:,:)  !<<<----wrong!!! can not use the k-diag trick!!
           !
           expHt_tmp = matmul(Qhop_dag(it,:,:),Rhop(it,:,:))  
           Ht(Ns+1:2*Ns,1:Ns) = expHt_tmp + conjg(transpose(sl_lgr(it,:,:)))
           !
           expHt_tmp = matmul(Qhop_dag(it,:,:),Qhop(it,:,:))
           Ht(Ns+1:2*Ns,Ns+1:2*Ns) = expHt_tmp
           !
           !
           expHt = expHt + Ht*tstep*0.5d0
           !
           !
           expHt_tmp = matmul(Rhop_dag(it-1,:,:),Rhop(it-1,:,:))
           Ht(1:Ns,1:Ns) = expHt_tmp
           !
           expHt_tmp = matmul(Rhop_dag(it-1,:,:),Qhop(it-1,:,:))
           Ht(1:Ns,Ns+1:2*Ns) = expHt_tmp + sl_lgr(it-1,:,:)
           !
           expHt_tmp = matmul(Qhop_dag(it-1,:,:),Rhop(it-1,:,:))
           Ht(Ns+1:2*Ns,1:Ns) = expHt_tmp + conjg(transpose(sl_lgr(it-1,:,:)))
           !
           expHt_tmp = matmul(Qhop_dag(it-1,:,:),Qhop(it-1,:,:))
           Ht(Ns+1:2*Ns,Ns+1:2*Ns) = expHt_tmp
           !
           expHt = expHt + Ht*tstep*0.5d0
           !
        end if
     end do
     !
   end subroutine get_hamiltonian_time_int_superc_
   !
   subroutine get_relative_time_FT(Ftt,Ftw,wre,deps_)
     real(8),dimension(:) :: wre
     complex(8),dimension(:,:) :: Ftw
     complex(8),dimension(Ntgf,Ntgf) :: Ftt
     complex(8),dimension(Ntgf)    :: Ft
     real(8),optional :: deps_
     integer :: Nw
     integer :: it,iit,iw,iti,jjt
     real(8) :: t,tt,deps
     
     deps=0.001d0
     if(present(deps_)) deps=deps_
     !
     Nw =size(wre)
     if(size(Ftw,1).ne.size(Ftt,1)) stop "wrong dimension 1 Ftw"
     if(size(Ftw,2).ne.Nw) stop "wrong dimension 2  Ftw"
     !
     do it=1,Ntgf
        Ft=Ftt(it,:)
        do iw=1,Nw
           Ftw(it,iw) = zero
           do iit=1,Ntgf-1
              iti = it + Nt0 -1
              jjt = iti + 1 - iit - 1
              t = t_grid(iti)
              tt = t_grid(jjt)
              Ftw(it,iw) = Ftw(it,iw) + exp(xi*(wre(iw)+xi*deps)*(t-tt))*Ft(iit)*tstep*0.5
              tt = t_grid(jjt+1)
              Ftw(it,iw) = Ftw(it,iw) + exp(xi*(wre(iw)+xi*deps)*(t-tt))*Ft(iit)*tstep*0.5
           end do
        end do
     end do
   end subroutine get_relative_time_FT
   !
 END MODULE GZ_neqGREENS_FUNCTIONS






 
   ! subroutine get_time_evolution_operators_superc_tt(ik,iti,jtj,Vt,Rhop,Qhop)
   !   implicit none
   !   complex(8),dimension(2*Ns,2*Ns) :: Vt
   !   complex(8),dimension(Nttgf,Ns,Ns) :: Rhop,Qhop
   !   complex(8),dimension(Ns,Ns) :: Rhop_dag,Qhop_dag
   !   integer :: ik,iti,jtj,is,js,iis,it
   !   complex(8),dimension(Ns,Ns) :: Hk,expHt_tmp
   !   complex(8),dimension(2*Ns,2*Ns) :: expHt,expHt_v,Ht,Vtmp
   !   real(8),dimension(2*Ns) :: expHt_e
   !   real(8) :: time
   !   !
   !   expHt=zero
   !   Ht=zero
   !   Vt=zero
   !   !
   !   do it=jtj+1,iti
   !      time = t_grid(it)          
   !      call get_Hk_t(Hk,ik,time)
   !      !
   !      Rhop_dag=conjg(transpose(Rhop(it,:,:)))
   !      Qhop_dag=conjg(transpose(Qhop(it,:,:)))
   !      !
   !      expHt_tmp = matmul(Hk,Rhop(it,:,:))
   !      expHt_tmp = matmul(Rhop_dag,expHt_tmp)
   !      Ht(1:Ns,1:Ns) = expHt_tmp
   !      expHt_tmp = matmul(Hk,Qhop(it,:,:))
   !      expHt_tmp = matmul(Rhop_dag,expHt_tmp)
   !      Ht(1:Ns,Ns+1:2*Ns) = expHt_tmp
   !      expHt_tmp = matmul(Hk,Rhop(it,:,:))
   !      expHt_tmp = matmul(Qhop_dag,expHt_tmp)  
   !      Ht(Ns+1:2*Ns,1:Ns) = expHt_tmp
   !      expHt_tmp = matmul(Hk,Qhop(it,:,:))
   !      expHt_tmp = matmul(Qhop_dag,expHt_tmp)
   !      Ht(Ns+1:2*Ns,Ns+1:2*Ns) = expHt_tmp
   !      !
   !      !
   !      expHt = expHt + Ht*tstep*0.5d0
   !      !
   !      !
   !      time = t_grid(it-1)          
   !      call get_Hk_t(Hk,ik,time)
   !      !
   !      Rhop_dag=conjg(transpose(Rhop(it-1,:,:)))
   !      Qhop_dag=conjg(transpose(Qhop(it-1,:,:)))
   !      !
   !      expHt_tmp = matmul(Hk,Rhop(it-1,:,:))
   !      expHt_tmp = matmul(Rhop_dag,expHt_tmp)
   !      Ht(1:Ns,1:Ns) = expHt_tmp
   !      expHt_tmp = matmul(Hk,Qhop(it-1,:,:))
   !      expHt_tmp = matmul(Rhop_dag,expHt_tmp)
   !      Ht(1:Ns,Ns+1:2*Ns) = expHt_tmp
   !      expHt_tmp = matmul(Hk,Rhop(it-1,:,:))
   !      expHt_tmp = matmul(Qhop_dag,expHt_tmp)  
   !      Ht(Ns+1:2*Ns,1:Ns) = expHt_tmp
   !      expHt_tmp = matmul(Hk,Qhop(it-1,:,:))
   !      expHt_tmp = matmul(Qhop_dag,expHt_tmp)
   !      Ht(Ns+1:2*Ns,Ns+1:2*Ns) = expHt_tmp
   !      !
   !      expHt = expHt + Ht*tstep*0.5d0
   !      !
   !   end do
   !   !
   !   expHt_v = expHt
   !   call matrix_diagonalize(expHt_v,expHt_e)     
   !   !
   !   do is=1,2*Ns
   !      do js=1,2*Ns
   !         Vt(is,js)=zero
   !         do iis=1,2*Ns
   !            Vt(is,js)=Vt(is,js)+expHt_v(is,iis)*conjg(expHt_v(js,iis))*exp(xi*expHt_e(iis))
   !         end do
   !      end do
   !   end do
   
   ! end subroutine get_time_evolution_operators_superc_tt




    
   ! subroutine get_time_evolution_operators_superc(ik,Vt,Vt_dag,Rhop,Qhop)
   !   implicit none
   !   complex(8),dimension(Nttgf,2*Ns,2*Ns) :: Vt,Vt_dag
   !   complex(8),dimension(Nttgf,Ns,Ns) :: Rhop,Qhop
   !   complex(8),dimension(Nttgf,Ns,Ns) :: Rhop_dag,Qhop_dag
   !   integer :: ik,is,js,iis,it,iit
   !   complex(8),dimension(Ns,Ns) :: Hk,expHt_tmp
   !   complex(8),dimension(2*Ns,2*Ns) :: expHt,expHt_v,Ht,Vtmp
   !   complex(8),dimension(2*Ns,2*Ns) :: expHt_,expHt_v_
   !   real(8),dimension(2*Ns) :: expHt_e,expHt_e_
   !   real(8) :: time
   !   !
   !   expHt=zero
   !   expHt_=zero
     
   !   Ht=zero
   !   do it=1,Nttgf
   !      !
   !      Rhop_dag(it,:,:) = conjg(transpose(Rhop(it,:,:)))
   !      Qhop_dag(it,:,:) = conjg(transpose(Qhop(it,:,:)))        
   !      !

   !      ! expHt=zero
   !      ! Ht=zero
   !      ! do iit=1,it-1
   !      !    !
   !      !    time = t_grid(iit)          
   !      !    call get_Hk_t(Hk,ik,time)
   !      !    !
   !      !    expHt_tmp = matmul(Hk,Rhop(iit,:,:))
   !      !    expHt_tmp = matmul(Rhop_dag(iit,:,:),expHt_tmp)
   !      !    Ht(1:Ns,1:Ns) = expHt_tmp
   !      !    !
   !      !    expHt_tmp = matmul(Hk,Qhop(iit,:,:))
   !      !    expHt_tmp = matmul(Rhop_dag(iit,:,:),expHt_tmp)
   !      !    Ht(1:Ns,Ns+1:2*Ns) = expHt_tmp
   !      !    !
   !      !    expHt_tmp = matmul(Hk,Rhop(iit,:,:))
   !      !    expHt_tmp = matmul(Qhop_dag(iit,:,:),expHt_tmp)  
   !      !    Ht(Ns+1:2*Ns,1:Ns) = expHt_tmp
   !      !    !
   !      !    expHt_tmp = matmul(Hk,Qhop(iit,:,:))
   !      !    expHt_tmp = matmul(Qhop_dag(iit,:,:),expHt_tmp)
   !      !    Ht(Ns+1:2*Ns,Ns+1:2*Ns) = expHt_tmp
   !      !    !
   !      !    expHt = expHt + Ht*tstep*0.5d0
   !      !    !
   !      !    time = t_grid(iit+1)          
   !      !    call get_Hk_t(Hk,ik,time)
   !      !    !
   !      !    expHt_tmp = matmul(Hk,Rhop(iit+1,:,:))
   !      !    expHt_tmp = matmul(Rhop_dag(iit+1,:,:),expHt_tmp)
   !      !    Ht(1:Ns,1:Ns) = expHt_tmp
   !      !    !
   !      !    expHt_tmp = matmul(Hk,Qhop(iit+1,:,:))
   !      !    expHt_tmp = matmul(Rhop_dag(iit+1,:,:),expHt_tmp)
   !      !    Ht(1:Ns,Ns+1:2*Ns) = expHt_tmp
   !      !    !
   !      !    expHt_tmp = matmul(Hk,Rhop(iit+1,:,:))
   !      !    expHt_tmp = matmul(Qhop_dag(iit+1,:,:),expHt_tmp)  
   !      !    Ht(Ns+1:2*Ns,1:Ns) = expHt_tmp
   !      !    !
   !      !    expHt_tmp = matmul(Hk,Qhop(iit+1,:,:))
   !      !    expHt_tmp = matmul(Qhop_dag(iit+1,:,:),expHt_tmp)
   !      !    Ht(Ns+1:2*Ns,Ns+1:2*Ns) = expHt_tmp

   !      !    expHt = expHt + Ht*tstep*0.5d0
   !      !    !
   !      ! end do


        
   !      if(it.gt.1) then
   !         time = t_grid(it)          
   !         call get_Hk_t(Hk,ik,time)
   !         !
   !         expHt_tmp = matmul(Hk,Rhop(it,:,:))
   !         expHt_tmp = matmul(Rhop_dag(it,:,:),expHt_tmp)
   !         Ht(1:Ns,1:Ns) = expHt_tmp
   !         !
   !         expHt_tmp = matmul(Hk,Qhop(it,:,:))
   !         expHt_tmp = matmul(Rhop_dag(it,:,:),expHt_tmp)
   !         Ht(1:Ns,Ns+1:2*Ns) = expHt_tmp
   !         !
   !         expHt_tmp = matmul(Hk,Rhop(it,:,:))
   !         expHt_tmp = matmul(Qhop_dag(it,:,:),expHt_tmp)  
   !         Ht(Ns+1:2*Ns,1:Ns) = expHt_tmp
   !         !
   !         expHt_tmp = matmul(Hk,Qhop(it,:,:))
   !         expHt_tmp = matmul(Qhop_dag(it,:,:),expHt_tmp)
   !         Ht(Ns+1:2*Ns,Ns+1:2*Ns) = expHt_tmp
   !         !
   !         !
   !         expHt = expHt + Ht*tstep*0.5d0
   !         expHt_ = expHt_ - Ht*tstep*0.5d0
   !         !
   !         !
   !         time = t_grid(it-1)
   !         call get_Hk_t(Hk,ik,time)
   !         !
   !         expHt_tmp = matmul(Hk,Rhop(it-1,:,:))
   !         expHt_tmp = matmul(Rhop_dag(it-1,:,:),expHt_tmp)
   !         Ht(1:Ns,1:Ns) = expHt_tmp
   !         !
   !         expHt_tmp = matmul(Hk,Qhop(it-1,:,:))
   !         expHt_tmp = matmul(Rhop_dag(it-1,:,:),expHt_tmp)
   !         Ht(1:Ns,Ns+1:2*Ns) = expHt_tmp
   !        !
   !         expHt_tmp = matmul(Hk,Rhop(it-1,:,:))
   !         expHt_tmp = matmul(Qhop_dag(it-1,:,:),expHt_tmp)
   !         Ht(Ns+1:2*Ns,1:Ns) = expHt_tmp
   !         !
   !         expHt_tmp = matmul(Hk,Qhop(it-1,:,:))
   !         expHt_tmp = matmul(Qhop_dag(it-1,:,:),expHt_tmp)
   !         Ht(Ns+1:2*Ns,Ns+1:2*Ns) = expHt_tmp
   !         !
   !         !
   !         expHt = expHt + Ht*tstep*0.5d0
   !         !
   !         if(ik.eq.10) write(*,*) it,it-1
   !         !
   !      end if
   !     !
   !     ! if(ik.eq.10) write(501,'(20F18.10)') t_grid(it),dreal(expHt(:,:))
   !     ! if(ik.eq.10) write(502,'(20F18.10)') t_grid(it),dimag(expHt(:,:))
   !     ! if(ik.eq.10) write(503,'(20F18.10)') t_grid(it),dreal(Vt_dag(it,:,:))
   !     ! if(ik.eq.10) write(504,'(20F18.10)') t_grid(it),dimag(Vt_dag(it,:,:))
   !     ! if(ik.eq.10) write(505,'(20F18.10)') t_grid(it),dreal(Vtmp(:,:))
   !     ! if(ik.eq.10) write(506,'(20F18.10)') t_grid(it),dimag(Vtmp(:,:))
   !     !
   !     ! if(ik.eq.10) write(501,'(20F18.10)') dreal(expHt)
   !     ! if(ik.eq.10) write(502,'(20F18.10)') dimag(expHt)
   !     !
   !     !
   !     expHt_v = expHt
   !     call matrix_diagonalize(expHt_v,expHt_e)
   !     !
   !     !
   !     do is=1,2*Ns
   !        do js=1,2*Ns
   !           Vt(it,is,js)=zero
   !           Vt_dag(it,is,js)=zero
   !           do iis=1,2*Ns
   !              Vt(it,is,js)=Vt(it,is,js)+expHt_v(is,iis)*conjg(expHt_v(js,iis))*exp(xi*expHt_e(iis))
   !              Vt_dag(it,is,js)=Vt_dag(it,is,js)+expHt_v(is,iis)*conjg(expHt_v(js,iis))*exp(-xi*expHt_e(iis))
   !           end do
   !         end do
   !      end do
   !      !
   !      !
   !      ! if(ik.eq.10) write(503,'(20F18.10)') t_grid(it),dreal(Vt_dag(it,:,:))
   !      ! if(ik.eq.10) write(504,'(20F18.10)') t_grid(it),dimag(Vt_dag(it,:,:))
   !      ! if(ik.eq.10) write(505,'(20F18.10)') t_grid(it),dreal(Vt(it,:,:))
   !      ! if(ik.eq.10) write(506,'(20F18.10)') t_grid(it),dimag(Vt(it,:,:))
   !      !
   !   end do
   !   !
   ! end subroutine get_time_evolution_operators_superc
   !
   !
