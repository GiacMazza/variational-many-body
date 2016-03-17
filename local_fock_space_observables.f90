! Local density !
function local_density(cc,ca) result(ni)
  real(8),dimension(Ns,nFock,nFock) :: cc,ca
  real(8),dimension(nFock,nFock) :: Id
  real(8),dimension(Ns,nFock,nFock) :: ni
  integer                        :: i,ispin,iorb,istate
  !
  do ispin=1,2
     do iorb=1,Norb
        istate=index(ispin,iorb)
        ni(istate,:,:)=matmul(cc(istate,:,:),ca(istate,:,:))
     end do
  end do
  !
end function local_density


function local_density_matrix(cc,ca) result(ni)
  real(8),dimension(Ns,nFock,nFock) :: cc,ca
  real(8),dimension(nFock,nFock) :: Id
  real(8),dimension(Ns,Ns,nFock,nFock) :: ni
  integer                        :: i,ispin,iorb,istate,jspin,jorb,jstate
  !
  do ispin=1,2
     do jspin=1,2
        do iorb=1,Norb
           do jorb=1,Norb
              istate=index(ispin,iorb)
              jstate=index(jspin,jorb)
              ni(istate,jstate,:,:)=matmul(cc(istate,:,:),ca(jstate,:,:))
           end do
        end do
     end do
  end do
  !
end function local_density_matrix


function local_doubly(cc,ca) result(di)
  real(8),dimension(Ns,nFock,nFock) :: cc,ca
  real(8),dimension(nFock,nFock) :: Id,tmp_up,tmp_dw
  real(8),dimension(Norb,nFock,nFock) :: di
  integer                        :: i,ispin,iorb,istate
  !
  do iorb=1,Norb
     !
     ispin=1
     istate=index(ispin,iorb)
     tmp_up = matmul(cc(istate,:,:),ca(istate,:,:))
     ispin=2
     istate=index(ispin,iorb)
     tmp_dw = matmul(cc(istate,:,:),ca(istate,:,:))
     !
     di(iorb,:,:) = matmul(tmp_up,tmp_dw)
     !
  end do
  !    
end function local_doubly


function local_density_density_orb(cc,ca) result(ddi)
  real(8),dimension(Ns,nFock,nFock) :: cc,ca
  real(8),dimension(nFock,nFock) :: Id,tmp_ni,tmp_nj
  real(8),dimension(Norb,Norb,nFock,nFock) :: ddi
  integer                        :: i,ispin,iorb,istate,jstate,jorb
  !
  do iorb=1,Norb
     do jorb=1,Norb
        ddi(iorb,jorb,:,:)=0.d0
        tmp_ni=0.d0
        tmp_nj=0.d0
        do ispin=1,2
           istate=index(ispin,iorb)
           jstate=index(ispin,jorb)
           tmp_ni=tmp_ni+matmul(cc(istate,:,:),ca(istate,:,:))
           tmp_nj=tmp_nj+matmul(cc(jstate,:,:),ca(jstate,:,:))             
        end do
        ddi(iorb,jorb,:,:)=matmul(tmp_ni,tmp_nj)
     end do
  end do
  !    
end function local_density_density_orb





function local_spin_flip(cc,ca) result(sfi)
  real(8),dimension(Ns,nFock,nFock) :: cc,ca
  real(8),dimension(nFock,nFock) :: Id,tmp
  real(8),dimension(Norb,Norb,nFock,nFock) :: sfi
  integer                        :: i,ispin,iorb,istate,jstate,jorb
  integer                        :: is_up,is_dn,js_up,js_dn
  !
  sfi=0.d0
  do iorb=1,Norb
     do jorb=1,Norb
        !
        is_up = index(1,iorb)
        is_dn = index(2,iorb)
        !
        js_up = index(1,jorb)
        js_dn = index(2,jorb)
        !             
        if(iorb.ne.jorb) then
           tmp=matmul(CC(js_dn,:,:),CA(js_up,:,:))
           tmp=matmul(CA(is_dn,:,:),tmp)
           tmp=matmul(CC(is_up,:,:),tmp)
           sfi(iorb,jorb,:,:) = tmp
        end if
     end do
  end do
  !    
end function local_spin_flip



function local_pair_hopping(cc,ca) result(phi)
  real(8),dimension(Ns,nFock,nFock) :: cc,ca
  real(8),dimension(nFock,nFock) :: Id,tmp
  real(8),dimension(Norb,Norb,nFock,nFock) :: phi
  integer                        :: i,ispin,iorb,istate,jstate,jorb
  integer                        :: is_up,is_dn,js_up,js_dn
  !
  phi=0.d0
  do iorb=1,Norb
     do jorb=1,Norb
        !
        is_up = index(1,iorb)
        is_dn = index(2,iorb)
        !
        js_up = index(1,jorb)
        js_dn = index(2,jorb)
        !             
        if(iorb.ne.jorb) then
           tmp=matmul(CA(js_dn,:,:),CA(js_up,:,:))
           tmp=matmul(CC(is_dn,:,:),tmp)
           tmp=matmul(CC(is_up,:,:),tmp)
           phi(iorb,jorb,:,:) = tmp
        end if
     end do
  end do
  !    
end function local_pair_hopping


function local_sc_order(cc,ca) result(phi)
  real(8),dimension(Ns,nFock,nFock) :: cc,ca
  real(8),dimension(nFock,nFock) :: Id,tmp
  real(8),dimension(Ns,Ns,nFock,nFock) :: phi
  integer                        :: i,ispin,iorb,istate,jstate,jorb
  integer                        :: is_up,is_dn,js_up,js_dn
  !
  phi=0.d0
  do istate=1,NS
     do jstate=1,NS
        phi(istate,jstate,:,:) = matmul(CC(istate,:,:),CC(jstate,:,:))
     end do
  end do
  !    
end function local_sc_order

