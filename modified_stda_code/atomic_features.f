      ! This file provides a implemation for atomic reactivity description
      ! The basic idea is to atomistic weighted averaging molecular orbitals
      ! energy according to the electronic population as well as bolzman
      ! sampled orbital energies.


      subroutine atomic_frontier_orbitals(ncent,nmo,nbf,foccu,
     .         eps,iaoat,c,ehoao,eluao)
      use, intrinsic :: ieee_arithmetic
      implicit none
      real*8,intent(in) :: foccu(nmo),c(nbf*nmo),eps(nmo)
      integer,intent(in):: nmo,nbf,ncent,iaoat(nbf)
      integer :: i,j,jj,k,kk,hindex,lindex,m
      real*8 occ,virt,ps,ehomo,elumo,fermi,weight,value,tmp,
     . tmp1,tmp2,tmp3,tmp4,bposum,bpvsum
      real*8,parameter :: ev=27.2113961318065d0 ! eh to ev
      real*8,allocatable :: S(:,:),po(:,:),pv(:,:),
     .weighted_sum(:),sumofweight(:),weighted_sum2(:),sumofweight2(:)
      real*8,intent(out) :: ehoao(ncent),eluao(ncent)



      allocate(S(nbf,nbf))
      ! get overlap matrix
      call ovmat(nbf,c,nmo,S)

      ! get frontional orbial
      ehomo=ieee_value(ehomo, ieee_quiet_nan)
      elumo=ieee_value(elumo, ieee_quiet_nan)
      do i=1,nmo
        occ=foccu(i)
        if(occ.gt.1d-5)then
          ehomo=eps(i) ! find homo
          hindex=i
        endif
        virt=1.000-occ
        if((virt.gt.1d-5).and.(ieee_is_nan(elumo)))then
          elumo=eps(i) ! find lumo
          lindex=i
        endif
      enddo

      if((ieee_is_nan(ehomo)).or.(ieee_is_nan(elumo)))then
      deallocate(S)
        if(.not.ieee_is_nan(ehomo))then
          do m=1,ncent
            ehoao(m)=ehomo*ev
          enddo
        else
          do m=1,ncent
            ehoao(m)=-1d10
          enddo
        endif
        if(.not.ieee_is_nan(elumo))then
          do m=1,ncent
            eluao(m)=elumo*ev
          enddo
        else
          do m=1,ncent
            eluao(m)=1d10
          enddo
        endif
      return
      endif

      ! print out fermi level
      fermi = 0.5000*(ehomo+elumo)
      write(*,'(A,f8.2,A)') "Fermi level = ",fermi*ev,"eV"
      write(*,'(A,f8.2,A)') "HL-gap = ",(elumo-ehomo)*ev,"eV"


      ! get Mulliken pop
      allocate(po(ncent,nmo),pv(ncent,nmo))
      po=0
      pv=0
      do i=1,nmo ! loop MO
        ! occ.
        occ=foccu(i)
        if(occ.gt.1d-5)then
          do j=1,nbf ! loop basis
            jj=iaoat(j) ! get atom index
            do k=1,j-1 ! non-diagonal elements
              kk=iaoat(k)
              ps=occ*c(nbf*(i-1)+j)*c(nbf*(i-1)+k)*S(k,j)
              po(kk,i)=po(kk,i)+ps
              po(jj,i)=po(jj,i)+ps
            enddo
            ps=S(j,j)*c(nbf*(i-1)+j)*c(nbf*(i-1)+j)*occ
            po(jj,i)=po(jj,i)+ps ! diagonal elements
          enddo
        endif
        ! virt.
        virt=1.000-occ
        if(virt.gt.1d-5)then
          do j=1,nbf
            jj=iaoat(j)
            do k=1,j-1 ! non-diagonal elements
              kk=iaoat(k)
              ps=virt*c(nbf*(i-1)+j)*c(nbf*(i-1)+k)*S(j,k)
              pv(kk,i)=pv(kk,i)+ps
              pv(jj,i)=pv(jj,i)+ps
            enddo
            ps=S(j,j)*c(nbf*(i-1)+j)*c(nbf*(i-1)+j)*virt
            pv(jj,i)=pv(jj,i)+ps ! diagonal elements
          enddo
        endif
      enddo
      deallocate(S)

      ! get average for each atom

      allocate(weighted_sum(ncent),sumofweight(ncent),
     .         weighted_sum2(ncent),sumofweight2(ncent))
      weighted_sum = 0
      sumofweight= 0
      weighted_sum2 = 0
      sumofweight2= 0
      do i=1,nmo
        ! occ.
        occ=foccu(i)
        if(occ.gt.1d-5)then
          value=(fermi-eps(i))*ev
          tmp=value**2
          do m=1,ncent
            weight=po(m,i)/(tmp*occ)
            if(weight.gt.0) then
              weighted_sum(m)=weighted_sum(m)+weight*value
              sumofweight(m) = sumofweight(m)+weight
            endif
          enddo
        endif
        ! virt.
        virt=1.0000-occ
        if(virt.gt.1d-5)then
          value=(eps(i)-fermi)*ev
          tmp=value**2
          do m=1,ncent
            weight=pv(m,i)/(tmp*virt)
            if(weight.gt.0) then
              weighted_sum2(m)=weighted_sum2(m)+weight*value
              sumofweight2(m) = sumofweight2(m)+weight
            endif
          enddo
        endif
      enddo



      ! get and print out res.
      write(*,*)
      write(*,*) "atom  HOAO (eV)    LUAO (eV)    %HOMO
     .      %HOMO-1     %LUMO      %LUMO+1"
      ! get res.
      ehoao=0
      eluao=0
      do m=1,ncent
        ehoao(m)= weighted_sum(m)/sumofweight(m)
        ehoao(m) = fermi*ev-ehoao(m)/sum(po(m,:))

        eluao(m)=weighted_sum2(m)/sumofweight2(m)
        eluao(m) = fermi*ev+eluao(m)/sum(pv(m,:))

        ! get ele % in FMO
        tmp1=po(m,hindex)
        tmp2=po(m,hindex-1)
        tmp3=pv(m,lindex)
        tmp4=pv(m,lindex+1)
        write(*,'(I3,6(F12.4))') m,ehoao(m),eluao(m),tmp1,tmp2,tmp3,tmp4
      enddo
      open(unit=26,file='charges')
      do m=1,ncent
      write(26,'(f12.4)') sum(po(m,:)*2)
      enddo
      close(26)

      open(unit=25,file='pop_fmo.csv')
      write(25,'(a)') "atom, %HOMO,
     .     %HOMO-1,   %LUMO,     %LUMO+1"
      do m=1,ncent
        tmp1=po(m,hindex)
        tmp2=po(m,hindex-1)
        tmp3=pv(m,lindex)
        tmp4=pv(m,lindex+1)
        write(25,'(I3,a)', advance="no") m,','
        write(25,'(f12.4,a)', advance="no") po(m,hindex),','
        write(25,'(f12.4,a)', advance="no") po(m,hindex-1),','
        write(25,'(f12.4,a)', advance="no") pv(m,lindex),','
        write(25,'(f12.4)') po(m,lindex+1)
      enddo
      close(25)

      deallocate(po,pv,weighted_sum,
     . sumofweight,weighted_sum2,sumofweight2)
      endsubroutine atomic_frontier_orbitals


      subroutine foccu(nmo,occ,imethod,ccspin,focca,foccb)
      ! get alpha and beta fractional occupation in each MO
      implicit none
      real*8,intent(in) :: occ(nmo)
      integer,intent(in) :: nmo,imethod,ccspin(nmo)
      integer i,j,nalpha
      real*8 nela,nelb
      real*8,intent(out) :: foccb(nmo),focca(nmo)

      write(*,*) "Getting fractional occupation in each MO..."

      focca=0.0d0
      foccb=0.0d0
      nela=0.0d0
      nelb=0.0d0
      nalpha=0
      if(imethod.ne.2) then ! For RKS: focca=fooccb=occ/2
        do i=1,nmo
          focca(i)=0.5000*occ(i)
          nela=nela+focca(i)
          foccb(i)=0.5000*occ(i)
          nelb=nelb+foccb(i)
        enddo
      else ! For UKS: read ccspin
        do i=1,nmo
          if (ccspin(i).eq.1) then ! alpha
            focca(i)=occ(i)
            nela=nela+focca(i)
          else
          if (nalpha.eq.0)then
          nalpha=i ! alpha and beta are stored separately
          endif
          foccb(i-nalpha+1)=occ(i)
          nelb=nelb+foccb(i)
          endif
        enddo
      endif

      write(*,'(A,f4.1)') "alpha electron number: ",nela
      write(*,'(A,f4.1)') "beta electron number: ",nelb
      write(*,'(A,f4.1)') "tot electron number: ",nela+nelb

      endsubroutine foccu

      subroutine ovmat(nbf,c,nmo,S)
      implicit none
      real*8,intent( in ) :: c(nbf*nmo)
      integer,intent(in) :: nbf,nmo
      real*8, allocatable :: ovvec(:)
      real*8,intent(out) :: S(nbf,nbf)

      allocate(ovvec(nbf*(nbf+1)/2))
      ! read overlap vector
      open(unit=40,file='sint',form='unformatted',status='old')
      read(40)ovvec
      close(40)

      ! blow up vector to matrix
      call blow(nbf,ovvec,S)

      deallocate(ovvec)
      endsubroutine ovmat


      subroutine eps_check(nmo,eps,occ,focca,foccb)
      implicit none
      integer, intent( in ) :: nmo
      real*8, intent( in ) :: eps(nmo),occ(nmo),
     .  focca(nmo),foccb(nmo)
      integer i
      call header("eps",0)
      write(*,*)"       MO     Occupation  alpha   beta
     .    Energy/Eh   Energy/eV"
      do i=1,nmo
        write(*,'(3x,i6,3x,f8.2,3x,i6,3x,i6,5x,f8.4,5x,f8.4)')
     .   i,occ(i),int(focca(i)),int(foccb(i)),
     . eps(i),eps(i)*27.2113961318065d0
      enddo
      endsubroutine

      subroutine dmat(nbf,focc,nmo,c,P)
      implicit none
      real*8, intent( in ) :: c(nbf*nmo),focc(nmo)
      integer, intent( in ) :: nmo,nbf
      real*8, intent( out ) :: P(nbf,nbf)
      integer i,j,k

      do i=1,nbf
        do j=1,i
          P(i,j)=0.0d0
          P(j,i)=0.0d0
          do k=1,nmo
            P(i,j)=P(i,j)+focc(k)*c(nbf*(k-1)+i)*c(nbf*(k-1)+j)
          enddo
          P(j,i)=P(i,j)
        enddo
      enddo
      endsubroutine dmat



      subroutine print_MLfeaturefile(out,ncent,co,
     .       ehoaoa,eluaoa,
     .       ehoaob,eluaob)
      implicit none
      integer, intent(in) :: out,ncent
      real*8,intent(in) :: ehoaoa(ncent),eluaoa(ncent),
     .ehoaob(ncent),eluaob(ncent),co(ncent,4)
      integer :: i
      write(out,'(a)') 'Atom,   '//
     .                 'HOAO_a (eV),   LUAO_a (eV),   '//
     .                 'HOAO_b (eV),   LUAO_b (eV)'
      do i = 1, ncent
        write(out,'(i2,a)', advance="no") int(co(i,4)),','
        write(out,'(f14.4,a)', advance="no") ehoaoa(i), ','
        write(out,'(f14.4,a)', advance="no") eluaoa(i), ','
        write(out,'(f14.4,a)', advance="no") ehoaob(i), ','
        write(out,'(f14.4)') eluaob(i)
      enddo
      write(*,*)
      endsubroutine

      subroutine get_embo(ncent,nbf,nmo,c,co,focca,foccb,
     . eps,iaoat,mbo,embo)
      implicit none
      real*8, intent( in ) :: c(nbf*nmo),co(ncent,4),
     . focca(nmo),foccb(nmo),eps(nmo)
      integer, intent( in ) :: nmo,nbf,ncent,iaoat(nbf)
      integer i,j,k,m,jj,kk
      real*8  tmp1,tmp2,rab,occ,ps,tmp,occe
      real*8, intent( out ) :: embo(ncent,ncent),mbo(ncent,ncent)
      real*8, allocatable :: psa(:,:),psb(:,:),
     . dmata(:,:),dmatb(:,:),S(:,:),alpha(:,:),po(:,:)
      allocate(dmata(nbf,nbf),dmatb(nbf,nbf),
     .psa(nbf,nbf),psb(nbf,nbf),S(nbf,nbf),alpha(ncent,ncent),
     .po(ncent,ncent))

      psa= 0.0d0
      psb =0.0d0
      ! construct density matrix
      call dmat(nbf,focca,nmo,c,dmata)
      call dmat(nbf,foccb,nmo,c,dmatb)
      ! overlap matrix
      call ovmat(nbf,c,nmo,S)

      !get alpha
      alpha=0
      do i=1,nmo ! loop MO
        ! occ.
        occ=focca(i)+foccb(i)
        if(occ.gt.1d-5)then
          do j=1,nbf ! loop basis
            jj=iaoat(j) ! get atom index
            do k=1,j-1 ! non-diagonal elements
              kk=iaoat(k)
              tmp=occ*c(nbf*(i-1)+j)*c(nbf*(i-1)+k)*S(k,j)
              alpha(kk,jj)=alpha(kk,jj)+tmp
              alpha(jj,kk)=alpha(jj,kk)+tmp
            enddo
          enddo
        endif
      enddo

      po = 0.0d0
      do i=1,nmo
        occ= focca(i) +foccb(i)
        if (occ.gt.0.0001d0)then
            occe = occ*eps(i)*27.2113961318065d0
            do j =1,nbf
                jj=iaoat(j)
                do k=1,j-1
                    kk=iaoat(k)
                    if (kk.ne.jj)then
                        tmp = c(nbf*(i-1)+j)*c(nbf*(i-1)+k)*occe*S(j,k)
                        po(kk,jj)=po(kk,jj)+tmp
                    endif
                enddo
            enddo
        endif
      enddo


      ! construct ps matrix
      call dsymm('l','l',nbf,nbf,1.d0,dmata,nbf,S,nbf,0.d0,psa,nbf)
      call dsymm('l','l',nbf,nbf,1.d0,dmatb,nbf,S,nbf,0.d0,psb,nbf)
      deallocate(dmata,dmatb,S)

      mbo = 0
      embo=0
      do i=1,ncent ! atom i
        do j=1,i ! atom j
          tmp1 = 0.0d0
          rab=sum((co(i,1:3) - co(j,1:3))**2)
          if(rab < 100.0d0)then
            do k=1,nbf
              do m=1,nbf
                if (iaoat(k).eq.i.and.iaoat(m).eq.j)then
                  tmp1 = tmp1 + psa(k,m)*psa(m,k)+psb(k,m)*psb(m,k)
                endif
              enddo
            enddo
          endif
          mbo(i,j) = 2 * tmp1
          mbo(j,i) = 2 * tmp1
          if (mbo(i,j).gt.0.4d0)then
          embo(i,j)= po(i,j)* alpha(i,j)
          embo(j,i)= po(j,i)* alpha(j,i)
          endif
        enddo
      enddo

      deallocate(psa,psb,alpha,po)
      write(*,*) "at1      at2     MBO    eMBO"
      do i=1,ncent
        do j = 1, i-1
            if(mbo(j,i).gt.0.2d0)then
            write(*,'(2(I3),2(F8.4))') j,i,mbo(j, i),-embo(j, i)
            endif
        enddo
      enddo
      endsubroutine get_embo

      subroutine print_embofile(out,ncent,mbo,embo)
      implicit none
      integer, intent(in) :: out,ncent
     .
      real*8,intent(in) :: mbo(ncent,ncent),embo(ncent,ncent)
      integer :: i, j
      do i = 1, ncent
        do j = 1, i-1
            if (mbo(j, i) > 0.2)then
                write(out, '(2(I3),2(F8.4))') j, i,mbo(j, i),-embo(j, i)
            endif
        end do
      end do
      endsubroutine