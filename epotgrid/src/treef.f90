! ***************************************
!      treef.f90 for epotgrid 
!  7 Dec. 2015    written by D.KAWATA
! ***************************************
! *** Definition of treeforce() ***
! * calculate acceraration *


subroutine treeforce(ngp,npt,ndmt)
      use gcdp_const
      use gcdp_system
      use gcdp_baryon
      use gcdp_btree
      use gcdp_kernel
      use gcdp_dm
      use gcdp_dmtree
      use grid_particle


      implicit none
      include 'mpif.h'

      integer,intent(in) :: ngp,npt,ndmt
      integer i,ip,level,is,np
      integer nval,isend,npj,srank,rrank,irecv
      integer maxntw
! * Number of Notfinished particle, temp *      
      integer nlist,tnlist
! * Node searching for each Particle *      
      integer nd,pn,pnj,ntf,flagp
      double precision r2,ir2,rij,ir1,ir3
! * for calculate gradients *
      double precision dwi,dwj
! * for test far away *            
      double precision theta2
      double precision l_r2
      double precision xij,yij,zij,tij
      double precision dphiir,dphijr,si,sj,hsi,hsj
      double precision phii,phij
! *** potential and correction adaptive softening ***
      double precision dphidrr
! *** information particles need communication ***
      integer iproc,snval,rnval,snival,npjs
      integer ncomp,ncompt,ndivcom,idivcom
      integer ncomptmax,npprocmax,ndivcompt
      integer ncompi
! *** for pseudo node 
      integer ndp,pnodess
! *** flag for tree ***
      integer flagtr
! *** for test file *** 
      character fileo*60
! for work for communication fixed size array
      integer pncompt(0:MNB-1),proccompt(0:MNB-1)
      integer ierr
! pn_nfp can be bigger than ng+ns in this proc, when receiving the particles
      integer pn_nfp(0:MNB-1)
! for work allocatable
      integer,allocatable :: npltf(:),flago(:)
      integer,allocatable :: pncomp(:),flagproc(:),flagcom(:)
      integer,allocatable :: npproc(:),isproc(:),ieproc(:),cproc(:)
      integer,allocatable :: nppnd(:),nextpnd(:),daupnd(:) &
       ,ndpnd(:),idpnd(:)
      integer,allocatable :: idisp(:),jjlen(:),ireqs(:),npjr(:)
      integer,allocatable :: list(:),node(:)
      integer,allocatable :: tivs(:),tivr(:)
      double precision,allocatable :: xijp(:),yijp(:),zijp(:) &
       ,r2p(:)
      double precision,allocatable :: xpnd(:),ypnd(:),zpnd(:) &
       ,lpnd(:),deltapnd(:),mpnd(:),cxpnd(:),cypnd(:),czpnd(:) &
       ,trbuf(:)
      double precision,allocatable :: tdvs(:),tdvr(:)

! *** Basic Value ***
      theta2 = THETA*THETA
      nval=8

! maximum number of tree walk
      maxntw=(npt+ndmt)*3+nprocs*MAXNODESEND

      if(myrank.eq.0) then
        write(6,*) ' treef(), ngp,npt,ndmt=',ngp,npt,ndmt
      endif

!allocate work space for MPI
      allocate(idisp(0:nprocs))
      allocate(jjlen(0:nprocs))
      allocate(ireqs(0:nprocs))
      allocate(npproc(0:nprocs))
      allocate(isproc(0:nprocs))
      allocate(ieproc(0:nprocs))
      allocate(cproc(0:nprocs))
      allocate(npjr(0:nprocs))

      allocate(flagproc(0:ngp-1))
      allocate(flagcom(0:ngp-1))
      allocate(pncomp(0:ngp-1))

      flagtr=0

   71 allocate(tdvr(0:nval*ngp-1))

      do i=0,ngp-1
        pn=i
        pn_nfp(i)=pn
        tdvr(i)=xp(pn)
        tdvr(i+ngp)=yp(pn)
        tdvr(i+ngp*2)=zp(pn)
        tdvr(i+ngp*3)=hp(pn)
        flagproc(i)=myrank
        flagcom(i)=-1
      enddo

! *** initialization for a?_p ***
      do i=ngp*4,nval*ngp-1
        tdvr(i)=0.0d0
      enddo

      npj=ngp
      ncomp=0
      ncompt=0
      iproc=0
      idivcom=0
      ndivcom=-1
      do i=0,nprocs-1
        npproc(i)=0
      enddo

      if(npt.eq.0) then
        flagtr=1
      endif

      if(flagtr.eq.0) then
        if(myrank.eq.0) then
          write(6,*) ' start treeforce from baryon.'
        endif
        goto 70
      else if(flagtr.eq.1.and.ndmt.gt.0) then
        if(myrank.eq.0) then
          write(6,*) ' start treeforce from DM.'
        endif
        if(flagtr.ne.1) then
          flagtr=1
        endif
        goto 91
      else
        goto 94
      endif

   70 if(np_tr(0).eq.0) then
        goto 93
      endif
      if(iproc.ne.0.and.proc_tr(0).ne.myrank) then
        goto 93
      endif
! *** initialization ***

      allocate(list(0:npj))
      allocate(node(0:npj))
      allocate(flago(0:npj))

      do i=0,npj-1
        list(i)=i
        node(i)=0 
      enddo
      nlist=npj
      level=0

! *** calculate contribution from SPH Particle ***
   77 if(nlist.le.0) then

        deallocate(list)
        deallocate(node)
        deallocate(flago)

        goto 93
      endif
      level = level+1
      if(level.gt.maxntw) then
        write(6,*) ' Error in treef(): failure in walking tree'
        write(6,*) 'myrank,iproc.idivcom,ndivcom=' &
         ,myrank,iproc,idivcom,ndivcom
        write(fileo,'(a6,i3.3)') 'errtfb',myrank
        open(60,file=fileo,status='unknown')
        do i=0,nlist-1
          pn=list(i)
          write(60,'(3(1pE13.5),I10)') tdvr(pn),tdvr(pn+npj),tdvr(pn+npj*2),pn
        enddo
        close(60)
        stop
      endif

      allocate(npltf(0:nlist))
      allocate(xijp(0:nlist))
      allocate(yijp(0:nlist))
      allocate(zijp(0:nlist))
      allocate(r2p(0:nlist))

      ntf=0
      do i=0,nlist-1
        pn = list(i)         
        nd = node(pn)
        flago(pn)=0
        if(np_tr(nd).gt.1) then
! *** check for particles within a node ***
          xij=tdvr(pn)-cx_tr(nd)
          yij=tdvr(pn+npj)-cy_tr(nd)		
          zij=tdvr(pn+npj*2)-cz_tr(nd)
          if(dabs(xij).lt.0.6d0*l_tr(nd)) then
          if(dabs(yij).lt.0.6d0*l_tr(nd)) then
          if(dabs(zij).lt.0.6d0*l_tr(nd)) then
! *** go to daughter node ***
            flago(pn)=1
          endif
          endif
          endif
        else if(np_tr(nd).eq.1.and.iproc.eq.0) then
! skip own particle
          if(proc_tr(nd).eq.myrank.and.pn_nfp(pn).eq.pn_tr(nd)) then
            flago(pn)=-1
            node(pn)=next_tr(nd)		              
          endif
        endif

        if(flago(pn).eq.0) then
          xij=tdvr(pn)-cmx_tr(nd)
          yij=tdvr(pn+npj)-cmy_tr(nd)		
          zij=tdvr(pn+npj*2)-cmz_tr(nd)
          r2=xij*xij+yij*yij+zij*zij
          if(r2.gt.0.0d0) then
! *** register the particles need distance check ***
            npltf(ntf)=pn
! *** store distance from the center of the mass ***
            xijp(ntf)=xij
            yijp(ntf)=yij
            zijp(ntf)=zij
            r2p(ntf)=r2
            ntf=ntf+1
          else if(np_tr(nd).eq.1) then
            node(pn)=next_tr(nd)		
          else
            flago(pn)=1
          endif
        endif
      enddo
! *** for step 0 or continue ***
      do i=0,ntf-1
        pn=npltf(i)
        nd=node(pn)
        if(np_tr(nd).eq.1) then
          if(proc_tr(nd).eq.myrank) then
            flago(pn)=2
          else
            flago(pn)=3
          endif
        else
          l_r2=(l_tr(nd)+delta_tr(nd))*(l_tr(nd)+delta_tr(nd))/r2p(i)
          if(l_r2.le.theta2) then
            flago(pn)=3
          else
! *** go to daughter ***
            flago(pn)=1
          endif
        endif
      enddo

! *** tree force calculation
      do i=0,ntf-1
        pn=npltf(i)
        nd=node(pn)
        if(flago(pn).eq.2) then
          if(iproc.eq.0.or.nd.gt.nodese_tr) then
            pnj = pn_tr(nd)
            rij=dsqrt(r2p(i))
            ir2=1.0d0/r2p(i)
            ir3=ir2/rij
! *** dphi/dr(hi) and phii ***
            hsi=tdvr(pn+npj*3)
            si=rij/hsi
            if(si.lt.1.0d0) then
              is=int(si*dnktab)
              if(is.lt.0) then
                is=0
              else if(is.ge.NKTAB) then
                is=NKTAB-1
              endif
! dphi/dr/r
              dphiir=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is)) &
               *(si-s_tb(is))*dnktab
              dphiir=dphiir/(hsi**3)
! phii
              phii=hphi_tb(is)+(hphi_tb(is+1)-hphi_tb(is)) &
               *(si-s_tb(is))*dnktab
              phii=phii/hsi
            else
              dphiir=ir3
              phii=-1.0d0/rij
            endif
! *** dphi/dr(hj) and phij ***
            hsj=h_p(pnj)
            sj=rij/hsj
            if(sj.lt.1.0d0) then
              is=int(sj*dnktab)
              if(is.lt.0) then
                is=0
              else if(is.ge.NKTAB) then
                is=NKTAB-1
              endif
! dphi/dr/r
              dphijr=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is)) &
               *(sj-s_tb(is))*dnktab
              dphijr=dphijr/(hsj**3)
! phij
              phij=hphi_tb(is)+(hphi_tb(is+1)-hphi_tb(is)) &
               *(sj-s_tb(is))*dnktab
              phij=phij/hsj
            else
              dphijr=ir3
              phij=-1.0d0/rij
            endif
! *** force calculation ***
! potential
            tdvr(pn+npj*4)=tdvr(pn+npj*4)+0.50d0*m_p(pnj)*G*(phii+phij)
! *** dvx_p,dvy_p,dvz_p ***
            dphidrr=m_p(pnj)*G*0.5d0*(dphiir+dphijr)
            tdvr(pn+npj*5)=tdvr(pn+npj*5)-xijp(i)*dphidrr
            tdvr(pn+npj*6)=tdvr(pn+npj*6)-yijp(i)*dphidrr
            tdvr(pn+npj*7)=tdvr(pn+npj*7)-zijp(i)*dphidrr
          endif
! * update node *
          node(pn)=next_tr(nd)
        else if(flago(pn).eq.3) then
          if(iproc.eq.0.or.nd.gt.nodese_tr) then
            ir2=1.0d0/r2p(i)
            r2=r2p(i)
            rij=dsqrt(r2)       
            ir1=1.0d0/rij
            ir3=ir2*ir1
! *** use softening of particle i
! *** dphi/dr(hi) ***
            hsi=tdvr(pn+npj*3)
            si=rij/hsi
            if(si.lt.1.0d0) then
              is=int(si*dnktab)
              if(is.lt.0) then
                is=0
              else if(is.ge.NKTAB) then
                is=NKTAB-1
              endif
! dphi/dr/r
              dphiir=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is)) &
               *(si-s_tb(is))*dnktab
              dphiir=dphiir/(hsi**3)
! phi
              phii=hphi_tb(is)+(hphi_tb(is+1)-hphi_tb(is))*(si-s_tb(is))*dnktab
              phii=phii/hsi
            else
              dphiir=ir3
              phii=-1.0d0/rij
            endif
! *** dphi/dr(hj) and phij ***
            hsj=hm_tr(nd)
            sj=rij/hsj
            if(sj.lt.1.0d0) then
              is=int(sj*dnktab)
              if(is.lt.0) then
                is=0
              else if(is.ge.NKTAB) then
                is=NKTAB-1
              endif
! dphi/dr/r
              dphijr=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is)) &
               *(sj-s_tb(is))*dnktab
              dphijr=dphijr/(hsj**3)
! phij
              phij=hphi_tb(is)+(hphi_tb(is+1)-hphi_tb(is))*(sj-s_tb(is))*dnktab
              phij=phij/hsj
            else
              dphijr=ir3
              phij=-1.0d0/rij
            endif

! *** force calculation ***
! potential
            tdvr(pn+npj*4)=tdvr(pn+npj*4)+0.50d0*mass_tr(nd)*G*(phii+phij)
! *** dvx_p,dvy_p,dvz_p ***
            dphidrr=mass_tr(nd)*G*0.5d0*(dphiir+dphijr)
            tdvr(pn+npj*5)=tdvr(pn+npj*5)-xijp(i)*dphidrr
            tdvr(pn+npj*6)=tdvr(pn+npj*6)-yijp(i)*dphidrr
            tdvr(pn+npj*7)=tdvr(pn+npj*7)-zijp(i)*dphidrr
          endif
! *** update node
          node(pn)=next_tr(nd)
        endif
      enddo
      do i=0,nlist-1
        pn=list(i)
        if(flago(pn).eq.1) then
          nd=node(pn)
! * check if the pseudo node or not *
          if(iproc.eq.0.and.daughter_tr(nd).eq.-1) then
            if(proc_tr(nd).ne.myrank.and.flagproc(pn).ne.proc_tr(nd)) then
! *** proc_tr will be the same for different domain within the same proc
              flagproc(pn)=proc_tr(nd)
! *** store pn ***
              if(flagcom(pn).lt.0) then
                flagcom(pn)=proc_tr(nd)
                pncomp(ncomp)=pn_nfp(pn)
                ncomp=ncomp+1
              endif
              if(ncompt.lt.MNB) then
                pncompt(ncompt)=pn_nfp(pn)
                proccompt(ncompt)=proc_tr(nd)
              endif
              ncompt=ncompt+1
              npproc(proc_tr(nd))=npproc(proc_tr(nd))+1
            endif
            node(pn)=next_tr(nd)
          else
            node(pn)=daughter_tr(nd)
          endif
        endif
      enddo

      deallocate(npltf)
      deallocate(xijp)
      deallocate(yijp)
      deallocate(zijp)
      deallocate(r2p)

! * update not-finished particle list *
      tnlist = nlist
      nlist = 0
      if(iproc.eq.0.and.next_tr(0).ne.0) then
        do i=0,tnlist-1
          if(node(list(i)).gt.0) then
            list(nlist)=list(i)
            nlist=nlist+1
          else if(node(list(i)).eq.0) then
            node(list(i))=next_tr(0)
            list(nlist)=list(i)
            nlist=nlist+1
          endif
        enddo
      else
        do i=0,tnlist-1
          if(node(list(i)).gt.0) then
            list(nlist)=list(i)
            nlist=nlist+1
          endif
        enddo
      endif
      goto 77      

! ***   Calculate force contributed from DM particle ***
! * Initialization *
   91 if(np_dmtr(0).eq.0) then
        goto 93
      endif
      if(iproc.ne.0.and.proc_dmtr(0).ne.myrank) then
        goto 93
      endif

      allocate(list(0:npj))
      allocate(node(0:npj))
      allocate(flago(0:npj))

      do i=0,npj-1
        list(i)=i
        node(i)=0 
      enddo
      nlist=npj
      level=0 
  777 if(nlist.le.0) then

        deallocate(list)
        deallocate(node)
        deallocate(flago)

        goto 93
      endif
      level = level+1
      if(level.gt.maxntw) then
        write(6,*) ' Error in treef(): failure in walking dmtree'
        write(6,*) 'myrank,iproc.idivcom,ndivcom=',myrank,iproc,idivcom,ndivcom
        write(fileo,'(a7,i3.3)') 'errtfdt',myrank
        open(60,file=fileo,status='unknown')
        do i=0,nlist-1
          pn=list(i)
          write(60,'(3(1pE13.5),I10)') tdvr(pn),tdvr(pn+npj),tdvr(pn+npj*2),pn
        enddo
        close(60)
        stop
      endif

      allocate(npltf(0:nlist))
      allocate(xijp(0:nlist))
      allocate(yijp(0:nlist))
      allocate(zijp(0:nlist))
      allocate(r2p(0:nlist))

      ntf=0
      do i=0,nlist-1
        pn = list(i)         
        nd = node(pn)
        flago(pn)=0
        if(np_dmtr(nd).gt.1) then
          xij=tdvr(pn)-cx_dmtr(nd)
          yij=tdvr(pn+npj)-cy_dmtr(nd)		
          zij=tdvr(pn+npj*2)-cz_dmtr(nd)
          if(dabs(xij).lt.0.6d0*l_dmtr(nd)) then
          if(dabs(yij).lt.0.6d0*l_dmtr(nd)) then
          if(dabs(zij).lt.0.6d0*l_dmtr(nd)) then
            flago(pn)=1
          endif
          endif
          endif
        endif
        if(flago(pn).eq.0) then
          xij=tdvr(pn)-cmx_dmtr(nd)
          yij=tdvr(pn+npj)-cmy_dmtr(nd)		
          zij=tdvr(pn+npj*2)-cmz_dmtr(nd)
          r2=xij*xij+yij*yij+zij*zij
          if(r2.gt.0.0d0) then
! *** register the particles need distance calculation ***
            npltf(ntf)=pn
! *** store distance from the center of the mass ***
            xijp(ntf)=xij
            yijp(ntf)=yij
            zijp(ntf)=zij
            r2p(ntf)=r2
            ntf=ntf+1
          else if(np_dmtr(nd).eq.1) then
            node(pn)=next_dmtr(nd)		
          else
            flago(pn)=1
          endif
        endif
      enddo
! *** for step 0 or continue ***
      do i=0,ntf-1
        pn=npltf(i)
        nd=node(pn)
        if(np_dmtr(nd).eq.1) then
          if(proc_dmtr(nd).eq.myrank) then
            flago(pn)=2
          else
            flago(pn)=3
          endif
        else
          l_r2=(l_dmtr(nd)+delta_dmtr(nd))*(l_dmtr(nd)+delta_dmtr(nd))/r2p(i)
          if(l_r2.le.theta2) then
            flago(pn)=3
          else
! *** go to daughter ***
            flago(pn)=1
          endif
        endif
      enddo
! *** tree force calculation
      do i=0,ntf-1
        pn=npltf(i)
        nd=node(pn)
        if(flago(pn).eq.2) then
          if(iproc.eq.0.or.nd.gt.nodese_dmtr) then
            pnj = pn_dmtr(nd)
            rij=dsqrt(r2p(i))
            ir2=1.0d0/r2p(i)
            ir3=ir2/rij
! *** dphi/dr(hi) and dwi ***
            hsi=tdvr(pn+npj*3)
            si=rij/hsi
            if(si.lt.1.0d0) then
              is=int(si*dnktab)
              if(is.lt.0) then
                is=0
              else if(is.ge.NKTAB) then
                is=NKTAB-1
              endif
! dphi/dr/r
              dphiir=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is)) &
               *(si-s_tb(is))*dnktab
              dphiir=dphiir/(hsi**3)
! phii
              phii=hphi_tb(is)+(hphi_tb(is+1)-hphi_tb(is)) &
               *(si-s_tb(is))*dnktab
              phii=phii/hsi
            else
              dphiir=ir3
              phii=-1.0d0/rij
            endif
! *** dphi/dr(hj) ***
            hsj=h_dm(pnj)
            sj=rij/hsj
            if(sj.lt.1.0d0) then
              is=int(sj*dnktab)
              if(is.lt.0) then
                is=0
              else if(is.ge.NKTAB) then
                is=NKTAB-1
              endif
! dphi/dr/r
              dphijr=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is)) &
               *(sj-s_tb(is))*dnktab
              dphijr=dphijr/(hsj**3)
! phij
              phij=hphi_tb(is)+(hphi_tb(is+1)-hphi_tb(is)) &
               *(sj-s_tb(is))*dnktab
              phij=phij/hsj
            else
              dphijr=ir3
              phij=-1.0d0/rij
            endif
! *** force calculation ***
            tdvr(pn+npj*4)=tdvr(pn+npj*4)+0.5d0*m_dm(pnj)*G*(phii+phij)
! *** dvx_p,dvy_p,dvz_p ***
            dphidrr=m_dm(pnj)*0.5d0*G*(dphiir+dphijr)
            tdvr(pn+npj*5)=tdvr(pn+npj*5)-xijp(i)*dphidrr
            tdvr(pn+npj*6)=tdvr(pn+npj*6)-yijp(i)*dphidrr
            tdvr(pn+npj*7)=tdvr(pn+npj*7)-zijp(i)*dphidrr
          endif
! * update node *
          node(pn)=next_dmtr(nd)
        else if(flago(pn).eq.3) then
          if(iproc.eq.0.or.nd.gt.nodese_dmtr) then
            ir2=1.0d0/r2p(i)
            r2=r2p(i)
            rij=dsqrt(r2)       
            ir1=1.0d0/rij
            ir3=ir2*ir1
! *** use softening of particle i
! *** dphi/dr(hi) and phii ***
            hsi=tdvr(pn+npj*3)
            si=rij/hsi
            if(si.lt.1.0d0) then
              is=int(si*dnktab)
              if(is.lt.0) then
                is=0
              else if(is.ge.NKTAB) then
                is=NKTAB-1
              endif
! dphi/dr/r
              dphiir=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is)) &
               *(si-s_tb(is))*dnktab
              dphiir=dphiir/(hsi**3)
! phii
              phii=hphi_tb(is)+(hphi_tb(is+1)-hphi_tb(is)) &
               *(si-s_tb(is))*dnktab
              phii=phii/hsi

            else
              dphiir=ir3
              phii=-1.0d0/rij
            endif
! *** dphi/dr(hj) and phij ***
            hsj=hm_dmtr(nd)
            sj=rij/hsj
            if(sj.lt.1.0d0) then
              is=int(sj*dnktab)
              if(is.lt.0) then
                is=0
              else if(is.ge.NKTAB) then
                is=NKTAB-1
              endif
! dphi/dr/r
              dphijr=dphidr_r_tb(is)+(dphidr_r_tb(is+1)-dphidr_r_tb(is)) &
               *(sj-s_tb(is))*dnktab
              dphijr=dphijr/(hsj**3)
! phij
              phij=hphi_tb(is)+(hphi_tb(is+1)-hphi_tb(is)) &
               *(sj-s_tb(is))*dnktab
              phij=phij/hsj
            else
              dphijr=ir3
              phij=-1.0d0/rij
            endif
! *** force calculation ***
            tdvr(pn+npj*4)=tdvr(pn+npj*4)+0.5d0*mass_dmtr(nd)*G*(phii+phij)
! *** dvx_p,dvy_p,dvz_p ***
            dphidrr=mass_dmtr(nd)*G*0.5d0*(dphiir+dphijr)
            tdvr(pn+npj*5)=tdvr(pn+npj*5)-xijp(i)*dphidrr
            tdvr(pn+npj*6)=tdvr(pn+npj*6)-yijp(i)*dphidrr
            tdvr(pn+npj*7)=tdvr(pn+npj*7)-zijp(i)*dphidrr
          endif
! * update node *
          node(pn)=next_dmtr(nd)
        endif
      enddo
      do i=0,nlist-1
        pn = list(i)
        if(flago(pn).eq.1) then
          nd = node(pn)		
! * check if the pseudo node or not *
! * check if the pseudo node or not *
          if(iproc.eq.0.and.daughter_dmtr(nd).eq.-1) then
            if(proc_dmtr(nd).ne.myrank &
              .and.flagproc(pn).ne.proc_dmtr(nd)) then
! *** proc_tr will be the same for different domain within the same proc
              flagproc(pn)=proc_dmtr(nd)
! *** store pn ***
              if(flagcom(pn).lt.0) then
                flagcom(pn)=proc_dmtr(nd)
                pncomp(ncomp)=pn_nfp(pn)
                ncomp=ncomp+1
              endif
              if(ncompt.lt.MNB) then
                pncompt(ncompt)=pn_nfp(pn)
                proccompt(ncompt)=proc_dmtr(nd)
              endif
              ncompt=ncompt+1
              npproc(proc_dmtr(nd))=npproc(proc_dmtr(nd))+1
            endif
            node(pn)=next_dmtr(nd)
          else
            node(pn) = daughter_dmtr(nd)
          endif
        endif
      enddo

      deallocate(npltf)
      deallocate(xijp)
      deallocate(yijp)
      deallocate(zijp)
      deallocate(r2p)

! * update not-finished particle list *
      tnlist = nlist
      nlist = 0
      if(iproc.eq.0.and.next_dmtr(0).ne.0) then
        do i=0,tnlist-1
          if(node(list(i)).gt.0) then
            list(nlist)=list(i)
            nlist=nlist+1
          else if(node(list(i)).eq.0) then
            node(list(i))=next_dmtr(0)
            list(nlist)=list(i)
            nlist=nlist+1
          endif
        enddo
      else
        do i=0,tnlist-1
          if(node(list(i)).gt.0) then
            list(nlist)=list(i)
            nlist=nlist+1
          endif
        enddo
      endif
      goto 777        

! *** end itteration within the proc ***
! *** update variables ***
   93 if(iproc.ne.0) then
! *** sending back the results ***
! *** new snval to send the data back
        snval=4

        allocate(tivs(0:nprocs))

        do i=0,nprocs-1
! *** keep npjr: data received from each proc -> now sending***
          tivs(i)=npjr(i)
! *** ireqs=jjlen number of particles sent -> now receiving.
          npjr(i)=ireqs(i)
        enddo
! *** send and receive double precision values ***
        isend=0
        do i=0,nprocs-1
          idisp(i)=isend 
! *** tivs is numbr of particles received from each proc -> now sending
          jjlen(i)=tivs(i)*snval
          isend=isend+jjlen(i)
        enddo

        deallocate(tivs)
        allocate(tdvs(0:npj*snval))

! *** setting the sending back data ***
        do i=0,npj-1
! *** ax, ay, az
          tdvs(snval*i)=tdvr(i+npj*4)
          tdvs(snval*i+1)=tdvr(i+npj*5)
          tdvs(snval*i+2)=tdvr(i+npj*6)
          tdvs(snval*i+3)=tdvr(i+npj*7)
        enddo
! *** sending and reveiving the data ***
        irecv=0
! *** number of particles sent is ncomp 
        npjs=npj
! *** ncompi: number of particles originally sent from iproc=0
        npj=ncompi

        deallocate(tdvr)
        allocate(tdvr(0:npj*nval))

        do ip=0,nprocs-1

          allocate(trbuf(0:npjr(ip)*snval))

          call MPI_SCATTERV(tdvs,jjlen,idisp,MPI_DOUBLE_PRECISION &
           ,trbuf,npjr(ip)*snval,MPI_DOUBLE_PRECISION &
           ,ip,MPI_COMM_WORLD,ierr)
! *** set the data to tdvr ***
          do i=0,npjr(ip)-1
            tdvr(irecv+npj*4)=trbuf(snval*i)
            tdvr(irecv+npj*5)=trbuf(snval*i+1)
            tdvr(irecv+npj*6)=trbuf(snval*i+2)
            tdvr(irecv+npj*7)=trbuf(snval*i+3)
            irecv=irecv+1
          enddo

          deallocate(trbuf)

        enddo

        deallocate(tdvs)    

      endif
! *** update the values ***
      do i=0,npj-1
        pn=pn_nfp(i)
        epotp(pn)=epotp(pn)+tdvr(i+npj*4)
        dvxp(pn)=dvxp(pn)+tdvr(i+npj*5)
        dvyp(pn)=dvyp(pn)+tdvr(i+npj*6)
        dvzp(pn)=dvzp(pn)+tdvr(i+npj*7)
      enddo

      deallocate(tdvr)

! *** end of neighbour search ***
      if(nprocs.le.1.or.(iproc.ne.0.and.idivcom.gt.ndivcom)) then
        flagtr=flagtr+1
        if(flagtr.eq.1) then
          goto 71
        else
          goto 94
        endif
      endif

      if(nprocs.gt.1.and.iproc.eq.0) then
        iproc=1
! *** check if need communication ***
        ncomptmax=0
        call MPI_ALLREDUCE(ncompt,ncomptmax,1,MPI_INTEGER &
        ,MPI_MAX,MPI_COMM_WORLD,ierr)
        if(ncomptmax.gt.0) then

          allocate(tivr(0:nprocs))

! *** get how many particles each proc receives
          do i=0,nprocs-1
            tivr(i)=0
          enddo
          call MPI_ALLREDUCE(npproc,tivr,nprocs,MPI_INTEGER &
           ,MPI_SUM,MPI_COMM_WORLD,ierr)
          npprocmax=0
          do i=0,nprocs-1
            if(tivr(i).gt.npprocmax) then
              npprocmax=tivr(i)
            endif
          enddo

          deallocate(tivr)

          if(npprocmax.gt.MNB) then
! *** number of times receiving and calculating for the other nodes
            ndivcom=int((npprocmax+nprocs*nprocs)/MNB)+1
          else
            ndivcom=1
          endif
          ndivcompt=0
          if(ncomptmax.gt.MNB) then
            ndivcompt=int((ncomptmax+nprocs*nprocs)/MNB)+1
            if(ndivcompt.gt.ndivcom) then
              ndivcom=ndivcompt
            endif
          endif
          idivcom=1
        else
! *** no communication required ***
          flagtr=flagtr+1
          if(flagtr.eq.1) then
            goto 71
          else
            goto 94
          endif
        endif
      endif
! *** do communication ndivcom times
      if(idivcom.le.ndivcom) then
        if(ndivcom.eq.1) then
! *** can use pncompt
! *** preparation for sending the data to the other procs ***
          isend=0          
! *** store particle list in the order of sending procs in list() ***
          do ip=0,nprocs-1
            idisp(ip)=isend           
            jjlen(ip)=0
            do i=0,ncompt-1
              if(proccompt(i).eq.ip) then
                pn_nfp(isend)=pncompt(i)
                jjlen(ip)=jjlen(ip)+1
                isend=isend+1
              endif
            enddo
! keep original jjlen 
            ireqs(ip)=jjlen(ip)
          enddo
        else
! *** set the range of number of particles for sending for each proc
          do i=0,nprocs-1
            call para_range(0,npproc(i)-1,ndivcom,idivcom-1 &
             ,isproc(i),ieproc(i))
            cproc(i)=0
          enddo
! *** preparation for sending the data to the other procs ***
          isend=0          
! *** store particle list in the order of sending procs in list() ***
          do ip=0,nprocs-1
            idisp(ip)=isend           
            jjlen(ip)=0
            if(npproc(ip).gt.0.and.ip.ne.myrank) then
! *** get the coordinate and etc. for pseudo node for the proc

              if(flagtr.eq.0) then

                if(proc_tr(0).eq.myrank) then
                  pnodess=next_tr(0)
                else
                  pnodess=0
                endif

                ndp=num_tr-pnodess-1
                allocate(idpnd(0:num_tr))
                allocate(nppnd(0:ndp))
                allocate(ndpnd(0:ndp))
                allocate(xpnd(0:ndp))
                allocate(ypnd(0:ndp))
                allocate(zpnd(0:ndp))
                allocate(lpnd(0:ndp))
                allocate(mpnd(0:ndp))
                allocate(cxpnd(0:ndp))
                allocate(cypnd(0:ndp))
                allocate(czpnd(0:ndp))
                allocate(deltapnd(0:ndp))
                allocate(nextpnd(0:ndp))
                allocate(daupnd(0:ndp))

                ndp=0
                do nd=pnodess,num_tr-1
                  if(proc_tr(nd).eq.ip) then
                    ndpnd(ndp)=nd
                    idpnd(nd)=ndp
                    nppnd(ndp)=np_tr(nd)
                    xpnd(ndp)=cmx_tr(nd)
                    ypnd(ndp)=cmy_tr(nd)
                    zpnd(ndp)=cmz_tr(nd)
                    lpnd(ndp)=l_tr(nd)
                    mpnd(ndp)=mass_tr(nd)
                    cxpnd(ndp)=cx_tr(nd)
                    cypnd(ndp)=cy_tr(nd)
                    czpnd(ndp)=cz_tr(nd)
                    deltapnd(ndp)=delta_tr(nd)
                    ndp=ndp+1
                  endif
                enddo
                do i=0,ndp-1
                  nd=ndpnd(i)                  
                  if(next_tr(nd).gt.pnodess &
                   .and.next_tr(nd).le.ndpnd(ndp-1)) then
                    nextpnd(i)=idpnd(next_tr(nd))      
                  else
                    nextpnd(i)=ndp
                  endif
                  if(daughter_tr(nd).ne.-1) then
                    daupnd(i)=idpnd(daughter_tr(nd))
                  else 
                    daupnd(i)=-1
                  endif
                enddo
              else if(flagtr.eq.1) then
                if(proc_dmtr(0).eq.myrank) then
                  pnodess=next_dmtr(0)
                else
                  pnodess=0
                endif

                ndp=num_dmtr-pnodess-1
                allocate(idpnd(0:num_dmtr))
                allocate(nppnd(0:ndp))
                allocate(ndpnd(0:ndp))
                allocate(xpnd(0:ndp))
                allocate(ypnd(0:ndp))
                allocate(zpnd(0:ndp))
                allocate(lpnd(0:ndp))
                allocate(mpnd(0:ndp))
                allocate(cxpnd(0:ndp))
                allocate(cypnd(0:ndp))
                allocate(czpnd(0:ndp))
                allocate(deltapnd(0:ndp))
                allocate(nextpnd(0:ndp))
                allocate(daupnd(0:ndp))

                ndp=0
                do nd=pnodess,num_dmtr-1
                  if(proc_dmtr(nd).eq.ip) then
                    ndpnd(ndp)=nd
                    idpnd(nd)=ndp
                    nppnd(ndp)=np_dmtr(nd)
                    xpnd(ndp)=cmx_dmtr(nd)
                    ypnd(ndp)=cmy_dmtr(nd)
                    zpnd(ndp)=cmz_dmtr(nd)
                    lpnd(ndp)=l_dmtr(nd)
                    mpnd(ndp)=mass_dmtr(nd)
                    cxpnd(ndp)=cx_dmtr(nd)
                    cypnd(ndp)=cy_dmtr(nd)
                    czpnd(ndp)=cz_dmtr(nd)
                    deltapnd(ndp)=delta_dmtr(nd)
                    ndp=ndp+1
                  endif
                enddo
                do i=0,ndp-1
                  nd=ndpnd(i)                  
                  if(next_dmtr(nd).gt.pnodess &
                   .and.next_dmtr(nd).le.ndpnd(ndp-1)) then
                    nextpnd(i)=idpnd(next_dmtr(nd))      
                  else
                    nextpnd(i)=ndp
                  endif
                  if(daughter_dmtr(nd).ne.-1) then
                    daupnd(i)=idpnd(daughter_dmtr(nd))
                  else 
                    daupnd(i)=-1
                  endif
                enddo
              endif
! *** search particles need communication with ip
              do i=0,ncomp-1
                pn=pncomp(i)
                nd=0
   72           flagp=0
! *** check for particles within a node ***
                if(nppnd(nd).gt.1) then
                  xij=xp(pn)-cxpnd(nd)
                  yij=yp(pn)-cypnd(nd)
                  zij=zp(pn)-czpnd(nd)
                  if(dabs(xij).lt.0.6d0*lpnd(nd)) then
                  if(dabs(yij).lt.0.6d0*lpnd(nd)) then
                  if(dabs(zij).lt.0.6d0*lpnd(nd)) then
! *** go to daughter node ***
                    flagp=1
                  endif
                  endif
                  endif
                  if(flagp.eq.0) then
                    xij=xp(pn)-xpnd(nd)
                    yij=yp(pn)-ypnd(nd)
                    zij=zp(pn)-zpnd(nd)
                    r2=xij*xij+yij*yij+zij*zij
                    if(r2.gt.0.0d0) then
                      l_r2=(lpnd(nd)+deltapnd(nd))*(lpnd(nd)+deltapnd(nd))/r2
                    else
                      l_r2 = 0.0d0
                    endif
                    if(l_r2.gt.theta2) then
                      flagp=1
                    endif
                  endif
                endif
                if(flagp.eq.1) then
                  if(daupnd(nd).eq.-1) then
                    if(cproc(ip).ge.isproc(ip).and.cproc(ip).le.ieproc(ip)) then
                      pn_nfp(isend)=pn
                      jjlen(ip)=jjlen(ip)+1
                      isend=isend+1
                    endif
                    cproc(ip)=cproc(ip)+1
                    goto 95
                  endif
                  nd=daupnd(nd)
                else
                  nd=nextpnd(nd)
                endif
                if(nd.ge.ndp) then
                  goto 95
                endif
                goto 72
   95         enddo

              if(cproc(ip).ne.npproc(ip)) then
                write(6,*) ' Error in treef():npproc,cproc' &
                 ,',myrank,idiv,ndiv,ip,flagtr=' &
                ,npproc(ip),cproc(ip),myrank,idivcom,ndivcom,ip,flagtr
                stop
              endif

              deallocate(idpnd)
              deallocate(nppnd)
              deallocate(ndpnd)
              deallocate(xpnd)
              deallocate(ypnd)
              deallocate(zpnd)
              deallocate(lpnd)
              deallocate(mpnd)
              deallocate(cxpnd)
              deallocate(cypnd)
              deallocate(czpnd)
              deallocate(deltapnd)
              deallocate(nextpnd)
              deallocate(daupnd)

            endif
! keep original jjlen 
            ireqs(ip)=jjlen(ip)
          enddo
        endif
! *** update idivcom
        idivcom=idivcom+1
! *** np for this communication
        ncompi=isend
! *** getting the total number of particles received at each proc ***
        npj=0
        do ip=0,nprocs-1
          irecv=0
          call MPI_SCATTER(jjlen,1,MPI_INTEGER &
           ,irecv,1,MPI_INTEGER,ip,MPI_COMM_WORLD,ierr)
          npjr(ip)=irecv
! *** update npj ***
          npj=npj+irecv
        enddo
        if(npj.gt.MNB) then
          write(6,*) ' Error in treef(): npj > MNB'
          write(6,*) ' npj,MNB=',npj,MNB
          call MPI_ABORT(MPI_COMM_WORLD,ierr)
          stop
        endif

        snval=4
        allocate(tdvs(0:ncompi*snval))

        do i = 0,ncompi-1
          pn=pn_nfp(i)
! *** pn_nfp is not used in the other proc, so keep this ***
          tdvs(snval*i)=xp(pn)
          tdvs(snval*i+1)=yp(pn)
          tdvs(snval*i+2)=zp(pn)             
          tdvs(snval*i+3)=hp(pn)
        enddo
! *** reset sending parameters ***
        do i=0,nprocs-1
          idisp(i)=idisp(i)*snval
          jjlen(i)=ireqs(i)*snval
        enddo

        allocate(tdvr(0:npj*nval))

! *** sending and reveiving the data ***
        irecv=0
        do ip=0,nprocs-1

          allocate(trbuf(0:npjr(ip)*snval))

          call MPI_SCATTERV(tdvs,jjlen,idisp,MPI_DOUBLE_PRECISION &
           ,trbuf,npjr(ip)*snval,MPI_DOUBLE_PRECISION,ip,MPI_COMM_WORLD,ierr)
! *** set the data to tdvr ***
          do i=0,npjr(ip)-1
              tdvr(irecv)=trbuf(snval*i)
              tdvr(irecv+npj)=trbuf(snval*i+1)
              tdvr(irecv+npj*2)=trbuf(snval*i+2)
              tdvr(irecv+npj*3)=trbuf(snval*i+3)
              irecv=irecv+1
          enddo

          deallocate(trbuf)

        enddo       

        deallocate(tdvs)

! *** check ***
        if(irecv.ne.npj) then
          write(6,*) ' Error in treef(): irecv,npj=',irecv,npj
          write(6,*) ' after sending the data to the other pe'
          stop
        endif
! *** initialization ***
! *** setting 0 for ax,ay,az ***
        do i=npj*4,npj*nval-1
          tdvr(i)=0.0d0
        enddo
        if(flagtr.eq.0) then
          goto 70
        else if(flagtr.eq.1) then
          goto 91
        else 
          if(myrank.eq.0) then
            write(*,*) ' Error in treef(): flagtr=',flagtr
          endif
          call  MPI_ABORT(MPI_COMM_WORLD,ierr)
          stop
        endif            
      endif

   94 deallocate(idisp)
      deallocate(jjlen)
      deallocate(ireqs)
      deallocate(npproc)
      deallocate(isproc)
      deallocate(ieproc)
      deallocate(cproc)
      deallocate(npjr)

      deallocate(flagproc)
      deallocate(flagcom)
      deallocate(pncomp)

!      write(fileo,'(a4,i3.3)') 'fptf',myrank
!      open(60,file=fileo,status='unknown')
!      do i=0,nag-1
!        if(x_p(i)**2+y_p(i)**2+z_p(i)**2.gt.0.0d0) then
!        write(60,'(9(1pE13.5))') x_p(i),y_p(i),z_p(i) &
!         ,ax_p(i),ay_p(i),az_p(i) &
!         ,dsqrt(x_p(i)**2+y_p(i)**2+z_p(i)**2) &
!         ,-(ax_p(i)*x_p(i)+ay_p(i)*y_p(i)+az_p(i)*z_p(i)) &
!         /dsqrt(x_p(i)**2+y_p(i)**2+z_p(i)**2)  &
!         ,h_p(i)
!        else
!        write(60,'(9(1pE13.5))') x_p(i),y_p(i),z_p(i) &
!         ,dvx_p(i),dvy_p(i),dvz_p(i) &
!         ,dsqrt(x_p(i)**2+y_p(i)**2+z_p(i)**2) &
!         ,-(ax_p(i)*x_p(i)+ay_p(i)*y_p(i)+az_p(i)*z_p(i)) &
!         ,h_p(i)
!        endif
!      enddo
!      close(60)

end subroutine

