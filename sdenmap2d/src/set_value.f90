
! ****************************************************
!    set_value.f90 for grid3d
!  14  Aug. 2014    written by D.KAWATA
! ****************************************************

! ***********************************************************
!    Definition of set_value()
!  This program set density,div(v),|rot(v)|
! ***********************************************************

subroutine set_value(np)      
      use gcdp_const
      use gcdp_system
      use particle
      use gcdp_gtree
      use gcdp_kernel

      implicit none
      include 'mpif.h'

      integer,intent(in) :: np
      integer npb
      integer i,j,nit,mnit,mnitb,npj,srank,rrank,ip,nval,isend &
       ,irecv,nival,is,iproc,npj0,snval,rnval,snival,npjs,snvalmax
      integer flagnb,nssl
      integer ntw,maxntw
      integer ierr
! * Number of Notfinished particle, temp *      
      integer nlist,tnlist
! * Particle Number in this node *      
      integer pni,pn,nd
      double precision r,s,rc,rh
! * Smoothing length *      
      double precision hsi,ihsi,wij
! * for calculate div(), rot() *
      double precision dwr,dwh,dphih
      double precision xij,yij,zij,tij
      double precision vxij,vyij,vzij
      double precision crit,ipi
! *** information for particles need communication ***
      integer ncomp,ncompt,ndivcom,idivcom
      integer ncomptmax,npprocmax,ndivcompt
      integer ncompi
      integer pncompt(0:MNB-1),proccompt(0:MNB-1)
! *** for pseudo node 
      integer ndp,pnodess
! for test output
      character fileo*60
! pn_nfp can be bigger than ng+ns in this proc, when receiving the particles
      integer pn_nfp(0:MNB-1)
! for work allocatable
      integer,allocatable :: flagi(:)
      integer,allocatable :: pncomp(:),flagproc(:),flagcom(:)
      integer,allocatable :: npproc(:),isproc(:),ieproc(:),cproc(:)
      double precision,allocatable :: hn0p(:)
! *** for pseudo node 
      integer,allocatable :: nextpnd(:),daupnd(:),ndpnd(:),idpnd(:)
      double precision,allocatable :: xpnd(:),ypnd(:),zpnd(:),lpnd(:)
      double precision,allocatable ::  trbuf(:)
      integer,allocatable :: tibuf(:)
! *** for work ***
      integer,allocatable :: idisp(:),jjlen(:),ireqs(:)
      integer,allocatable :: nd_nfp(:),list(:),node(:),tivr(:) &
        ,tivs(:),npjr(:),talist(:),nalist(:),slgtrlist(:)
      double precision,allocatable :: tx(:),ty(:),tz(:),du0(:),tdvr(:) &
        ,tdvs(:),fhx(:),fx(:),hu(:),lu(:),av(:),adv(:)
      double precision,allocatable :: omghp(:),rhop(:)
      integer,allocatable :: nnbp(:)

      allocate(nnbp(0:np-1))
      allocate(rhop(0:np-1))
      allocate(omghp(0:np-1))

! *** allocate work space for MPI ***
      allocate(idisp(0:nprocs))
      allocate(jjlen(0:nprocs))
      allocate(ireqs(0:nprocs))
      allocate(npproc(0:nprocs))
      allocate(isproc(0:nprocs))
      allocate(ieproc(0:nprocs))
      allocate(cproc(0:nprocs))
      allocate(npjr(0:nprocs))

! **** Basic Value ***
      ipi = 1.0d0/M_PI
      crit=dsqrt(3.0d0)*0.5d0
      nval = 14
      nival = 1

! *** initialization ***
      snval=1
      snival=1
! * Calculate Density & div(v) & rot(v) *	  
! * Initialization *
! *** maximum iteration to switch to bi-section ***
      mnitb = 10
! *** maximum iteration ***
      mnit = 100
      nit = 0
! *** for communication ***
      npj=np
! *** initialisation ***
      idivcom=0
      ndivcom=-1
      do i=0,nprocs-1
        npproc(i)=0
      enddo
! *** number of baryon particle: size of list_ap, therefore pn_nfp ***
      npb=np
! *** allocate work space for each baryon particle ***

      allocate(flagi(0:npb))
      allocate(fhx(0:npb))
      allocate(fx(0:npb))
      allocate(hu(0:npb))
      allocate(lu(0:npb))
      allocate(du0(0:npb))
      allocate(tx(0:npb))
      allocate(ty(0:npb))
      allocate(tz(0:npb))
      allocate(av(0:npb))
      allocate(adv(0:npb))
      allocate(nd_nfp(0:npb))
      allocate(hn0p(0:npb))

      do i=0,np-1
        pn_nfp(i)=i
      enddo
      do i=0,np-1
! *** flag for iteration 1: ready for bi-section ***
        pni=pn_nfp(i)
        flagi(pni)=0
        fhx(pni)=INF
        fx(pni)=-INF
        hu(pni)=-INF
        lu(pni)=-INF
! store original h since 8/7 2008, try but too strict condition
        du0(pni)=hp(pni)
! largest h when nnb_p=0
        hn0p(pni)=0.0d0
      enddo

! *** allocate work space for active particles ***

      allocate(tdvr(0:np*nval))
      allocate(tivr(0:np*nival))
      allocate(flagproc(0:np))
      allocate(flagcom(0:np))
      allocate(pncomp(0:np))
      allocate(talist(0:np))

      do i=0,np-1
        pni=pn_nfp(i)
        tdvr(i)=xp(pni)
        tdvr(i+np)=yp(pni)
        tdvr(i+np*2)=zp(pni)
        tdvr(i+np*3)=vxp(pni)
        tdvr(i+np*4)=vyp(pni)
        tdvr(i+np*5)=vzp(pni)
        tdvr(i+np*6)=hp(pni)
        talist(i)=pni
        flagproc(i)=myrank
        flagcom(i)=-1
        tivr(i)=0
      enddo
! rho_p,omgh,div
      do i=np*7,np*nval-1
        tdvr(i)=0.0d0
      enddo 
! *** tree walk for all the proc ***
      npj = np
! *** initialization ***
   71 do i=0,npj-1
        pni=pn_nfp(i)
        omghp(pni)=0.0d0
        nnbp(pni)=0
      enddo

      do i=0,npj-1
        pni=pn_nfp(i)
        rhop(pni)=0.0d0
      enddo
! *** number of particles need communication ***
      ncomp=0
      ncompt=0
! *** iproc=0: in proc, iproc=1: tree walk for particles in the other node ***
      iproc=0 
      maxntw=num_gtr*10

   70 if(np_gtr(0).eq.0) then
        goto 99
      endif
      if(iproc.ne.0.and.proc_gtr(0).ne.myrank) then
        goto 99
      endif

      allocate(list(0:npj))
      allocate(node(0:npj))

      do i=0,npj-1
        list(i)=i
! *** node has to start from 0
        node(i)=0
      enddo 
      nlist=npj


! **********    start tree walk   *********
      ntw=0
   77 if(nlist.eq.0) then

        deallocate(list)
        deallocate(node)

        goto 99
      endif
      ntw=ntw+1

      if(ntw.gt.maxntw) then
        write(6,*) ' Error set_value(): ntw>',maxntw
        write(6,*) ' myrank,iproc,npj,nlist,numgtr,ntw=' &
         ,myrank,iproc,npj,nlist,num_gtr,ntw
!        if(myrank.eq.0) then
        write(fileo,'(a5,i3.3,i8.8)') 'errsv',myrank,ntw
        open(60,file=fileo,status='unknown')
        do i=0,nlist-1
          pni=list(i)
          write(60,'(4(1pE13.5),3I10)') tdvr(pni),tdvr(pni+npj) &
           ,tdvr(pni+npj*2),tdvr(pni+npj*6),pni,node(pni)
        enddo                  
        close(60)
        write(fileo,'(a5,i3.3,a1,i3.3)') 'gtree',myrank
        open(60,file=fileo,status='unknown')
        do i=0,num_gtr-1
          write(60,'(6I10,5(1pE13.5))') i,np_gtr(i),pn_gtr(i) &
           ,next_gtr(i),daughter_gtr(i),proc_gtr(i),l_gtr(i) &
           ,cx_gtr(i),cy_gtr(i),cz_gtr(i),hm_gtr(i)
        enddo
        close(60)
!        endif
!        if(ntw.gt.maxntw+10) then
          stop
!        endif
      endif

      do i=0,nlist-1
        pni=list(i)
        nd=node(pni)         
        xij=tdvr(pni)-cx_gtr(nd)
        yij=tdvr(pni+npj)-cy_gtr(nd)
        zij=tdvr(pni+npj*2)-cz_gtr(nd)
        rc=xij*xij+yij*yij+zij*zij
! *** only take into account gather condition ***
        rh=tdvr(pni+npj*6)+crit*l_gtr(nd)
        rh=rh*rh
        if(np_gtr(nd).eq.1) then
        if(rc.lt.rh) then
        if(proc_gtr(nd).eq.myrank) then
          pn = pn_gtr(nd)
          r=rc
          if(r.ne.0.0d0) then
            r=dsqrt(r)
          endif            
! *** calculate density for both gas and star,
!     in case if some neighbour are changed to stars ****
          hsi=tdvr(pni+npj*6)
          ihsi = 1.0d0/hsi
          s = r*ihsi
          if(s.lt.1.0d0) then
! * satisfy the criterion for neigbour particles *
! * calculate density if flag!=1 *
            is=int(s*dnktab)
            if(is.lt.0) then
              is=0
            else if(is.ge.NKTAB) then
              is=NKTAB-1
            endif
            wij=w_tb(is)+(w_tb(is+1)-w_tb(is))*(s-s_tb(is))*dnktab
            wij=wij*(ihsi**3)
! *** update density ***
            tdvr(pni+npj*7)=tdvr(pni+npj*7)+massp(pn)*wij
! *** calculation of omega ***
            dwr=dwds_s_tb(is)+(dwds_s_tb(is+1)-dwds_s_tb(is)) &
             *(s-s_tb(is))*dnktab
! *** dW/ds/r/h
            dwr=dwr*(ihsi**5)
! *** calculating omgh_p ***
! *** dW/dh ***
            dwh=-s*r*dwr-3.0d0*wij*ihsi
! *** for omgh_p ***
            tdvr(pni+npj*8)=tdvr(pni+npj*8)+massp(pn)*dwh
            if(iproc.ne.0.or.talist(pni).ne.pn) then
! *** count neighbours, not including particle itself ***
              tivr(pni)=tivr(pni)+1
! *** velocity difference ***
              vxij=tdvr(pni+npj*3)-vxp(pn)
              vyij=tdvr(pni+npj*4)-vyp(pn)
              vzij=tdvr(pni+npj*5)-vzp(pn)
! *** calculation for del_v and rot_v ***
            endif
          endif
        endif
        endif
        endif
! * update node *		
        if((np_gtr(nd).eq.1.and.proc_gtr(nd).eq.myrank).or.rc.ge.rh) then
          node(pni)=next_gtr(nd)
        else
! * check if the pseudo node or not *
          if(iproc.eq.0.and.daughter_gtr(nd).eq.-1) then
            if(proc_gtr(nd).ne.myrank.and.flagproc(pni).ne.proc_gtr(nd)) then
! proc_gtr will be the same for different domain within the same proc
              flagproc(pni)=proc_gtr(nd)
! *** store pn ***
              if(flagcom(pni).lt.0) then
                flagcom(pni)=proc_gtr(nd)
                pncomp(ncomp)=pn_nfp(pni)
                ncomp=ncomp+1
              endif
              if(ncompt.lt.MNB) then
                pncompt(ncompt)=pn_nfp(pni)
                proccompt(ncompt)=proc_gtr(nd)
              endif
              ncompt=ncompt+1
              npproc(proc_gtr(nd))=npproc(proc_gtr(nd))+1
            endif
            node(pni)=next_gtr(nd)
          else
! *** information is in the other node ***
            node(pni)=daughter_gtr(nd)
          endif
        endif          
      enddo
! * update not-finished particle list *
      tnlist = nlist
      nlist = 0
      if(iproc.eq.0.and.next_gtr(0).ne.0) then
        do i=0,tnlist-1
          if(node(list(i)).gt.0) then
            list(nlist)=list(i)
            nlist=nlist+1
          else if(node(list(i)).eq.0) then
            node(list(i))=next_gtr(0)
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
! *** end itteration within the proc ***

! *** the results whitn the proc ***
! *** update variables ***
   99 if(iproc.ne.0) then
! *** setting numbers for sending and receiving
! *** snval -> rnval: received nval
        rnval=snval
! *** new nval to send the results
        snval=7
        if(allocated(tivs)) then
          deallocate(tivs)
        endif
        allocate(tivs(0:nprocs))
        do i=0,nprocs-1
! *** keep npjr: data received from each proc ***
          tivs(i)=npjr(i)
! *** ireqs(i)=jjlen(i)number of particles sent -> now receiving.
          npjr(i)=ireqs(i)
        enddo

! *** for double precision ***
        isend=0
        do i=0,nprocs-1
          idisp(i)=isend 
! *** tivs is numbr of particles received from each proc
          jjlen(i)=tivs(i)*snval
          isend=isend+jjlen(i)          
        enddo

! *** setting the sending back data ***
        if(allocated(tdvs)) then
          deallocate(tdvs)
        endif
        allocate(tdvs(0:npj*snval))
        do i=0,npj-1
          tdvs(snval*i)=tdvr(i+npj*7)
          tdvs(snval*i+1)=tdvr(i+npj*8)
          tdvs(snval*i+2)=tdvr(i+npj*9)
          tdvs(snval*i+3)=tdvr(i+npj*10)
          tdvs(snval*i+4)=tdvr(i+npj*11)
          tdvs(snval*i+5)=tdvr(i+npj*12)
          tdvs(snval*i+6)=tdvr(i+npj*13)
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

          if(allocated(trbuf)) then
            deallocate(trbuf)
          endif
          allocate(trbuf(0:npjr(ip)*snval))

          call MPI_SCATTERV(tdvs,jjlen,idisp,MPI_DOUBLE_PRECISION &
          ,trbuf,npjr(ip)*snval,MPI_DOUBLE_PRECISION,ip,MPI_COMM_WORLD,ierr)
! *** set the data to tdvr ***
          do i=0,npjr(ip)-1
            tdvr(irecv+npj*7)=trbuf(snval*i)
            tdvr(irecv+npj*8)=trbuf(snval*i+1)
            tdvr(irecv+npj*9)=trbuf(snval*i+2)
            tdvr(irecv+npj*10)=trbuf(snval*i+3)
            tdvr(irecv+npj*11)=trbuf(snval*i+4)
            tdvr(irecv+npj*12)=trbuf(snval*i+5)
            tdvr(irecv+npj*13)=trbuf(snval*i+6)
            irecv=irecv+1
          enddo
        enddo    

        deallocate(tdvs)
        deallocate(trbuf)

        if(irecv.ne.npj) then
          write(6,*) ' Error in set_value(): irecv,npj=',irecv,npj
          write(6,*) ' when sending back the results.'
          stop
        endif
! *** for integer values ***
        isend=0
        snival=1
        do i=0,nprocs-1
          idisp(i)=isend 
! *** tivs is numbr of particles received from each proc
          jjlen(i)=tivs(i)*snival
          isend=isend+jjlen(i)          
        enddo
! *** setting the sending back data, use npjs ***

        if(allocated(tivs)) then
          deallocate(tivs)
        endif
        allocate(tivs(0:npjs))

        do i=0,npjs-1
          tivs(snival*i)=tivr(i)
        enddo
! *** sending and reveiving the data ***
        irecv=0
! *** ncompi: number of particles originally sent from iproc=0
        npj=ncompi

        deallocate(tivr)
        allocate(tivr(0:npj*snival))

        do ip=0,nprocs-1

          allocate(tibuf(0:npjr(ip)*snival))

          call MPI_SCATTERV(tivs,jjlen,idisp,MPI_INTEGER,tibuf &
           ,npjr(ip)*snival,MPI_INTEGER,ip,MPI_COMM_WORLD,ierr)
! *** set the data to tdvr ***
          do i=0,npjr(ip)-1
            tivr(irecv)=tibuf(snival*i)
            irecv=irecv+1
          enddo

          deallocate(tibuf)

        enddo          

        deallocate(tivs)

      endif

      do i=0,npj-1
        pni = pn_nfp(i)
        rhop(pni)=rhop(pni)+tdvr(i+npj*7)
        nnbp(pni)=nnbp(pni)+tivr(i)
      enddo
      do i=0,npj-1
        pni = pn_nfp(i)
        omghp(pni)=omghp(pni)+tdvr(i+npj*8)
      enddo   

      deallocate(tdvr)
      deallocate(tivr)

! *** This has to be here, in case no communication required, then come to 72
      if(iproc.ne.0.and.idivcom.gt.ndivcom) then
        do i=0,npj0-1
          pn_nfp(i)=nd_nfp(i)
        enddo   
        npj=npj0
      endif

! *** restore pn_nfp ***
! *** set omgh_p and etc. ***
   72 if(nprocs.eq.1.or.(iproc.ne.0.and.idivcom.gt.ndivcom)) then

        do i=0,npj-1
          pni=pn_nfp(i)
          if(rhop(pni).gt.0.0d0) then
            omghp(pni)=1.0d0+THIRD*(hp(pni)/rhop(pni))*omghp(pni)
          endif
        enddo
      endif

      if(nprocs.gt.1.and.iproc.eq.0) then
        iproc=1
! *** check if need communication ***
        ncomptmax=0
        call MPI_ALLREDUCE(ncompt,ncomptmax,1,MPI_INTEGER &
         ,MPI_MAX,MPI_COMM_WORLD,ierr)

        if(ncomptmax.gt.0) then
! *** get how many particles each proc receives

          allocate(tivr(0:nprocs))

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
! *** store original pn_nfp in nd_nfp ***
          npj0=npj
          do i=0,npj0-1
            nd_nfp(i)=pn_nfp(i)
          enddo
        else
! *** no communication required ***
          goto 72
        endif
      endif

! *** do communication ndivcom times
      if(idivcom.le.ndivcom) then

        if(ndivcom.eq.1) then
!        if(2.eq.1) then
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
          if(isend.ne.ncompt) then    
            write(6,*) ' Error in set_value(): isend.ne.ncomp'
            write(6,*) ' when counting N particles need communication'
            write(6,*) ' myrank,isend,ncomp=',myrank,isend,ncomp
            stop
          endif
        else
! *** set the range of number of particles for sending for each proc
          do i=0,nprocs-1
            call para_range(0,npproc(i)-1,ndivcom,idivcom-1 &
             ,isproc(i),ieproc(i))
! initialise for counting particles for each proc
            cproc(i)=0
          enddo
! *** preparation for sending the data to the other procs ***
          isend=0          
! gtree id for starting pseudo node
          if(proc_gtr(0).eq.myrank) then
            pnodess=next_gtr(0)
          else
! in case, if there is no local tree, but pseudo nodes
            pnodess=0
          endif
          ndp=num_gtr-pnodess

          allocate(idpnd(0:num_gtr))
          allocate(ndpnd(0:ndp))
          allocate(xpnd(0:ndp))
          allocate(ypnd(0:ndp))
          allocate(zpnd(0:ndp))
          allocate(lpnd(0:ndp))
          allocate(nextpnd(0:ndp))
          allocate(daupnd(0:ndp))

! *** store particle list in the order of sending procs in list() ***
! *** get the coordinate and etc. for pseudo node for the proc
          do ip=0,nprocs-1
            idisp(ip)=isend           
            jjlen(ip)=0
            if(npproc(ip).gt.0.and.ip.ne.myrank) then
              ndp=0
              do nd=pnodess,num_gtr-1
                if(proc_gtr(nd).eq.ip) then
                  ndpnd(ndp)=nd
                  idpnd(nd)=ndp
                  xpnd(ndp)=cx_gtr(nd)
                  ypnd(ndp)=cy_gtr(nd)
                  zpnd(ndp)=cz_gtr(nd)
                  lpnd(ndp)=l_gtr(nd)
                  ndp=ndp+1
                endif
              enddo
              do i=0,ndp-1
                nd=ndpnd(i)                  
                if(next_gtr(nd).gt.pnodess &
                 .and.next_gtr(nd).le.ndpnd(ndp-1)) then
                  nextpnd(i)=idpnd(next_gtr(nd))      
                else
                  nextpnd(i)=ndp
                endif
                if(daughter_gtr(nd).ne.-1) then
                  daupnd(i)=idpnd(daughter_gtr(nd))
                else 
                  daupnd(i)=-1
                endif
              enddo
! *** search particles need communication with ip
              do i=0,ncomp-1
                pn=pncomp(i)
                nd=0
   73           xij=xp(pn)-xpnd(nd)
                yij=yp(pn)-ypnd(nd)
                zij=zp(pn)-zpnd(nd)
                rc=xij*xij+yij*yij+zij*zij
! *** only take into account gather condition ***
                rh=hp(pn)+crit*lpnd(nd)
                rh=rh*rh
                if(rc.lt.rh) then
                  if(daupnd(nd).eq.-1) then
                    if(cproc(ip).ge.isproc(ip).and.cproc(ip).le.ieproc(ip)) then
                      pn_nfp(isend)=pn
                      jjlen(ip)=jjlen(ip)+1
                      isend=isend+1
                    endif
                    cproc(ip)=cproc(ip)+1
                    goto 90
                  endif
                  nd=daupnd(nd)
                else
                  nd=nextpnd(nd)
                endif
                if(nd.ge.ndp) then
                  goto 90
                endif
                goto 73
   90         enddo
              if(cproc(ip).ne.npproc(ip)) then
                write(6,*) ' Error in set_value():npproc,cproc,' &
                 ,'myrank,ip,ndp,nit,idiv,ndivcom,ncomp=' &
                 ,npproc(ip),cproc(ip),myrank,ip,ndp,nit,idivcom &
                 ,ndivcom,ncomp
                 write(fileo,'(a8,i3.3)') 'setverrt',myrank
                 open(60,file=fileo,status='unknown')
                 do i=0,ndp-1
                   nd=ndpnd(i)
                   write(60,'(9I10,5(1pE13.5))') i,ndpnd(i),idpnd(nd) &
                    ,nextpnd(i),daupnd(i),np_gtr(nd),next_gtr(nd) &
                    ,daughter_gtr(nd),proc_gtr(nd),xpnd(i),ypnd(i),zpnd(i) &
                    ,lpnd(i),cy_gtr(nd)
                 enddo
                 close(60)
                 write(fileo,'(a8,i3.3)') 'setverrp',myrank
                 open(60,file=fileo,status='unknown')
                 do i=0,ncomp-1
                   pn=pncomp(i)
                   write(60,'(4(1pE13.5))') xp(pn),yp(pn),zp(pn),hp(pn)
                enddo
                close(60)
                call MPI_BARRIER(MPI_COMM_WORLD,ierr)
                call MPI_FINALIZE(ierr)
                stop
              endif
            endif

! keep original jjlen 
            ireqs(ip)=jjlen(ip)
          enddo
          deallocate(idpnd)
          deallocate(ndpnd)
          deallocate(xpnd)
          deallocate(ypnd)
          deallocate(zpnd)
          deallocate(lpnd)
          deallocate(nextpnd)
          deallocate(daupnd)
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
          npj=npj+irecv
        enddo
        if(npj.gt.MNB) then
          write(6,*) ' Error in set_value(): npj receiving >MNB'
          write(6,*) ' npj,MNB,ndiv,idiv=',npj,MNB,ndivcom,idivcom
          call MPI_ABORT(MPI_COMM_WORLD,ierr)
          stop
        endif

        snval=7
! ncompi = number of particles to be sent

        allocate(tdvs(0:ncompi*snval))

        do i=0,ncompi-1
          pni=pn_nfp(i)
          tdvs(snval*i)=xp(pni)
          tdvs(snval*i+1)=yp(pni)
          tdvs(snval*i+2)=zp(pni)
          tdvs(snval*i+3)=vxp(pni)
          tdvs(snval*i+4)=vyp(pni)
          tdvs(snval*i+5)=vzp(pni)
          tdvs(snval*i+6)=hp(pni)
        enddo
! *** sending proc ***
        do i=0,nprocs-1
          idisp(i)=idisp(i)*snval
          jjlen(i)=ireqs(i)*snval
        enddo

! reallocate tdvr and tivr to receive npj particles
        allocate(tdvr(0:npj*nval))
        allocate(tivr(0:npj*nival))

! *** sending and reveiving the data ***
        irecv=0
        do ip=0,nprocs-1

          allocate(trbuf(0:npjr(ip)*snval))

          call MPI_SCATTERV(tdvs,jjlen,idisp,MPI_DOUBLE_PRECISION &
           ,trbuf,npjr(ip)*snval,MPI_DOUBLE_PRECISION,ip,MPI_COMM_WORLD,ierr)
! *** sethe data to tdvr ***
          do i=0,npjr(ip)-1
            tdvr(irecv)=trbuf(snval*i)
            tdvr(irecv+npj)=trbuf(snval*i+1)
            tdvr(irecv+npj*2)=trbuf(snval*i+2)
            tdvr(irecv+npj*3)=trbuf(snval*i+3)
            tdvr(irecv+npj*4)=trbuf(snval*i+4)
            tdvr(irecv+npj*5)=trbuf(snval*i+5)
            tdvr(irecv+npj*6)=trbuf(snval*i+6)
            irecv=irecv+1
          enddo

          deallocate(trbuf)

        enddo       

        deallocate(tdvs)

        if(irecv.ne.npj) then
          write(6,*) ' Error in set_value(): irecv,npj=',irecv,npj
          write(6,*) ' after sending the data to the other pe'
          stop
        endif
! *** initialization ***
        do i=0,npj*nival-1
          tivr(i)=0
        enddo
        do i=npj*snval,npj*nval-1
          tdvr(i)=0.0d0
        enddo

        goto 70
      endif

      deallocate(flagproc)
      deallocate(flagcom)
      deallocate(talist)

! * Check smoothing length: tdvr(i+np*6)=hp(pni)
! * criterion h=etah (m/rho)^1/3

      allocate(list(0:npj))

      tnlist=npj
      do i=0,tnlist-1
        list(i)=pn_nfp(i)
      enddo
      nlist = 0

! *** nalist strores the particle whose density=0 ***
      nssl=0
      if(nit.lt.mnitb) then

        allocate(nalist(0:tnlist))

        do i = 0,tnlist-1
          pn = list(i)
          if(nnbp(pn).le.0) then
! *** density is zero. can happen for star ***
            nalist(nssl)=pn 
            nssl=nssl+1
            pn_nfp(nlist)=pn
            nlist = nlist+1
            flagi(pn)=0
            if(hp(pn).gt.hn0p(pn)) then
              hn0p(pn)=hp(pn)
            endif
          else
! *** target density ***
            av(pn)=massp(pn)*((ETAH/hp(pn))**3)
! *** difference between target and real density ***
            adv(pn)=av(pn)-rhop(pn)
! *** new h ***
            av(pn)=hp(pn)-adv(pn)/(-3.0d0*rhop(pn)*omghp(pn)/hp(pn))
! tried below, but convergence is slower. 8/7 2008
!              av(pn)=hp(pn)-adv(pn)/(-3.0d0*av(pn)*omghp(pn)/hp(pn))
            if(dabs(hp(pn)-av(pn))/hp(pn).gt.ERRH) then
! tried below, but too strict condtion 8/7 2008
!              if(dabs(hp(pn)-av(pn))/du0(pn).gt.ERRH) then
! *** store the data for bi-section
              if(adv(pn).gt.0.0d0.and.hp(pn).gt.hu(pn)) then
! largest hu
                hu(pn)=hp(pn)
                fhx(pn)=adv(pn)
              else if(adv(pn).lt.0.0d0.and.hp(pn).lt.lu(pn)) then
! smallest lu
                lu(pn)=hp(pn)
                fx(pn)=adv(pn)
                endif
! *** update ***
              if(av(pn).gt.DHFLIM*hp(pn)) then
                hp(pn)=DHFLIM*hp(pn)
              else if(av(pn).gt.0.0d0) then
                 hp(pn)=av(pn)
              else
                hp(pn)=0.5d0*hp(pn)
              endif
              pn_nfp(nlist)=pn
              nlist = nlist+1
            endif
          endif
        enddo
        if(nssl.gt.0) then

!            allocate(slgtrlist(0:nssl))
!
!            do i=0,nssl-1
!              slgtrlist(i)=nalist(i)
!            enddo
! *** set h for particles whose density is zero, using local _gtr ***
!            call setslgtr(nssl,slgtrlist,npb,0)
!
!            deallocate(slgtrlist)

          do i=0,nssl-1
            pn=nalist(i)
            hp(pn)=1.2d0*hp(pn)
          enddo
        endif

        deallocate(nalist)

      else
        do i = 0,tnlist-1
          pn = list(i)
! *** target density ***
          av(pn)=massp(pn)*((ETAH/hp(pn))**3)
! *** difference between target and real density ***
          adv(pn)=av(pn)-rhop(pn)
        enddo
        if(nit.eq.mnitb) then
          if(nssl.gt.0) then
            write(6,*) ' Error in set_value(): density is zero.after iteration.'
            write(6,*) ' myrank,nssl=',myrank,nssl
            call MPI_ABORT(MPI_COMM_WORLD,ierr)
            stop
          endif
          do i=0,tnlist-1
            pn=list(i)
            if(hu(pn).gt.0.0d0.and.lu(pn).gt.0.0d0) then
              flagi(pn)=0
            else if(lu(pn).lt.0.0d0) then
! this includes lu(pn).lt.0.0d0
              flagi(pn)=-1
            else if(hu(pn).lt.0.0d0) then     
!                hp(pn)=lu(pn)
              flagi(pn)=1
              if(lu(pn).lt.0.0d0) then
                write(6,*) ' Error in set_value(): rank,id=',myrank,idp(pn)
                write(6,*) ' no guess for initial h'
                write(6,*) ' h,rho,flagi=',hp(pn),rhop(pn),flagi(pn)
                write(6,*) ' lu,hu,av,adv=',lu(pn),hu(pn),av(pn),adv(pn)
                stop
              endif
            else
              write(6,*) ' Error in set_value(): rank,id=',myrank,idp(pn)
              write(6,*) ' no guess for initial h'
              write(6,*) ' h,rho=',hp(pn),rhop(pn)
              call MPI_ABORT(MPI_COMM_WORLD,ierr)
              stop
            endif
          enddo
        endif
! *** bi-section ***
        do i=0,tnlist-1
          pn=list(i)
          if(flagi(pn).eq.0) then
            if(adv(pn)*fhx(pn).lt.0.0d0) then
              lu(pn)=hp(pn)
            else
              hu(pn)=hp(pn)
              fhx(pn)=adv(pn)
            endif
            av(pn)=0.5d0*(lu(pn)+hu(pn))
          else if(flagi(pn).eq.1) then
! *** need to find f(hu)>0 ***
            if(adv(pn).gt.0.0d0.and.nnbp(pn).gt.0) then
              fhx(pn)=adv(pn)
              hu(pn)=hp(pn)
              av(pn)=0.5d0*(lu(pn)+hu(pn))
              flagi(pn)=0
            else if(nnbp(pn).le.0) then
              if(hp(pn).gt.hn0p(pn)) then
                hn0p(pn)=hp(pn)
              endif 
              av(pn)=0.5d0*(lu(pn)+hn0p(pn))
            else
! *** believe that decreasing h would increase f(h) ***
              if(hp(pn).lt.lu(pn)) then
                lu(pn)=hp(pn)
                fx(pn)=adv(pn)
                av(pn)=0.5d0*(hp(pn)+hn0p(pn))
              else
                av(pn)=0.8d0*lu(pn)
                if(av(pn).lt.hn0p(pn)) then
                 av(pn)=0.5d0*(hp(pn)+hn0p(pn))
                endif
              endif
            endif
          else if(flagi(pn).eq.-1) then
! **** Note: hu < lu ***
            if(adv(pn).lt.0.0d0.and.nnbp(pn).gt.0) then
              lu(pn)=hp(pn)
              if(hu(pn).gt.0.0d0) then
                av(pn)=0.5d0*(lu(pn)+hu(pn))
                flagi(pn)=0
              else
                fx(pn)=adv(pn)
! *** now turned to flagi(pn) and decreasing h to find hu
                flagi(pn)=1
              endif
            else
              av(pn)=1.2d0*hp(pn)
            endif
          endif
          if(dabs(hp(pn)-av(pn))/hp(pn).gt.ERRH) then
! tried below, but too strict condtion 8/7 2008
!            if(dabs(hp(pn)-av(pn))/du0(pn).gt.ERRH
! *** update ***
            hp(pn)=av(pn)
            pn_nfp(nlist)=pn
            nlist = nlist+1
         else if(flagi(pn).eq.1.and. &
          (nnbp(pn).le.0.or.(lu(pn)-hn0p(pn))/lu(pn).gt.ERRH)) then
! *** update ***
            hp(pn)=av(pn)
            pn_nfp(nlist)=pn
            nlist = nlist+1
         else if(flagi(pn).ne.0.and.flagi(pn).ne.1) then
! *** update ***
            hp(pn)=av(pn)
            pn_nfp(nlist)=pn
            nlist = nlist+1
          endif
        enddo
      endif

      deallocate(list)

      flagnb=0
      call MPI_ALLREDUCE(nlist,flagnb,1,MPI_INTEGER &
       ,MPI_SUM,MPI_COMM_WORLD,ierr)

      if(flagnb.gt.0) then
        if(nit.gt.mnit) then
          if(myrank.eq.0) then
            write(6,*) myrank &
             ,' Warning: iteration does not work in set_value' &
             ,nlist,flagnb
          endif
          if(nlist.gt.0) then
            write(fileo,'(a5,i3.3)') 'error',myrank
            open(60,file=fileo,status='unknown')
            do i=0,nlist-1
              pn=pn_nfp(i)
              write(60,'(15(1pE13.5),5I10)') xp(pn),yp(pn),zp(pn),hp(pn) &
               ,rhop(pn),ETAH*((massp(pn)/rhop(pn))**THIRD) &                
               ,av(pn),adv(pn),hu(pn) &
               ,lu(pn),fhx(pn),fx(pn),du0(pn),hn0p(pn),omghp(pn) &
               ,nnbp(pn),pn,idp(pn),flagi(pn),flagi(pn)
            enddo
            write(6,*) ' myrank,np,nlist=',myrank,np,nlist
          endif

!           write(fileo,'(a4,i3.3)') 'allp',myrank
!            open(60,file=fileo,status='unknown')
!            do i=0,ng-1
!              pn=list_ap(i)
!              write(60,'(5(1pE13.5),3I10)') xp(pn),yp(pn),zp(pn),hp(pn) &
!               ,rhop(pn),nnbp(pn),idp(pn),pn
!            enddo
!            close(60)

          call MPI_BARRIER(MPI_COMM_WORLD,ierr)
          call MPI_ABORT(MPI_COMM_WORLD,ierr)
          stop
        endif
! *** since gas h not changed, neighbour search is needed only for 
! *** stars or gas which changed h.
        npj=nlist
! reallocate tdvr and tivr to receive npj particles
        allocate(tdvr(0:npj*nval))
        allocate(tivr(0:npj*nival))
        allocate(flagproc(0:npj))
        allocate(flagcom(0:npj))
        allocate(talist(0:npj))
        do i=0,npj-1
          pni = pn_nfp(i)

! for test
!            if(idp(pni).eq.7935) then
!              write(6,'(a2,2I6,9(1pE13.5),2I6)') 'h=' &
!              ,nit,idp(pni),hp(pni) &
!              ,hnmaxinp(pni),hnmaxoutp(pni) &
!              ,hu(pni),lu(pni),adv(pni),av(pni),fx(pni),fhx(pni) &
!              ,nnbp(pni),flagi(pni)
!            endif


          tdvr(i)=xp(pni)
          tdvr(i+npj)=yp(pni)
          tdvr(i+npj*2)=zp(pni)
          tdvr(i+npj*3)=vxp(pni)
          tdvr(i+npj*4)=vyp(pni)
          tdvr(i+npj*5)=vzp(pni)
          tdvr(i+npj*6)=hp(pni)
          talist(i)=pni
          flagproc(i)=myrank
          flagcom(i)=-1
          tivr(i)=0
        enddo
        do i=npj*7,npj*nval-1
          tdvr(i)=0.0d0
        enddo 
        nit = nit+1
! *** initialze proc iteration
        iproc=0
! *** initialisation ***
        idivcom=0
        ndivcom=-1
        do i=0,nprocs-1
          npproc(i)=0
        enddo
        goto 71
      endif
! * define div(v) & |rot(v)| *
  999 do i=0,np-1 
        pn=i
      enddo

!      write(fileo,'(a5,i3.3)') 'setvp',myrank
!      open(60,file=fileo,status='unknown')
!      do i=0,np-1
!        pni=i
!        write(60,'(7(1pE13.5),I10)') xp(pni),yp(pni),zp(pni) &
!         ,rhop(pni),hp(pni),omghp(pni),massp(pni),nnbp(pni)
!      enddo
!      close(60)
!      stop

! deallocate
      deallocate(idisp)
      deallocate(jjlen)
      deallocate(ireqs)
      deallocate(npproc)
      deallocate(isproc)
      deallocate(ieproc)
      deallocate(cproc)
      deallocate(npjr)

      deallocate(flagi)
      deallocate(fhx)
      deallocate(fx)
      deallocate(hu)
      deallocate(lu)
      deallocate(du0)
      deallocate(tx)
      deallocate(ty)
      deallocate(tz)
      deallocate(av)
      deallocate(adv)
      deallocate(nd_nfp)
      deallocate(hn0p)

      if(allocated(list)) then
        deallocate(list)
        deallocate(node)
      endif

      deallocate(nnbp)
      deallocate(rhop)
      deallocate(omghp)

end subroutine
