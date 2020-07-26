! /****************************************************
!   output.f for bin2asc
!  4 Oct. 2017    written by D.KAWATA
! ****************************************************/

subroutine output(ngt,ng,ndmt,ndm,ndm1t,ndm1,nst,ns &
       ,flagr,step,nskip,asc,flagbo,flago,flagselp,rrange,zrange)
      use gcdp_const
      use gcdp_system
      use gcdp_baryon
      use gcdp_dm
 
      implicit none
      include 'mpif.h'
      
      integer,intent(in) ::  ngt,ng,ndmt,ndm,ndm1t,ndm1,nst,ns,flagr,nskip &
        ,step,flagbo,flago,flagselp
      double precision,intent(in) :: asc,rrange(0:1),zrange(0:1)
      double precision asci
      integer nvals,nivals
      integer i,ip,ip2,npj,pn,np,ngo,ndmo,nso,nc,nselp
      double precision nhp,temp,rp,gam
! for test particle info
      double precision rmp,eccp
! for work
      integer ierr
      integer,allocatable :: istatus(:)
      integer,allocatable :: npjr(:),ireqs(:),ireqr(:)
      integer,allocatable :: tivs(:),tivr(:)
      double precision,allocatable :: tdvs(:),tdvr(:)
! /**** List *****/
      character fileo*60

      allocate(istatus(MPI_STATUS_SIZE))

! *** for binary output ***
! allocate array for MPI
      allocate(npjr(0:nprocs-1))
      allocate(ireqs(0:nprocs-1))
      allocate(ireqr(0:nprocs-1))

! 1/a
      asci=1.0d0/asc

      gam=5.0d0/3.0d0

! number of selected particles
      nselp=0
      if(flagselp.ne.0) then
        if(flagselp.eq.1) then
          do i=0,ng-1
            pn=list_ap(i)
            rp=dsqrt(x_p(pn)**2+y_p(pn)**2+z_p(pn)**2)
            if(rp.gt.rrange(0).and.rp.lt.rrange(1)) then
              if(z_p(pn).gt.zrange(0).and.z_p(pn).lt.zrange(1)) then   
                nselp=nselp+1
              endif
            endif
          enddo 
       else if(flagselp.eq.2) then
          do i=0,ndm-1
            rp=dsqrt(x_dm(i)**2+y_dm(i)**2+z_dm(i)**2)
            if(rp.gt.rrange(0).and.rp.lt.rrange(1)) then
              if(z_dm(i).gt.zrange(0).and.z_dm(i).lt.zrange(1)) then
                nselp=nselp+1
              endif
            endif
          enddo 
        else if(flagselp.eq.2) then
          do i=ng,ng+ns-1
            pn=list_ap(i)
            rp=dsqrt(x_p(pn)**2+y_p(pn)**2+z_p(pn)**2)
            if(rp.gt.rrange(0).and.rp.lt.rrange(1)) then
              if(z_p(pn).gt.zrange(0).and.z_p(pn).lt.zrange(1)) then   
                nselp=nselp+1
              endif
            endif
          enddo 
        endif
        allocate(tivs(0:0))
        allocate(tivr(0:0))
        tivs(0)=nselp
        tivr(0)=0
        call MPI_ALLREDUCE(tivs,tivr,1,MPI_INTEGER &
             ,MPI_MAX,MPI_COMM_WORLD,ierr)
        nselp=tivr(0)
        deallocate(tivs)        
        deallocate(tivr)
        if(myrank.eq.0) then
          write(6,*) ' N selected particle=',nselp
        endif
      endif
      
      ngo=0
      ndmo=0
      nso=0
! *** output gas particles ***
      if(ngt.gt.0) then
        if(flagr.eq.0) then
          nvals=9
          nivals=1
        else
          nvals=20
          nivals=3
        endif
        np=0
        do i=0,ng-1,nskip
          np=np+1
        enddo

        allocate(tivs(0:ng*nivals-1))
        allocate(tivr(0:ng*nivals-1))
        allocate(tdvs(0:ng*nvals-1))
        allocate(tdvr(0:ng*nvals-1))

        nc=0
        do i=0,ng-1,nskip
          pn=list_ap(i)
          tdvs(nc)=x_p(pn)*asci
          tdvs(nc+np)=y_p(pn)*asci
          tdvs(nc+np*2)=z_p(pn)*asci
          tdvs(nc+np*3)=vx_p(pn)
          tdvs(nc+np*4)=vy_p(pn)
          tdvs(nc+np*5)=vz_p(pn)
          tdvs(nc+np*6)=m_p(pn)
          tdvs(nc+np*7)=rho_p(pn)
          tdvs(nc+np*8)=u_p(pn)
          tivs(nc)=id_p(pn)
          nc=nc+1
        enddo
        if(flagr.gt.0) then
          nc=0
          do i=0,ng-1,nskip
            pn=list_ap(i)
            tdvs(nc+np*9)=mzHe_p(pn)
            tdvs(nc+np*10)=mzC_p(pn)
            tdvs(nc+np*11)=mzN_p(pn)
            tdvs(nc+np*12)=mzO_p(pn)
            tdvs(nc+np*13)=mzNe_p(pn)
            tdvs(nc+np*14)=mzMg_p(pn)
            tdvs(nc+np*15)=mzSi_p(pn)
            tdvs(nc+np*16)=mzFe_p(pn)
            tdvs(nc+np*17)=mzZ_p(pn)
            tdvs(nc+np*18)=h_p(pn)
            tdvs(nc+np*19)=myu_p(pn)
            tivs(nc+np)=flagfd_p(pn)
            tivs(nc+np*2)=nnb_p(pn)
            nc=nc+1
          enddo
        endif
! *** set number of receiving data 
        if(myrank.eq.0) then
          do ip=0,nprocs-1
            if(ip.eq.myrank) then
              npjr(ip)=np
            else
              npjr(ip)=0
              call MPI_IRECV(npjr(ip),1,MPI_INTEGER,ip,ip &
               ,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
            endif
          enddo
        else
          ip=myrank
          call MPI_ISEND(np,1,MPI_INTEGER,0,myrank &
            ,MPI_COMM_WORLD,ireqs(ip),ierr)
          call MPI_WAIT(ireqs(ip),istatus,ierr)
        endif    
! *** receiving id_p and double data and write
        if(myrank.eq.0) then
          write(fileo,'(a8,i6.6)') 'output/g',step
          open(61,file=fileo,status='unknown')
          if(flagselp.eq.1) then
            write(fileo,'(a11,i6.6,a4)') 'output/selg',step,'.bin'
            open(63,file=fileo,status='unknown',form='unformatted')
            write(63) nselp
          endif             
          do ip=0,nprocs-1
            if(ip.eq.myrank) then
              do i=0,npjr(myrank)*nivals-1
                tivr(i)=tivs(i)
              enddo
              do i=0,npjr(myrank)*nvals-1
                tdvr(i)=tdvs(i)
              enddo
            else
              call MPI_IRECV(tivr,npjr(ip)*nivals,MPI_INTEGER,ip,ip &
                ,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
              ip2=ip*2
              call MPI_IRECV(tdvr,npjr(ip)*nvals,MPI_DOUBLE_PRECISION &
               ,ip,ip2,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
            endif
! *** write the data
            npj=npjr(ip)
            ngo=ngo+npj
            if(flagselp.eq.1) then
              do i=0,npj-1
                rp=dsqrt(tdvr(i)**2+tdvr(i+npj)**2+tdvr(i+npj*2)**2)
                if(rp.gt.rrange(0).and.rp.lt.rrange(1)) then
                  if(tdvr(i+npj*2).gt.zrange(0) &
                    .and.tdvr(i+npj*2).lt.zrange(1)) then
                    write(63) tdvr(i),tdvr(i+npj),tdvr(i+npj*2) &
                      ,tdvr(i+npj*3),tdvr(i+npj*4),tdvr(i+npj*5)
                  endif
                endif
              enddo
            endif  
            if(flagr.eq.0) then
              do i=0,npj-1
                write(61,161) tdvr(i),tdvr(i+npj),tdvr(i+npj*2) &
               ,tdvr(i+npj*3),tdvr(i+npj*4),tdvr(i+npj*5) &
               ,tdvr(i+npj*6),tdvr(i+npj*7),tdvr(i+npj*8) &
               ,dsqrt(tdvr(i)**2+tdvr(i+npj)**2+tdvr(i+npj*2)**2) &
               ,tivr(i)
 161            format(10(1pE13.5),I10)
              enddo
            else
              do i=0,npj-1
                rp=dsqrt(tdvr(i)**2+tdvr(i+npj)**2+tdvr(i+npj*2)**2)
                nhp=((tdvr(i+npj*6)-((tdvr(i+npj*9)+tdvr(i+npj*17)) &
                  /MUSM))/tdvr(i+npj*6))*tdvr(i+npj*7)*(DU/MP)
                p_p(i)=(gam-1.0d0)*tdvr(i+npj*7)*tdvr(i+npj*8)
                temp=(p_p(i)*tdvr(i+npj*19)/(tdvr(i+npj*7)*TPRHO*MYU))*1.0e4
                cs_p(i)=dsqrt(gam*p_p(i)/tdvr(i+npj*7))
                write(61,162) tdvr(i),tdvr(i+npj),tdvr(i+npj*2) &
                  ,tdvr(i+npj*3),tdvr(i+npj*4),tdvr(i+npj*5) &
! *** m,rho,u
                  ,tdvr(i+npj*6),tdvr(i+npj*7),tdvr(i+npj*8) &
! *** mzHe,C,N,O,Ne,Mg,Si,Fe,Z
                  ,tdvr(i+npj*9),tdvr(i+npj*10),tdvr(i+npj*11) &
                  ,tdvr(i+npj*12),tdvr(i+npj*13),tdvr(i+npj*14) &
                  ,tdvr(i+npj*15),tdvr(i+npj*16),tdvr(i+npj*17) &
! *** id,flagfd,h,myu,nhp,temp,rp,nnb
                  ,tivr(i),tivr(i+npj),tdvr(i+npj*18),tdvr(i+npj*19) &
                  ,nhp,temp,rp,tivr(i+npj*2)
 162            format(18(1pE13.5),2I10,5(1pE13.5),I10)
              enddo
            endif
          enddo
          close(61)
          if(flagselp.eq.1) then
            close(63)
          endif  
        else
! *** sending the data to rank 0 ***
          ip=myrank
          call MPI_ISEND(tivs,np*nivals,MPI_INTEGER,0,myrank &
            ,MPI_COMM_WORLD,ireqs(ip),ierr)
          call MPI_WAIT(ireqs(ip),istatus,ierr)
          ip2=myrank*2
          call MPI_ISEND(tdvs,np*nvals,MPI_DOUBLE_PRECISION,0,ip2 &
           ,MPI_COMM_WORLD,ireqs(ip),ierr)
          call MPI_WAIT(ireqs(ip),istatus,ierr)
        endif    

        deallocate(tivs)
        deallocate(tivr)
        deallocate(tdvs)
        deallocate(tdvr)
      endif

      if(ndmt.gt.0) then
! *** output DM
        nvals=13
        nivals=1

        allocate(tivs(0:ndm*nivals-1))
        allocate(tivr(0:ndm*nivals-1))
        allocate(tdvs(0:ndm*nvals-1))
        allocate(tdvr(0:ndm*nvals-1))

        np=0
        do i=0,ndm-1,nskip
          np=np+1
        enddo
        nc=0
        do i=0,ndm-1,nskip
          pn=list_adm(i)
          tdvs(nc)=x_dm(pn)*asci
          tdvs(nc+np)=y_dm(pn)*asci
          tdvs(nc+np*2)=z_dm(pn)*asci
          tdvs(nc+np*3)=vx_dm(pn)
          tdvs(nc+np*4)=vy_dm(pn)
          tdvs(nc+np*5)=vz_dm(pn)
          tdvs(nc+np*6)=m_dm(pn)
          tdvs(nc+np*7)=rho_dm(pn)
          tdvs(nc+np*8)=h_dm(pn)
          tdvs(nc+np*9)=tadd_dm(pn)
          tdvs(nc+np*10)=rperi_dm(pn)
          tdvs(nc+np*11)=rapo_dm(pn)
          tdvs(nc+np*12)=zmax_dm(pn)
          tivs(nc)=id_dm(pn)
          nc=nc+1
        enddo
! *** set number of receiving data 
        if(myrank.eq.0) then
          do ip=0,nprocs-1
            if(ip.eq.myrank) then
              npjr(ip)=np
            else
              npjr(ip)=0
              call MPI_IRECV(npjr(ip),1,MPI_INTEGER,ip,ip &
                ,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
            endif
          enddo
        else
          ip=myrank
          call MPI_ISEND(np,1,MPI_INTEGER,0,myrank &
            ,MPI_COMM_WORLD,ireqs(ip),ierr)
          call MPI_WAIT(ireqs(ip),istatus,ierr)
        endif    
! *** receiving id_p and double data and write
        if(myrank.eq.0) then
          if(flagbo.eq.0) then
            write(fileo,'(a8,i6.6,a4)') 'output/d',step,'.bin'
            open(61,file=fileo,status='unknown',form='unformatted')
            write(61) ndm1t
          else
            write(fileo,'(a8,i6.6)') 'output/d',step
            open(61,file=fileo,status='unknown')
            write(fileo,'(a9,i6.6)') 'output/ld',step
            open(62,file=fileo,status='unknown')
          endif
          if(flagselp.eq.2) then
            write(fileo,'(a11,i6.6,a4)') 'output/seld',step,'.bin'
!            open(63,file=fileo,status='unknown')
!            write(63,*) '# ',nselp                        
            open(63,file=fileo,status='unknown',form='unformatted')
            write(63) nselp
          endif             
          do ip=0,nprocs-1
            if(ip.eq.myrank) then
              do i=0,npjr(myrank)-1
                tivr(i)=tivs(i)
              enddo
              do i=0,npjr(myrank)*nvals-1
                tdvr(i)=tdvs(i)
              enddo
            else
              call MPI_IRECV(tivr,npjr(ip),MPI_INTEGER,ip,ip &
                ,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
              ip2=ip*2
              call MPI_IRECV(tdvr,npjr(ip)*nvals,MPI_DOUBLE_PRECISION &
               ,ip,ip2,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
            endif
! *** write the data
            npj=npjr(ip)
            ndmo=ndmo+npj
            if(flagselp.eq.2) then
              do i=0,npj-1
                rp=dsqrt(tdvr(i)**2+tdvr(i+npj)**2+tdvr(i+npj*2)**2)
                if(rp.gt.rrange(0).and.rp.lt.rrange(1)) then
                  if(tdvr(i+npj*2).gt.zrange(0) &
                    .and.tdvr(i+npj*2).lt.zrange(1)) then
                    write(63) tdvr(i),tdvr(i+npj),tdvr(i+npj*2) &
                         ,tdvr(i+npj*3),tdvr(i+npj*4),tdvr(i+npj*5)
!                    write(63,*) tdvr(i),tdvr(i+npj),tdvr(i+npj*2) &
!                      ,tdvr(i+npj*3),tdvr(i+npj*4),tdvr(i+npj*5)
                  endif
                endif
              enddo
            endif  
            if(flagbo.eq.0) then
              do i=0,npj-1
                if(tivr(i).lt.ndm1t) then
                  write(61) tdvr(i),tdvr(i+npj),tdvr(i+npj*2),tivr(i)
                endif
              enddo
            else
              if(flago.eq.1) then
! test particle output
                do i=0,npj-1
                  if(tdvr(i+npj*11).gt.0.0d0) then
! mean radius
                    rmp=0.5d0*(tdvr(i+npj*10)+tdvr(i+npj*11))
! eccentricity
                    eccp=(tdvr(i+npj*11)-tdvr(i+npj*10)) &
                      /(tdvr(i+npj*11)+tdvr(i+npj*10))
                  else
                    rmp=0.0d0
                    eccp=0.0d0
                  endif
                  write(61,'(12(1pE13.5),I10)') tdvr(i),tdvr(i+npj) &
                   ,tdvr(i+npj*2),tdvr(i+npj*3),tdvr(i+npj*4) &
                   ,tdvr(i+npj*5),tdvr(i+npj*9),tdvr(i+npj*10) &
                   ,tdvr(i+npj*11),tdvr(i+npj*12),rmp,eccp,tivr(i)
                enddo
              else
                do i=0,npj-1
                  if(tivr(i).lt.ndm1t) then
                    write(61,163) tdvr(i),tdvr(i+npj),tdvr(i+npj*2) &
                      ,tdvr(i+npj*3),tdvr(i+npj*4),tdvr(i+npj*5) &
                      ,tdvr(i+npj*6),tdvr(i+npj*7),tdvr(i+npj*8) &
                      ,tdvr(i+npj*9),tivr(i)
                  else
                    write(62,163) tdvr(i),tdvr(i+npj),tdvr(i+npj*2) &
                      ,tdvr(i+npj*3),tdvr(i+npj*4),tdvr(i+npj*5) &
                      ,tdvr(i+npj*6),tdvr(i+npj*7),tdvr(i+npj*8) &
                      ,tdvr(i+npj*9),tivr(i)
                  endif
 163              format(10(1pE13.5),I10)
                enddo
              endif
            endif
          enddo
          close(61)
          if(flagbo.ne.0) then
            close(62)
          endif
          if(flagselp.eq.2) then
            close(63)
          endif
        else
          ip=myrank
          call MPI_ISEND(tivs,np,MPI_INTEGER,0,myrank &
            ,MPI_COMM_WORLD,ireqs(ip),ierr)
          call MPI_WAIT(ireqs(ip),istatus,ierr)
          ip2=myrank*2
          call MPI_ISEND(tdvs,np*nvals,MPI_DOUBLE_PRECISION,0,ip2 &
            ,MPI_COMM_WORLD,ireqs(ip),ierr)
          call MPI_WAIT(ireqs(ip),istatus,ierr)
        endif    

        deallocate(tivs)
        deallocate(tivr)
        deallocate(tdvs)
        deallocate(tdvr)

      endif

! *** output star
      if(nst.gt.0) then
        if(flagr.eq.0) then
          nvals=7
          nivals=1
        else
          nvals=20
          nivals=3
        endif

        allocate(tivs(0:ns*nivals-1))
        allocate(tivr(0:ns*nivals-1))
        allocate(tdvs(0:ns*nvals-1))
        allocate(tdvr(0:ns*nvals-1))

        np=0
        do i=0,ns-1,nskip
          np=np+1
        enddo
        nc=0
        do i=0,ns-1,nskip
          pn=list_ap(i+ng)
          tdvs(nc)=x_p(pn)*asci
          tdvs(nc+np)=y_p(pn)*asci
          tdvs(nc+np*2)=z_p(pn)*asci
          tdvs(nc+np*3)=vx_p(pn)
          tdvs(nc+np*4)=vy_p(pn)
          tdvs(nc+np*5)=vz_p(pn)
          tdvs(nc+np*6)=m_p(pn)
          tivs(nc)=id_p(pn)
          nc=nc+1
        enddo
        if(flagr.gt.0) then
          nc=0
          do i=0,ns-1,nskip
            pn=list_ap(i+ng)
            tdvs(nc+np*9)=mzHe_p(pn)
            tdvs(nc+np*10)=mzC_p(pn)
            tdvs(nc+np*11)=mzN_p(pn)
            tdvs(nc+np*12)=mzO_p(pn)
            tdvs(nc+np*13)=mzNe_p(pn)
            tdvs(nc+np*14)=mzMg_p(pn)
            tdvs(nc+np*15)=mzSi_p(pn)
            tdvs(nc+np*16)=mzFe_p(pn)
            tdvs(nc+np*17)=mzZ_p(pn)
            tdvs(nc+np*18)=h_p(pn)
            tdvs(nc+np*19)=ts_p(pn)
            tivs(nc+np)=flagfd_p(pn)
            tivs(nc+np*2)=nnb_p(pn)
            nc=nc+1
          enddo
        endif
! *** set number of receiving data 
        if(myrank.eq.0) then
          do ip=0,nprocs-1
            if(ip.eq.myrank) then
              npjr(ip)=np
            else
              npjr(ip)=0
              call MPI_IRECV(npjr(ip),1,MPI_INTEGER,ip,ip &
                ,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
            endif
          enddo
        else
          ip=myrank
          call MPI_ISEND(np,1,MPI_INTEGER,0,myrank &
            ,MPI_COMM_WORLD,ireqs(ip),ierr)
          call MPI_WAIT(ireqs(ip),istatus,ierr)
        endif    
! *** receiving id_p and double data and write
        if(myrank.eq.0) then
          write(fileo,'(a8,i6.6)') 'output/s',step
          open(61,file=fileo,status='unknown')
          if(flagselp.eq.1) then
            write(fileo,'(a11,i6.6,a4)') 'output/sels',step,'.bin'
            open(63,file=fileo,status='unknown',form='unformatted')
            write(63) nselp
          endif             
          do ip=0,nprocs-1
            if(ip.eq.myrank) then
              do i=0,npjr(myrank)*nivals-1
                tivr(i)=tivs(i)
              enddo
              do i=0,npjr(myrank)*nvals-1
                tdvr(i)=tdvs(i)
              enddo
            else
!            write(6,*) ' star recv nivals,nvails,npjr(ip),ip=' &
!              ,nivals,nvals,npjr(ip),ip
              call MPI_IRECV(tivr,npjr(ip)*nivals,MPI_INTEGER,ip,ip &
                ,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
              ip2=ip*2
              call MPI_IRECV(tdvr,npjr(ip)*nvals,MPI_DOUBLE_PRECISION &
                ,ip,ip2,MPI_COMM_WORLD,ireqr(ip),ierr)
              call MPI_WAIT(ireqr(ip),istatus,ierr)
            endif
! *** write the data
            npj=npjr(ip)
            nso=nso+npj
            if(flagselp.eq.3) then
              do i=0,npj-1
                rp=dsqrt(tdvr(i)**2+tdvr(i+npj)**2+tdvr(i+npj*2)**2)
                if(rp.gt.rrange(0).and.rp.lt.rrange(1)) then
                  if(tdvr(i+npj*2).gt.zrange(0) &
                    .and.tdvr(i+npj*2).lt.zrange(1)) then
                    write(63) tdvr(i),tdvr(i+npj),tdvr(i+npj*2) &
                      ,tdvr(i+npj*3),tdvr(i+npj*4),tdvr(i+npj*5)
                  endif
                endif
              enddo
            endif  
            if(flagr.eq.0) then
              do i=0,npj-1
                write(61,164) tdvr(i),tdvr(i+npj),tdvr(i+npj*2) &
                 ,tdvr(i+npj*3),tdvr(i+npj*4),tdvr(i+npj*5) &
                 ,tdvr(i+npj*6),tivr(i) &
                 ,dsqrt(tdvr(i)**2+tdvr(i+npj)**2+tdvr(i+npj*2)**2)
 164           format(7(1pE13.5),I10,1pE13.5)
              enddo
            else
              do i=0,npj-1
                rp=dsqrt(tdvr(i)**2+tdvr(i+npj)**2+tdvr(i+npj*2)**2) 
                write(61,165) tdvr(i),tdvr(i+npj),tdvr(i+npj*2) &
                 ,tdvr(i+npj*3),tdvr(i+npj*4),tdvr(i+npj*5),tdvr(i+npj*6) &
! *** mzHe,C,N,O,Ne,Mg,Si,Fe,Z
                 ,tdvr(i+npj*9),tdvr(i+npj*10),tdvr(i+npj*11) &
                 ,tdvr(i+npj*12),tdvr(i+npj*13),tdvr(i+npj*14) &
                 ,tdvr(i+npj*15),tdvr(i+npj*16),tdvr(i+npj*17) &
! *** ts,id,flagfd,h,rp,nnb
                 ,tdvr(i+npj*19),tivr(i),tivr(i+npj) &
                 ,tdvr(i+npj*18),rp,tivr(i+npj*2)
 165             format(17(1pE13.5),2I10,2(1pE13.5),I10)
              enddo
            endif
          enddo
          close(61)
          if(flagselp.eq.3) then
            close(63)
          endif   
        else
          ip=myrank
          call MPI_ISEND(tivs,np*nivals,MPI_INTEGER,0,myrank &
            ,MPI_COMM_WORLD,ireqs(ip),ierr)
          call MPI_WAIT(ireqs(ip),istatus,ierr)
          ip2=myrank*2
          call MPI_ISEND(tdvs,np*nvals,MPI_DOUBLE_PRECISION,0,ip2 &
           ,MPI_COMM_WORLD,ireqs(ip),ierr)
          call MPI_WAIT(ireqs(ip),istatus,ierr)
        endif    

        if(myrank.eq.0) then
          write(6,*) ' output ng,ndm,ns=',ngo,ndmo,nso
        endif
      endif

      deallocate(istatus)
      deallocate(npjr)
      deallocate(ireqs)
      deallocate(ireqr)

 90   return
      end
