!
!   main-vmetmap.f90 
! ver.4 18  Jan.  2015      written by D. Kawata
!

program vmetmap

      implicit none

      integer pgopen
      integer,parameter :: nxp=4
      integer,parameter :: nyp=1
! *** for work ***
      integer nix,nstep,step,nval
      integer i,j,k,ctype,ip,flagc,istep
      integer maxi,maxj,mini,minj
! *** input parameter ***
! *** grid data to read ***
      integer nx,ny
      double precision xrange(2),yrange(2),tu
      integer flaginc
      double precision iang,angc
! for plot
      integer lwdef
      real xl,yl,xh,yh
      real chdef
! *** for mesh data ***
      integer,allocatable :: idr(:,:)
      double precision txp,typ,tzp,tfdp,rtmp(7),vroty,vzivrot
      double precision,allocatable :: xm(:,:),ym(:,:),denxy(:,:),denxz(:,:) &
       ,vrotm(:,:),vradm(:,:),vzm(:,:),metm(:,:),rpix(:,:),ddenm(:,:)
      double precision,allocatable :: fdval(:,:)
      real,allocatable :: fp(:,:),fdenp(:,:)
! for radial mean value
      integer ndr
      integer,allocatable :: nprm(:)
      double precision drrm,rmaxrm,rminrm,vrotrmp,vradrmp,denrmp,metrmp
      double precision,allocatable :: mrm(:),vrotrm(:),vradrm(:) &
       ,rrm(:),metrm(:),msrm(:)
! for viewpoint
      real rxvp(2),ryvp(2),dxvp,dyvp
      real xlvp(nxp),xhvp(nxp),ylvp(nyp),yhvp(nyp)
! *** for plot ***
      real fmax,fmin,bright,contra,fdenmin,fdenmax
      real fmaxv(5),fminv(5)
      integer ncon
      real fmaxc,fminc,alev,dlev
      real tr(6)      
! plot circle
      integer npcirc,nradc,ic
      real ddeg,radc,xcirc(100),ycirc(100),degp
! ***  for filename ***
      character filen*60,mlab*60
! solar metallicity
      double precision,parameter :: ZSOL=0.019d0
      double precision,parameter :: METMIN=-2.0d0
      double precision,parameter :: M_PI=3.141592653589793238d0

      chdef=2.0
      lwdef=3
      
! *** open input file ****
      open(50,file='ini/inputp-vmetmap.dat',status='old')
      read(50,*) nstep,step
      read(50,*) flagc
      read(50,*) fminv(1),fmaxv(1)
      read(50,*) fminv(2),fmaxv(2)
      read(50,*) fminv(3),fmaxv(3)
      read(50,*) fminv(4),fmaxv(4)
      read(50,*) fminv(5),fmaxv(5)
      close(50)
      nval=4

      if(flagc.eq.0) then      
        write(6,*) ' plot for gas'
      else if(flagc.eq.1) then
        write(6,*) ' plot for DM'
      else
        write(6,*) ' plot for star'
      endif
      write(6,*) ' rad velocity map range (km/s) =',fmaxv(1),fminv(1)
      write(6,*) ' rot velocity map range (km/s) =',fmaxv(2),fminv(2)
      write(6,*) ' [Z] map range (dex) =',fmaxv(3),fminv(3)
      write(6,*) ' image density range =',fmaxv(4),fminv(4)
      write(6,*) ' density contrast range =',fmaxv(5),fminv(5)

      ctype=2

      if(nstep.gt.1) then
        open(50,file='ini/file.dat',status='old') 
      endif
      do istep=1,nstep
        if(nstep.gt.1) then
          read(50,*) step
        endif
        if(flagc.eq.0) then
          write(filen,'(a13,i6.6)') 'output/gden2d',step
        else if(flagc.eq.1) then
          write(filen,'(a13,i6.6)') 'output/dden2d',step
        else
          write(filen,'(a13,i6.6)') 'output/sden2d',step
        endif
        write(6,*) ' reading ',filen
! ****    set filename for data file   *****
        open(51,file=filen,status='old',form='unformatted')
        read(51) nx,ny
        read(51) xrange(1),xrange(2)
        read(51) yrange(1),yrange(2)
        read(51) tu
        read(51) flaginc
        read(51) iang
        write(6,*) ' nx,ny=',nx,ny
        write(6,*) ' x range =',xrange(1),xrange(2)
        write(6,*) ' y range =',yrange(1),yrange(2)
        write(6,*) ' time =',tu
        if(flaginc.ne.0) then
          write(6,*) 'inclination angle=',iang
        endif
! *** print output file ***
        allocate(xm(nx,ny))
        allocate(ym(nx,ny))
        allocate(denxy(nx,ny))
        allocate(denxz(nx,ny))
        allocate(vrotm(nx,ny))
        allocate(vradm(nx,ny))
        allocate(vzm(nx,ny))
        allocate(metm(nx,ny))
        allocate(fdval(nx,ny))
        allocate(fp(nx,ny))
        allocate(fdenp(nx,ny))
        allocate(rpix(nx,ny))
        allocate(ddenm(nx,ny))
        allocate(idr(nx,ny))

        do j=1,ny
          do i=1,nx
            read(51) xm(i,j),ym(i,j),denxy(i,j),denxz(i,j),vrotm(i,j) &
             ,vradm(i,j),vzm(i,j),metm(i,j)
! metallicity -> [Z]
            if(metm(i,j)/ZSOL.gt.10.0d0**METMIN) then
              metm(i,j)=dlog10(metm(i,j)/ZSOL)
            else
              write(6,*) ' [Z] is smaller than METMIN:',i,j,metm(i,j)
              metm(i,j)=METMIN
            endif
          enddo
        enddo

! set ddenm, total stellar density
        if(flagc.eq.2) then
          do j=1,ny
            do i=1,nx
              ddenm(i,j)=denxy(i,j)
            enddo
          enddo
        else
! read stellar data
          write(filen,'(a13,i6.6)') 'output/sden2d',step
! ****    set filename for data file   *****
          open(51,file=filen,status='old',form='unformatted',err=90)
          goto 70
  90      write(6,*) ' Error: no stellar data. Please run sden2d for star'
          stop
  70      do i=1,6
            read(51)
          enddo
          do j=1,ny
            do i=1,nx
              read(51) rtmp(1),rtmp(2),ddenm(i,j),rtmp(3),rtmp(4),rtmp(5) &  
                ,rtmp(6),rtmp(7)
            enddo
          enddo
        endif
! get mean value 
        rminrm=0.0d0
!        rmaxrm=dsqrt(2.0d0)*0.5d0*(xrange(2)-xrange(1))
        rmaxrm=0.5d0*(xrange(2)-xrange(1))
        write(6,*) ' R range for radial mean value=',rminrm,rmaxrm
        ndr=int(nx/6)
        drrm=(rmaxrm-rminrm)/dble(ndr)
        write(6,*) ' Ndr,drrm,dpix=',ndr,drrm,(xrange(2)-xrange(1))/dble(nx)

        allocate(nprm(ndr))
        allocate(mrm(ndr))
        allocate(msrm(ndr))
        allocate(vradrm(ndr))
        allocate(vrotrm(ndr))
        allocate(metrm(ndr))
        allocate(rrm(ndr))

        do i=1,ndr
          nprm(i)=0
          mrm(i)=0.0d0
          msrm(i)=0.0d0
          vradrm(i)=0.0d0
          vrotrm(i)=0.0d0
          metrm(i)=0.0d0
          rrm(i)=rminrm+(dble(i)-0.5d0)*drrm
        enddo

! set radius for each pixel
        if(flaginc.eq.0) then
          do j=1,ny
            do i=1,nx
              rpix(i,j)=dsqrt(xm(i,j)**2+ym(i,j)**2)
            enddo
          enddo       
        else  
! deg->rad
          iang=iang*M_PI/180.0d0    
          angc=1.0d0/dcos(iang)
          do j=1,ny
            do i=1,nx
              typ=ym(i,j)*angc
! radius after the inclination correction
              rpix(i,j)=dsqrt(xm(i,j)**2+typ**2)
            enddo
          enddo
        endif

!        open(60,file='testpix.dat',status='unknown')
        do j=1,ny
          do i=1,nx
            if(rpix(i,j).lt.rmaxrm) then
              idr(i,j)=(rpix(i,j)-rminrm)/drrm+1
!            write(60,'(3(1pE13.5),I10)') xm(i,j),ym(i,j),rpix(i,j),idr(i,j)
              if(idr(i,j).gt.0.and.idr(i,j).lt.ndr+1) then
                nprm(idr(i,j))=nprm(idr(i,j))+1
                mrm(idr(i,j))=mrm(idr(i,j))+denxy(i,j)
                vradrm(idr(i,j))=vradrm(idr(i,j))+denxy(i,j)*vradm(i,j)
                vrotrm(idr(i,j))=vrotrm(idr(i,j))+denxy(i,j)*vrotm(i,j)
                metrm(idr(i,j))=metrm(idr(i,j))+metm(i,j)
                msrm(idr(i,j))=msrm(idr(i,j))+ddenm(i,j)
              else
                write(6,*) ' cannot find idr at x,y,r,idr=',xm(i,j),ym(i,j) &
                 ,rpix(i,j),idr(i,j)
              endif
            endif
          enddo
        enddo
! getting mean and output
        if(flagc.eq.0) then
          write(filen,'(a14,i6.6)') 'output/gvmetrm',step
        else
          write(filen,'(a14,i6.6)') 'output/svmetrm',step
        endif
        open(60,file=filen,status='unknown')
        write(60,'(a9,I10,(1pE13.5))') '# ndr, dr',ndr,drrm
        do i=1,ndr
          if(nprm(i).gt.0) then
            vradrm(i)=vradrm(i)/mrm(i)
            vrotrm(i)=vrotrm(i)/mrm(i)
            mrm(i)=mrm(i)/dble(nprm(i))
            metrm(i)=metrm(i)/dble(nprm(i))
            msrm(i)=msrm(i)/dble(nprm(i))
          endif
          write(60,'(6(1pE13.5),I10)') rrm(i),mrm(i),vradrm(i),vrotrm(i) &
           ,metrm(i),msrm(i),nprm(i)
        enddo        
        close(60)
! subtract the mean
        open(60,file='testp.dat',status='unknown')
        do j=1,ny
          do i=1,nx
            if(rpix(i,j).lt.rmaxrm) then
              if(idr(i,j).ge.1.and.idr(i,j).le.ndr) then
                if(rpix(i,j).gt.rrm(idr(i,j)).or.idr(i,j).eq.1) then
                  vradrmp=vradrm(idr(i,j))+(rpix(i,j)-rrm(idr(i,j))) &
                   *(vradrm(idr(i,j)+1)-vradrm(idr(i,j)))/drrm
                  vrotrmp=vrotrm(idr(i,j))+(rpix(i,j)-rrm(idr(i,j))) &
                   *(vrotrm(idr(i,j)+1)-vrotrm(idr(i,j)))/drrm
                  denrmp=msrm(idr(i,j))+(rpix(i,j)-rrm(idr(i,j))) &
                   *(msrm(idr(i,j)+1)-msrm(idr(i,j)))/drrm
                  metrmp=metrm(idr(i,j))+(rpix(i,j)-rrm(idr(i,j))) &
                   *(metrm(idr(i,j)+1)-metrm(idr(i,j)))/drrm
                else
                  vradrmp=vradrm(idr(i,j)-1)+(rpix(i,j)-rrm(idr(i,j)-1)) &
                   *(vradrm(idr(i,j))-vradrm(idr(i,j)-1))/drrm
                  vrotrmp=vrotrm(idr(i,j)-1)+(rpix(i,j)-rrm(idr(i,j)-1)) &
                   *(vrotrm(idr(i,j))-vrotrm(idr(i,j)-1))/drrm
                  denrmp=msrm(idr(i,j)-1)+(rpix(i,j)-rrm(idr(i,j)-1)) &
                   *(msrm(idr(i,j))-msrm(idr(i,j)-1))/drrm
                  metrmp=metrm(idr(i,j)-1)+(rpix(i,j)-rrm(idr(i,j)-1)) &
                   *(metrm(idr(i,j))-metrm(idr(i,j)-1))/drrm
                endif
              else
! use the last point
                vradrmp=vradrm(ndr)
                vrotrmp=vrotrm(ndr)
                denrmp=msrm(ndr)
                metrmp=metrm(ndr)
              endif
              write(60,'(5(1pE13.5))') rpix(i,j),vradrmp,vrotrmp,denrmp,metrmp
! subtract the mean value
              vradm(i,j)=(vradm(i,j)-vradrmp)
              vrotm(i,j)=(vrotm(i,j)-vrotrmp)
              ddenm(i,j)=(ddenm(i,j)-denrmp)/denrmp
              metm(i,j)=metm(i,j)-metrmp
              vzivrot=0.0d0
              if(flaginc.ne.0.and.rpix(i,j).gt.0.0d0) then
! vrotrmp contribution to vz
! original y poisiotn
                typ=ym(i,j)*angc
! original y-direction velocity from vrotrmp
                vroty=-vrotrmp*xm(i,j)/rpix(i,j)
! vz direction correction after inclination
                vzivrot=vroty*dsin(iang)
                vzm(i,j)=vzm(i,j)-vzivrot
              endif
              write(60,'(8(1pE13.5))') rpix(i,j),vradrmp,vrotrmp,denrmp,metrmp &
               ,xm(i,j),ym(i,j),vzivrot
            else
! set minimum value
              vradm(i,j)=0.0d0
              vrotm(i,j)=0.0d0
              vzm(i,j)=0.0d0
              metm(i,j)=fminv(3)
              ddenm(i,j)=fminv(4)
            endif
          enddo
        enddo
        close(60)


        deallocate(nprm)
        deallocate(mrm)
        deallocate(vradrm)
        deallocate(vrotrm)
        deallocate(rrm)

! *** Input Parameter ***
        xl = real(xrange(1))
        yl = real(yrange(1))
        xh = real(xrange(2))
        yh = real(yrange(2))

        write(*,*) ' xl,yh,yl,yh     = ',xl,xh,yl,yh
! *** set transfer array ***
        tr(2) = (xh-xl)/real(nx)
        tr(1) = xl-0.5*tr(2)
        tr(3) = 0.0
        tr(5) = 0.0      
        tr(6) = (yh-yl)/real(ny)
        tr(4) = yl-0.5*tr(6)

! set view point
        rxvp(1) = 0.05
        rxvp(2) = 0.975
        dxvp = (rxvp(2)-rxvp(1))/real(nxp)
        do i=1,nxp
          xlvp(i)=rxvp(1)+real(i-1)*dxvp
          xhvp(i)=rxvp(1)+real(i)*dxvp
        enddo

        ryvp(1) = 0.2
        ryvp(2) = 0.9
        dyvp = (ryvp(2)-ryvp(1))/real(nyp)
        do i=1,nyp
          ylvp(i)=ryvp(2)-real(i)*dyvp
          yhvp(i)=ryvp(2)-real(i-1)*dyvp
        enddo

! set output file name
        if(flagc.eq.0) then
          write(filen,'(a13,i6.6,a8)') 'output/pgvmet',step,'.gif/gif'
        else if(flagc.eq.1) then
          write(filen,'(a13,i6.6,a8)') 'output/pdvmet',step,'.gif/gif'
        else
          write(filen,'(a13,i6.6,a8)') 'output/psvmet',step,'.gif/gif'
        endif

! *** plot 2d image at ix=nix ***
! *** open device ****
        if (pgopen(filen) .ne.1) then
          write(*,*) ' Cannot open device!'
          stop
        endif
        call pgpap(12.,(dxvp/dyvp)*(yh-yl)/(xh-xl))
        call pgsch(chdef)
        call pgslw(lwdef)

        call pgscr(1,0.0,0.0,0.0)
        call pgscr(0,1.0,1.0,1.0)
        call pgeras

! *** write header ***
        call pgsvp(0.0,1.0,ryvp(2),1.0)
        call pgswin(0.0,1.0,0.0,1.0)
        write(6,*) ' step, t=',step,tu
        if(flaginc.eq.0) then
          write(mlab,'(a3,f7.3,a6)') ' t=',tu,' (Gyr)'
        else
          iang=iang*180.0d0/M_PI
          write(mlab,'(a3,f7.3,a9,f5.1,a6)') ' t=',tu,' (Gyr) i=',iang,' (deg)'
          iang=iang*M_PI/180.0
        endif
        call pgptxt(0.5,0.5,0.,0.5,mlab)

        do ip=1,nxp
! *** Set up window and viewport ***
          call pgsvp(xlvp(ip),xhvp(ip),ylvp(1),yhvp(1))
          call pgswin(xl,xh,yl,yh)

! set fdval
          if(ip.eq.1) then
            do j=1,ny
              do i=1,nx
                if(flaginc.eq.0) then
                  fdval(i,j)=vradm(i,j)
                else
                  fdval(i,j)=vzm(i,j)
                endif
              enddo
            enddo
          else if(ip.eq.2) then
            do j=1,ny
              do i=1,nx
                fdval(i,j)=vrotm(i,j)
              enddo
            enddo
          else if(ip.eq.3) then
            do j=1,ny
              do i=1,nx
                fdval(i,j)=metm(i,j)
              enddo
            enddo
          else
            do j=1,ny
              do i=1,nx
                fdval(i,j)=denxy(i,j)
              enddo
            enddo
          endif

! ***** set plotting file *****
          fmin = 1.0e12
          fmax = -1.0e12
          fdenmin=1.0e12
          fdenmax=-1.0e12
          if(ip.ge.4) then
            do j = 1,ny
              do i = 1,nx
                if(fdval(i,j).gt.0.0d0) then
                  fp(i,j)=real(dlog10(fdval(i,j)))
                  if(fp(i,j).gt.fmax) then
                    fmax = fp(i,j)
                    maxi = i
                    maxj = j                        
                  endif
                  if(fp(i,j).lt.fmin) then
                    fmin = fp(i,j)
                    mini = i
                    minj = j            
                  endif
                else
                  fp(i,j)=-999.9
                endif
              enddo
            enddo
          else
            do j = 1,ny
              do i = 1,nx
                fp(i,j)=real(fdval(i,j))
                if(fp(i,j).gt.fmax) then
                  fmax = fp(i,j)
                  maxi = i
                  maxj = j                        
                endif
                if(fp(i,j).lt.fmin) then
                  fmin = fp(i,j)
                  mini = i
                  minj = j            
                endif
! density contrast
                fdenp(i,j)=real(ddenm(i,j))
                if(fdenp(i,j).gt.fdenmax) then
                  fdenmax = fdenp(i,j)
                endif
                if(fdenp(i,j).lt.fdenmin) then
                  fdenmin = fdenp(i,j)
                endif
              enddo
            enddo
          endif
          write(*,*) ip,' fmax, fmin      = ',fmax,fmin
          write(6,*) ' colour image density range=',fmaxv(ip),fminv(ip)
          if(ip.eq.1) then
            write(*,*) ' density contrast max,min =',fdenmax,fdenmin
          endif
          fmax=fmaxv(ip)
          fmin=fminv(ip)

          call pgbox('bcts',0.,0,'bcts',0.,0)
          if(ip.eq.1) then
            call pgbox('bcts',0.,0,'bncts',0.,0)
            call pgmtxt('L',2.0,0.5,0.5,'(kpc)')
!          else
!            call pgbox('bncts',0.,0,'bcts',0.,0)
!            call pgmtxt('B',2.0,0.5,0.5,'(kpc)')
          endif
! *** Set up the color map ***
          bright = 0.5
          contra = 1.0
          call palett(ctype,contra,bright)
! *** Draw the map with pgimg ***
          call pgimag(fp(1,1),nx,ny,1,nx,1,ny,fmin,fmax,tr)

! density contour 
          ncon=5
          dlev=(fmaxv(5)-fminv(5))/real(ncon)
          do i=1,ncon+1
            alev=fminv(5)+dlev*real(i-1)
            call pgcont(fdenp,nx,ny,1,nx,1,nx,alev,-1,tr)
          enddo
! circle
          call pgsci(2)
          call pgsfs(2)
          call pgsls(2)
          npcirc=100
          ddeg=2.0*real(M_PI)/real(npcirc)
          nradc=3
          do ic=1,nradc  
            if(ic.eq.1) then
              radc=4.5
            else if(ic.eq.2) then
              radc=6.3
            else
              radc=8.0
            endif
!            radc=4.0+real(ic-1)*2.0
!            write(6,*) ' draw circle R=',radc
            do i=1,npcirc
              xcirc(i)=radc*cos(ddeg*real(i))
              ycirc(i)=radc*sin(ddeg*real(i))
              if(flaginc.ne.0) then
                ycirc(i)=ycirc(i)*cos(real(iang))
              endif
            enddo
            call pgline(npcirc,xcirc,ycirc)
          enddo
          call pgsls(1)
          call pgsci(1)

! *** Draw a Wedge ***
          if(ip.eq.1) then
            if(flaginc.eq.0) then
              call pgwedg('BI',0.,4.0,fminv(ip),fmaxv(ip),'\gDV\drad\u (km/s)') 
            else
              call pgwedg('BI',0.,4.0,fminv(ip),fmaxv(ip),'\gDV\dLOS\u (km/s)') 
            endif
          else if(ip.eq.2) then
            call pgwedg('BI',0.,4.0,fminv(ip),fmaxv(ip),'\gDV\drot\u (km/s)') 
          else if(ip.eq.3) then
            call pgwedg('BI',0.,4.0,fminv(ip),fmaxv(ip),'\gD[Z]') 
          endif

        enddo

        deallocate(xm)
        deallocate(ym)
        deallocate(fdval)
        deallocate(denxy)
        deallocate(denxz)
        deallocate(vrotm)
        deallocate(vradm)
        deallocate(vzm)
        deallocate(metm)
        deallocate(fp)
        deallocate(fdenp)
        deallocate(rpix)
        deallocate(ddenm)
        deallocate(idr)

! *** close the device and exit ***
        call pgend

      enddo
      close(51)


      stop
end program vmetmap

!***** from pgdemo4.f *****
      
SUBROUTINE PALETT(TYPE, CONTRA, BRIGHT)
!-----------------------------------------------------------------------
! Set a "palette" of colors in the range of color indices used by
! PGIMAG.
!C-----------------------------------------------------------------------
      INTEGER TYPE
      REAL CONTRA, BRIGHT
!
      REAL GCL(2), GCR(2), GCG(2), GCB(2)      
      REAL GL(2), GR(2), GG(2), GB(2)
      REAL RL(9), RR(9), RG(9), RB(9)
      REAL HL(5), HR(5), HG(5), HB(5)
      REAL WL(10), WR(10), WG(10), WB(10)
      REAL AL(20), AR(20), AG(20), AB(20)

      DATA GL /0.0, 1.0/
      DATA GR /0.0, 1.0/
      DATA GG /0.0, 1.0/
      DATA GB /0.0, 1.0/

      DATA RL /-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7/
      DATA RR / 0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0/
      DATA RG / 0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0/
      DATA RB / 0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0/

      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
      DATA HR /0.0, 0.5, 1.0, 1.0, 1.0/
      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
      DATA HB /0.0, 0.0, 0.0, 0.3, 1.0/

!c      DATA HL /0.0, 0.2, 0.4, 0.6, 1.0/
!c      DATA HB /0.0, 0.5, 1.0, 1.0, 1.0/
!c      DATA HG /0.0, 0.0, 0.5, 1.0, 1.0/
!c      DATA HR /0.0, 0.0, 0.0, 0.3, 1.0/

      DATA WL /0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0/
      DATA WR /0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0/
      DATA WG /0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0/
      DATA WB /0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0/

      DATA AL /0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, &
              0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0/
      DATA AR /0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, &
              0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/
      DATA AG /0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8, &
              0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0/
      DATA AB /0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9, &
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/

      IF (TYPE.EQ.1) THEN
!        -- gray scale
         CALL PGCTAB(GL, GR, GG, GB, 2, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.2) THEN
!        -- rainbow
         CALL PGCTAB(RL, RR, RG, RB, 9, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.3) THEN
!        -- heat
         CALL PGCTAB(HL, HR, HG, HB, 5, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.4) THEN
!        -- weird IRAF
         CALL PGCTAB(WL, WR, WG, WB, 10, CONTRA, BRIGHT)
      ELSE IF (TYPE.EQ.5) THEN
!        -- AIPS
         CALL PGCTAB(AL, AR, AG, AB, 20, CONTRA, BRIGHT)
      END IF
END subroutine
