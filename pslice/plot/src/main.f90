!
!   main.f90 for plot_slval
! 14 Aug. 2017      written by D. Kawata
!

program plot_slval

      implicit none

      integer pgopen
      integer,parameter :: nxp=1
      integer,parameter :: nyp=1
      double precision,parameter :: VUKMS=207.4d0
! *** for work ***
      integer i,j,k,ctype,ip,ixp,iyp,flaglog,flagv,ipval
      integer step
      integer maxi,maxj,mini,minj
      integer flagsp
      integer nval,iv
      integer nskipa
! *** input parameter ***
! *** grid data to read ***
      integer nx,ny,flagid,ids,ide,flag3d,flagrot
      double precision xld,xhd,yld,yhd,zpd,tpd,zhd,zld,hmp
      double precision dx,dy,ds
      real vdx,vrot,rxyp
      real xar(2),yar(2)
! for plot
      integer lwdef
      real xl,yl,xh,yh
      real chdef
! *** for mesh data ***
      double precision txp,typ,tzp,tfdp
      double precision,allocatable :: xm(:,:),ym(:,:)
      double precision,allocatable :: fdval(:,:,:)
      double precision drtmp(4)
      real alev(1)
      real,allocatable :: fp(:,:)
! for star data
      integer nsp
      real,allocatable :: xsp(:),ysp(:),zsp(:),rxysp(:),vrotsp(:),agesp(:)
! for viewpoint
      real rxvp(2),ryvp(2),dxvp,dyvp
      real xlvp(nxp),xhvp(nxp),ylvp(nyp),yhvp(nyp)
! *** for plot ***
      real fmax,fmin,bright,contra,fcmin,fcmax
      real fmaxv,fminv
      integer ncon
      real fmaxc,fminc
      real tr(6)      
! ***  for filename ***
      character filen*60,mlab*60,rchr*5,wlab*60

      chdef=1.5
      lwdef=3
!  plot type, 1: gray scale, 2: rainbow
      ctype=2
      
! *** open input file ****
      open(50,file='ini/input-plot_slval.dat',status='old')
      read(50,*) step
      read(50,*) flagsp
      read(50,*) fminv,fmaxv
      read(50,*) flaglog
      read(50,*) flagv
      read(50,*) ipval
      read(50,'(a60)') wlab
      close(50)

      write(6,*) ' step=',step
      if(flagsp.eq.0) then
        write(6,*) ' for gas'
      else if(flagsp.eq.1) then
        write(6,*) ' for DM'
      else if(flagsp.eq.2) then
        write(6,*) ' for star'
      else 
        write(6,*) ' flagsp should be between 0 and 2'
        stop
      endif
      write(6,*) ' ipval=',ipval
      write(6,*) ' 0:m, 1:vx, 2:vy, 3:vz, 4:p, 5:T, 6:mH, 7:Z'
      if(flaglog.ne.0) then
        write(6,*) ' contour map in log'
      endif
      write(6,*) ' image max and min range =',fmaxv,fminv
      if(flagv.ne.0) then
        write(6,*) ' with velocity arrows'
      endif
      write(6,*) ' label=',wlab

! set view point
      rxvp(1) = 0.1
      rxvp(2) = 0.9
      dxvp = (rxvp(2)-rxvp(1))/real(nxp)
      do i=1,nxp
        xlvp(i)=rxvp(1)+real(i-1)*dxvp
        xhvp(i)=rxvp(1)+real(i)*dxvp
      enddo

      ryvp(1) = 0.1
      ryvp(2) = 0.9
      dyvp = (ryvp(2)-ryvp(1))/real(nyp)
      do i=1,nyp
        ylvp(i)=ryvp(2)-real(i)*dyvp
        yhvp(i)=ryvp(2)-real(i-1)*dyvp
      enddo

      ip=0
      iyp=1
      do ixp=1,nxp
        ip=ip+1
!                                     123456789012345678901234
        if(flagsp.eq.0) then
          write(filen,'(a16,i6.6,a4)') '../output/slmapg',step,'.dat'
        else
          write(filen,'(a15,i6.6,a4)') '../output/slmap',step,'.dat'
        endif
        write(6,*) ' reading ',filen
! ****    set filename for data file   *****
        open(51,file=filen,status='old',form='unformatted')
        read(51) nx,ny
        read(51) flag3d,flagrot
        read(51) xld,xhd
        read(51) yld,yhd
        read(51) zld,zhd
        read(51) nval
        read(51) zpd
        read(51) tpd   
        read(51) hmp
        read(51) flagid
        read(51) ids,ide

        write(6,*) ' nx,ny=',nx,ny
        write(6,*) ' flag3d,flagrot=',flag3d,flagrot
        write(6,*) ' x range =',xld,xhd
        write(6,*) ' y range =',yld,yhd
        write(6,*) ' z range =',zld,zhd
        write(6,*) ' nval=',nval
        write(6,*) ' zpd=',zpd
        write(6,*) ' time=',tpd
        write(6,*) ' hmp=',hmp
        write(6,*) ' flagid,ids,ide=',flagid,ids,ide

        allocate(xm(nx,ny))
        allocate(ym(nx,ny))
        allocate(fdval(nx,ny,nval))
        allocate(fp(nx,ny))

! read the data
        do iv=1,nval
          do j=1,ny
            read(51) (fdval(i,j,iv),i=1,nx)
          enddo
        enddo

! *** Input Parameter ***
        xl = real(xld)
        yl = real(yld)
        xh = real(xhd)
        yh = real(yhd)

! set area
        dx = dble((xh-xl)/real(nx))
        dy = dble((yh-yl)/real(ny))
        ds = dx*dy

! *** set transfer array ***
        tr(2) = (xh-xl)/real(nx)
        tr(1) = xl-0.5*tr(2)
        tr(3) = 0.0
        tr(5) = 0.0      
        tr(6) = (yh-yl)/real(ny)
        tr(4) = yl-0.5*tr(6)

        if(ip.eq.1) then
! *** plot 2d image at ix=nix ***
! *** open device ****
          if (pgopen('?') .ne.1) then
            write(*,*) ' Cannot open device!'
            stop
          endif
          call pgpap(0.,1.0)
          call pgsch(chdef)
          call pgslw(lwdef)

          call pgscr(1,0.0,0.0,0.0)
          call pgscr(0,1.0,1.0,1.0)
          call pgeras

        endif

! *** Set up window and viewport ***
        call pgsvp(xlvp(ixp),xhvp(ixp),ylvp(iyp),yhvp(iyp))
        call pgswin(xl,xh,yl,yh)

        iv=ipval+1
        fmin = 1.0e12
        fmax = -1.0e12
        if(flaglog.ne.0) then
          do j = 1,ny
            do i = 1,nx
              if(fdval(i,j,iv).gt.0.0d0) then
                fp(i,j)=real(dlog10(fdval(i,j,iv)))
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
              fp(i,j)=real(fdval(i,j,iv))
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
            enddo
          enddo
        endif
        write(*,*) ip,' fmax, fmin      = ',fmax,fmin
        fmax=fmaxv
        fmin=fminv

! *** Set up the color map ***
        bright = 0.5
        contra = 1.0
        call palett(ctype,contra,bright)
! *** Draw the map with pgimg ***
        call pgimag(fp(1,1),nx,ny,1,nx,1,ny,fmin,fmax,tr)

! *** dwarf the velocity field ***
        if(flagv.ne.0) then
          nskipa = 5
! work for +-2 kpc
          vdx = real(dx)*3.0
          write(6,*) ' vdx=',vdx
!          call pgsah(2,45.0,0.2)
          call pgsah(1,45.0,0.3)
!          call pgsch(0.2)
          call pgsch(0.5)
          do i=nskipa,nx-1,nskipa
            do j=nskipa,ny-1,nskipa
              xar(1) = xl+real(dx)*(real(i)-0.5)
              yar(1) = yl+real(dy)*(real(j)-0.5)
! 5 pixel / 100 km/s
              xar(2) = xar(1)+real(fdval(i,j,2)*VUKMS/50.0d0)*vdx
              yar(2) = yar(1)+real(fdval(i,j,3)*VUKMS/50.0d0)*vdx
! vrot
              rxyp=sqrt(xar(1)**2+yar(1)**2)
              call pgsci(0)
              call pgslw(lwdef)
              if(rxyp.gt.0.0) then
                vrot=((real(fdval(i,j,2)*yar(1)-real(fdval(i,j,3))*xar(1))) &
                  /rxyp)*VUKMS
                if(vrot.lt.0.0) then
                  write(6,*) ' negative vrot=',vrot,' at ',xar(1),yar(1)
                  call pgsci(1)
                  call pgslw(5)
                endif
              else
                vrot=0.0
              endif
              call pgarro(xar(1),yar(1),xar(2),yar(2))
              call pgsci(1)
            enddo
          enddo
          call pgsch(chdef)
        endif

        call pgsci(0)
        call pgbox('bncts',0.,0,'bcts',0.,0)
        call pgsci(1)
        call pgbox('nts ',0.,0,'',0.,0)
        call pgmtxt('B',2.0,0.5,0.5,'(kpc)')

        if(ixp.eq.1) then
          call pgbox(' ',0.,0,'nts',0.,0)
          call pgmtxt('L',2.0,0.5,0.5,'(kpc)')
        endif

! *** Draw a Wedge ***
!      call pgwedg('RI',0.,4.,fmin,fmax,'N\s\u') 
       call pgwedg('RI',0.,3.,fmin,fmax,wlab) 

        deallocate(xm)
        deallocate(ym)
        deallocate(fdval)
        deallocate(fp)

      enddo
      close(51)

! *** close the device and exit ***
      call pgend

end program 

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
