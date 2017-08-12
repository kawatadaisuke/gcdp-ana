!
!   main.f90 for plot_rvrot
! 12 Aug. 2017      written by D. Kawata
!

program plot_rvrot

      implicit none

      integer pgopen
      integer,parameter :: nxp=1
      integer,parameter :: nyp=1
! *** for work ***
      integer i,j,k,ctype,ip,ixp,iyp,flaglog
      integer maxi,maxj,mini,minj
      integer flagplstar
! *** input parameter ***
! *** grid data to read ***
      integer nx,ny
      double precision xld,xhd,yld,yhd
! for plot
      integer lwdef
      real xl,yl,xh,yh
      real chdef
! *** for mesh data ***
      double precision txp,typ,tzp,tfdp
      double precision,allocatable :: xm(:,:),ym(:,:)
      double precision,allocatable :: fdval(:,:)
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
      character filen*60,mlab*60,rchr*5

      chdef=2.5
      lwdef=3
!  plot type, 1: gray scale, 2: rainbow
      ctype=2
      
! *** open input file ****
      open(50,file='ini/input-plotrvrot.dat',status='old')
      read(50,*) flaglog
      read(50,*) fminv,fmaxv
      read(50,*) flagplstar
      close(50)

      if(flaglog.eq.0) then
        write(6,*) ' contour map in log'
      endif
      write(6,*) ' image density range =',fmaxv,fminv
      if(flagplstar.ne.0) then
        write(6,*) ' overplot stars from output/starR_Vrot.asc'
      endif

! set view point
      rxvp(1) = 0.15
      rxvp(2) = 0.95
      dxvp = (rxvp(2)-rxvp(1))/real(nxp)
      do i=1,nxp
        xlvp(i)=rxvp(1)+real(i-1)*dxvp
        xhvp(i)=rxvp(1)+real(i)*dxvp
      enddo

      ryvp(1) = 0.15
      ryvp(2) = 0.975
      dyvp = (ryvp(2)-ryvp(1))/real(nyp)
      do i=1,nyp
        ylvp(i)=ryvp(2)-real(i)*dyvp
        yhvp(i)=ryvp(2)-real(i-1)*dyvp
      enddo

      ip=0
      iyp=1
      do ixp=1,nxp
        ip=ip+1
!                             123456789012345678901234
        write(filen,'(a24)') '../output/smapR_Vrot.dat'
        write(6,*) ' reading ',filen
! ****    set filename for data file   *****
        open(51,file=filen,status='old',form='unformatted')
        read(51) nx,ny
        read(51) xld,xhd
        read(51) yld,yhd
        write(6,*) ' nx,ny=',nx,ny
        write(6,*) ' x range =',xld,xhd
        write(6,*) ' y range =',yld,yhd

        allocate(xm(nx,ny))
        allocate(ym(nx,ny))
        allocate(fdval(nx,ny))
        allocate(fp(nx,ny))

        do j=1,ny
          do i=1,nx
            read(51) xm(i,j),ym(i,j),fdval(i,j)
          enddo
        enddo
        close(51)

        if(flagplstar.ne.0) then
          open(50,file='../output/starR_Vrot.asc',status='old')
          read(50,'(a5,I10)') rchr,nsp
          write(6,*) ' N star data=',nsp
          read(50,*) 
       
          allocate(xsp(0:nsp-1))
          allocate(ysp(0:nsp-1))
          allocate(zsp(0:nsp-1))
          allocate(rxysp(0:nsp-1))
          allocate(vrotsp(0:nsp-1))
          allocate(agesp(0:nsp-1))
  
          do i=0,nsp-1
            read(50,'(6(1pE13.5))') xsp(i),ysp(i),zsp(i),rxysp(i) &
              ,vrotsp(i),agesp(i)
          enddo
          close(50)
        endif

! *** Input Parameter ***
        xl = real(xld)
        yl = real(yld)
        xh = real(xhd)
        yh = real(yhd)

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
          call pgpap(0.,0.7)
          call pgsch(chdef)
          call pgslw(lwdef)

          call pgscr(1,0.0,0.0,0.0)
          call pgscr(0,1.0,1.0,1.0)
          call pgeras

        endif

! *** Set up window and viewport ***
        call pgsvp(xlvp(ixp),xhvp(ixp),ylvp(iyp),yhvp(iyp))
        call pgswin(xl,xh,yl,yh)

        fmin = 1.0e12
        fmax = -1.0e12
        if(flaglog.eq.0) then
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

        if(flagplstar.ne.0) then
          call pgsci(0)
          call pgsch(0.25)
          call pgpt(nsp,rxysp,vrotsp,17)
          call pgsch(chdef)
          call pgsci(1)
        endif

        call pgsci(0)
        call pgbox('bncts',0.,0,'bcts',0.,0)
        call pgsci(1)
        call pgbox('nts ',0.,0,'',0.,0)
        call pgmtxt('B',2.0,0.5,0.5,'R (kpc)')

        if(ixp.eq.1) then
          call pgbox(' ',0.,0,'nts',0.,0)
          call pgmtxt('L',2.0,0.5,0.5,'V\drot\u')
        endif

! *** Draw a Wedge ***
!      call pgwedg('RI',0.,4.,fmin,fmax,'N\s\u') 

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
