
!   main.f90 for plot for 2D distribution of Vc from ../epotsqgrid
! 01 Aug.  2017      written by D. Kawata
!

program plot2dvc

      implicit none

      integer pgopen
      integer,parameter :: nxp=1
      integer,parameter :: nyp=1
      double precision,parameter :: VUKMS=207.4
      double precision,parameter :: LUKPC=100.0
! *** for work ***
      integer nix,nstep,step,nval
      integer i,j,k,ctype,ip,istep
      integer maxi,maxj,mini,minj
! for vc
      double precision rm,frm
! *** input parameter ***
! *** grid data to read ***
      integer nx,ny,flaginc,nz
      double precision xrange(2),yrange(2),tu,iang,zrange(2)
! for plot
      real xl,yl,xh,yh
      real chdef
! *** for mesh data ***
      double precision txp,typ,tzp,tfdp
      double precision,allocatable :: xm(:,:),ym(:,:),zm(:,:)
      double precision,allocatable :: epotm(:,:)
      double precision,allocatable :: dvxm(:,:),dvym(:,:),dvzm(:,:)
      double precision,allocatable :: vcm(:,:)
      double precision drtmp(4)
      real,allocatable :: fp(:,:)
! for viewpoint
      real rxvp(2),ryvp(2),dxvp,dyvp
      real xlvp(nxp),xhvp(nxp),ylvp(nyp),yhvp(nyp)
! *** for plot ***
      real fmax,fmin,bright,contra,fcmin,fcmax
      real fmaxv,fminv
      real fmaxc,fminc
      real tr(6)      
! for contour
      integer ncon,flagcon
      real dlev,alev
! ***  for filename ***
      character filen*60,mlab*60

      chdef=1.5
      
! *** open input file ****
      open(50,file='ini/inputp.dat',status='old')
      read(50,*) nstep,step
      read(50,*) fminv,fmaxv
      close(50)
      nval=4

      write(6,*) ' image Vc range =',fmaxv,fminv
   
 777  write(*,*) ' plot type, 1: gray scale, 2: rainbow '
      read(5,*) ctype
      if (ctype.lt.1.or.ctype.gt.5) then
        goto 777
      endif        

      if(nstep.gt.1) then
        open(50,file='ini/file.dat',status='old') 
      endif
      do istep=1,nstep
        if(nstep.gt.1) then
          read(50,*) step
        endif
        write(filen,'(a17,i6.6,a4)') 'output/epotsqgrid',step,'.dat'
        write(6,*) ' reading ',filen
! ****    set filename for data file   *****
        open(51,file=filen,status='old',form='unformatted')
        read(51) nx,ny,nz
        read(51) xrange(1),xrange(2)
        read(51) yrange(1),yrange(2)
        read(51) zrange(1),zrange(2)


        write(6,*) ' nx,ny,nz=',nx,ny,nz
        write(6,*) ' x range =',xrange(1),xrange(2)
        write(6,*) ' y range =',yrange(1),yrange(2)
! *** print output file ***
        allocate(xm(nx,ny))
        allocate(ym(nx,ny))
        allocate(zm(nx,ny))
        allocate(epotm(nx,ny))       
        allocate(dvxm(nx,ny))
        allocate(dvym(nx,ny))
        allocate(dvzm(nx,ny))
        allocate(vcm(nx,ny))
        allocate(fp(nx,ny))

        open(60,file='vcmap.dat',status='unknown')
        do j=1,ny
          do i=1,nx
            read(51) xm(i,j),ym(i,j),zm(i,j),epotm(i,j),dvxm(i,j) &
             ,dvym(i,j),dvzm(i,j)
            rm=dsqrt(xm(i,j)**2+ym(i,j)**2)
            if(rm.gt.0.0011d0) then
              frm=(dvxm(i,j)*xm(i,j)+dvym(i,j)*ym(i,j))/rm
            else
              frm=0.0d0
            endif
            vcm(i,j)=dsqrt(-frm*rm)*VUKMS
            write(60,'(5(1pE13.5))') xm(i,j),ym(i,j),vcm(i,j),rm,frm
          enddo
        enddo
        close(60)

! *** Input Parameter ***
        xl = real(xrange(1)*LUKPC)
        yl = real(yrange(1)*LUKPC)
        xh = real(xrange(2)*LUKPC)
        yh = real(yrange(2)*LUKPC)

        write(*,*) ' xl,yh,yl,yh     = ',xl,xh,yl,yh
! *** set transfer array ***
        tr(2) = (xh-xl)/real(nx)
        tr(1) = xl-0.5*tr(2)
        tr(3) = 0.0
        tr(5) = 0.0      
        tr(6) = (yh-yl)/real(ny)
        tr(4) = yl-0.5*tr(6)

! set view point
        rxvp(1) = 0.1
        rxvp(2) = 0.85
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

! set output file name
        write(filen,'(a12,i6.6,a8)') 'output/vcmap',step,'.gif/gif'

! *** plot 2d image at ix=nix ***
! *** open device ****
        if (pgopen(filen) .ne.1) then
          write(*,*) ' Cannot open device!'
          stop
        endif
        call pgpap(0.,(dxvp/dyvp)*(yh-yl)/(xh-xl))
        call pgsch(chdef)

! *** write header ***
        call pgsvp(0.0,1.0,ryvp(2),1.0)
        call pgswin(0.0,1.0,0.0,1.0)
        write(6,*) ' step, t=',step,tu
!        write(mlab,'(a3,f6.3,a5)') ' t=',tu,'(Gyr)'
!        call pgptxt(0.5,0.5,0.,0.5,mlab)

!        open(60,file='test.dat',status='unknown') 
        do ip=1,2
! *** Set up window and viewport ***
          call pgsvp(xlvp(ip),xhvp(ip),ylvp(1),yhvp(1))
          call pgswin(xl,xh,yl,yh)
   
! ***** set plotting file *****
          fmin = 1.0e12
          fmax = -1.0e12
          do j = 1,ny
            do i = 1,nx
              fp(i,j)=real(vcm(i,j))
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
          write(*,*) ip,' fmax, fmin      = ',fmax,fmin
          write(6,*) ' colour image density range=',fmaxv,fminv
          fmax=fmaxv
          fmin=fminv

          call pgbox('bncts',0.,0,'bncts',0.,0)
          call pgmtxt('B',2.0,0.5,0.5,'x (kpc)')
          call pgmtxt('L',2.0,0.5,0.5,'y (kpc)')

! *** Set up the color map ***
          bright = 0.5
          contra = 1.0
          call palett(ctype,contra,bright)
! *** Draw the map with pgimg ***
          call pgimag(fp(1,1),nx,ny,1,nx,1,ny,fmin,fmax,tr)

! *** Draw a Wedge ***
          call pgwedg('RI',0.,4.,fmin,fmax,'Vcirc') 

        enddo

        deallocate(xm)
        deallocate(ym)
        deallocate(zm)
        deallocate(epotm)
        deallocate(dvxm)
        deallocate(dvym)
        deallocate(dvzm)
        deallocate(vcm)
        deallocate(fp)

! *** close the device and exit ***
        call pgend

      enddo
      close(51)


      stop
end program plot2dvc

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
