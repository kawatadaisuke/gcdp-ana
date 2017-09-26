!
!   main.f90 for plot for sval2d
! 03 Feb.  2016      written by D. Kawata
!

program plot2d

      implicit none

      integer pgopen
      integer,parameter :: nxp=2
      integer,parameter :: nyp=1
! *** for work ***
      integer nix,nstep,step,nval
      integer i,j,k,ctype,ip,flagc,flaglog,istep
      integer maxi,maxj,mini,minj
! *** input parameter ***
! *** grid data to read ***
      integer nx,ny,flaginc
      double precision xrange(2),yrange(2),tu,iang
! for plot
      real xl,yl,xh,yh
      real chdef
! *** for mesh data ***
      double precision txp,typ,tzp,tfdp
      double precision,allocatable :: xm(:,:),ym(:,:),denxy(:,:),denxz(:,:)
      double precision,allocatable :: fdval(:,:)
      double precision drtmp(4)
      real,allocatable :: fp(:,:)
! for viewpoint
      real rxvp(2),ryvp(2),dxvp,dyvp
      real xlvp(nxp),xhvp(nxp),ylvp(nyp),yhvp(nyp)
! *** for plot ***
      real fmax,fmin,bright,contra,fcmin,fcmax
      real fmaxv(2),fminv(2)
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
      read(50,*) flagc
      read(50,*) flaglog
      read(50,*) fminv(1),fmaxv(1)
      read(50,*) fminv(2),fmaxv(2)
      read(50,*) flagcon
      close(50)
      nval=4

      if(flagc.eq.0) then      
        write(6,*) ' plot for gas'
      else if(flagc.eq.1) then
        write(6,*) ' plot for DM'
      else
        write(6,*) ' plot for star'
      endif
      if(flaglog.eq.0) then
        write(6,*) ' plot in log'
      endif
      write(6,*) ' image density range (x-y) =',fmaxv(1),fminv(1)
      write(6,*) ' image density range (x-z) =',fmaxv(2),fminv(2)
      if(flagcon.ne.0) then
        write(6,*) ' draw contours'
      endif
   
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
          write(6,*) ' inclination angle=',iang
        endif
! *** print output file ***
        allocate(xm(nx,ny))
        allocate(ym(nx,ny))
        allocate(denxy(nx,ny))
        allocate(denxz(nx,ny))
        allocate(fdval(nx,ny))
        allocate(fp(nx,ny))

        do j=1,ny
          do i=1,nx
            read(51) xm(i,j),ym(i,j),denxy(i,j),denxz(i,j),drtmp(1),drtmp(2) &
             ,drtmp(3),drtmp(4)
          enddo
        enddo


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
        if(flagc.eq.0) then
          write(filen,'(a14,i6.6,a8)') 'output/pgden2d',step,'.gif/gif'
        else if(flagc.eq.1) then
          write(filen,'(a14,i6.6,a8)') 'output/pdden2d',step,'.gif/gif'
        else
          write(filen,'(a14,i6.6,a8)') 'output/psden2d',step,'.gif/gif'
        endif

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
        write(mlab,'(a3,f6.3,a5)') ' t=',tu,'(Gyr)'
        call pgptxt(0.5,0.5,0.,0.5,mlab)



!        open(60,file='test.dat',status='unknown') 
        do ip=1,2
! *** Set up window and viewport ***
          call pgsvp(xlvp(ip),xhvp(ip),ylvp(1),yhvp(1))
          call pgswin(xl,xh,yl,yh)

! set fdval
          if(ip.eq.1) then
            do j=1,ny
              do i=1,nx
                fdval(i,j)=denxy(i,j)
              enddo
            enddo
          else
            do j=1,ny
              do i=1,nx
                fdval(i,j)=denxz(i,j)
              enddo
            enddo
          endif

!          if(ip.eq.1) then
!            do j=1,ny
!              do i=1,nx
!                write(60,'(3(1pE13.5))') xm(i,j),ym(i,j),fdval(i,j)
!              enddo
!            enddo
!          endif
   
! ***** set plotting file *****
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
          write(6,*) ' colour image density range=',fmaxv(ip),fminv(ip)
          fmax=fmaxv(ip)
          fmin=fminv(ip)

          if(ip.eq.1) then
            call pgbox('bncts',0.,0,'bncts',0.,0)
            call pgmtxt('B',2.0,0.5,0.5,'x-y (kpc)')
          else
            call pgbox('bncts',0.,0,'bcts',0.,0)
            call pgmtxt('B',2.0,0.5,0.5,'x-z (kpc)')
          endif
! *** Set up the color map ***
          bright = 0.5
          contra = 1.0
          call palett(ctype,contra,bright)
! *** Draw the map with pgimg ***
          call pgimag(fp(1,1),nx,ny,1,nx,1,ny,fmin,fmax,tr)

! density concour
          if(flagcon.ne.0) then
            ncon=10
            dlev=(fmax-fmin)/real(ncon)
            do i=1,ncon+1
              alev=fmin+dlev*real(i-1)
              call pgcont(fp(1,1),nx,ny,1,nx,1,nx,alev,-1,tr)
            enddo
          endif

! *** Draw a Wedge ***
!      call pgwedg('RI',0.,4.,fmin,fmax,'N\s\u') 

        enddo

        deallocate(xm)
        deallocate(ym)
        deallocate(fdval)
        deallocate(denxy)
        deallocate(denxz)
        deallocate(fp)

! *** close the device and exit ***
        call pgend

      enddo
      close(51)


      stop
end program plot2d

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
