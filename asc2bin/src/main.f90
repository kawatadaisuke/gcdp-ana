! /***************************************
!   asc2bin
!  30 Jan. 2018  written by D.Kawata
! ****************************************/

program asc2bin
      use gcdp_const
      use gcdp_system
      use gcdp_baryon
      use gcdp_dm
   
      implicit none

      integer i,is,nval
      integer npt,ngt,nst,ndmt,ndm1t
      integer ng,ndm,ndm1,ns,np
      integer nstep,step
      integer ierr
      integer nivalb,ndbval,ndbhyd,ndbmet,ndbsf
! *** Time, redshift, back ground density ***
      double precision tu,zu,ai
! ***  for filename ***
      character fileg*60,filedm*60,files*60,fileo*60

! *****   Open Input File   ******
      open(50,file='./ini/input.dat',status='old')
      read(50,*) ngt
      read(50,'(a60)') fileg
      read(50,*) ndmt
      read(50,'(a60)') filedm
      read(50,*) nst
      read(50,'(a60)') files

      read(50,*) tu
      read(50,*) step
      close(50)

      tu=2.62
      ai=1.0
      nivalb=3
      ndbval=9
      ndbhyd=0
      ndbmet=0
      ndbsf=1

      write(6,*) ' time=',tu
      write(6,*) ' step name=',step
      write(6,*) ' ng,ndm,ns=',ngt,ndmt,nst
      if(ngt.gt.0) then
        write(6,*) ' read gas data =',fileg
      endif
      if(nst.gt.0) then
        write(6,*) ' read star data =',files
      endif

! *** this version only for gas and star particles ***
      ndmt=0

      ng=ngt
      ndm=ndmt
      ns=nst
      npt=ngt+nst
      np=npt

! allocate
      call allocate_baryon(np) 

! *** read the data ***
      np=0
      if(ngt.gt.0) then
        open(50,file=fileg,status='old')
        ! if there is a header
        do i=0,ngt-1
          read(50,*) x_p(np),y_p(np),z_p(np),vx_p(np),vy_p(np),vz_p(np) &
            ,mzZ_p(np),m_p(np)
          ! convert unit
          x_p(np)=x_p(np)/LUKPC
          y_p(np)=y_p(np)/LUKPC
          z_p(np)=z_p(np)/LUKPC
          vx_p(np)=vx_p(np)/VUKMS
          vy_p(np)=vy_p(np)/VUKMS
          vz_p(np)=vz_p(np)/VUKMS
          mzZ_p(np)=mzZ_p(np)*m_p(np)
          ts_p(np)=0.0d0
          m_p(np)=m_p(np)/MUSM
          id_p(np)=np
          list_ap(np)=np
          flagc_p(np)=-1
          rho_p(np)=1.0d0
          u_p(np)=1.0d0
          np=np+1
        enddo
        close(50)
      endif
      if(nst.gt.0) then
        open(50,file=files,status='old')
        do i=0,ns-1
          read(50,*) x_p(np),y_p(np),z_p(np),vx_p(np),vy_p(np),vz_p(np) &
            ,mzZ_p(np),ts_p(np),m_p(np)
          ! convert unit
          x_p(np)=x_p(np)/LUKPC
          y_p(np)=y_p(np)/LUKPC
          z_p(np)=z_p(np)/LUKPC
          vx_p(np)=vx_p(np)/VUKMS
          vy_p(np)=vy_p(np)/VUKMS
          vz_p(np)=vz_p(np)/VUKMS
          mzZ_p(np)=mzZ_p(np)*m_p(np)
          ts_p(np)=(tu-ts_p(np))/TMUGYR
          m_p(np)=m_p(np)/MUSM
          id_p(np)=np
          list_ap(np)=np
          flagc_p(np)=1
          rho_p(np)=1.0d0
          u_p(np)=1.0d0
          np=np+1
        enddo
        close(50)
      endif

      if(np.ne.npt) then
        write(6,*) ' Error in reading baryon data. npt and np are inconsistent.'
        write(6,*) ' npt, np=',npt,np
        stop
      endif

      ! convert tu
      tu=tu/TMUGYR

! output the data
      myrank=0
      write(fileo,'(a15,i6.6,a1,i4.4)') &
        './output/bbvals',step,'n',myrank

      open(61,file=fileo,status='unknown',form='unformatted')
      write(61) npt,ndmt,ndmt,ai,tu
      write(61) 1,1,nivalb,ndbval
      write(61) ng,ndm,ns,ndm,ng,ndm,ns,np
      write(61) id_p
      write(61) flagc_p
      write(61) list_ap
      write(61) x_p
      write(61) y_p
      write(61) z_p
      write(61) vx_p
      write(61) vy_p
      write(61) vz_p
      write(61) m_p
      write(61) rho_p
      write(61) u_p
      close(61)

      write(fileo,'(a15,i6.6,a1,i4.4)') &
        './output/bbhyds',step,'n',myrank
      open(61,file=fileo,status='unknown',form='unformatted')
      write(61) npt,ndmt,ndmt,ai,tu
      write(61) 1,1,ndbhyd
      close(61)       
      write(fileo,'(a15,i6.6,a1,i4.4)') &
        './output/bbmets',step,'n',myrank
      open(61,file=fileo,status='unknown',form='unformatted')
      write(61) npt,ndmt,ndmt,ai,tu
      write(61) 1,1,ndbmet
      close(61)       

      write(fileo,'(a15,i6.6,a1,i4.4)') &
        './output/bbsfis',step,'n',myrank
      open(61,file=fileo,status='unknown',form='unformatted')
      write(61) npt,ndmt,ndmt,ai,tu
      write(61) 1,1,ndbsf
      write(61) ts_p
      close(61)       

end program
