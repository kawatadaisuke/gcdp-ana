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
      if(nst.gt.0) then
        write(6,*) ' read star data =',files
      endif

! *** this version only for star particles ***
      ngt=0
      ndmt=0

      ng=ngt
      ndm=ndmt
      ns=nst
      npt=ngt+nst
      np=npt

! allocate
      call allocate_baryon(np) 

! *** read the data ***
      if(nst.gt.0) then
        open(50,file=files,status='old')
        ! if there is a header
        read(50,*)
        do i=0,ns-1
          read(50,*) x_p(i),y_p(i),z_p(i),vx_p(i),vy_p(i),vz_p(i) &
            ,mzZ_p(i),ts_p(i),m_p(i)
          ! convert unit
          x_p(i)=x_p(i)/LUKPC
          y_p(i)=y_p(i)/LUKPC
          z_p(i)=z_p(i)/LUKPC
          vx_p(i)=vx_p(i)/VUKMS
          vy_p(i)=vy_p(i)/VUKMS
          vz_p(i)=vz_p(i)/VUKMS
          mzZ_p(i)=mzZ_p(i)*m_p(i)
          ts_p(i)=(tu-ts_p(i))/TMUGYR
          m_p(i)=m_p(i)/MUSM
          id_p(i)=i
          list_ap(i)=i
          flagc_p(i)=1
          rho_p(i)=1.0d0
          u_p(i)=1.0d0
        enddo
        close(50)
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
