! /***************************************
!   lbsels
! ****************************************/

program lbsels
      use gcdp_const
      use gcdp_system
      use gcdp_baryon
      use gcdp_dm
   
      implicit none
      include 'mpif.h'

      integer id,i,j,k,pn,flag,ns,ng,np,istep
      integer npt,ngt,nst,ndmt,ndm1t,ndm,ndm1
      integer step,is,ic
      integer nval
      integer flagc,flagcom,flagr,flagtargetf,flagdup
      integer idr(0:1)
      integer ierr
      character chread*9
! *** Time, redshift, back ground density ***
      double precision tu,zu,ai
! for output
      integer nskip,nsel,nselout,nselfin
! solar position and velocity
      double precision rzasun(0:2),xsun(0:1),vsun(0:2)
! l, b and distance range and age range 
      double precision glonran(0:1),glatran(0:1),dran(0:1),ageran(0:1)
      double precision zran(0:1)
! particle data
      integer nagep,nsearch
      integer,allocatable :: listagep(:),flagtaken(:),listsel(:) &
        ,listselfin(:)
      double precision modp,vradxyp
      double precision,allocatable :: glonp(:),glatp(:),dxyp(:) &
       ,d3dp(:),vglonp(:),vglatp(:),vlosp(:),agep(:),vradgalp(:),vrotgalp(:) &
       ,rxygalp(:),phigalp(:),fehp(:)
      double precision vrmselp,vtmselp,vzmselp,vr2mselp,vt2mselp,vz2mselp &
       ,mtotselp
      double precision vrsig,vtsig,vzsig
! for random pick up
      integer jsel
      double precision rsel,xsel,ysel,dminsel,thsel,dsel
! target star data from axsymdiskm-fit_sels.asc
! made with /Users/dkawata/work/obs/projs/Cepheids-kinematics/py/axsymdiskm-fit.py
      integer ntargs,pnts
      double precision degrange,dang,danglim,distlimmax,distlimmin
      double precision dispts,disptsmin,distxyts
      double precision,allocatable :: glons(:),glats(:),distxys(:),hrvs(:) &
       ,vlons(:),errhrvs(:),errvlons(:),mods(:),errmods(:) &
       ,dists(:),distmins(:),distmaxs(:),xts(:),yts(:),zts(:) &
       ,ras(:),decs(:),pmras(:),pmdecs(:),errpmras(:),errpmdecs(:) &
       ,pmradec_corrs(:),logps(:),loggs(:),teffs(:),vmags(:),vicols(:),avs(:)
      character(len=20),allocatable :: names(:)
! to sample
      integer ncanp
      integer,allocatable :: listcanp(:)
! for work
      integer,allocatable :: tivr(:),tivs(:)
      double precision,allocatable :: tdvr(:),tdvs(:)
      double precision,allocatable :: tx(:),ty(:)
! external
      integer idum
      real ran1
      external ran1

! ***** MPI Initialization *****
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      if(myrank.eq.0) then
        write(6,*) ' lbsels'
        print *, "Process ",myrank, " of ", nprocs, " is alive"
      endif
      if(nprocs.gt.NCPU) then
        if(myrank.eq.0) then 
          write(6,*) ' Error: increase NCPU=',NCPU
        endif
        call MPI_FINALIZE()
        stop
      endif

      idum=-194020

! *****   Open Input File   ******
      if(myrank.eq.0) then
        open(50,file='./ini/input.dat',status='old')
        read(50,*) step
        read(50,*) nskip,nselout
! solar position R, z, angle for Position of Sun (degree) 
!   from (x,y) = (-rzasun(0),0) anti-clockwise ***/
        read(50,*) rzasun(0),rzasun(1),rzasun(2)
! *** Sun's radial (outward+), rotation, vertical velocity (clockwise) ***/
        read(50,*) vsun(0),vsun(1),vsun(2)
        read(50,*) glonran(0),glonran(1)
        read(50,*) glatran(0),glatran(1)
        read(50,*) dran(0),dran(1)
        read(50,*) zran(0),zran(1)
        read(50,*) ageran(0),ageran(1)
        read(50,*) flagtargetf
        read(50,*) degrange
        read(50,*) flagdup
        close(50)
        write(6,*) ' step=',step 
        write(6,*) ' output nskip for lb*.dat=',nskip
        if(nselout.gt.0) then
          write(6,*) ' only output selected particles np<',nselout
        endif
        write(6,*) ' Solar position angle from (-',rzasun(0),',',rzasun(1) &
          ,')=',rzasun(2)
        write(6,*) ' Solar motion Vrad,V,W=',vsun(0),vsun(1),vsun(2)
        write(6,*) ' selection l range min,max =',glonran(0),glonran(1)
        write(6,*) '           b range min,max =',glatran(0),glatran(1)
        write(6,*) '           z range min,max =',zran(0),zran(1)
        write(6,*) '    distance range min,max =',dran(0),dran(1)
        write(6,*) ' for star age range=',ageran(0),ageran(1)
        if(flagtargetf.ne.0) then
          if(flagtargetf.eq.1) then
            write(6,*) ' read ini/axsymdiskm-fit_sels.asc, target position file'
          else
            write(6,*) ' read ini/sels_rv.asc, target position file'
          endif
          write(6,*) ' search degree range=',degrange
          if(flagdup.ne.0) then
            write(6,*) ' allow duplicates'
          endif
        endif
      endif

      nval=5
      allocate(tivr(0:nval-1))
      if(myrank.eq.0) then
        tivr(0)=step
        tivr(1)=nskip
        tivr(2)=nselout
        tivr(3)=flagtargetf
        tivr(4)=flagdup
      endif
      call MPI_BCAST(tivr,nval,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      step=tivr(0)
      nskip=tivr(1)
      nselout=tivr(2)
      flagtargetf=tivr(3)
      flagdup=tivr(4)
      deallocate(tivr)

      nval=17
      allocate(tdvr(0:nval-1))
      if(myrank.eq.0) then
        tdvr(0)=rzasun(0)
        tdvr(1)=rzasun(1)
        tdvr(2)=rzasun(2)
        tdvr(3)=vsun(0)
        tdvr(4)=vsun(1)
        tdvr(5)=vsun(2)
        tdvr(6)=glonran(0)
        tdvr(7)=glonran(1)
        tdvr(8)=glatran(0)
        tdvr(9)=glatran(1)
        tdvr(10)=dran(0)
        tdvr(11)=dran(1)
        tdvr(12)=ageran(0)
        tdvr(13)=ageran(1)
        tdvr(14)=degrange
        tdvr(15)=zran(0)
        tdvr(16)=zran(1)
      endif
      call MPI_BCAST(tdvr,nval,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      rzasun(0)=tdvr(0)
      rzasun(1)=tdvr(1)
      rzasun(2)=tdvr(2)
      vsun(0)=tdvr(3)
      vsun(1)=tdvr(4)
      vsun(2)=tdvr(5)
      glonran(0)=tdvr(6)
      glonran(1)=tdvr(7)
      glatran(0)=tdvr(8)
      glatran(1)=tdvr(9)
      dran(0)=tdvr(10)
      dran(1)=tdvr(11)
      ageran(0)=tdvr(12)
      ageran(1)=tdvr(13)
      degrange=tdvr(14)
      zran(0)=tdvr(15)
      zran(1)=tdvr(16)
      deallocate(tdvr)

! solar angle deg -> radian
      rzasun(2)=rzasun(2)*M_PI/180.0d0

! read full data
      flagr=1

! *** read the data ***
      call rdata(step,ngt,ng,ndmt,ndm,ndm1t,ndm1,nst,ns,ai,tu,flagr)
      if(myrank.eq.0) then
        write(6,*) ' step, t(GYR)=',step,tu*TMUGYR
        write(6,*) ' ngt,ndmt,ndm1t,nst=',ngt,ndmt,ndm1t,nst
      endif

! transfer to heliocentric position and velocity
      npt=ngt+nst
      np=ng+ns
      if(np.gt.0) then
! for baryon
        allocate(tx(0:np-1))
        allocate(ty(0:np-1))
        allocate(vradgalp(0:np-1))
        allocate(vrotgalp(0:np-1))
        allocate(rxygalp(0:np-1))
        allocate(phigalp(0:np-1))

        do i=0,np-1
! kpc and km/s unit
          x_p(i)=x_p(i)*LUKPC
          y_p(i)=y_p(i)*LUKPC
          z_p(i)=z_p(i)*LUKPC 
          vx_p(i)=vx_p(i)*VUKMS
          vy_p(i)=vy_p(i)*VUKMS
          vz_p(i)=vz_p(i)*VUKMS
! Galactic rotation and radial velocity
          rxygalp(i)=dsqrt(x_p(i)**2+y_p(i)**2)
          vradgalp(i)=(vx_p(i)*x_p(i)+vy_p(i)*y_p(i))/rxygalp(i)
          vrotgalp(i)=(vx_p(i)*y_p(i)-vy_p(i)*x_p(i))/rxygalp(i)
! rotation
          tx(i)=x_p(i)*dcos(-rzasun(2))-y_p(i)*dsin(-rzasun(2))
          ty(i)=x_p(i)*dsin(-rzasun(2))+y_p(i)*dcos(-rzasun(2))     
! add R and z
          x_p(i)=tx(i)+rzasun(0)
          y_p(i)=ty(i)
          z_p(i)=z_p(i)-rzasun(1)
! velocity 
          tx(i)=vx_p(i)*dcos(-rzasun(2))-vy_p(i)*dsin(-rzasun(2))
          ty(i)=vx_p(i)*dsin(-rzasun(2))+vy_p(i)*dcos(-rzasun(2))     
! add solar velocity
          vx_p(i)=tx(i)+vsun(0)
          vy_p(i)=ty(i)-vsun(1)
          vz_p(i)=vz_p(i)-vsun(2)
        enddo      
        deallocate(tx)
        deallocate(ty)

! l-b coordinate
        allocate(d3dp(0:np-1))
        allocate(dxyp(0:np-1))
        allocate(glonp(0:np-1))
        allocate(glatp(0:np-1))
        allocate(vglonp(0:np-1))
        allocate(vglatp(0:np-1))
        allocate(vlosp(0:np-1))
        allocate(agep(0:np-1))
        allocate(fehp(0:np-1))

        do i=0,np-1
          d3dp(i)=dsqrt(x_p(i)**2+y_p(i)**2+z_p(i)**2)
          dxyp(i)=dsqrt(x_p(i)**2+y_p(i)**2)
          glatp(i)=dasin(z_p(i)/d3dp(i))
          glonp(i)=dacos(x_p(i)/dxyp(i))
! angle between star and sun from the center
          phigalp(i)=M_PI-dacos((x_p(i)-rzasun(0))/rxygalp(i))
          if(y_p(i).lt.0.0d0) then
            phigalp(i)=-phigalp(i)
          endif
          if(y_p(i).lt.0.0d0) then
            glonp(i)=2.0d0*M_PI-glonp(i)
          endif
! velocity        
          vlosp(i)=(vx_p(i)*x_p(i)+vy_p(i)*y_p(i)+vz_p(i)*z_p(i))/d3dp(i)
          vglonp(i)=(-vx_p(i)*y_p(i)+vy_p(i)*x_p(i))/dxyp(i)
          vradxyp=(vx_p(i)*x_p(i)+vy_p(i)*y_p(i))/dxyp(i)
          vglatp(i)=(-vradxyp*z_p(i)+vz_p(i)*dxyp(i))/d3dp(i)
! rad to deg
          glonp(i)=glonp(i)*180.0d0/M_PI
          glatp(i)=glatp(i)*180.0d0/M_PI
! age and [Fe/H](=0.0 for now)  for star
          agep(i)=(tu-ts_p(i))*TMUGYR
          fehp(i)=0.0d0
        enddo
! if the glon range is setted with negative values
        if(glonran(0).lt.0.0d0) then
          do i=0,np-1
            if(glonp(i).gt.180.0d0) then
              glonp(i)=glonp(i)-360.0d0
            endif
          enddo
        endif

! output *** not parallelised ***
! for gas
        if(ng.gt.0) then
          open(60,file='output/lbg.dat',status='unknown')
          do i=0,ng-1,nskip
            pn=list_ap(i)
            write(60,160) x_p(pn),y_p(pn),z_p(pn),vx_p(pn),vy_p(pn),vz_p(pn) &
              ,glonp(pn),glatp(pn),d3dp(pn),vglonp(pn),vlosp(pn) &
              ,rxygalp(pn),vradgalp(pn),vrotgalp(pn),phigalp(pn)
 160       format(15(1pE13.5))
          enddo
          close(60)

          allocate(listsel(0:ng-1))
          allocate(listselfin(0:ng-1))
          nsel=0
          ! selecting the particle in the region
          do i=0,ng-1
            pn=list_ap(i)
            if(glonp(pn).ge.glonran(0).and.glonp(pn).le.glonran(1)) then
            if(glatp(pn).ge.glatran(0).and.glatp(pn).le.glatran(1)) then
            if(d3dp(pn).ge.dran(0).and.d3dp(pn).le.dran(1)) then
              listsel(nsel)=pn
              nsel=nsel+1
            endif
            endif
            endif    
          enddo
          write(6,*) ' lbselg nsel=',nsel
          if(nselout.eq.0.or.nselout.gt.nsel) then
            nselfin=nsel
            do i=0,nsel-1
              pn=listsel(i)
              listselfin(i)=pn
            enddo
          else
            nselfin=nselout
            allocate(flagtaken(0:nsel-1))
            do i=0,nsel-1
              flagtaken(i)=0
            enddo
            ! open(62,file='selpos.asc',status='unknown')
            do i=0,nselout-1
              rsel=dsqrt(dble(ran1(idum)))*dran(1)
              thsel=2.0*M_PI*dble(ran1(idum))
              xsel=rsel*dcos(thsel)
              ysel=rsel*dsin(thsel)
              ! write(62,'(4(1pE13.5))') xsel,ysel,rsel,thsel
              dminsel=INF
              jsel=-1
              do j=0,nsel-1
                pn=listsel(j)
                if(flagtaken(j).eq.0) then
                  if(dabs(x_p(pn)-xsel).lt.dminsel &
                    .and.dabs(y_p(pn)-ysel).lt.dminsel) then
                    dsel=dsqrt((x_p(pn)-xsel)**2+(y_p(pn)-ysel)**2)
                    if(dsel.lt.dminsel) then
                      dminsel=dsel
                      jsel=j
                    endif
                  endif
                endif
              enddo 
              if(jsel.lt.0) then
                write(6,*) ' Error: no gas particle found for random selection'
                write(6,*) i,' around x,y,dminsel=',xsel,ysel,dminsel
                stop
              endif
              flagtaken(jsel)=1
              pn=listsel(jsel)
              listselfin(i)=pn
            enddo
            deallocate(flagtaken)
            ! close(62)
          endif
          !output
          open(60,file='output/lbselg.dat',status='unknown')
          do i=0,nselfin-1
            pn=listselfin(i)
            write(60,160) x_p(pn),y_p(pn),z_p(pn),vx_p(pn),vy_p(pn) & 
             ,vz_p(pn),glonp(pn),glatp(pn),d3dp(pn),vglonp(pn),vlosp(pn) &
             ,rxygalp(pn),vradgalp(pn),vrotgalp(pn),phigalp(pn)
          enddo
          close(60)
          deallocate(listsel)  
          deallocate(listselfin)  
        endif
! for star
        if(ns.gt.0) then
          open(60,file='output/lbs.dat',status='unknown')
          do i=ng,np-1,nskip
            pn=list_ap(i)
            write(60,161) x_p(pn),y_p(pn),z_p(pn),vx_p(pn),vy_p(pn),vz_p(pn) & 
              ,glonp(pn),glatp(pn),d3dp(pn),vglonp(pn),vlosp(pn),agep(pn) &
              ,rxygalp(pn),vradgalp(pn),vrotgalp(pn),phigalp(pn)
 161       format(16(1pE13.5))
          enddo
          close(60)

          allocate(listsel(0:ns-1))
          allocate(listselfin(0:ns-1))
          nsel=0
          ! selecting the particle in the region
          do i=ng,np-1
            pn=list_ap(i)
            if(glonp(pn).ge.glonran(0).and.glonp(pn).le.glonran(1)) then
            if(glatp(pn).ge.glatran(0).and.glatp(pn).le.glatran(1)) then
            if(d3dp(pn).ge.dran(0).and.d3dp(pn).le.dran(1)) then
            if(z_p(pn).ge.zran(0).and.z_p(pn).le.zran(1)) then
            if(agep(pn).ge.ageran(0).and.agep(pn).le.ageran(1)) then
              listsel(nsel)=pn
              nsel=nsel+1
            endif
            endif
            endif
            endif
            endif    
          enddo
          write(6,*) ' lbsels nsel=',nsel
          if(nselout.eq.0.or.nselout.gt.nsel) then
            nselfin=nsel
            do i=0,nsel-1
              pn=listsel(i)
              listselfin(i)=pn
            enddo
          else
            nselfin=nselout
            allocate(flagtaken(0:nsel-1))
            do i=0,nsel-1
              flagtaken(i)=0
            enddo
            ! open(62,file='selpos.asc',status='unknown')
            do i=0,nselout-1
              rsel=dsqrt(dble(ran1(idum)))*dran(1)
              thsel=2.0*M_PI*dble(ran1(idum))
              xsel=rsel*dcos(thsel)
              ysel=rsel*dsin(thsel)
              ! write(62,'(4(1pE13.5))') xsel,ysel,rsel,thsel
              dminsel=INF
              jsel=-1
              do j=0,nsel-1
                pn=listsel(j)
                if(flagtaken(j).eq.0) then
                  if(dabs(x_p(pn)-xsel).lt.dminsel &
                    .and.dabs(y_p(pn)-ysel).lt.dminsel) then
                    dsel=dsqrt((x_p(pn)-xsel)**2+(y_p(pn)-ysel)**2)
                    if(dsel.lt.dminsel) then
                      dminsel=dsel
                      jsel=j
                    endif
                  endif
                endif
              enddo 
              if(jsel.lt.0) then
                write(6,*) ' Error: no star particle found for random selection'
                write(6,*) i,' around x,y,dminsel=',xsel,ysel,dminsel
                stop
              endif
              flagtaken(jsel)=1
              pn=listsel(jsel)
              listselfin(i)=pn
            enddo
            deallocate(flagtaken)
            ! close(62)
          endif
          !output
          mtotselp=0.0d0
          vrmselp=0.0d0
          vtmselp=0.0d0
          vzmselp=0.0d0
          vr2mselp=0.0d0
          vt2mselp=0.0d0
          vz2mselp=0.0d0
          open(60,file='output/lbsels.dat',status='unknown')
          write(60,'(a7,I10)') '# step=',step 
          write(60,'(a30,1pE13.5,a1,1pE13.5,a2,1pE13.5)') &
!            123456789012345678901234567890
            '# Solar position angle from (-',rzasun(0),',',rzasun(1) &
            ,')=',rzasun(2)*180.0d0/M_PI
!                                       123456789012345678901234567890
          write(60,'(a21,3(1pE13.5))') '# Solar motion U,V,W=' &
             ,vsun(0),vsun(1),vsun(2)
!                                       123456789012345678901234567890
          write(60,'(a30,2(1pE13.5))') '# selection l range min,max =' &
            ,glonran(0),glonran(1)
          write(60,'(a30,2(1pE13.5))') '#           b range min,max =' &
            ,glatran(0),glatran(1)
          write(60,'(a30,2(1pE13.5))') '#    distance range min,max =' &
            ,dran(0),dran(1)
          write(60,'(a30,2(1pE13.5))') '#        for star age range =' &
            ,ageran(0),ageran(1)
          write(60,'(a63)') &
!            123456789012345678901234567890123456789012345678901234567890123
            '# x y z vx vy vz glon glat d3d vglon vlos age Rxy vrad vrot phi'
          do i=0,nselfin-1
            pn=listselfin(i)
            write(60,161) x_p(pn),y_p(pn),z_p(pn),vx_p(pn),vy_p(pn),vz_p(pn) &
             ,glonp(pn),glatp(pn),d3dp(pn),vglonp(pn),vlosp(pn),agep(pn) &
             ,rxygalp(pn),vradgalp(pn),vrotgalp(pn),phigalp(pn)
            mtotselp=mtotselp+m_p(pn)
            vrmselp=vrmselp+vradgalp(pn)*m_p(pn)
            vtmselp=vtmselp+vrotgalp(pn)*m_p(pn)
            vzmselp=vzmselp+vz_p(pn)*m_p(pn)
            vr2mselp=vr2mselp+(vradgalp(pn)**2)*m_p(pn)
            vt2mselp=vt2mselp+(vrotgalp(pn)**2)*m_p(pn)
            vz2mselp=vz2mselp+(vz_p(pn)**2)*m_p(pn)
          enddo
          close(60)
          if(mtotselp.gt.0.0d0) then
            vrmselp=vrmselp/mtotselp
            vtmselp=vtmselp/mtotselp
            vzmselp=vzmselp/mtotselp           
            vr2mselp=vr2mselp/mtotselp
            vt2mselp=vt2mselp/mtotselp
            vz2mselp=vz2mselp/mtotselp           
          endif
          if(myrank.eq.0) then
            open(61,file='output/lbsels-vmvsig.asc',status='unknown')     
            write(6,*) ' Mean velocity rad, th, z=',vrmselp,vtmselp,vzmselp
            vrsig=dsqrt(vr2mselp-vrmselp**2)
            vtsig=dsqrt(vt2mselp-vtmselp**2)
            vzsig=dsqrt(vz2mselp-vzmselp**2) 
            write(6,*) ' vsig rad, th, z=',vrsig,vtsig,vzsig
!                              12345678901234567890123456789012345
            write(61,'(a35)') '# Vradm Vthm Vzm vrsig vthsig vzsig'
            write(61,'(6(1pE13.5))') vrmselp,vtmselp,vzmselp &
             ,vrsig,vtsig,vzsig
            close(61)
          endif

! read target star files
          if(flagtargetf.ne.0) then
            if(myrank.eq.0) then
              if(flagtargetf.eq.1) then
                open(50,file='ini/axsymdiskm-fit_sels.asc',status='unknown')
              else
                open(50,file='ini/sels_rv.asc',status='unknown')
                read(50,'(A9,I10)') chread,ntargs
              endif
              write(6,*) ' N target stars=',ntargs
            endif
            call MPI_BCAST(ntargs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! allocate
            allocate(glons(0:ntargs-1))
            allocate(glats(0:ntargs-1))
            allocate(distxys(0:ntargs-1))
            allocate(hrvs(0:ntargs-1))
            allocate(vlons(0:ntargs-1))
            allocate(errhrvs(0:ntargs-1))
            allocate(errvlons(0:ntargs-1))
            allocate(mods(0:ntargs-1))
            allocate(errmods(0:ntargs-1))
            allocate(dists(0:ntargs-1))
            allocate(distmins(0:ntargs-1))
            allocate(distmaxs(0:ntargs-1))
            allocate(xts(0:ntargs-1))
            allocate(yts(0:ntargs-1))
            allocate(zts(0:ntargs-1))
            allocate(ras(0:ntargs-1))   
            allocate(decs(0:ntargs-1))   
            allocate(pmras(0:ntargs-1))   
            allocate(pmdecs(0:ntargs-1))   
            allocate(errpmras(0:ntargs-1))   
            allocate(errpmdecs(0:ntargs-1))   
            allocate(pmradec_corrs(0:ntargs-1))   
            allocate(logps(0:ntargs-1))   
            allocate(loggs(0:ntargs-1))   
            allocate(teffs(0:ntargs-1))   
            allocate(names(0:ntargs-1))
            allocate(vmags(0:ntargs-1))
            allocate(vicols(0:ntargs-1))
            allocate(avs(0:ntargs-1))
    
            if(myrank.eq.0) then
               if(flagtargetf.eq.1) then
                 do i=0,ntargs-1
                   read(50,'(17(1pE13.5),1x,a20)') glons(i),glats(i) &
                     ,distxys(i),hrvs(i),vlons(i),errhrvs(i),errvlons(i) &
                     ,mods(i),errmods(i),ras(i),decs(i),pmras(i) &
                     ,pmdecs(i),errpmras(i),errpmdecs(i) &
                     ,pmradec_corrs(i),logps(i),names(i)
                 enddo
               else
                 do i=0,ntargs-1
                   read(50,'(15(1pE13.5))') xts(i),yts(i),zts(i) &
                     ,hrvs(i),errhrvs(i),vmags(i),vicols(i),ras(i),decs(i) &
                     ,glons(i),glats(i),dists(i),teffs(i),loggs(i),avs(i)
                 enddo
               endif
            endif

            call MPI_BCAST(glons,ntargs,MPI_DOUBLE_PRECISION,0 &
              ,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(glats,ntargs,MPI_DOUBLE_PRECISION,0 &
              ,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(hrvs,ntargs,MPI_DOUBLE_PRECISION,0 &
              ,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(errhrvs,ntargs,MPI_DOUBLE_PRECISION,0 &
              ,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(ras,ntargs,MPI_DOUBLE_PRECISION,0 &
              ,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(decs,ntargs,MPI_DOUBLE_PRECISION,0 &
              ,MPI_COMM_WORLD,ierr)

            if(flagtargetf.eq.0) then
              ! info only from axsymdiskm-fit_sels.asc
              call MPI_BCAST(distxys,ntargs,MPI_DOUBLE_PRECISION,0 &
                ,MPI_COMM_WORLD,ierr)
              call MPI_BCAST(vlons,ntargs,MPI_DOUBLE_PRECISION,0 &
                ,MPI_COMM_WORLD,ierr)
              call MPI_BCAST(errvlons,ntargs,MPI_DOUBLE_PRECISION,0 &
                ,MPI_COMM_WORLD,ierr)
              call MPI_BCAST(mods,ntargs,MPI_DOUBLE_PRECISION,0 &
                ,MPI_COMM_WORLD,ierr)
              call MPI_BCAST(errmods,ntargs,MPI_DOUBLE_PRECISION,0 &
                ,MPI_COMM_WORLD,ierr)
              call MPI_BCAST(pmras,ntargs,MPI_DOUBLE_PRECISION,0 &
                ,MPI_COMM_WORLD,ierr)
              call MPI_BCAST(pmdecs,ntargs,MPI_DOUBLE_PRECISION,0 &
                ,MPI_COMM_WORLD,ierr)
              call MPI_BCAST(errpmras,ntargs,MPI_DOUBLE_PRECISION,0 &
                ,MPI_COMM_WORLD,ierr)
              call MPI_BCAST(errpmdecs,ntargs,MPI_DOUBLE_PRECISION,0 &
                ,MPI_COMM_WORLD,ierr)
              call MPI_BCAST(pmradec_corrs,ntargs,MPI_DOUBLE_PRECISION,0 &
                ,MPI_COMM_WORLD,ierr)
              call MPI_BCAST(logps,ntargs,MPI_DOUBLE_PRECISION,0 &
               ,MPI_COMM_WORLD,ierr)
              call MPI_BCAST(names,ntargs,MPI_CHARACTER,0 &
                ,MPI_COMM_WORLD,ierr)
            else
              ! info from sels_rv.asc
              call MPI_BCAST(xts,ntargs,MPI_DOUBLE_PRECISION,0 &
                ,MPI_COMM_WORLD,ierr)
              call MPI_BCAST(yts,ntargs,MPI_DOUBLE_PRECISION,0 &
                ,MPI_COMM_WORLD,ierr)
              call MPI_BCAST(zts,ntargs,MPI_DOUBLE_PRECISION,0 &
                ,MPI_COMM_WORLD,ierr)
              call MPI_BCAST(vmags,ntargs,MPI_DOUBLE_PRECISION,0 &
                ,MPI_COMM_WORLD,ierr)
              call MPI_BCAST(vicols,ntargs,MPI_DOUBLE_PRECISION,0 &
                ,MPI_COMM_WORLD,ierr)
              call MPI_BCAST(teffs,ntargs,MPI_DOUBLE_PRECISION,0 &
                ,MPI_COMM_WORLD,ierr)
              call MPI_BCAST(loggs,ntargs,MPI_DOUBLE_PRECISION,0 &
                ,MPI_COMM_WORLD,ierr)
              call MPI_BCAST(avs,ntargs,MPI_DOUBLE_PRECISION,0 &
                ,MPI_COMM_WORLD,ierr)
            endif
            if(flagtargetf.eq.0) then
! distance range
              do i=0,ntargs-1
                dists(i)=(10.0d0**((mods(i)+5.0d0)/5.0d0))*0.001d0
                distmins(i)=(10.0d0**((mods(i)-errmods(i)+5.0d0)/5.0d0))*0.001d0
                distmaxs(i)=(10.0d0**((mods(i)+errmods(i)+5.0d0)/5.0d0))*0.001d0
! x,y,z position
                distxyts=dists(i)*dcos(glats(i))
                xts(i)=distxyts*dcos(glons(i))
                yts(i)=distxyts*dsin(glons(i))
                zts(i)=dists(i)*dsin(glats(i))
! rad -> deg
                glons(i)=glons(i)*180.0d0/M_PI
                glats(i)=glats(i)*180.0d0/M_PI
              enddo
            else
              do i=0,ntargs-1
                ! assuming 15% uncertainty
                distmins(i)=dists(i)-0.15*dists(i)
                distmaxs(i)=dists(i)+0.15*dists(i)
                distxys(i)=dists(i)*cos(glats(i)*M_PI/180.0d0)
                vlons(i)=0.0d0
                errvlons(i)=0.0d0
                mods(i)=0.0d0
                errmods(i)=0.0d0
                pmras(i)=0.0d0
                errpmras(i)=0.0d0
                pmdecs(i)=0.0d0
                errpmdecs(i)=0.0d0
                pmradec_corrs(i)=0.0d0
                logps(i)=0.0d0
                names(i)=' '
              enddo
            endif
    
! target selection
            allocate(listagep(0:np-1))
            allocate(flagtaken(0:np-1))
            do i=0,np-1
              flagtaken(i)=0
            enddo
            nagep=0
            do i=ng,np-1
              pn=list_ap(i)
              if(agep(pn).ge.ageran(0).and.agep(pn).le.ageran(1)) then
                listagep(nagep)=pn
                nagep=nagep+1
              endif
            enddo
! find particles 
            allocate(listcanp(0:ns-1))
            if(flagtargetf.eq.0) then
              open(60,file='output/lbsels_targets.dat',status='unknown')
              write(60,'(a62,a62,a34)') &
!            12345678901234567890123456789012345678901234567890123456789012'
            '# Glon Glat Dxy HRV Vlon e_HRV e_Vlon Mod e_Mod Glon_t Glat_t ' &
!             12345678901234567890123456789012345678901234567890123456789012'
            ,'Dxy_t HRV_t Mod_t Dmin_tp xp yp zp xt yt zt e_PMRA e_PMDEC PMR' &
!             12345678901234567890123456789012345678901234567890123456789012'
            ,'ADEC_corr logPer Vlat Vx Vy Vz D3d'
            else
              open(60,file='output/lbsels_targets_rv.dat',status='unknown')
              write(60,'(a46,a49)') &
!               12345678901234567890123456789012345678901234567890123456789012'
               '# Glon Glat D3D xp yp zp vxp vyp vzp HRV Vlon ' &
              ,'FeH Age V VIcol xt yt zt HRV_t e_HRV Teff logg Av'
!               12345678901234567890123456789012345678901234567890123456789012'
            endif
            do i=0,ntargs-1
              ncanp=0
              danglim=degrange
              distlimmax=distmaxs(i)
              distlimmin=distmins(i)
              nsearch=0
 70           do j=0,nagep-1
                pn=listagep(j)
                if(flagtaken(pn).eq.0) then
                  if(d3dp(pn).gt.distlimmin.and.d3dp(pn).lt.distlimmax) then
                  if(dabs(glonp(pn)-glons(i)).lt.danglim.and. &
                     dabs(glatp(pn)-glats(i)).lt.danglim) then
                    dang=dsqrt((glonp(pn)-glons(i))**2 &
                              +(glatp(pn)-glats(i))**2)
                    if(dang.lt.danglim) then
                      listcanp(ncanp)=pn
                      ncanp=ncanp+1
                    endif
                  endif
                  endif
                endif
              enddo
              if(ncanp.eq.0) then
                distlimmin=dists(i)-(dists(i)-distlimmin)*2.0
                if(distlimmin.lt.0.0d0) then
                  distlimmin=0.0d0 
                endif
                distlimmax=dists(i)+(distlimmax-dists(i))*2.0
                danglim=danglim*1.2d0
                nsearch=nsearch+1
                goto 70
              endif 
              if(nsearch.gt.3) then
                write(6,*) i,'star, Nsearch=',nsearch &
                          ,' final ang,dist min max=',danglim &
                  ,distlimmin,distlimmax
              endif
              disptsmin=INF
              pn=-1
              do j=0,ncanp-1
                pn=listcanp(j)
                dispts=dsqrt((xts(i)-x_p(pn))**2+(yts(i)-y_p(pn))**2 &
                            +(zts(i)-z_p(pn))**2)
                if(dispts.lt.disptsmin) then
                  disptsmin=dispts
                  pnts=pn
                endif
              enddo
              if(pn.eq.-1) then
                write(6,*) i,' star could not find particle'
                stop
              else if(flagdup.eq.0) then
                flagtaken(pnts)=1
              endif
              modp=5.0d0*dlog10(d3dp(pnts)*1000.0d0)-5.0d0
              if(flagtargetf.eq.0) then
                write(60,'(31(1pE13.5))') &
                 glonp(pnts),glatp(pnts),dxyp(pnts),vlosp(pnts),vglonp(pnts) &
                ,errhrvs(i),errvlons(i),modp,errmods(i),glons(i) &
                ,glats(i),distxys(i),hrvs(i),vlons(i),mods(i) &
                ,disptsmin,x_p(pnts),y_p(pnts),z_p(pnts),xts(i) &
                ,yts(i),zts(i),errpmras(i),errpmdecs(i),pmradec_corrs(i) &
                ,logps(i),vglatp(pnts),vx_p(pnts),vy_p(pnts),vz_p(pnts) &
                ,d3dp(pnts)
              else
                write(60,'(23(1pE13.5))')  &
                  glonp(pnts),glatp(pnts),d3dp(pnts),x_p(pnts),y_p(pnts) &
                 ,z_p(pnts),vx_p(pnts),vy_p(pnts),vz_p(pnts),vlosp(pnts) &
                 ,vglonp(pnts),fehp(pnts),agep(pnts),vmags(i),vicols(i) &
                 ,xts(i),yts(i),zts(i),hrvs(i),errhrvs(i) &
                 ,teffs(i),loggs(i),avs(i)
              endif
            enddo
! deallocate
            deallocate(glons)
            deallocate(glats)
            deallocate(distxys)
            deallocate(hrvs)
            deallocate(vlons)
            deallocate(errhrvs)
            deallocate(errvlons)
            deallocate(mods)
            deallocate(errmods)
            deallocate(dists)
            deallocate(distmins)
            deallocate(distmaxs)
            deallocate(xts)
            deallocate(yts)
            deallocate(zts)
            deallocate(listagep)
            deallocate(flagtaken)
          endif
        endif
        deallocate(glonp)    
        deallocate(glatp)    
        deallocate(d3dp)    
        deallocate(dxyp)    
        deallocate(vglonp)    
        deallocate(vglatp)    
        deallocate(vlosp)
        deallocate(agep)
        deallocate(vradgalp)
        deallocate(vrotgalp)
        deallocate(rxygalp)
        deallocate(phigalp)
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)

end program
