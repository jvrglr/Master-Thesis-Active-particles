      program softcore_active_2D
        !Compute movement of particles which interact according to a
        !soft core potential-->GEM-alpha potential
        !Active matter, meaning self-propelled particles
        !We integrate the SDE according to Euler method
        !(Aditive noise--> Euler method=Milhstein method)
        !---------Notation----------------
        !-->PBC==Periodic Boundary Conditions
        implicit none
        integer*4 N,L,GEM,Npairs,lmax
        double precision D,R,e,alph,h,U,pi,Dr,bin,L2
        parameter (Dr=0.1,pi=4.0d0*datan(1.0d0),
     c      N=2000,L=1,D=0.0001d0,e=0.0333d0,R=0.1d0,GEM=3,U=3.0d0,
     c  bin=0.01,Npairs=int(0.5*N*(N-1)),
     c  lmax=int(L*dsqrt(2.0d0)/(2.0*bin)),L2=dble(L)/2.0d0)
        double precision x(N),y(N),xnew(N),ynew(N),frcex(N),frcey(N),
     c  angle(N),sn(N),cs(N),Eulerx(N),Eulery(N),g(0:lmax)
        double precision dran_g,updateF,normalize,v_mean
        double precision noise,noiser,xj,yj,fx,fy,dist
        integer i,j,k,tmax,cont(N),dummy,label
        real start,finish !To compute CPU time
        character(len=20) str !Character function
        write(*,*) "--------------------------"
        write(*,*) "There are",N,"particles"
        write(*,*) "Dimensionless parameters:"
        write(*,*) U/dsqrt(Dr*D),e*Dr/U**2.0d0,Dr*R/U
        write(*,*) "--------------------------"
        call cpu_time(start)
        call dran_ini(1994)
        call initial_condition(x,y,N,L)
        call initial_orientation(angle,N,pi)

        tmax=1 !Thermalization time
        h=0.001d0 !Integration time step
        noise=dsqrt(2.0d0*D*h)
        noiser=dsqrt(2.0d0*Dr*h)
        sn=0.0d0
        cs=0.0d0
        frcex=0.0d0
        frcey=0.0d0
        dummy=0
        do k=1,int(tmax/h)
        !do k=1,300
c          dummy=dummy+1
c          if (dummy.eq.1) then
c            dummy=0
C          open(unit=2, file='data.'//trim(str(k))//'.dat',
C     c            status="unknown")
          open(unit=4, file='force.'//trim(str(k))//'.dat',
     c            status="unknown")
c          open(unit=5, file='self_U.'//trim(str(k/1))//'.dat',
c     c            status="unknown")
          do i=1,N
c           write(2,*) x(i)/R,y(i)/R
           write(4,*) x(i)/R,y(i)/R,frcex(i),frcey(i),
     c       U*cos(angle(i)),U*sin(angle(i))
c           write(5,*) x(i)/R,y(i)/R,U*cos(angle(i)),U*sin(angle(i))
          enddo
c          close(2)
          close(4)
c         close(5)
c          endif

!----------------------------------------------------Euler-Muruyama method starts
          call deterministic_term_plus_radial !updates frcex,y and angle of displacement
     c    (frcex,frcey,angle,x,y,N,L,GEM,R,e,noiser,g,lmax,bin)
          do i=1,N
            Eulerx(i)=U*dcos(angle(i))+frcex(i) !I need to save it for Heun method
            Eulery(i)=U*dsin(angle(i))+frcey(i)
            !UPDATE POSITION: deterministic force+self propelled movement+brownian motion
            xnew(i)=x(i)+h*Eulerx(i)+noise*dran_g()
            ynew(i)=y(i)+h*Eulery(i)+noise*dran_g()
          enddo
!-----------------------------------------------------Euler-Muruyama method ends
!-----------------------------------------------------Heun method starts
          call deterministic_term !updates frcex,y and angle of displacement
     c    (frcex,frcey,angle,xnew,ynew,N,L,GEM,R,e,noiser)
           do i=1,N
             Eulerx(i)=Eulerx(i)+U*dcos(angle(i))+frcex(i) !I need to save it for Heun method
             Eulery(i)=Eulery(i)+U*dsin(angle(i))+frcey(i)
             !UPDATE POSITION: deterministic force+self propelled movement+brownian motion
             xnew(i)=x(i)+0.5d0*h*Eulerx(i)+noise*dran_g()
             ynew(i)=y(i)+0.5d0*h*Eulery(i)+noise*dran_g()
           enddo
!------------------------------------------------------Heun method ends
!------------------------------------------------------Print radial distribution
          open(unit=1,file="radial."//trim(str(k))//'.dat',
     c            status="unknown")
          do label=0,int(L2/bin)
            write(1,*) label*bin/R,normalize(g(label),Npairs,pi,bin)
          enddo
          close(1)
!-------------------------------------------------------Boundary conditions
          do i=1,N
            if (xnew(i).gt.L) then
              xnew(i)=mod(xnew(i),dble(L))
            endif
            if (ynew(i).gt.L) then
              ynew(i)=mod(ynew(i),dble(L))
            endif

            if (xnew(i).lt.0) then
              xnew(i)=L-mod(-xnew(i),dble(L))
            endif
            if (ynew(i).lt.0) then
              ynew(i)=L-mod(-ynew(i),dble(L))
            endif
            x(i)=xnew(i)
            y(i)=ynew(i)
          enddo
        enddo
      open(unit=1, file="Vicsek.dat",status="unknown")
        write(1,*) "Vicsek order parameter",v_mean(N,angle),"Rv=2.0"
        write(*,*) "Vicsek order parameter",v_mean(N,angle),"Rv=2.0"
        close(1)

!----------------------------------------------------------Save last state
        open(unit=2, file="data.dat0", status="unknown")
        open(unit=4, file="force.dat0", status="unknown")
        open(unit=5, file="self_U.dat0", status="unknown")
      write(2,*) "#L=",L,"N=",N,"D=",D,"epsil=",e,"GEM=",GEM,"Dr=",Dr
      write(2,*) "#x(i)/R y(i)/R"
      write(4,*) "#x(i)/R y(i)/R","Fx","Fy"
      write(2,*) "#x(i)/R y(i)/R","Ux","Uy"
        do i=1,N
          write(2,*) x(i)/R,y(i)/R
          write(4,*) x(i)/R,y(i)/R,frcex(i),frcey(i)
          write(5,*) x(i)/R,y(i)/R,U*cos(angle(i)),U*sin(angle(i))
        enddo
        close(2)
        close(4)
        close(5)

        call cpu_time(finish)
        write(*,*) "CPU time",finish-start !Time in seconds
        !call coarse_grain(x,y,N,L,R,h,noise,GEM,e)

      end
