      program softcore_active_2D
        !Compute movement of particles which interact according to a
        !soft core potential-->GEM-alpha potential
        !Active matter, meaning self-propelled particles
        !We integrate the SDE according to Euler method
        !(Aditive noise--> Euler method=Milhstein method)
        !---------Notation----------------
        !-->PBC==Periodic Boundary Conditions
        implicit none
        integer*4 N,L,GEM
        double precision D,R,e,alph,h,U,pi,Dr
        parameter (Dr=0.1,pi=4.0d0*datan(1.0d0),
     c      N=2000,L=1,D=0.0001d0,e=0.0333d0,R=0.1d0,GEM=3,U=3.0d0)
        double precision x(N),y(N),xnew(N),ynew(N),frcex(N),frcey(N),
     c  angle(N)
        double precision dran_g,updateF
        double precision noise,noiser,xj,yj,fx,fy,dist,aux
        integer i,j,k,tmax
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

        do k=1,int(tmax/h)
          frcex=0.0d0
          frcey=0.0d0
          do i=1,N !compute forces over each particle and new positions
            do j=i+1,N ! compute all interactions
              !Having into account that we are leading with to body-interactions
              !(Third newton's law F(i-->j)=-F(j-->i))
                call distance(x(i),y(i),x(j),y(j),xj,yj,L,dist) !return xj,yj,dist
                !dist=distance between particle i and j having into account PBC
                call compute_force
     c      (x(i),y(i),x(j),y(j),xj,yj,L,dist,GEM,R,e,frcex,frcey,N,i,j)
            enddo
            !Angular diffusion
            angle(i)=angle(i)+noiser*dran_g()
            !UPDATE POSITION: deterministic force+self propelled movement+brownian motion
            xnew(i)=x(i)+h*frcex(i)+noise*dran_g()+U*h*cos(angle(i))
            ynew(i)=y(i)+h*frcey(i)+noise*dran_g()+U*h*sin(angle(i))
          enddo

          do i=1,N !Boundary conditions
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

        open(unit=2, file="data.dat", status="unknown")
        open(unit=4, file="force.dat", status="unknown")
        open(unit=5, file="self_U.dat", status="unknown")
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
