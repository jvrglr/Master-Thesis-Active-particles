      program softcore_active_2D
        !Compute movement of particles which interact according to a
        !soft core potential-->GEM-alpha potential
        !Active matter, meaning self-propelled particles
        !We integrate the SDE according to Euler method
        !(Aditive noise--> Euler method=Milhstein method)
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
          do i=1,N !compute force over each particle and new position
            do j=i+1,N ! compute all interactions
                call distance(x(i),y(i),x(j),y(j),xj,yj,L,dist) !return xj,yj,dist
                call compute_force
     c      (x(i),y(i),x(j),y(j),xj,yj,L,dist,GEM,R,e,frcex,frcey,N,i,j)
            enddo
            angle(i)=angle(i)+noiser*dran_g()
            !-----------------------------UPDATE POSITION
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

      subroutine compute_force
     c      (x,y,x1,y1,xj,yj,L,dist,GEM,R,e,frcex,frcey,N,i,j)
        !x=x(i),x1=x(j)
        double precision x,y,x1,y1,xj,yj,dist,R,e,fx,fy
        double precision updateF
        integer*4 L,GEM,N,i,j
        double precision frcex(N),frcey(N)
        call distance(x,y,x1,y1,xj,yj,L,dist) !return xj,yj,dist
        !PROGRAM COULD FAIL IF  y=yj
        aux=updateF(GEM,R,e,dist,x,xj)
        fx=aux*(x-xj)
        fy=aux*(y-yj)
        frcex(i)=frcex(i)+fx
        frcex(j)=frcex(j)-fx
        frcey(i)=frcey(i)+fy
        frcey(j)=frcey(j)-fy
        return
      end subroutine

      subroutine coarse_grain(x,y,N,L,R,h,noise,GEM,e)
        !Compute coarse grained density, store it in file
        integer*4 N,L,Ncel,iter,GEM
        double precision R,w,h,noise,dist,xj,yj,fx,fy,dran_g,updateF
        double precision e
        double precision x(N),y(N),xnew(N),ynew(N),frcex(N),frcey(N)
        parameter(w=0.017d0*0.1,Ncel=int(1/w)+1)
         !w=!width of cell
         !Ncel=number of cells in one axis
        integer*4 C(Ncel*Ncel)
        integer*4 i,j,k

        iter=500 !number of iterations to compute averages
        C=0

        do k=1,iter
          frcex=0.0d0
          frcey=0.0d0
          do i=1,N !compute force over each particle and new position
            do j=i+1,N ! compute all interactions
                call distance(x(i),y(i),x(j),y(j),xj,yj,L,dist) !return xj,yj,dist
                call compute_force
     c      (x(i),y(i),x(j),y(j),xj,yj,L,dist,GEM,R,e,frcex,frcey,N,i,j)
            enddo
            xnew(i)=x(i)+h*frcex(i)+noise*dran_g()
            ynew(i)=y(i)+h*frcey(i)+noise*dran_g()
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

          do i=1,N
            j=int(x(i)/w)+int(y(i)/w)*Ncel !label of cell
            C(j)=C(j)+1
          enddo
        enddo
        C=C/(dble(N)*w**2.0d0*dble(iter)) !normalization

        open(unit=3, file="coarse.dat",status="unknown")

        do i=0,Ncel-1
          do k=0,Ncel-1
            j=i+k*Ncel !label of cell
            write(3,*) dble(i)*w/R,dble(k)*w/R,C(j)
           enddo
         enddo
        close(3)

        open(unit=3, file="line_coarse.dat",status="unknown")
        do i=0,Ncel-1
          j=int(3.1d0*R/W)+i*Ncel
          write(3,*) dble(i)*w/R,C(j)
        enddo
        close(3)

        return
      end subroutine


      subroutine distance(x1,y1,x2,y2,xj,yj,L,d)
        !Compute Euclidean distance (d) on torus between particles (x1,y1) and (x2,y2)
        !Saves reflexion of (x2,y2) on (xj,yj)
        !d=distance[(x1,y1),(xj,yj)]
        !d on squared box can't be larger than L/2
        integer*4 L
        double precision x1,y1,xj,yj,x2,y2
        double precision half,delta,d

        half=dble(L)/2.0d0
        xj=x2
        yj=y2

        !boundary conditions on x axis
        delta=x1-xj
        if (delta.gt.half) then
          xj=xj+dble(L)
        else if ((-delta).gt.half) then
          xj=xj-dble(L)
        endif

        !boundary conditions on y axis
        delta=y1-yj
        if (delta.gt.half) then
          yj=yj+dble(L)
        else if ((-delta).gt.half) then
          yj=yj-dble(L)
        endif

        d=dsqrt((x1-xj)**2.0d0+(y1-yj)**2.0d0)

        return
      end subroutine

      double precision function updateF(GEM,R,e,dist,xi,xj)
        !component x of Force that particle j exert on particle i
        integer GEM
        double precision R,e,dist,xi,xj,aux

        if (xi.eq.xj) then
          !soft core, particles can be in same point
          !without this if dist=0.0, computer does not solve
          !properly indetermination giving NaN (important when noise=0)
          updateF=0.0d0
        else
          aux=(R)**dble(GEM)
          updateF=(e*dble(GEM)/aux)*
     c      dexp(-dist**dble(GEM)/aux)*
     c      dist**dble(GEM-2)
        endif

        return
      end
      subroutine initial_orientation(angle,N,pi)
        !Random initial orientation for self-propelled movement
        integer N,i
        double precision angle(N), pi
        double precision dran_u

        do i=1,N
          angle(i)=-pi+dran_u()*2.0d0*pi
        enddo
        return
      end subroutine

      subroutine initial_condition(x,y,N,L)
        !Random initial position
        integer N
        double precision x(N),y(N)
        double precision dran_u
        integer i

        do i=1,N
          x(i)=dran_u()*dble(L)
          y(i)=dran_u()*dble(L)
        enddo

        return
      end

      character(len=20) function str(k)
      !   "Convert an integer to string."
          integer, intent(in) :: k
          write (str, *) k !write to a string
          str = adjustl(str)
      end function str
