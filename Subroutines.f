      subroutine deterministic_term_plus_radial
        !This sunroutine is just equal to "deterministic_term", but it also saves
        !the histogram for radial frequencies.

        !TO INCREASE EFFICIENCY, I COULD COMPUTE THE COS AND SIN JUST ONCE AND SAVE IT
        !I WOULD NEED ANOTHER TWO VARIABLES (MY ACTUAL SN AND CS) TO SAVE AVERAGE
     c    (fx,fy,angle,x,y,N,L,GEM,R,e,noiser,g,lmax,bin)
       integer*4 N,i,j,GEM,L,label
       double precision R,e,noiser,xj,yj,dist,dran_g,bin
       double precision fx(N),fy(N),angle(N),x(N),y(N),cs(N),sn(N),
     c  g(0:lmax)
       integer*4 cont(N)

       cont=1 !counts of number of "Vicsek-neighbours" of each bird
       !each bird is one of its own "Vicsek-neighbours"
       sn=0.0d0
       cs=0.0d0
       fx=0.0d0
       fy=0.0d0
       g=0.0d0
       do i=1,N
         cs(i)=cs(i)+dcos(angle(i))! we need them to compute angular averages properly
         sn(i)=sn(i)+dsin(angle(i))
         do j=i+1,N ! compute all interactions
           !Having into account that we are dealing with two body-interactions
           !(Third newton's law F(i-->j)=-F(j-->i))
           call distance(x(i),y(i),x(j),y(j),xj,yj,L,dist) !returns xj,yj,dist
           !dist=distance between particle i and j having into account PBC
           !xj,yj are reflexion of x(j),y(j) on boundaries (needed to compute sense of force)
           call compute_force (x(i),y(i),x(j),y(j), !Repulsive soft core force
     c      xj,yj,L,dist,GEM,R,e,fx,fy,N,i,j)
           call vicsek_interaction(dist,sn,cs,angle,i,j,cont,N) !Vicsek-like alignment, updates sn,cs,cont
           label=int(dist/bin)
           g(label)=g(label)+1.0d0/dble(label)
         enddo
         !Angular diffusion + Vicsek force-->update self-propelled angle
         cs(i)=cs(i)/dble(cont(i))
         sn(i)=sn(i)/dble(cont(i))
         angle(i)=atan2(sn(i),cs(i))+noiser*dran_g()
       enddo
       return
      end subroutine

      subroutine deterministic_term
        !TO INCREASE EFFICIENCY, I COULD COMPUTE THE COS AND SIN JUST ONCE AND SAVE IT
        !I WOULD NEED ANOTHER TWO VARIABLES (MY ACTUAL SN AND CS) TO SAVE AVERAGE
     c    (fx,fy,angle,x,y,N,L,GEM,R,e,noiser)
       integer*4 N,i,j,GEM,L
       double precision R,e,noiser,xj,yj,dist,dran_g
       double precision fx(N),fy(N),angle(N),x(N),y(N),cs(N),sn(N)
       integer*4 cont(N)

       cont=1 !counts of number of "Vicsek-neighbours" of each bird
       !each bird is one of its own "Vicsek-neighbours"
       sn=0.0d0
       cs=0.0d0
       fx=0.0d0
       fy=0.0d0
       do i=1,N
         cs(i)=cs(i)+dcos(angle(i))! we need them to compute angular averages properly
         sn(i)=sn(i)+dsin(angle(i))
         do j=i+1,N ! compute all interactions
           !Having into account that we are dealing with two body-interactions
           !(Third newton's law F(i-->j)=-F(j-->i))
           call distance(x(i),y(i),x(j),y(j),xj,yj,L,dist) !returns xj,yj,dist
           !dist=distance between particle i and j having into account PBC
           !xj,yj are reflexion of x(j),y(j) on boundaries (needed to compute sense of force)
           call compute_force (x(i),y(i),x(j),y(j), !Repulsive soft core force
     c      xj,yj,L,dist,GEM,R,e,fx,fy,N,i,j)
           call vicsek_interaction(dist,sn,cs,angle,i,j,cont,N) !Vicsek-like alignment, updates sn,cs,cont
         enddo
         !Angular diffusion + Vicsek force-->update self-propelled angle
         cs(i)=cs(i)/dble(cont(i))
         sn(i)=sn(i)/dble(cont(i))
         angle(i)=atan2(sn(i),cs(i))+noiser*dran_g()
       enddo
       return
      end subroutine

      subroutine vicsek_interaction
     c   (dist,sn,cs,angle,i,j,cont,N)
!Vicsek-like alingment interaction
c        REFERENCES:
C       ---------------------------------------------------------------------------
C        -->Vicsek, T., Czirók, A., Ben-Jacob, E., Cohen, I., & Shochet, O. (1995).
C        Novel type of phase transition in a system of self-driven particles.
C        Physical review letters, 75(6), 1226.
C       ---------------------------------------------------------------------------
        integer N,i,j
        double precision angle(N),sn(N),cs(N)
        integer cont(N)
        double precision dran_g,Vic_R,dist
        parameter (Vic_R=0.2d0)

        if (dist.le.Vic_R) then
          sn(i)=sn(i)+dsin(angle(j)) !it will be more efficient to store sin(angle(i))
          sn(j)=sn(j)+dsin(angle(i))

          cs(i)=cs(i)+dcos(angle(j))
          cs(j)=cs(j)+dcos(angle(i))

          cont(i)=cont(i)+1
          cont(j)=cont(j)+1

        endif

        return
      end

      double precision function v_mean(N,angle)
C       Computes order parameter for Viksec model
C       If "birds" are aligned the outcome will be one
C       If final state is isotrope, the outcome will be zero
c        REFERENCES:
C       ---------------------------------------------------------------------------
C        -->Vicsek, T., Czirók, A., Ben-Jacob, E., Cohen, I., & Shochet, O. (1995).
C        Novel type of phase transition in a system of self-driven particles.
C        Physical review letters, 75(6), 1226.
C       ---------------------------------------------------------------------------
        integer*4 N,i
        double precision angle(N)
        double precision sumx,sumy

        sumx=0.0d0
        sumy=0.0d0
        do i=1,N
          sumx=sumx+cos(angle(i))
          sumy=sumy+sin(angle(i))
        enddo
        v_mean=dsqrt(sumx**2.0d0+sumy**2.0d0)/dble(N)
        return
      end

      double precision function normalize(x,N,pi,dr)
        !Normalize histogram to have a density probability function
        !x-->frequency in particular bin
        !N-->total number of samples
        !output--> in [0,1], is the prob of finding a particle in the bin
c        REFERENCES:
C       ---------------------------------------------------------------------------
c        Delfau, J. B., López, C., & Hernández-García, E. (2017).
c        Active cluster crystals. New Journal of Physics, 19(9), 095001.
c       -----------------------------------------------------------------------------
        integer*4 N
        double precision pi,dr,x

        normalize=x/(2.0d0*pi*(dr**2.0d0)*dble(N))
        return
      end

      subroutine compute_force
     c      (x,y,x1,y1,xj,yj,L,dist,GEM,R,e,frcex,frcey,N,i,j)
        !x=x(i),x1=x(j)
        double precision x,y,x1,y1,xj,yj,dist,R,e,fx,fy,aux
        double precision updateF
        integer*4 L,GEM,N,i,j
        double precision frcex(N),frcey(N)
        !call distance(x,y,x1,y1,xj,yj,L,dist) !return xj,yj,dist
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

      subroutine distance(x1,y1,x2,y2,xj,yj,L,d)
        !this routine will  return xj,yj and d
        !--------> xj,yj
        !reflexion of (x2,y2) is saved on (xj,yj)
        ! (xj,yj) is the same point that (x2,y2) in the torus
        ! it is needed to save (xj,yj) in order to compute direcction/sence of forces
        !--------> d
        !Compute Euclidean distance (d) on torus between particles (x1,y1) and (x2,y2)
        !d=distance[(x1,y1),(xj,yj)]
        !dx=sqrt((x1-xj)^2) on torus can't be larger than L/2
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
        integer N,L
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
