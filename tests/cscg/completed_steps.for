      program completed_steps
      implicit none
      include 'dimnls.inc'
      include 'eprmat.inc'
      double complex b(MXDIM),x(MXDIM),al(MXSTEP),bl(MXSTEP)
      double complex a11,a22,a12,z,det,q11,q12,phase
      double precision tol,error,value(1),expected
      integer preset,mode,limit,ndone,deriv
      logical complete

c     Complex symmetric 2-by-2 matrix with starting vector (1,0).
      a11=(2d0,0.3d0)
      a22=(3d0,-0.2d0)
      a12=(0.7d0,0.15d0)
      zdiag(:,1)=(/dble(a11),dimag(a11)/)
      zdiag(:,2)=(/dble(a22),dimag(a22)/)
      jzmat(1:3)=(/1,2,2/)
      kzmat(1:3)=(/1,2,2/)
      izmat(1)=2
      izmat(MXEL)=2
      zmat(1)=dimag(a12)
      zmat(MXEL)=dble(a12)
      b(1:2)=(/ (1d0,0d0),(0d0,0d0) /)
      z=(0.12d0,0.85d0)
      phase=exp(dcmplx(0d0,28d0*acos(-1d0)/180d0))

      do preset=0,1
         do mode=0,3
            al(1:3)=preset*(11d0,-4d0)
            limit=2
            select case(mode)
            case(0)
               tol=2d0
            case(1)
               tol=0.45d0
            case(2)
               tol=1d-12
            case(3)
               tol=1d-12
               limit=3
            end select
            call cscg(b,2,limit,tol,(0.2d0,0.4d0),x,al,bl,ndone,error)
            select case(mode)
            case(0)
               if (ndone.ne.0) stop 1
               call cfs(al,bl,0,dimag(z),-0.1d0,1,dble(z),
     #                  0,28d0,value,complete)
               if (complete) stop 2
               cycle
            case(1)
               if (ndone.ne.1) stop 3
               q11=1d0/(a11+z)
               q12=(0d0,0d0)
            case(2,3)
               if (mode.eq.2.and.ndone.ne.-2) stop 4
               if (mode.eq.3.and.ndone.ne.2) stop 5
               det=(a11+z)*(a22+z)-a12*a12
               q11=(a22+z)/det
               q12=-a12/det
            end select
            do deriv=0,1
               call cfs(al,bl,abs(ndone),dimag(z),-0.1d0,1,dble(z),
     #                  deriv,28d0,value,complete)
               expected=dble(phase*q11)
               if (deriv.eq.1)
     #            expected=dble(phase*(0d0,1d0)*(q11*q11+q12*q12))
               if (.not.complete) stop 6
               if (.not.(abs(value(1)-expected).lt.1d-12)) stop 7
            enddo
         enddo
      enddo
      end
