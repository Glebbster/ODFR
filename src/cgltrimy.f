c        NLS VERSION 20 May 1992, modified AXV
c**********************************************************************
c
c                        SUBROUTINE : CGLTRI
c                        -------------------
c
c       This subroutine will generate the "Lanczos" tridiagonal 
c       from some of the quantities which arise naturally in the
c       course of the conjugate gradients algorithm (in the absence
c       of preconditioning).  This is accomplished by multiplying 
c       out the factors of the Cholesky decomposition of the Lanczos 
c       tridiagonal matrix which are implicitly constructed during 
c       the conjugate gradients procedure.  See pp. 370-1 in "Matrix
c       Computations", G. Golub and C. Van Loan, Johns Hopkins Univ.
c       Press, 1983 and pp.  in " ", J. Cullum and R. Willoughby,
c       Springer-Verlag, 1985 for more information.
c
c       written by DJS 20-SEP-87
c
c**********************************************************************
c
      subroutine cgltri(ndone,al,bl,shiftr,shifti)
c
      include 'dimnls.inc'
      include 'rndoff.inc'
c
      
      integer ndone
      double precision shiftr,shifti
cAXV      double precision a,b
cAXV      dimension a(2,MXSTEP),b(2,MXSTEP)
      complex*16 shiftc
      complex*16 al(MXSTEP),bl(MXSTEP)
      complex*16 c(MXSTEP),tchalf
c
      integer i,j
      double precision amp,phase,tr,ti,trhalf,tihalf      
cAXV      double precision c
cAXV      dimension c(2,MXSTEP)
      parameter (CZERO=(0.0D0,0.0D0))

c
c######################################################################
c
      if (ndone .lt. 1) return

      if (ndone .gt. MXSTEP) then
         write(*,*) 'CGLTRI error: ndone exceeds MXSTEP'
         return
      endif
      do 10 i=1,ndone
        c(i)=al(i)
 10   continue
c
      do 20 i=2,ndone
        j=i-1
        tr=dreal(bl(i))
        ti=dimag(bl(i))
c
        amp=(tr*tr+ti*ti)**0.25D0
        if (amp.gt.rndoff) then
          phase=0.5D0*datan2(ti,tr)
          trhalf=amp*cos(phase)
          tihalf=amp*sin(phase)
          if (abs(trhalf).lt.rndoff) trhalf=0.0D0
          if (abs(tihalf).lt.rndoff) tihalf=0.0D0
          tchalf=dcmplx(trhalf,tihalf)          
        else
          tchalf=CZERO
        end if
c
        bl(j)=c(j)*tchalf
        if (dreal(bl(j)).lt.0.0D0) then
          bl(j)=-bl(j)
        end if
        al(i)=c(i)+c(j)*bl(i)
 20   continue
c
      shiftc=dcmplx(shiftr,shifti)
      do 30 i=1,ndone
        al(i)=al(i)-shiftc
 30   continue
c
      bl(ndone)=CZERO
c
c----------------------------------------------------------------------
c     return to caller
c----------------------------------------------------------------------
c
      return 
      end
