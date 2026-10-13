c NLSL Version 1.5.1 beta 11/25/95, modified AXV
c*********************************************************************
c                       ==================
c                        subroutine:EPRLS
c                       ==================
c
c       Subroutine version of EPRLL family of programs by D. Schneider.
c       This routine is intended for use with nonlinear least-squares
c       applications. The routine calculates the tridiagonal matrix
c       for the parameters passed in the fparm and iparm arrays using
c       the Conjugate Gradients version of the Lanczos algorithm.
c
c       Calculation parameters are input in the arrays fparm and iparm.
c       They are expected in the same order as that found in common 
c       /eprprm/. The input flag icalc is specified as a three-digit
c       number ijk, where
c
c       j.ne.0  => Perform matrix and CG calculations (matrll/cscg) 
c       k.ne.0  => Perform starting vector and CG calculations (stvect/cscg)
c
c       The continued-fractions calculation is always performed.
c
c          Subroutine arguments are output as follows:
c
c             al(ndone)    Diagonal of tridiagonal matrix for spectrum
c             be(ndone)    Off-diagonal of tridiagonal matrix for spectrum
c             ndone        Number of CG steps taken. Zero for Lanczos
c                              calculations, and set negative if CG did
c                              not converge
c             ierr         Error flag. Meaning of values assigned to ierr
c                          are defined in 'errmsg.inc'.
c
c       Written by DEB, 26-Sep-91 based on programs from DJS 
c
c       Uses :
cAXV  pmatrl.f and pstvec.f are not used as they do not work with given bss matrix
c               cscg.f       Conjugate gradients tridiagonalization
c               cfs.f        Continued-fraction spectral calculation
c
c*********************************************************************
c
      subroutine eprls(icalc,jsite,al,be,spectr,nft,nstp,ndone,ierr)
      implicit none
c
      include 'dimnls.inc'
      include 'eprprm.inc'
      include 'eprmat.inc'
      include 'prmeqv.inc'
      include 'errmsg.inc'
      include 'pidef.inc'
c
      double precision EIGHT,ONE,THIRD,TWO,ZERO,EPS
      double complex CZERO
      parameter (ZERO=0.0D0,ONE=1.0D0,THIRD=0.333333333333333D0,
     #           TWO=2.0D0,EIGHT=8.0D0,EPS=1.0D-3,CZERO=(0.0D0,0.0D0))
c
      integer ndone,icalc,ierr,nft,jsite	
c
      double precision spectr((MXSPT))
      double complex al(MXDIM),be(MXDIM)
      double complex stv(MXDIM),y(MXDIM)

      logical matrix,start
c
      double precision cgerr,wline,gib
      integer i,nstp,iadd
c
      double precision dblint
      external dblint
cAXV:
     	double precision fldiShifted,AaddTemp(MaxAddNucl),sin2psi
c
c#####################################################################
c
      matrix = mod(icalc,100)/10 .ne. 0
      start  = mod(icalc,10) .ne. 0
c
      if (matrix) then
         call matrll(ierr)
c
         if (ierr.gt.FATAL) then
            write(*,*)eprerr(ierr)
            stop
         endif
      end if
c
c----------------------------------------------------------------------
c        Generate starting vector
c----------------------------------------------------------------------
      if (start) then
            call stvect(stv,ierr)
      end if
c
      if (start .or. matrix) then
c
c---------------------------------------------------------------------
c       Call conjugate gradients tridiagonalization routine
c---------------------------------------------------------------------
c
         do i=1,ndim
            y(i)=CZERO
         end do
c
c       -----------------------------------------------------
c        Set the number of steps for CG tridiagonalization
c        If nstep is unspecified, use the matrix dimension or
c        the number of steps previously taken, whichever is
c        smaller
c       -----------------------------------------------------
c
         call cscg( stv,ndim,nstp,cgtol,dcmplx(shiftr,shifti),
     #              y,al,be,ndone,cgerr )
c
c         Check for user input or halt before proceeding
c
c         *** CG did not converge
         if (ndone.lt.0) then
            ndone=-ndone
            ierr=NOCONVRG
         end if
c
      end if	 !start or matrix
c
c     ---------------------------------------------
c     Continued-fraction calculation of spectrum
c     ---------------------------------------------
cAXV Shifting of the calculated spectrum
      fldiShifted=b0-fldi+SpcShift+fldShftSite(jsite)
      call cfs( al,be,ndone,fldiShifted,-dfld,nfld,w0+lbl,ideriv,phase,
     #          spectr )
c
c     ------------------------------------
c     Convolution with Gaussian lineshape
c     ------------------------------------
      sin2psi=dsin( RADIAN*psi )**2
     	if (IFLGGIB.eq.0) then
         gib=gib0+gib2*sin2psi
     	   else
         gib=gib0+(gibx-gib0)*sin2psi
     	   endif           
      wline=dsqrt(gib*gib+lbg*lbg)
c
      do iadd=1,NAddN
        AaddTemp(iadd)=Aadd(2*iadd-1)+
     #	  (Aadd(2*iadd)-Aadd(2*iadd-1))*sin2psi
      enddo
c
      call gconvl(spectr,wline,dfld,nfld,nft,AaddTemp,in2Add,NAddN)
c
c     ------------------
c     Normalize spectrum (PI comes from Lorentzian shape)
c     ------------------
      do i=1,nfld
         spectr(i)=spectr(i)/PI
      end do
c
      return
c
      end
