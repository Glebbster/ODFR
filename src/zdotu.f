c Version 1.4 10/10/94, modified AXV
c*********************************************************************
c
c       This complex double precision function returns the dot
c       product of the two input vectors without complex 
c       conjugation.
c
c       Calling Parameters:
c       ------------------
c
c       vectors:
c               x    : input vector
c               y    : input vector
c
c       scalars:
c               ndim : dimension of the vectors x and y
c
c
c       Local Variables:
c       ---------------
c               ir   : row counter
c               accr : real part of dot product
c               acci : imaginary part of dot product
c
c       Notes:
c       -----
c               It is assumed that the input vectors are stored in a
c               double precision complex array or as a real double 
c               precision array dimensioned as 2 X mxdim in the 
c               calling routine.
c
c
c
c       written by DJS 26-AUG-87
c
c*********************************************************************
c
      function zdotu(x,y,ndim)
      implicit none
c
      double complex zdotu
c
       include 'dimnls.inc'
      include 'rndoff.inc'
c
      integer ndim
cAXV      double precision x,y
cAXV      dimension x(2,MXDIM),y(2,MXDIM)
      complex*16 x(MXDIM),y(MXDIM)
      double precision xr,xi,yr,yi      
      
c
      integer ir
      double precision accx,accy,accr,acci,scale
c
      double precision ZERO
      parameter (ZERO=0.0D0)
c
c######################################################################
c
      accx=ZERO
      accy=ZERO
      accr=ZERO
      acci=ZERO
c
      do ir=1,ndim
           xr=dreal(x(ir))
           xi=dimag(x(ir))
           yr=dreal(y(ir))
           yi=dimag(y(ir))                               
         accx=accx+xr*xr+xi*xi
         accy=accy+yr*yr+yi*yi 
         accr=accr+xr*yr-xi*yi
         acci=acci+xr*yi+xi*yr
      end do
c
      scale=sqrt(accx*accy)
c
      if (dabs(accr)/scale.lt.RNDOFF) accr=ZERO
      if (dabs(acci)/scale.lt.RNDOFF) acci=ZERO
c
      zdotu=dcmplx(accr,acci)
c
      return
      end
