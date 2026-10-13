c NLSL Version 1.4 10/10/94 
c----------------------------------------------------------------------
c                    ======================
c                      subroutine REALFT
c                    ======================
c     From Numerical Recipes by Press et al.
c----------------------------------------------------------------------
      subroutine realft(dat,n,isign)
      implicit none
      integer i,i1,i2,i3,i4,isign,n,n2p3
      double precision c1,c2,h1r,h1i,h2r,h2i,wrs,wis,wr,wi,wpr,wpi,
     #                 wtemp,theta
      double precision dat(2*n)
c
      theta=6.28318530717959d0/2.0d0/dble(n)
      c1=0.5d0
      if (isign.eq.1) then
        c2=-0.5d0
        call four1(dat,n,+1)
      else
        c2=0.5d0
        theta=-theta
      endif
      wpr=-2.0d0*dsin(0.5d0*theta)**2
      wpi=dsin(theta)
      wr=1.d0+wpr
      wi=wpi
      n2p3=2*n+3
      do 11 i=2,n/2+1
        i1=2*i-1
        i2=i1+1
        i3=n2p3-i2
        i4=i3+1
        wrs=sngl(wr)
        wis=sngl(wi)
        h1r=c1*(dat(i1)+dat(i3))
        h1i=c1*(dat(i2)-dat(i4))
        h2r=-c2*(dat(i2)+dat(i4))
        h2i=c2*(dat(i1)-dat(i3))
        dat(i1)=h1r+wrs*h2r-wis*h2i
        dat(i2)=h1i+wrs*h2i+wis*h2r
        dat(i3)=h1r-wrs*h2r+wis*h2i
        dat(i4)=-h1i+wrs*h2i+wis*h2r
        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
11    continue
      if (isign.eq.1) then
        h1r=dat(1)
        dat(1)=h1r+dat(2)
        dat(2)=h1r-dat(2)
      else
        h1r=dat(1)
        dat(1)=c1*(h1r+dat(2))
        dat(2)=c1*(h1r-dat(2))
        call four1(dat,n,-1)
      endif
      return
      end

c----------------------------------------------------------------------
c                    ======================
c                      subroutine FOUR1
c                    ======================
c     From Numerical Recipes by Press et al.
c----------------------------------------------------------------------
      subroutine four1(dat,nn,isign)
      implicit none
      integer i,istep,j,m,mmax,n,nn,isign
      double precision wr,wi,wpr,wpi,wtemp,tempi,tempr,theta
      double precision dat(2*nn)
c
      n=2*nn
      j=1
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=dat(j)
          tempi=dat(j+1)
          dat(j)=dat(i)
          dat(j+1)=dat(i+1)
          dat(i)=tempr
          dat(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        go to 1
        endif
        j=j+m
11    continue
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*dsin(0.5d0*theta)**2
        wpi=dsin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
          do 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*dat(j)-sngl(wi)*dat(j+1)
            tempi=sngl(wr)*dat(j+1)+sngl(wi)*dat(j)
            dat(j)=dat(i)-tempr
            dat(j+1)=dat(i+1)-tempi
            dat(i)=dat(i)+tempr
            dat(i+1)=dat(i+1)+tempi
12        continue
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
        mmax=istep
      go to 2
      endif
      return
      end
