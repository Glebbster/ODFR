C for ODFR 10/04/2012
c*********************************************************************
c
c                  Wigner reduces matrix elements d(j,ms,m,bet)
c                   -----------------------------
c*********************************************************************
c
      double precision function Wig(j,ms,m,bet)
c
      integer j,m,ms
      double precision bet
c
      external FACmy
      double precision FACmy
      integer nu,imin,imax,n1,n2
      double precision temp1,temp2,temp3,bet_half
c      double precision tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
c
      intrinsic dsqrt
c
c#####################################################################
c
c      if ((m.lt.0).or.(m.gt.j).or.(abs(z).gt.1.0D0)) then
c        go to 9999
c      end if
c      if ((ms.lt.0).or.(ms.gt.j)) then
c        go to 9999
c      end if
       temp1=FACmy(j+m)
       temp1=temp1*FACmy(j-m)
       temp1=temp1*FACmy(j+ms)
       temp1=temp1*FACmy(j-ms)
c
       temp2=dsqrt(temp1)
       imax=MIN(j-ms,j+m)
       imin=MAX(0,m-ms)
       bet_half=bet/2.d0
       Wig=0.0d0
       do nu=imin,imax
          n1=2*j + m - ms - 2*nu
          n2=ms - m + 2*nu
          temp3=FACmy(j-ms-nu)
          temp3=temp3*FACmy(j+m-nu)
          temp3=temp3*FACmy(nu+ms-m)
          temp3=temp3*FACmy(nu)
c
          Wig=Wig+(-1)**nu*dcos(bet_half)**n1
     # *(-dsin(bet_half))**n2/temp3
       enddo
       Wig=Wig*temp2
c
          return
c
c---------------------------------------------------------------------
c     exit from program if improper arguments detected
c---------------------------------------------------------------------
c
c 9999     write (*.*) 'Improper parameters for d(j,m,ms,bet) !'
c          stop
          end

c###############################################################
       DOUBLE PRECISION FUNCTION FACmy(I)
C
        FACmy=1.d0 
	  II=I
        DO while (II.gt.1)
          FACmy=FACmy*II
          II=II-1
        enddo
       RETURN
       END
c###############################################################
