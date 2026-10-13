c For ODFR, Aleksey Bogdanov
	subroutine ODFint(kNet,anegODF)
c           
      include 'dimnls.inc'
      include 'pidef.inc'
      include 'dimodfl.inc'

      integer jODF
      double precision anegODF  !,result, absresult
      double precision, dimension (2,2) ::  Arres
      double precision,dimension (2) :: ttnet 
      double precision,dimension (2,2)::XxL

      common /A_jmn/A_jmn,jODF,numODFs
      external ODFL_Ro_for_Int

c/////////////////////////////////////////////////////
 	
   	ttnet(1)=kNet
      ttnet(2)=kNet 
   
      XxL(1,1)=0.D0
      XxL(1,2)=2.d0*PI

      XxL(2,1)=0
      XxL(2,2)=PI

      CALL INKOR(2,2,ttnet,1,XxL,ODFL_Ro_for_Int,Arres,0.0d0)
c
      anegODF=dabs((Arres(2,1)-Arres(1,1))/2/Arres(2,1))
c      write(*,*)Arres(2,1),Arres(1,1)
c
      end subroutine

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine ODFL_Ro_for_Int(NN, KK, XX, FF)
c  This subroutine interfaces ODFL_Ro with INKOR integration routine
      double precision, dimension (KK) :: XX
      double precision, dimension (NN) :: FF

      external ODFL_RoReFast
      double precision ODFL_RoReFast
c/////////////////////////////////////////////////////
c
      FF(1)=ODFL_RoReFast(XX(1),XX(2),0.d0)*dsin(XX(2))
      FF(2)=dabs(FF(1))
c
      end subroutine

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine OrdParPot(kNet,NINT,fNrm,ordTEMP)
c           
      include 'dimnls.inc'
      include 'pidef.inc'

      integer knet,NINT,ii
      double precision fNrm !potetials
      double precision, dimension (Jmax) ::ordTEMP ! order parameters
      double precision, dimension (Jmax,2) ::  Arres
      double precision,dimension (2) :: ttnet 
      double precision,dimension (2,2)::XxL

      external FunOrdParPot

c/////////////////////////////////////////////////////
 	
   	ttnet(1)=kNet
      ttnet(2)=kNet 
   
      XxL(1,1)=0.D0
      XxL(1,2)=2.d0*PI

      XxL(2,1)=0
      XxL(2,2)=PI

      do ii=1,Jmax
         ordTEMP(ii)=0.d0
      enddo 
	     
      CALL INKOR(Nint,2,ttnet,1,XxL,FunOrdParPot,Arres,0.0d0)

      do ii=1,Nint
         ordTEMP(ii)=Arres(ii,1)/4/PI/fNrm
      enddo 
     
      end subroutine

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      subroutine FunOrdParPot(NN, KK, XX, FF)
       IMPLICIT NONE
c  This subroutine interfaces ODFL_Ro with INKOR integration routine
	integer KK,NN
      double precision, dimension (KK) :: XX
      double precision, dimension (NN) :: FF
      integer	ii,Jtemp,mtemp
      double precision CPtl,pocexp,expm
      common /FunOrdPar/CPtl(5),Jtemp 
      external CompWigner
      double complex CompWigner
c/////////////////////////////////////////////////////
      pocexp= real(CPtl(1)*CompWigner(2,0,0,XX(1),XX(2),0.0d0)
     #            +2*CPtl(2)*CompWigner(2,2,0,XX(1),XX(2),0.0d0)
     #            +CPtl(3)*CompWigner(4,0,0,XX(1),XX(2),0.0d0)
     #            +2*CPtl(4)*CompWigner(4,2,0,XX(1),XX(2),0.0d0)
     #            +2*CPtl(5)*CompWigner(4,4,0,XX(1),XX(2),0.0d0))
      expm=exp(pocexp)     

      do ii=1,NN
	  mtemp=2*(ii-1)
        FF(ii)=dreal(CompWigner(Jtemp,mtemp,0,XX(1),XX(2),0.0d0))
     #        *expm*dsin(XX(2))
      enddo
cdeb	write(*,*) CompWigner(Jtemp,mtemp,0,0.0,XX(2),XX(1))

      end subroutine
