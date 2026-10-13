 
      SUBROUTINE librat(keylib,parlib,magin,magres)
cAXV This subprogram calculates the magnetic parameters (components of g- and A tensors)
c    overaged by high-frequent rotational vibrations around equilibrium orientation (librations)
      IMPLICIT NONE
       integer keylib
       double precision parlib(6),magin(6),magres(6)
c
	    if (keylib.eq.1)then
c low angular  oscillations around g-axes
     		call librSin3(parlib,magin,magres)
	    endif
          if (keylib.eq.2) then
c unlimited oscillation around g-axes (analog Griffits)
c uniform distribution into interval -alfa - +alfa
 		call Grif3(parlib,magin,magres)
	    endif
          if (keylib.eq.3) then
c unlimited oscillation around arbitrary axes(analog Timofeev)
 		call Tim3(parlib,magin,magres)
	    endif
          if (keylib.eq.4) then
c Our formulas as in PolymerPhys2009_47_107 and Nitroxides
 		call LibrFormul(parlib,magin,magres)
	    endif	    
	return
	end
c##############################################################################
      SUBROUTINE librSin3(par,magin,magres)
c   
       IMPLICIT NONE
       include 'pidef.inc'       
       
	 double precision par(6),magin(6),magres(6)
       double precision DHin(3),Slibx,Sliby,Slibz,gxxin,gyyin,gzzin

	Slibx=par(1)*RADIAN
	Sliby=par(2)*RADIAN
	Slibz=par(3)*RADIAN
	gxxin=magin(1)
	gyyin=magin(2)
	gzzin=magin(3)
	        DHin(1)=magin(4)
              DHin(2)=magin(5)
              DHin(3)=magin(6)
c
	magres(1)=gxxin+(gzzin-gxxin)*Sliby+(gyyin-gxxin)*Slibz
	magres(2)=gyyin+(gzzin-gyyin)*Slibx+(gxxin-gyyin)*Slibz
	magres(3)=gzzin+(gyyin-gzzin)*Slibx+(gxxin-gzzin)*Sliby	
c
	magres(4)=DHin(1)+(DHin(3)-DHin(1))*Sliby
     #        +(DHin(2)-DHin(1))*Slibz
	magres(5)=DHin(2)+(DHin(3)-DHin(2))*Slibx
     #        +(DHin(1)-DHin(2))*Slibz
	magres(6)=DHin(3)+(DHin(2)-DHin(3))*Slibx
     #        +(DHin(1)-DHin(3))*Sliby	
c
	 return
	end
C ##################################################
        SUBROUTINE Tim3(parlib,magin,magres)
	  IMPLICIT DOUBLE PRECISION (A-H),DOUBLE PRECISION (O-Z)
	 double precision parlib(6),magin(6),magres(6),prom(6)
       dimension U(3,3)
       include 'pidef.inc'
c unlimited oscillation around g-axes (analog Griffits)
c uniform distribution into interval -alfa - +alfa
c matrix is similar to ESRmod2
c fi and psi are exchanged !!!
        fi=parlib(4)*RADIAN
        te=parlib(5)*RADIAN
        psi=parlib(6)*RADIAN
c  as Gidomirov's book, p.71, fi and psi are exchanged, as 
c     notebook "Photoselection I", p. 0-1
       u(1,1) = dcos(psi)*dcos(fi)-dcos(te)*dsin(psi)*dsin(fi)
       u(1,2) = -dcos(psi)*dsin(fi)-dcos(te)*dsin(psi)*dcos(fi)
       u(1,3) = dsin(te)*dsin(psi)
       u(2,1) = dsin(psi)*dcos(fi)+dcos(te)*dcos(psi)*dsin(fi)
       u(2,2) = -dsin(psi)*dsin(fi)+dcos(te)*dcos(psi)*dcos(fi)
       u(2,3) = -dsin(te)*dcos(psi)
       u(3,1) = dsin(te)*dsin(fi)
       u(3,2) = dsin(te)*dcos(fi)
       u(3,3) = dcos(te) 
c
c  around X axis of librational sistem
	ALib2Alf=parlib(1)*RADIAN
	if (ALib2Alf.ne.0.d0)then
c equations (19) from Timofeev's articles
	  Q1=dsin(2.d0*ALib2Alf)/2.d0/ALib2Alf
  	  Q1= ( 1.d0 + Q1) /2.d0
  	  Q2= DSin(ALib2Alf) / ALib2Alf - Q1
        Q3= 1.d0 - DSin(ALib2Alf) /ALib2Alf - Q2
        Q4= DSin( 2.d0 *ALib2Alf) / 2.d0 /ALib2Alf
        Q4= ( 1.d0 - Q4) / 2.d0
c equations (14) from Timofeev's articles
       SXR=(3.d0*u(1,1)*u(1,1)-1.d0)/2.d0       
	 SYR=(3.d0*u(1,2)*u(1,2)-1.d0)/2.d0       
       SZR=(3.d0*u(1,3)*u(1,3)-1.d0)/2.d0  
       call tim1(magin(1),magin(2),magin(3),
     *           prom(1),prom(2),prom(3),
     *           SXR,SYR,SZR,Q1,Q2,Q3,Q4)      
       call tim1(magin(4),magin(5),magin(6),
     *           prom(4),prom(5),prom(6),
     *               SXR,SYR,SZR,Q1,Q2,Q3,Q4)      
	else
         prom(1)=magin(1)
         prom(2)=magin(2)
         prom(3)=magin(3)
           prom(4)=magin(4)
           prom(5)=magin(5)
           prom(6)=magin(6)
	endif

       magres(1)=prom(1)
       magres(2)=prom(2)
       magres(3)=prom(3)
         magres(4)=prom(4)
         magres(5)=prom(5)
         magres(6)=prom(6)

c  around Y axis of librational sistem
	ALib2Alf=parlib(2)*RADIAN
	if (ALib2Alf.ne.0.d0)then
	  Q1=dsin(2.d0*ALib2Alf)/2.d0/ALib2Alf
  	  Q1= ( 1.d0 + Q1) /2.d0
  	  Q2= DSin(ALib2Alf) / ALib2Alf - Q1
        Q3= 1.d0 - DSin(ALib2Alf) /ALib2Alf - Q2
        Q4= DSin( 2.d0 *ALib2Alf) / 2.d0 /ALib2Alf
        Q4= ( 1.d0 - Q4) / 2.d0
       SXR=(3.d0*u(2,1)*u(2,1)-1.d0)/2.d0       
	 SYR=(3.d0*u(2,2)*u(2,2)-1.d0)/2.d0       
       SZR=(3.d0*u(2,3)*u(2,3)-1.d0)/2.d0  
       call tim1(prom(1),prom(2),prom(3),
     *           magres(1),magres(2),magres(3),
     *           SXR,SYR,SZR,Q1,Q2,Q3,Q4)
       call tim1(prom(4),prom(5),prom(6),
     *           magres(4),magres(5),magres(6),
     *           SXR,SYR,SZR,Q1,Q2,Q3,Q4)
       prom(1)=magres(1)
       prom(2)=magres(2)
       prom(3)=magres(3)
         prom(4)=magres(4)
         prom(5)=magres(5)
         prom(6)=magres(6)
	endif
c  around Z axis of librational sistem
	ALib2Alf=parlib(3)*RADIAN
	if (ALib2Alf.ne.0.d0)then
	  Q1=dsin(2.d0*ALib2Alf)/2.d0/ALib2Alf
  	  Q1= ( 1.d0 + Q1) /2.d0
  	  Q2= DSin(ALib2Alf) / ALib2Alf - Q1
        Q3= 1.d0 - DSin(ALib2Alf) /ALib2Alf - Q2
        Q4= DSin( 2.d0 *ALib2Alf) / 2.d0 /ALib2Alf
        Q4= ( 1.d0 - Q4) / 2.d0
       SXR=(3.d0*u(3,1)*u(3,1)-1.d0)/2.d0       
	 SYR=(3.d0*u(3,2)*u(3,2)-1.d0)/2.d0       
       SZR=(3.d0*u(3,3)*u(3,3)-1.d0)/2.d0  
       call tim1(prom(1),prom(2),prom(3),
     *           magres(1),magres(2),magres(3),
     *           SXR,SYR,SZR,Q1,Q2,Q3,Q4)
       call tim1(prom(4),prom(5),prom(6),
     *           magres(4),magres(5),magres(6),
     *           SXR,SYR,SZR,Q1,Q2,Q3,Q4)
	endif
	return
	end
C ###################################################
           SUBROUTINE tim1(valin1,valin2,valin3,
     *                       valres1,valres2,valres3,
     *                       SXR,SYR,SZR,Q1,Q2,Q3,Q4)  
          IMPLICIT DOUBLE PRECISION (A-H),DOUBLE PRECISION (O-Z)
c equations for oscillation around arbitrary axis 
c equations (16'), (17') and (18') from Timofeev's articles
	 parameter (dvetret=2.d0/3.d0)
  	QQ3=(valin1+valin2+valin3)/3.d0
	QQ3=QQ3+dvetret*(valin1*SXR+valin2*SYR+valin3*SZR)
	QQ3=QQ3*Q3/3.d0
	sl1=(Q1+dvetret*Q2*(2.d0*SXR+1.d0))*valin1
	sl2=QQ3*(2.d0*SXR+1.d0)
	sl3=Q4*((2.d0*SYR+1.d0)*valin3+(2.d0*SZR+1.d0)*valin2)/3.d0
	valres1=sl1+sl2+sl3
c
	sl1=(Q1+dvetret*Q2*(2.d0*SYR+1.d0))*valin2
	sl2=QQ3*(2.d0*SYR+1.d0)
	sl3=Q4*((2.d0*SXR+1.d0)*valin3+(2.d0*SZR+1.d0)*valin1)/3.d0
	valres2=sl1+sl2+sl3
c
	sl1=(Q1+dvetret*Q2*(2.d0*SZR+1.d0))*valin3
	sl2=QQ3*(2.d0*SZR+1.d0)
	sl3=Q4*((2.d0*SYR+1.d0)*valin1+(2.d0*SXR+1.d0)*valin2)/3.d0
	valres3=sl1+sl2+sl3
	return
	end
C ##################################################
      SUBROUTINE grif3(parlib,magin,magres)   
       IMPLICIT NONE
       include 'pidef.inc'
	 double precision parlib(6),magin(6),magres(6),prom(6)
       double precision ALib2Alf,P1
c
	ALib2Alf=parlib(1)*RADIAN*2.d0
	if (ALib2Alf.ne.0.d0)then
	   P1=dsin(ALib2Alf)/ALib2Alf
	   call grif1(P1,magin(2),magin(3),prom(2),prom(3))
	   prom(1)=magin(1)
	   call grif1(P1,magin(5),magin(6),prom(5),prom(6))
	   prom(4)=magin(4)
	else
         prom(1)=magin(1)
         prom(2)=magin(2)
         prom(3)=magin(3)
           prom(4)=magin(4)
           prom(5)=magin(5)
           prom(6)=magin(6)
	endif

       magres(1)=prom(1)
       magres(2)=prom(2)
       magres(3)=prom(3)
         magres(4)=prom(4)
         magres(5)=prom(5)
         magres(6)=prom(6)

	ALib2Alf=parlib(2)*RADIAN*2.d0
	if (ALib2Alf.ne.0.d0)then
	   P1=dsin(ALib2Alf)/ALib2Alf
	   call grif1(P1,prom(3),prom(1),magres(3),magres(1))
         call grif1(P1,prom(6),prom(4),magres(6),magres(4))
         prom(1)=magres(1)
         prom(3)=magres(3)
           prom(4)=magres(4)
           prom(6)=magres(6)
	endif
	ALib2Alf=parlib(3)*RADIAN*2.d0
	if (ALib2Alf.ne.0.d0)then
	   P1=dsin(ALib2Alf)/ALib2Alf
	   call grif1(P1,prom(1),prom(2),magres(1),magres(2))
	   call grif1(P1,prom(4),prom(5),magres(4),magres(5))
	endif
      return
	end
C ##################################################
           SUBROUTINE grif1(P1,valin2,valin3,valres2,valres3)  
          IMPLICIT DOUBLE PRECISION (A-H),DOUBLE PRECISION (O-Z)
c formula for oscillation around first axis 
c analog the Griffits equations  
	valres2=(valin2*(1.d0+P1)+valin3*(1.d0-P1))/2.D0
	valres3=(valin2*(1.d0-P1)+valin3*(1.d0+P1))/2.D0
c	write(*,*)valres2
	return
	end
c################################################################
        DOUBLE PRECISION FUNCTION ANGFUN(ANGLE)
       IMPLICIT DOUBLE PRECISION (A-H),DOUBLE PRECISION (O-Z)
	   ANGFUN=ANGLE
	   if (ANGLE.gt.360)THEN
	  PROM=DABS(ANGLE/360.D0)
	  ANGFUN=360.D0*(PROM-Aint(PROM))
	  endif
	  if (ANGLE.lt.-360)THEN
	  PROM=DABS(ANGLE)/360.D0
	  ANGFUN=-360.D0*(PROM-Aint(PROM))
	  endif
	  return
	  end	
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE LibrFormul(parlib,magin,magres)
c   
       IMPLICIT NONE
       include 'rndoff.inc'
       include 'pidef.inc'       
       
	 double precision parlib(6),magin(6),magres(6)
       double precision P(3),DHin(3),gxxin,gyyin,gzzin,Alib2
       integer ii

       do ii=1,3 
         ALib2=parlib(ii)*RADIAN*2.d0
         if (dabs(ALib2).gt.rndoff) then
            P(ii)=0.5d0*(1.d0-dsin(ALib2)/ALib2)
         else 
            P(ii)=0.d0
         endif
       enddo
                                     	   
	gxxin=magin(1)
	gyyin=magin(2)
	gzzin=magin(3)
	      DHin(1)=magin(4)
            DHin(2)=magin(5)
            DHin(3)=magin(6)
c
	magres(1)=gxxin+(gzzin-gxxin)*P(2)+(gyyin-gxxin)*P(3)
	magres(2)=gyyin+(gzzin-gyyin)*P(1)+(gxxin-gyyin)*P(3)
	magres(3)=gzzin+(gyyin-gzzin)*P(1)+(gxxin-gzzin)*P(2)	
c
	magres(4)=DHin(1)+(DHin(3)-DHin(1))*P(2)
     #        +(DHin(2)-DHin(1))*P(3)
	magres(5)=DHin(2)+(DHin(3)-DHin(2))*P(1)
     #        +(DHin(1)-DHin(2))*P(3)
	magres(6)=DHin(3)+(DHin(2)-DHin(3))*P(1)
     #        +(DHin(1)-DHin(3))*P(2)	
c
	 return
	end