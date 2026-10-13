cAXV, for ODFR4 
        SUBROUTINE choiceExp(NSET,Ntoch,FieldExp,SpecExp,
     #  ArrField,chExp)
c       This subroutine choose the NSET points from experimental spctrum
c       using linear interpolation.
c       This set of points wil be used for fitting
       IMPLICIT DOUBLE PRECISION (A-H),DOUBLE PRECISION (O-Z)
C
       include 'dimnls.inc'
      integer Naver,keyiparm 
       common /Naver/Naver,keyiparm(MXSPC)

       DIMENSION ArrField(MXSPT), chExp(MXSPT)
	 DIMENSION FieldExp(MXINP), SpecExp(MXINP)
c  take  the values 0 for points, where experimental
c  spectrum is absent
	     Ifirst=1
        DO WHILE (ArrField(Ifirst).LT.FieldExp(Naver))
	  chExp(Ifirst)=0.0d0
	  Ifirst=Ifirst+1
	  ENDDO
	     Ilast=NSET
        DO WHILE (ArrField(Ilast).GT.FieldExp(Ntoch-Naver))
	  chExp(Ilast)=0.0d0
	  Ilast=Ilast-1
	  ENDDO
c take the linear interpolated values as the experimental 
c points at the need fields
        DO 10 I=Ifirst,Ilast
        FIELD=ArrField(I)   
          NPLEFT=0
          NPRIGHT=0
          J=1
            DO WHILE (FieldExp(J).LT.FIELD)
                  NPLEFT=J
                  J=J+1
            END DO
            NPRIGHT=J 
            chExp(i)=0.D0 
            SUMX=0.D0
            SUMY=0.D0
            SUMX2=0.D0
            SUMXY=0.D0
            DO J=NPLEFT-Naver+1,NPRIGHT+Naver-1
               SUMX=SUMX+FieldExp(J)
               SUMY=SUMY+SpecExp(J)
               SUMX2=SUMX2+FieldExp(J)**2
               SUMXY=SUMXY+FieldExp(J)*SpecExp(J)
            END DO
            DIV=2.D0*Naver*SUMX2-SUMX*SUMX
            SLOPE=(2.D0*Naver*SUMXY-SUMX*SUMY)/DIV
            YINTR=(SUMY*SUMX2-SUMX*SUMXY)/DIV
            chExp(i)=YINTR+SLOPE*FIELD 
c            write(*,*)field,spe,FieldExp(npleft),FieldExp(npright)
c            write(*,*)specexp(npleft),specexp(npright)
c            pause
  10  CONTINUE
       RETURN
       END
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
