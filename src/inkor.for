c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c                    
               SUBROUTINE INKOR(N,K,T,L,XL,EXTERN,R,EPS)

       IMPLICIT DOUBLE PRECISION (A-H),DOUBLE PRECISION (O-Z)
	 
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C         THE INSTRUCTION TO THE PROGRAM OF CONSTANT LIMITS MULTIPLE
C         INTEGRATION BY KOROBOV'S METHODE - 'INKOR'
C
C
C         N      - A NUMBER OF SIMULTANEOUSLY COMPUTING INTEGRALS
C
C         K      - AN INTEGRATION MULTIPLICITY
C
C         T(2)   - INITIAL AND FINAL NUMBERS OF INTEGRATION NET POINTS
C       (IF T(1)=0, OR T(2)=0 - THE PROGRAM FIXES THEM ITSELF)
C
C         L      - AN ORDER OF A PERIODIZATION
C       (IT IS RECOMMENDED TO GIVE L=1,2,3)
C
C         XL(K,2) - AN ARRAY OF LOWER AND UPPER LIMITS OF INTEGRATION
C
C         EXTERN - A NAME OF AN OUTER ARITHMETICAL BLOCK OF USER
C       IT MUST BE DESCRIBED IN THE OPERATOR 'EXTERNAL'
C       OF A CALLING PROGRAM AND HAVE THE FORM:
C       'SUBROUTINE EXTERN(N,K,X,F)',
C       WHERE X(K) IS AN ARRAY OF INTEGRATION VARIABLES       
C       AND F(N) IS AN ARRAY OF INTEGRANDS VALUES  
C  
C         R(N,2) - AN ARRAY OF RESULTS  
C  
C         EPS    - ACCURACY. IT SHOWS CONDITIONS OF WORK OF THE PROGRAM:
C          EPS=0 - INTEGRATION ON A FIXED NET (ON THE TWO FORMULAS)     
C       R(J,1) AND R(J,2) ARE IN THIS CASE VALUES OF THE J-TH 
C       INTEGRAL, EVALUATED WITH THE FIRST AND THE SECOND     
C       FORMULAS RESPECTIVELY
C          EPS>0 - INTEGRATION WITH THE ABSOLUTE ACCURACY 'EPS'         
C       R(J,1) AND R(J,2) ARE VALUES OF THE J-TH INTEGRAL THE 
C       FIRST AND ITS ABSOLUTE ACCURACY THE SECOND 
C          EPS<0 - INTEGRATION WITH THE RELATIVE ACCURACY 'ABS(EPS)'    
C       R(J,1) AND R(J,2) ARE VALUES OF THE J-TH INTEGRAL AND 
C       ITS RELATIVE ACCURACY
C  
C       **********
C        THE INFORMATION ABOUT ERRORS IN INPUT DATA (TOO HIGH MULTIPLICI
C        OR ACCURACY) WILL BE WRITTEN BY THIS PROGRAM TO FILE WITH LUN=3
C        IT'S CONVENIENT TO ASSIGN LUN=3 TO YOUR OUTPUT FILE.
C       **********
C
C  
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C  
      COMMON /FLAG/FLAG  
      COMMON /INDAT/PT(30),AT(30),LT
      DIMENSION T(2),XL(K,2)        
      DIMENSION S(20),S1(20),CSIG(11),CTAU(11) 
      DOUBLE PRECISION X(5000),F(5000),R(N,2)     
      DATA LMAX,RFLAG,CP3,CMAX/10,1.E6,10.,20./

      IPRINT=0
      CALL DATKOR(K,IPRINT)   
      KR=K    
      LP=L    
      P3=3.**KR          
      PMAX=T(2)          
      IF(PMAX) 1,2,300   
    1 PMAX=CMAX*P3       
      GO TO 300          
    2 PMAX=PT(LT)        
  300 CONTINUE
      Q=T(1)  
      IF(Q) 7,9,9        
    7 CONTINUE
      IF(KR-10) 8,8,9    
    8 CONTINUE
      Q=P3    
    9 CONTINUE
      P3=CP3*P3          
      IF(Q-PT(LT)) 4,4,3 
    3 Q=PT(LT)
    4 IF(LP-LMAX) 6,6,5  
    5 LP=LMAX 
    6 DO 999 KT=1,LT
      FLAG=0.       
      P=PT(KT)      
      A=AT(KT)      
      IF(Q-P) 10,10,999         
   10 Z=2.*P        
      LFLAG=2       
      IF(LP) 14,19,19           
   14 LP=0          
      IF(Q-P3) 16,15,15         
   15 LP=1          
   16 IF(KR-6) 17,17,19         
   17 LP=LP+1       
   19 CONTINUE      
      IF(LP) 20,20,11           
   11 LFLAG=1       
C     PERIOD        
      APR=1.        
      ALPHA=APR*FLOAT(LP)/P**(1./FLOAT(KR)) 
      IF(ALPHA-.3) 1004,1004,1003           
 1003 ALPHA=.3      
 1004 CONTINUE      
      ALPHA1=1./ALPHA           
      ALPHA2=ALPHA/2.           
      BETA=1.-ALPHA 
      BETA1=1./BETA 
      L1=LP+1       
      AK=0.         
      AL=L1         
      A2L=LP+L1     
      B2L=A2L*AL/(AL-1.)        
      CL=1.         
      C2L=1.        
      DO 12 I=1,L1  
      CTAU(I)=CL/A2L
      CSIG(I)=CTAU(I)/(A2L+1.)  
      AK=AK+1.      
      AL=AL-1.      
      A2L=A2L-1.    
      CL=CL*AL/AK   
      C2L=C2L*A2L/AK
   12 CONTINUE      
      B2L=B2L*C2L*BETA1         
      DO 13 I=1,L1  
      CTAU(I)=B2L*CTAU(I)       
      CSIG(I)=B2L*ALPHA*CSIG(I) 
   13 CONTINUE      
   20 A1=1.         
      DO 23 JK=1,KR 
      IF(dMOD(A1,2.d0)) 21,21,22  
   21 A1=P-A1       
   22 S1(JK)=A1     
      S(JK)=S1(JK)*2.           
      A1=dMOD(A1*A,P)           
   23 CONTINUE      
      DO 24 IJ=1,N  
      R(IJ,1)=0.    
   24 R(IJ,2)=0.    
      M=P           
      DO 100 I=1,2  
      DO 90 J=1,M   
      FACT=1.       
      DO 70 JK=1,KR 
      Y1=S1(JK)/Z   
      S1(JK)=S1(JK)+S(JK)-Z     
      IF(S1(JK)) 31,32,32       
   31 S1(JK)=S1(JK)+Z           
   32 GO TO (33,61),LFLAG       
   33 W=Y1          
      Y=W*ALPHA1    
      IF(W-ALPHA) 50,40,40      
   40 IF(W-BETA) 41,41,42       
   41 Y1=(W-ALPHA2)*BETA1       
      TAU=BETA1     
      GO TO 60      
   42 Y=ALPHA1-Y    
   50 Y1=-Y         
      Y2=Y          
      SIG=CSIG(1)   
      TAU=CTAU(1)   
      DO 51 IL=2,L1 
      Y2=Y2*Y       
      SIG=SIG*Y1+CSIG(IL)     
      TAU=TAU*Y1+CTAU(IL)     
   51 CONTINUE      
      Y1=SIG*Y2*Y   
      IF(Y1) 53,54,54         
   53 CONTINUE      
      Y1=0.         
   54 CONTINUE      
      TAU=TAU*Y2    
      IF(W-ALPHA) 60,52,52    
   52 Y1=1.d0-Y1      
   60 FACT=FACT*TAU 
   61 CONTINUE      
      Y=XL(JK,2)-XL(JK,1)     
      X(JK)=Y*Y1+XL(JK,1)     
      FACT=FACT*Y   
   70 CONTINUE   
   80 CALL EXTERN(N,K,X,F)    
      IF(FLAG) 777,81,999     
   81 DO 90 IJ=1,N  
      R(IJ,I)=R(IJ,I)+FACT*F(IJ)        
   90 CONTINUE      
      DO 91 JK=1,KR 
      S1(JK)=S(JK)  
   91 CONTINUE      
      DO 100 IJ=1,N 
      R(IJ,I)=R(IJ,I)/P       
  100 CONTINUE      
      KFLAG=1       
      EPS1=EPS      
      DO 200 I=1,N  
      IF(EPS1) 110,200,110    
  110 CONTINUE      
      R(I,1)=0.5*(R(I,1)+R(I,2))        
      R1=DABS(R(I,1)-R(I,2))  
      R(I,2)=R1     
      IF(EPS1) 120,140,140    
  120 R(I,2)=RFLAG  
      IF(R(I,1)) 130,140,130  
  130 R(I,2)=DABS(R1/R(I,1))  
  140 IF(R(I,2)-DABS(EPS1)) 200,160,160  
  160 KFLAG=2       
  200 CONTINUE      
      FLAG=P        
      GO TO (777,201),KFLAG   
  201   WRITE(*,*)80
        IF(P-PMAX) 80,900,900   
  999 CONTINUE      
  900 FLAG=0.
  777 Continue       
c  777 IF(FLAG) 800,778,800    
c  778 WRITE(5,780) EPS 
c  780 FORMAT(//' THE ACCURACY DEMANDED EPS=',E12.3,'HAS NOT GOTTEN'/)   
c      WRITE(5,781) P
c  781 FORMAT(/' MAKE T(2)> ' E12.1 ' OR CHANGE L')
c      WRITE(5,782) (R(i,2),i=1,N)  
c  782 FORMAT(/' THE ACTUAL VALUES OF ACCURACIES ARE: '//10E12.3//)
  800 RETURN        
      END 
      SUBROUTINE DATKOR(K,IPRINT)        
C         
       IMPLICIT DOUBLE PRECISION (A-H),DOUBLE PRECISION (O-Z)
      COMMON /INDAT/PT(30),AT(30),LT    
      DIMENSION KP(27),KP1(26),KA1(181),KA2(68),KA(249)     
      double precision KP,KP1,KA1,KA2,KA
      EQUIVALENCE (KA1(1),KA(1)),(KA2(1),KA(182)) 
C         
C         
C         
      DATA KP       
     P /23.,53.,101.,151.,307.,523.,829.,1259.,2129.,3001.,4001.,5003.,
     P 6007.,8191.,
     P  10007.,13001.,20011.,30011.,40009.,50021.,75011.,100003.,
     P 200003.,500009.,
     P  1000003.,2000003.,5000011./        
C         
      DATA KP1      
     P /3.,5.,8.,13.,21.,34.,55.,89.,144.,233.,377.,610.,987.,1597.,
     P 2584.,4181.,6765.,
     P  10946.,17711.,28657.,46368.,75025.,121393.,196418.,317811.,
     P 514229./
C         
      DATA KA1      
     3 /5.,17.,53.,20.,75.,78.,116.,433.,359.,276.,722.,1476.,
     3  592.,739.,544.,2135.,2759.,2023.,2979.,3483.,7623.,6697.,
     3  19477.,55279.,90347.,191421.,613511.,
     4  5.,31.,15.,22.,42.,178.,45.,483.,766.,174.,113.,792.,1351.,
     4 2488.,1206.,778.,
     4  2891.,2243.,2775.,3353.,7601.,10101.,19293.,55251.,
     4 122713.,191447.,492615.,
     5  21.,13.,127.,64.,146.,302.,618.,408.,1534.,840.,509.,1386.,
     5 198.,830.,2803.,
     5 2311.,2945.,3581.,7639.,10019.,19777.,61357.,122713.,
     5 191447.,492615.,
     6  109.,195.,118.,109.,200.,41.,233.,1751.,2037.,312.,1632.,
     6 2240.,1940.,3169.,
     6  2347.,2875.,3431.,7647.,10049.,19277.,33363.,66681.,
     6 191401.,613509.,
     7  27.,90.,227.,596.,170.,115.,229.,721.,2107.,2304.,2116.,
     7 1773.,2253.,2965.,
     7  3401.,7555.,6703.,24561.,33339.,66671.,191443.,613509., 
     8  32.,392.,86.,993.,1850.,2362.,1253.,3787.,436.,1415.,
     8 3295.,3267.,3093.,
     8  3455.,7809.,10029.,19827.,61365.,90343.,191437.,492261.,
     9  636.,108.,1579.,526.,1597.,2139.,470.,1437.,2533.,2077.,
     9 3119.,3469.,7813.,
     9  10099.,19461.,61365.,90343.,191415.,613507.,
     *  431.,1262.,1751.,564.,1415.,1715.,2065.,2901.,3413.,
     * 7677.,6699.,19153.,
     *  33341.,107369.,191435.,613515./         
C       
      DATA KA2      
     1 /586.,1409.,2205.,3237.,4325.,5419.,8063.,10151.,24573.,
     1  55269.,107393.,191389.,613515.,
     2  2191.,3319.,4319.,3407.,8099.,10151.,24573.,61365.,107377.,
     2  191389.,613507.,
     3  5407.,8061.,10761.,24573.,61371.,107375.,191433.,478471.,
     4  5407.,8099.,10143.,19363.,61371.,106245.,191433.,478471.,
     5  5417.,8061.,10171.,20015.,61371.,106265.,191405.,478471.,
     6  61373.,106257.,191405.,492613.,
     7  61353.,106245.,191421.,492613.,
     8  61353.,106257.,191421.,478471.,
     9  61353.,106245.,191421.,478471.,
     *  61353.,106233.,191421.,478471./         
C       
      DO 100 I=1,30 
      PT(I)=0.      
  100 AT(I)=0.      
C       
      IF(K-20) 2,2,1
    1 WRITE(3,101)  
  101 FORMAT(//'  BOTH K AND N MUST BE NOT GREATER THEN 20')        
      STOP          
C       
    2 IF(K-3) 3,5,5 
    3 DO 4 I=1,25   
      J=I+1         
      PT(I)=KP1(J)  
      AT(I)=KP1(I)  
      LT=25
    4 CONTINUE      
      GO TO 21      
    5 L=2           
      I=0           
      GO TO 9       
    6 LT=0          
    7 I=I+1         
      LT=LT+1       
      IF(KA(I)-400000.)7,8,8     
    8 L=L+1         
    9 IF(L-K) 6,10,10           
   10 I=I-LT        
      J=27-LT       
      DO 20 IJ=1,LT 
      I=I+1         
      J=J+1         
      PT(IJ)=KP(J)  
      AT(IJ)=KA(I)  
   20 CONTINUE      
C       
   21 IF(IPRINT) 22,23,22        
   22 WRITE(3,102) K ,          
     *(PT(I),I=1,15),(AT(I),I=1,15),(PT(I),I=16,30),(AT(I),I=16,30) 
  102 FORMAT(//40X,'THE PARAMETERS OF KOROBOV''S NETS   K=',        
     *I2/2(/1X,3H P=,15F8.0/1X,3H A=,15F8.0)////)       
C       
   23 RETURN        
      END           

C ##################################################
