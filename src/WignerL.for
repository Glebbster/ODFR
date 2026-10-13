c For ODFR, Aleksey Bogdanov: 
c         This module contains subroutines and declarations
c         related to decomposition of the ODF by Wigner
c         matrices.
c
      subroutine ODFLInit(jsite,jODF,A_jmn)
      implicit none
c  This subroutine initializes matrix A_jmn and 
cAXV      fill A_jmn from array paramODF
      include 'dimnls.inc'
      include 'dimodfl.inc'
      include 'parcom.inc'
      include 'parODF.inc'
c     	 
      integer j,m,n,isi
      integer jord,ntemp,ntemp1,ii,jsite,jODF
      logical A_null
      common /WigODFIsNull/A_null
      dimension A_null(0:Jmax,-Jmax:Jmax,-Jmax:Jmax)
c/////////////////////////////////////////////////////
c
      isi=1
      jODF=0
      jord=0
      ntemp=0
      ntemp1=0
      do j=0,Jmax
        do m=-Jmax,Jmax
          do n=-Jmax,Jmax
            A_jmn(j,m,n)=(0.0d0, 0.0d0)
            A_null(j,m,n)=.true.
          end do
        end do
      end do
      if ((jsite.lt.1).or.(jsite.gt.MXSITE)) return
      select case (keyODF(jsite))
      case (:-1)
        return
      case (0)      
        isi=jsite
      case (1:)
        isi=1
      end select 
c                          
        jODF=Jmx(isi)
        if ((jODF.lt.0).or.(jODF.gt.Jmax)) return
cAXV:
      if (ODFtype(1).eq.2) then
        ntemp1=0
        do jord=0,jODF,2
           do ii=0,jord,2
              ntemp=1+ii/2
              if (ntemp1+ntemp.gt.nprmODF(isi)) return
              A_jmn(jord, ii,0)=paramODF(ntemp1+ntemp,isi)
           enddo		
           ntemp1=ntemp1+1+jord/2
        enddo !jord
      endif !(ODFtype(1).eq.2)
      if (ODFtype(1).eq.3) then
        ntemp1=0
        do jord=0,jODF,2
           do ii=0,jord
              ntemp=1+ii
              if (ntemp1+ntemp.gt.nprmODF(isi)) return
              A_jmn(jord, ii,0)=paramODF(ntemp1+ntemp,isi)
           enddo		
           ntemp1=ntemp1+1+jord
        enddo !jord
      endif !(ODFtype(1).eq.3)
c             
        call MakeODFLReal(A_jmn, jODF)  ! make ODF real
        call ODFL_FindNull(A_jmn, jODF)	
c	     
	return 
 	end subroutine
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	double complex function CompWigner(j, m, n, alpha, beta, gamma)
      implicit none
c  This function evaluates complex Wigner function D_jmn(alpha, beta, gamma)
c
	  integer j, m, n
        double precision alpha, beta, gamma
        double complex result
        external Wig
        double precision Wig
c/////////////////////////////////////////////////////
c
	  result=(0.0d0, -1.0d0)*(real(m,8)*alpha+real(n,8)*gamma)
	  CompWigner=cdexp(result)*Wig(j,m,n,beta)
c	write(*,*) result
c	write(*,*)'stop2'
c	  
	end function
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	        subroutine ODFL_FindNull(A_jmn, jODF)
              implicit none
	        include 'dimnls.inc'
	        include 'dimodfl.inc'
	        include 'rndoff.inc'
c
               integer jODF,j,m,n
	         common /WigODFIsNull/A_null
	         logical,dimension(0:Jmax,-Jmax:Jmax,-Jmax:Jmax)::A_null
c        
         do j=0,jODF,2
           do m=-j,j
             do n=-j,j
               if (cdabs(A_jmn(j,m,n))<RNDOFF) then
                  A_null(j,m,n)=.TRUE.
               else
                  A_null(j,m,n)=.FALSE.
               end if
             end do !n
           end do !m
         end do !j
c         
        end subroutine
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	double precision function 
     *               ODFL_RoReFast(alpha, beta, gamma)
c
c         This function returns the real part of ODF specified by
c         coefficients A_jmn at Euler angles alpha, beta, gamma
c         The summation is proceeds only by Wigner function with non-null
c         order parameters (A_null/=0), this enhances the speed of calculation
c
        implicit none
        include 'dimnls.inc'
	  include 'dimodfl.inc'
c
	  external CompWigner
	  double complex CompWigner
	  double precision alpha, beta, gamma
	  integer jODF
	  integer j,m,n,numODFs
	  double precision RoTemp
c
       common /A_jmn/A_jmn,jODF,numODFs
	 common /WigODFIsNull/A_null
	 logical, dimension(0:Jmax,-Jmax:Jmax,-Jmax:Jmax)::A_null 
c
c       double precision Wig
c	  external Wig
c/////////////////////////////////////////////////////
c
c 	  write (*,*) A_jmn(2,0,0)
c
		  ODFL_RoReFast=0.0d0
		  RoTemp=0.0d0
	  do j=0,jODF,2
	 	do m=-j, j
	      do n=-j, j
	       if(.NOT.A_null(j,m,n)) then
		    RoTemp=RoTemp+dreal(A_jmn(j,m,n)*
     *                 CompWigner(j,m,n,alpha,beta,gamma))
	       end if
	      end do
		end do
	  end do
        ODFL_RoReFast=RoTemp
	  return
	end function
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	double precision function 
     *       ODFL_RoReFast_Integrated
     *            (alphaSF, betaSF, gammaSF, psiFr)
c
c         This function returns the real part of ODF specified by
c         coefficients A_jmn, the sample is rotated in the resonator according to
c         Euler angles alphaSF, betaSF, gammaSF. The number density of particles
c         aligned at the angle psiFr to the magnetic field, are evaluated
c         The summation is proceeds only by Wigner function with non-null
c         order parameters (A_null/=0), this enhances the speed of calculation
c
        implicit none
        include 'dimnls.inc'
	  include 'dimodfl.inc'
c
	  external CompWigner
	  double complex CompWigner
	  double precision alphaSF, betaSF, gammaSF, psiFr
	  integer jODF
		  integer j,p,numODFs
cAXV	  integer m,n
c
		  double precision RoTemp
c
	 common /WigODFIsNull/A_null
	 logical, dimension(0:Jmax,-Jmax:Jmax,-Jmax:Jmax)::A_null 
cAXV
       common /A_jmn/A_jmn,jODF,numODFs
       double precision Wig
       external Wig
c/////////////////////////////////////////////////////
c
c 	  write (*,*) 'A200=',A_jmn(2,0,0),jODF
c
		  ODFL_RoReFast_Integrated=0.0d0
		  RoTemp=0.0d0
	  do j=0,jODF,2
c	    do m=-j, j
	      do p=-j, j
	       if(.NOT.A_null(j,p,0)) then
		    RoTemp=RoTemp+dreal(A_jmn(j,p,0)*
     *                 CompWigner(j,p,0,alphaSF,betaSF,gammaSF))*
     *                 Wig(j,0,0,psiFr)
		       end if
	      end do
c		end do
		  end do
	        ODFL_RoReFast_Integrated=RoTemp
c
		  return
		end function
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
      subroutine MakeODFLReal(A_jmn, jODF)
c
c       This function makes changes the A_jmn matrix so that the corresponding
c       ODF would be real. To escape conflicts, the user should specify only
c	  coefficients A_jmn with non-negative m, then make the function real.
       include 'dimnls.inc'
       include 'dimodfl.inc'
c
c/////////////////////////////////////////////////////
	 do j=0,jODF
	   do m=-j,0
	     do n=-j,j
	       if(mod(m-n, 2)==0) then
	          A_jmn(j,m,n)=dconjg(A_jmn(j,-m,-n))
	       else
	          A_jmn(j,m,n)=dconjg(-A_jmn(j,-m,-n))
	       end if
	       if ((m==0) .and. (n==0)) then
	          A_jmn(j,m,n)=dcmplx(dreal(A_jmn(j,m,n)), 0.0d0)
	       end if
	     end do
	   end do
	 end do
c
	end subroutine
c
c%%%%%%%%%%%%%%%%% END MODULE WignerL %%%%%%%%%%%%%%%%%%%
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
