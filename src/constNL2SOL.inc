C for ODFR3 18/01/2016
c*********************************************************************
c	  
      parameter   ( NPD   =  MXPT )    ! Max. length of the vector function.
      parameter   ( NVD   =  MXVAR )   ! Max. number of variables.
      parameter   ( MVD   =  93 + NPD*NVD + 3*NPD + NVD*(3*NVD + 33) )
      parameter   ( MIVD  =  60 + NVD )
      DOUBLE PRECISION   uRparm(1),  ! user real parameters passed through NL2SOL.
     &                   v(MVD)      ! NL2SOL working variable.                  
      integer            iv(MIVD),   ! NL2SOL working variable.
     &                   uIparm(1)   ! user integer parameters passed through NL2SOL.

c*********************************************************************


