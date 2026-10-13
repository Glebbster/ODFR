C for ODFR 01/12/2014
c
	SUBROUTINE  ODFWrite(FileName)
c    
       IMPLICIT none
      include 'dimnls.inc'
      include 'eprprm.inc'
      include 'parcom.inc'
      include 'rndoff.inc'
      include 'parODF.inc'
      include 'mpi_state.inc'
c
      integer mord
      integer isi,ntemp,ntemp1,jord,mtemp,i,j
      double precision ordP1(Jmax,Jmax),ordP2(Jmax,Jmax)
      character*30 FileName, Astring 
      character*1  NamTmp      
c
      DOUBLE PRECISION fparms(NFPRM,MXSITE,MXSPC),Cptl,fNrm
      common /FunOrdPar/CPtl(5),Jord 
c
      integer iparms(NIPRM,MXSITE,MXSPC)
      COMMON /COMPARMS/ iparms,fparms
c ###################################################################
c
      if (mpi_rank.ne.0) return
	OPEN(UNIT=3,FILE=FileName)
	write(3,114) ODFcontrol  
 114  FORMAT (2x,I2,2x, '  1 - axial ODF, 2 -frozen ODF,
     # 3 - orthorhombic to nonrhombic, 4 -reverse ')
      write(3,101)PenFLG,PenLev,knet
 101  FORMAT (12x,I2,2x, G12.4,2x, i5,' PenFlg, Penalty coef.,',1x,
     # 'Net for ODF integration')
	Astring='ODF with neg.part   '
      do isi=1,nnsite(1)
c 
        write(3,*)Astring, anegODF(isi)
        write(3,*)Jmx(isi),keyODF(isi), '   (order, keyODF)'
c orthorhombic      
        if (ODFtype(1).le.2)  write(3,103)
 103  FORMAT(1x,'An0',10x,'Fl',1x,'An2',10x,'Fl',1x,'An4',10x,'Fl',1x,
     #'An6',10x,'Fl',1x,'An8',10x,'Fl')
c nonorthorhombic      
        if (ODFtype(1).eq.3) write(3,113)
 113  FORMAT(1x,'An0',10x,'Fl',1x,'An1',10x,'Fl',1x,'An2',10x,'Fl',1x,
     #'An3',10x,'Fl',1x,'An4',10x,'Fl')    
        ntemp1=0
        do jord=0,Jmx(isi),2
             if (ODFtype(1).le.2)then
               ntemp=1+jord/2
             else
               ntemp=1+jord
             endif  
             write(3,104) (paramODF(ntemp1+i,isi), 
     *        iflagODF(ntemp1+i,isi),i=1,ntemp)
 104  FORMAT(5(1x,G12.5,1x,I2))
             ntemp1=ntemp1+ntemp
        enddo !jord
      enddo	 !isi
      write (3,1062)' '
1062  Format(A1)
c     --------------------------------------------------------
c     Calculate order parameters for macroscopic ODF
c     --------------------------------------------------------
      do isi=1,nnsite(1)
c   For macro ODF 
  	  if ((ODFtype(1).ge.2).and.(keyODF(isi).ge.0)) then 
           write (3,1061) isi
 1061 format(1x,'Order parameters for macro ODF,site',i2,':')
           if (ODFtype(1).eq.2)then
               write (3,1067) ODFtype(1)
 1067 format(1x,'ODFtype: ',i2,' j,m - even')  
           else 
               write (3,1068) ODFtype(1)
 1068 format(1x,'ODFtype: ',i2,' j - even, m - even and odd' )
           endif
           ntemp1=0
           do jord=0,Jmx(isi),2
             if (ODFtype(1).le.2)then
               ntemp=1+jord/2
             else               
               ntemp=1+jord
             endif  
             do i=1,ntemp  	     
               ordP2(i,ntemp)=paramODF(ntemp1+i,isi)/(2.d0*jord+1.d0)
             enddo
             write(3,105) jord,(ordP2(i,ntemp),i=1,ntemp)
 105  FORMAT(1x,i2,1x,6(G12.5,1x))	
             ntemp1=ntemp1+ntemp
           enddo !jord
        endif	!ODFtype
c     --------------------------------------------------------
c     Calculate order parameters for potential
c     --------------------------------------------------------
        mord=10
  	  if (iparms(IIPT,isi,1) .gt. 0) then 
            do j=1,5
              cptl(j)=fparms(IC20-1+j,isi,1)
            enddo
            write (3,1060) isi
 1060 format(1x,'Order parameters for potential',i2,':')
            write (3,1070) 
 1070 format(1x,' j,n - even') 
            Jord=0
		    ntemp=1		   		  	  
            call OrdParPot(kNet,ntemp,1.d0,ordP1(1,1))
            fNrm=ordP1(1,1)
            ordP1(1,1)= ordP1(1,1)/fNrm		    
            write(3,105) 0,ordP1(1,1)
            mord=max(iparms(ILEMX,isi,1),Jmx(isi))
            do jord=2,mord,2
              ntemp=1+jord/2		   		  	  
              call OrdParPot(kNet,ntemp,fNrm,ordP1(1,ntemp))
              write(3,1051) jord,(ordP1(i,ntemp),i=1,ntemp)
 1051  FORMAT(1x,i2,1x,12(G12.5,1x))              
            enddo !jord  
        endif 	!if ((ODFtype(1).le.2)then
c     --------------------------------------------------------
c     Calculate order parameters for macroscopic ODF 
c     defined by dlam potential
c     --------------------------------------------------------        
  	  if ((ODFtype(1).eq.1).or.(fparms(IDC20,isi,1) .gt. 0)) then 
              cptl(1)=fparms(IDC20,isi,1)
            do j=2,5
              cptl(j)=0
            enddo
            write (3,1063) isi
 1063 format(1x,'Order par-s for macroscopic potential dlam',i2,':')
                write (3,1069) 
 1069 format(1x,' j,m - even') 

            Jord=0
		    ntemp=1		   		  	  
            call OrdParPot(kNet,ntemp,1.d0,ordP2(1,1))
            fNrm=ordP2(1,1)
            do jord=2,mord,2
              ntemp=1+jord/2		   		  	  
              call OrdParPot(kNet,ntemp,fNrm,ordP2(1,ntemp))
              write(3,105) jord,(ordP2(i,ntemp),i=1,ntemp)
            enddo !jord
            ordP2(1,1)= ordP2(1,1)/fNrm		    
        endif 	!ODFtype.eq.1
c     --------------------------------------------------------
c     Output global order parameters for paramagnetic centers
c     Ajmn=(Aext_jm0)*(Aint_j0n)=ordP2(m)*ordP1(n)
c     --------------------------------------------------------        
  	  if ((ODFtype(1).ge.1).or.((iparms(IIPT,isi,1) .gt. 0)
     #                    .and.(fparms(IDC20,isi,1) .gt. 0))) then 
          write (3,1064) isi
 1064 format(1x,'Order par-s for paramagnetic centers',i2,':')
c
        if (ODFtype(1).eq.2)then 
          write (3,*)'   macro - Pjm0=(Pext_jm0)*(Pint_j00)' 
          write (3,*)'   micro - Pj0n=(Pext_j00)*(Pint_j0n)' 
          write (3,*)'   mixed - Pjmn=(Pext_jm0)*(Pint_j0n)'
          write (3,1065) ODFtype(1)
 1065 format(1x,'ODFtype: ',i2,' j,m,n - even')
c           
          write(3,107) 0,ordP2(1,1)*ordP1(1,1)        
          do jord=2,Jmx(isi),2
              ntemp=1+jord/2		   		  	  
            if (iparms(IIPT,isi,1) .gt. 1) then
         write(3,107) jord,(ordP2(i,ntemp)*ordP1(1,ntemp),i=1,ntemp)
  107  FORMAT(1x,'macro ',i2,10(G12.5,1x))
         write(3,108) jord,(ordP2(1,ntemp)*ordP1(i,ntemp),i=2,ntemp)  	
  108  FORMAT(1x,'micro ',i2,10(G12.5,1x))
               do mtemp=2,ntemp
         write(3,109) jord,(ordP2(mtemp,ntemp)*ordP1(i,ntemp),i=2,ntemp)
  109  FORMAT(1x,'mixed ',i2,10(G12.5,1x))
               enddo ! mtemp
            else
c macro:            
         write(3,107) jord,(ordP1(1,ntemp)*ordP2(i,ntemp),i=1,ntemp)
            endif   !iparms(IIPT,isi,1)               
          enddo !jord 
        endif   !ODFtype(1).eq.2
c        
        if (ODFtype(1).eq.3)then 
          write (3,*)'   macro - Pjm0=(Pext_jm0)*(Pint_j00)' 
          write (3,*)'   micro - Pj0n=(Pext_j00)*(Pint_j0n)' 
          write (3,*)'   mixed - Pjmn=(Pext_jm0)*(Pint_j0n)'
          write (3,1066) ODFtype(1)
 1066 format(1x,'ODFtype: ',i2,' j,n - even, m - even and odd')
c           
          write(3,107) 0,ordP2(1,1)*ordP1(1,1)        
          do jord=2,Jmx(isi),2
              ntemp=1+jord/2
              ntemp1=1+jord		   		  	  
            if (iparms(IIPT,isi,1) .gt. 1) then
c macro:                
        write(3,107) jord,(ordP2(i,ntemp1)*ordP1(1,ntemp),i=1,ntemp1)
c micro:         
         write(3,108) jord,(ordP2(1,ntemp1)*ordP1(i,ntemp),i=2,ntemp)
c mixed:           	
               do mtemp=2,ntemp1
         write(3,109) jord,(ordP2(mtemp,ntemp1)*ordP1(i,ntemp),
     #                     i=1,ntemp)
               enddo ! mtemp
            else
c macro            
         write(3,107) jord,(ordP2(i,ntemp1)*ordP1(1,ntemp),i=1,ntemp1)
            endif   !iparms(IIPT,isi,1)               
          enddo !jord 
        endif   !ODFtype(1).eq.3 
c  
        endif 	!ODFtype.ge.1
c
c ODF for graph  
  	if ((ODFtype(1).ge.1).or.((iparms(IIPT,isi,1) .gt. 0)
     #                    .and.(fparms(IDC20,isi,1) .gt. 0))) then 
c Macro
        if (ODFtype(1).le.2)then 
  	    if (keyODF(isi).eq.0) then 
            write(NamTmp,106)isi
106         FORMAT (I1)
            OPEN(UNIT=4,FILE='G'//NamTmp//FileName)
            write(4,110) ordP2(1,1)*ordP1(1,1)            
              do jord=2,Jmx(isi),2
                   ntemp=1+jord/2
                write(4,110) (ordP2(i,ntemp)*ordP1(1,ntemp),i=1,ntemp)
  110  FORMAT(1x,10(G12.5,1x))
              enddo !jord   
c Potential
            if (iparms(IIPT,isi,1) .gt. 1) then
              write(4,111) mord
 111          FORMAT (1x, 'For potential, ord = ',I3)
              write(4,110) ordP1(1,1)            
              do jord=2,mord,2
                ntemp=1+jord/2		   		  	  
                write(4,110) (ordP1(i,ntemp),i=1,ntemp)  	
              enddo !jord                          
c Micro
              write(4,110) ' Micro '
              write(4,110) ordP1(1,1)*ordP2(1,1)            
              do jord=2,Jmx(isi),2
                ntemp=1+jord/2		   		  	  
                write(4,110) (ordP2(1,ntemp)*ordP1(i,ntemp),i=1,ntemp) 
              enddo !jord 
c Mixed 
              write(4,110) ' Mixed '             
              do jord=2,Jmx(isi),2
                ntemp=1+jord/2		   		  	  
                do mtemp=2,ntemp
               write(4,110)(ordP2(mtemp,ntemp)*ordP1(i,ntemp),i=2,ntemp)
                enddo ! mtemp
              enddo !jord                            
           endif   !iparms(IIPT,isi,1)                            
            CLOSE (UNIT=4)
          endif !keyODF(isi)  	
      endif    !ODFtype(1).le.2  
      if (ODFtype(1).eq.3)then
c Macro      
  	    if (keyODF(isi).eq.0) then 
            write(NamTmp,106)isi
            OPEN(UNIT=4,FILE='G'//NamTmp//FileName)
            write(4,110) ordP2(1,1)*ordP1(1,1) 
            do jord=2,Jmx(isi),2
                   ntemp=1+jord/2
                   ntemp1=1+jord	
        write(4,110) (ordP2(i,ntemp1)*ordP1(1,ntemp),i=1,ntemp1)
            enddo !jord   
c Potential
            if (iparms(IIPT,isi,1) .gt. 1) then
              write(4,111) mord
              write(4,110) ordP1(1,1)            
              do jord=2,mord,2
                ntemp=1+jord/2		   		  	  
                write(4,110) (ordP1(i,ntemp),i=1,ntemp)  	
              enddo !jord                          
c Micro
              write(4,110) ' Micro '
              write(4,110) ordP1(1,1)*ordP2(1,1)            
              do jord=2,Jmx(isi),2
                ntemp=1+jord/2
                ntemp1=1+jord		   		  	  
        write(4,110)(ordP2(1,ntemp1)*ordP1(i,ntemp),i=1,ntemp)        
              enddo !jord 
c Mixed 
              write(4,110) ' Mixed '             
              do jord=2,Jmx(isi),2
                ntemp=1+jord/2
                ntemp1=1+jord               
               do mtemp=2,ntemp1
         write(4,110) (ordP2(mtemp,ntemp1)*ordP1(i,ntemp),
     #                     i=1,ntemp)
               enddo ! mtemp                
              enddo !jord                            
           endif   !iparms(IIPT,isi,1)                            
            CLOSE (UNIT=4)
          endif !keyODF(isi)  	
      endif    !ODFtype(1).le.3         
      endif    !ODFtype(1).ge.1)....
      enddo  !isi
  	CLOSE (UNIT=3)          	
c
	return
	end


		





