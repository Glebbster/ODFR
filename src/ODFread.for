C for ODFR 01/12/2014
c
      SUBROUTINE  ODFread(FileNAME)
      implicit none
      include 'dimnls.inc'
      include 'parcom.inc'
      include 'parODF.inc'
c
      integer JTMP,ICON,js,ntemp,ntemp1,jord,jpar,ISTRING
      double precision tempparODF(Jmax)
      integer tempflODF(Jmax)
      character*30 FileName, Astring 
c
      ISTRING=0
      ICON=0
      js=0
      OPEN(UNIT=3,FILE=FileNAME,status='old',IOSTAT=ICON)
      if (ICON.NE.0) go to 36
      ISTRING=ISTRING+1
      read(3,*)ODFcontrol
      ISTRING=ISTRING+1
      read(3,*)PenFLG,PenLev,knet   
c      
c ODF parameters for sites
c orthorhombic and nonorthorhombic ODF:
      if ((ODFtype(1).eq.2).or.(ODFtype(1).eq.3)) then
        if (ODFcontrol.le.2) then
          do js=1,nnsite(1)
            ISTRING=ISTRING+1
            READ(3,*,ERR=36,IOSTAT=ICON) Astring
            ISTRING=ISTRING+1
            READ(3,*,ERR=36,IOSTAT=ICON)Jmx(js), keyODF(js)
            if ((Jmx(js).lt.0).or.(Jmx(js).gt.Jmax)) then
              ICON=1
              go to 36
            endif
            ISTRING=ISTRING+1
            READ(3,*,ERR=36,IOSTAT=ICON) Astring
            nprmODF(js)=0
            do jord=0,Jmx(js),2
c orthorhombic ODF (ODFtype(1).eq.2):
              if (ODFtype(1).eq.2) ntemp=1+jord/2
c nonorthorhombic ODF (ODFtype(1).eq.3):
              if (ODFtype(1).eq.3) ntemp=1+jord
              if (nprmODF(js)+ntemp.gt.MaxParODF) then
                ICON=1
                go to 36
              endif
              ISTRING=ISTRING+1
             READ(3,*,ERR=36,IOSTAT=ICON)(paramODF(nprmODF(js)+jpar,js),
     *                      iflagODF(nprmODF(js)+jpar,js),jpar=1,ntemp) 
              if (ODFcontrol.eq.1) then
                do jpar=2,ntemp
                  paramODF(nprmODF(js)+jpar,js)=0.d0
                  iflagODF(nprmODF(js)+jpar,js)=0
                enddo              
              endif      
              nprmODF(js)=nprmODF(js)+ntemp
            enddo !jord
            if (keyODF(js).gt.0)then
              Jmx(js)=Jmx(1)
              nprmODF(js)=nprmODF(1)
              do jpar=1,nprmODF(1)
                paramODF(jpar,js)=paramODF(jpar,1)
                iflagODF(jpar,js)=0
              enddo !jpar
            endif
          enddo	 !js
        endif !(ODFcontrol.le.2)
c Change ortorhombic ODF in nonorthorombic (ODFcontrol.eq.3)
c    and opposite change (ODFcontrol.eq.4)
        if ((ODFcontrol.eq.3).or.(ODFcontrol.eq.4)) then 
          do js=1,nnsite(1)
            ISTRING=ISTRING+1
            READ(3,*,ERR=36,IOSTAT=ICON) Astring
            ISTRING=ISTRING+1
            READ(3,*,ERR=36,IOSTAT=ICON)Jmx(js), keyODF(js)
            if ((Jmx(js).lt.0).or.(Jmx(js).gt.Jmax)) then
              ICON=1
              go to 36
            endif
            ISTRING=ISTRING+1
            READ(3,*,ERR=36,IOSTAT=ICON) Astring
            nprmODF(js)=0
            do jord=0,Jmx(js),2
              ntemp=1+jord/2
              ntemp1=1+jord
              if (((ODFcontrol.eq.3).and.(nprmODF(js)+ntemp1.gt.MaxParODF))
     #        .or.((ODFcontrol.eq.4).and.(nprmODF(js)+ntemp.gt.MaxParODF))) then
                ICON=1
                go to 36
              endif
              ISTRING=ISTRING+1
              if (ODFcontrol.eq.3) then              
                READ(3,*,ERR=36,IOSTAT=ICON)(tempparODF(jpar),
     *                      tempflODF(jpar),jpar=1,ntemp)
                paramODF(nprmODF(js)+1,js)=tempparODF(1)  
                iflagODF(nprmODF(js)+1,js)=tempflODF(1)     
                do jpar=2,ntemp
                  paramODF(nprmODF(js)+(jpar-1)*2,js)=0.0d0  
                  iflagODF(nprmODF(js)+(jpar-1)*2,js)=0
                  paramODF(nprmODF(js)+(jpar-1)*2+1,js)=tempparODF(jpar)
                  iflagODF(nprmODF(js)+(jpar-1)*2+1,js)=tempflODF(jpar)
                enddo
                nprmODF(js)=nprmODF(js)+ntemp1                
              endif !(ODFcontrol.eq.3)                 
              if (ODFcontrol.eq.4) then              
                READ(3,*,ERR=36,IOSTAT=ICON)(tempparODF(jpar),
     *                      tempflODF(jpar),jpar=1,ntemp1)
                paramODF(nprmODF(js)+1,js)=tempparODF(1)  
                iflagODF(nprmODF(js)+1,js)=tempflODF(1)     
                do jpar=2,ntemp
                  paramODF(nprmODF(js)+jpar,js)=tempparODF((jpar-1)*2+1)
                  iflagODF(nprmODF(js)+jpar,js)=tempflODF((jpar-1)*2+1)
                enddo
                nprmODF(js)=nprmODF(js)+ntemp
               endif !(ODFcontrol.eq.4)                      
c
            enddo !jord
            if (keyODF(js).gt.0)then
              Jmx(js)=Jmx(1)
              nprmODF(js)=nprmODF(1)
              do jpar=1,nprmODF(1)
                paramODF(jpar,js)=paramODF(jpar,1)
                iflagODF(jpar,js)=0
              enddo !jpar
            endif
          enddo	 !js
          ODFcontrol=0
        endif !(ODFcontrol.eq.3)or 4)                 
      endif  !ODFtype(1).eq.2 or 3
c
  36    if (ICON.NE.0) then
           write(*,*)'File ',FileNAME
           write(*,*)' - Error in line No',ISTRING,', js=',js  
           write(2,*)'File ',FileNAME
           write(2,*)' - Error in line No',ISTRING,', js=',js         
           if ((js.ge.1).and.(js.le.MXSITE)) then
             if ((Jmx(js).lt.0).or.(Jmx(js).gt.Jmax)) then
                write(*,*)' Jmx out of range [0,',Jmax,']'
                write(2,*)' Jmx out of range [0,',Jmax,']'
             endif
           endif
        end if
c
        CLOSE (UNIT=3)
c
        nprm4=0
        if ( ODFcontrol.lt.2) then 
          do js=1,nnsite(1)
            if (keyODF(js).eq.0) then     
              do jtmp=1,nprmODF(js)
                if (iflagODF(jtmp,js).ge.1) nprm4=nprm4+1
              end do
            end if            
          enddo
        endif
c
      return
      end

