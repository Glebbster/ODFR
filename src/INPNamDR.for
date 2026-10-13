C for ODFR 01/12/2014

      SUBROUTINE INPNamDR(nameNames,numSim,NSpec) 
c      SUBROUTINE read the names of files to be used 
c      NamSpecInp - experimental spectra
c      NamSpecRes - results of simulation
c      ParInpEPR - initial parameters
c      ParResEPR - result parameters 
c      ParInpODF - initial parameters of ODF
c      ParResODF - result parameters of ODF
*********************************************************************
c       numsim - number of consecutive simulations 
c      NSpec - array, numbers of spectra to be simulated
*********************************************************************
       include 'dimnls.inc'
       include 'name.inc'
       include 'parODF.inc'
       include 'mpi_state.inc'
c
        character*30 NameNames
        integer numSim, NSpec(MaxCALC)

       CHARACTER*3 tempStr
      
c          write(*,*) ' INPUT FILENAME FOR SPECTRA ' 
          OPEN ( UNIT=1,FILE=nameNames,status='old')         
       numSim=0
         READ(1,*)tempStr
       DO 38 WHILE (((LGT(tempStr,'end')).OR.
     *(LLT(tempStr,'end'))).AND.(numSim.LE.MaxCALC))
         numSim=numSim+1
         if ( numSim.gt.MaxCALC) then 
           if (mpi_rank.eq.0) OPEN ( UNIT=2,FILE='history',
     #                               POSITION='APPEND')
           if (mpi_rank.eq.0) WRITE (2,*)
     #                ' Number of simulations is limited by MaxCALC = ',
     #                MaxCALC
           if (mpi_rank.eq.0) WRITE (*,*)
     #                ' Number of simulations is limited by MaxCALC = ',
     #                MaxCALC
           if (mpi_rank.eq.0) CLOSE ( UNIT=2)
           goto 38
         endif                            
c
         READ(1,*,ERR=36)ParInpODF(numSim),ParResODF(numSim),
     *  NSpec(numSim)
                if (NSpec(numSim).gt.MXSPC) then
                  if (mpi_rank.eq.0) OPEN ( UNIT=2,FILE='history',
     #                                      POSITION='APPEND')
                  if (mpi_rank.eq.0) WRITE (2,*)' Simulation ', numSim,
     #                       ', Number of spectra ',NSpec(numSim),
     #                       'exceeds MXSPC =', MXSPC
                  if (mpi_rank.eq.0) WRITE (*,*)' Simulation ', numSim,
     #                       ', Number of spectra ',NSpec(numSim),
     #                       'exceeds MXSPC =', MXSPC
                  if (mpi_rank.eq.0) CLOSE ( UNIT=2)
                  stop
                endif          	            

         if (mpi_rank.eq.0) WRITE(*,100) numSim,ParInpODF(numSim),
     #                       ParResODF(numSim)

         DO i=1,NSpec(numSim)
           READ(1,*,ERR=37) ParInpEPR(i,numSim),ParResEPR(i,numSim), 
     *        NamSpecInp(i,numSim),NamSpecRes(i,numSim),ODFtype(i)
           ODFtypeAll(i,numSim)=ODFtype(i)
           if (mpi_rank.eq.0) WRITE(*,102)ParInpEPR(i,numSim),
     *        ParResEPR(i,numSim), 
     *        NamSpecInp(i,numSim),NamSpecRes(i,numSim),ODFtype(i)
         ENDDO

 100         format(1x,'Task ',1x,I2,':',1x,A20,1x,A20)
c 101         format(10x,A12,1x,A12)
 102         format(1x,'files : ',4(A20,1x),1x,i1)

         READ(1,*,ERR=37)tempStr
         goto 38
  36   continue !error handling
       if (mpi_rank.eq.0) write(*,*) 'first line of task ', numsim,
     #     ' is improper'
       if (mpi_rank.eq.0) write(2,*) 'first line of task ', numsim,
     #     ' is improper'
       stop
  37       continue !error handling
       if (mpi_rank.eq.0) write(*,*) ' task ', numsim, ' is improper'
       if (mpi_rank.eq.0) write(2,*) ' task ', numsim, ' is improper'
       stop
c        
  38   enddo !while
         CLOSE (UNIT=1)
       
       return
       end
