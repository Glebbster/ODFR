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
       use, intrinsic :: iso_fortran_env, only: error_unit
       include 'mpif.h'
       include 'dimnls.inc'
       include 'name.inc'
       include 'parODF.inc'
       include 'mpi_state.inc'
c
        character(*) NameNames
        integer numSim, NSpec(MaxCALC)
c Read into larger fields to detect overflow before storing paths.
        character(32768) paths(4)

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
         paths(:2)=''
         READ(1,*,ERR=36)paths(:2),NSpec(numSim)
         if (any(len_trim(paths(:2)).gt.MaxPath)) goto 39
         ParInpODF(numSim)=paths(1)
         ParResODF(numSim)=paths(2)
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

         if (mpi_rank.eq.0) WRITE(*,100) numSim,
     #      ParInpODF(numSim)(:max(20,len_trim(ParInpODF(numSim)))),
     #      ParResODF(numSim)(:max(20,len_trim(ParResODF(numSim))))

         DO i=1,NSpec(numSim)
           paths=''
           READ(1,*,ERR=37) paths,ODFtype(i)
           if (any(len_trim(paths).gt.MaxPath)) goto 39
           ParInpEPR(i,numSim)=paths(1)
           ParResEPR(i,numSim)=paths(2)
           NamSpecInp(i,numSim)=paths(3)
           NamSpecRes(i,numSim)=paths(4)
           ODFtypeAll(i,numSim)=ODFtype(i)
           if (mpi_rank.eq.0) WRITE(*,102)
     *    ParInpEPR(i,numSim)(:max(20,len_trim(ParInpEPR(i,numSim)))),
     *    ParResEPR(i,numSim)(:max(20,len_trim(ParResEPR(i,numSim)))),
     *    NamSpecInp(i,numSim)(:max(20,len_trim(NamSpecInp(i,numSim)))),
     *    NamSpecRes(i,numSim)(:max(20,len_trim(NamSpecRes(i,numSim)))),
     *    ODFtype(i)
         ENDDO

 100         format(1x,'Task ',1x,I2,':',1x,A,1x,A)
 102         format(1x,'files : ',4(A,1x),1x,i1)

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
  39   write(error_unit,'(A,A,I0,A,I0,A)') trim(nameNames),
     #      ': path in task ',numSim,' exceeds ',MaxPath,' characters'
       flush(error_unit)
       call MPI_ABORT(MPI_COMM_WORLD,1,mpi_ierr)
       stop 1
       end
