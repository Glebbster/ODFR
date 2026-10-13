c
c    This file is a part of ODFR4 software for EPR Spectra simulation
c    Copyright (C) 2018-2019  Andrey Vorobiev
c
c    The details and updates of this program can be found on its website:
c    http://.../
c
c    Please cite this software as:
c    A. Kh. Vorobiev, A. V. Bogdanov, T. S. Yankova, N. A. Chumakova, 
c    "Spin probe determination of molecular orientation distribution 
c     and rotational mobility in liquid crystals. Model-free approach",
c     Phys. Chem. Chem. Phys., submitted, 2018
c
c
c    This program is free software: you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation, version 3.
c
c    This program is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License
c    along with this program.  If not, see <http://www.gnu.org/licenses/>.
c

c Based on  program NLSL,Version 1.5.1b
c----------------------------------------------------------------------
c
      program ODFR4
c
      include 'mpif.h'
      include 'mpi_state.inc'
      include 'dimnls.inc'
      include 'rndoff.inc'
      include 'eprprm.inc'
      include 'expdat.inc'
      include 'parcom.inc'
      include 'lmcom.inc'
      include 'name.inc'
      include 'parODF.inc'
      include 'version.inc'
c
      integer isp,jpnt,ndattmp	
c
        integer numSim,jSim
        integer NSpec(MaxCALC)
      	 Logical FileErr
c
        integer iparms(NIPRM,MXSITE,MXSPC)
        DOUBLE PRECISION fparms(NFPRM,MXSITE,MXSPC)
        COMMON /COMPARMS/ iparms,fparms
        integer NptsExp
        double precision FieldExp,SpecExp
	 common/iexp/NptsExp(MXSPC)
 	 common/exp/FieldExp(MXINP,MXSPC),SpecExp(MXINP,MXSPC)
       double precision Field(MXSPT,MXSPC),chexp(MXSPT)
       common /Field/ Field
c
      double precision delb 
      common /resSpn/resSpn(MXSPT)     
c	
       character*30, PARAMETER :: programName='ODFR5' 
       character*30, PARAMETER :: nameNames='namODFR.nam'
c	  
       integer   JVAR,Kcontrol
       common /control/Kcontrol
      logical flfirst
      common/flfst/flfirst,jsim
      double precision   DndatNprm     
      common /DndatNprm/ DndatNprm         
c
      double precision time,cpu
      external      lfun     
c nspc - number of spectra          	 
c For debugging:
c       double precision TempField(MXINP), TempSpec(MXINP)
c #####################################################################      
      call MPI_INIT(mpi_ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, mpi_ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, mpi_ierr)
      mpi_enabled = mpi_size .gt. 1

      if (mpi_rank.eq.0) then
         if (mpi_enabled) then
            write(*,*) 'MPI enabled, ranks=',mpi_size
         else
            write(*,*) 'Serial execution, ranks=',mpi_size
         endif
      endif

      call INPNamDR(nameNames,numSim,NSpec)
      jSim=0
      if (mpi_rank.eq.0) then
         OPEN ( UNIT=2,FILE='history', POSITION='APPEND')
         WRITE(2,103) trim(programName),trim(ODFR_VERSION)
      endif
c
 3     jSim=jSim+1
       call cpu_time(cpu)
c 
       call nlsinit
       nspc= NSpec(jSim)
       do isp=1,nspc
          ODFtype(isp)=ODFtypeAll(isp,jSim)
       enddo
c
       if (mpi_rank.eq.0) WRITE(2,102) jSim
       if (mpi_rank.eq.0) WRITE(*,102) jSim
 102   format(1x,'Simulation ',1x,I2,':' )
 103   format(1x,A,' version ',A)
c
       if (mpi_rank.eq.0) WRITE(2,100) ParInpODF(jSim),
     #                    ParResODF(jSim)
 100       format(1x,A12,1x,A12)
       DO isp=1,nspc
         if (mpi_rank.eq.0) WRITE(2,101) ParInpEPR(isp,jSim),
     * ParResEPR(isp,jSim),
     * NamSpecInp(isp,jSim),NamSpecRes(isp,jSim)
       END DO	!isp
 101       format(1x,'files : ',4(A20,1x),F6.2)
c
       FileErr=.False.
       if (mpi_rank.eq.0) CLOSE ( UNIT=2)
       if (mpi_rank.eq.0) OPEN ( UNIT=2,FILE='history',
     #                           POSITION='APPEND')
       do isp=1,nspc
            call parreadN(FileErr,ParInpEPR(isp,jSim),isp)
       if (FileErr) goto 2
       enddo !isp
c
       call ODFread(ParInpODF(jSim))
       nprm=nprm+nprm4
c
	 call xpack(1,JVAR,xfit)
c
       if (mpi_rank.eq.0) then
          OPEN ( UNIT=1,FILE='temp 1.par')
     	       write(1,*) Kcontrol
          CLOSE (  UNIT=1 )
       endif
       if (mpi_enabled) call MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
          SELECT CASE (Kcontrol)
          CASE (2)
		   nprm=0 !single calculation of spectra
          CASE (1)
             goto 9999
          CASE default	        	     
            if (JVAR.ne.nprm) then ! NUMBER OF PARAMETERS IS WRONG   
                   if (mpi_rank.eq.0) WRITE(*,*)
     #              '!!  NUMBER OF PARAMETERS IS WRONG !!'
                   if (mpi_rank.eq.0) WRITE(*,*)
     #              '       really      declared       '
                   if (mpi_rank.eq.0) write(*,*)jvar,nprm
              goto 9999
            else
              do ix=1,nprm
                xfdstp(ix)=dabs(gtol*xfit(ix))
                if(xfdstp(ix).lt.RNDOFF) xfdstp(ix)=gtol
              enddo
            endif
          END SELECT !kcontrol
c temporary - for debugging:
c        isp=1
c       do isi=1,MXSITE
c         do jpar=1,NFPRM
c           fparm(jpar,isi)=fparms(jpar,isi,isp)
c         enddo !jpar
c         do jpar=1,NIPRM
c           iparm(jpar,isi)=iparms(jpar,isi,isp)
c         enddo !jpar
c       enddo !isi
c        -------------------------------------------------------------
c         Find smallest power of 2 greater than or equal to the number
c         of data points (for Fourier-Transform applications)
c        -------------------------------------------------------------
       do isp=1,nspc
         npts(isp)=iparms(INFLD,1,isp)
         ndatot=ndatot+npts(isp)
         if (rmsn(nspc).le.RNDOFF) rmsn(nspc)=1.0d0
         nft(isp)=1
 9       nft(isp)=nft(isp)*2
         if (nft(isp).lt.npts(isp)) go to 9
       enddo !isp
c preparation of field array
        do isp=1,nspc
          delb=fparms(IDFLD,1,isp)
          DO jpnt=1,npts(isp)
	       Field(jpnt,isp)=fparms(IFLDI,1,isp)+delb*(jpnt-1)
          END DO
        enddo
c
c get experimental data
	    call specread(FileErr,jSim,nspc,keySp)
	    if (FileErr) goto 44
          ndattmp=0
          do isp=1,nspc
            call choiceExp(iparms(INFLD,1,isp),NptsExp(isp),
     #           FieldExp(1,isp),SpecExp(1,isp),Field(1,isp),chexp)
              ixsp(isp)=ndattmp+1
            DO jpnt=1,npts(isp)
                   fitdata(ndattmp+jpnt)=chexp(jpnt)
            END DO
            ndattmp=ndattmp+npts(isp)
c
c For debugging:
c          do jpnt=1,npts(isp)
c             jtmp=ixsp(isp)+jpnt-1
c	       TempSpec(jpnt)=fitdata(jtmp)
c	       TempField(jpnt)=Field(jpnt,1)
c          enddo
c          call promout(npts(isp),TempField,TempSpec)
c
          enddo !isp
c 
      DndatNprm=dfloat(ndatot-nprm)
      info=0
      if (nprm.eq.0) then
          call lfun(0,0,xfit,fvec,fvecPure,fjac,0,1) !fresh reporting pass
          call lfun(0,0,xfit,fvec,fvecPure,fjac,0,3)
          if (mpi_rank.eq.0) write(2,*) ' Finished'
      else !fitting
         if (Kcontrol.eq.3) then 
           call fitNL2 
         else      
           call fitl(info)
         endif        
         call lfun(0,0,xfit,fvec,fvecPure,fjac,0,1) !fresh reporting pass
         call lfun(0,0,xfit,fvec,fvecPure,fjac,0,3)

      endif !if (nprm.eq.0)
c   output:
        if (Kcontrol.eq.1) Kcontrol=3	         
        if (mpi_rank.eq.0) then
          do isp=1,nspc
            call SpecOUT(programName,NamSpecRes(isp,jSim),isp,
     #         NamSpecInp(isp,jSim),ParInpEPR(isp,jSim),
     #         ParInpODF(jSim),info)
            call PARWrite(ParResEPR(isp,jSim),isp)
          enddo ! isp
          call ODFwrite(ParResODF(jSim))
        endif
c
 2    continue
       call cpu_time(time)
       time=time-cpu
       if (mpi_rank.eq.0) write(*,*)'  CPU time: ',	time
       if (mpi_rank.eq.0) write(2,*)'  CPU time: ',	time
 44    if (numSim-jSim) 1, 1, 3
 1    continue
c
 9999 continue
       if (mpi_rank.eq.0) CLOSE ( UNIT=2)
       call MPI_FINALIZE(mpi_ierr)
      stop
c
      end
c**********************************************************************
c       SUBROUTINE promOUT(flagname,iort,NTOCH,spec1,spec2)
c       IMPLICIT DOUBLE PRECISION (A-H),DOUBLE PRECISION (O-Z)
c       dimension spec1(*),spec2(*)
c       CHARACTER*30 AME,tmp
c       Character*2 tmp1,tmp2
c       integer flagname,iort 
C
c      if (flagname.eq.0) then
c         write(*,*) ' INPUT FILENAME FOR promout '
c         read(*,15) ame
c 15      format(A12)
c      else
c       write(tmp1,101)flagname
c       write(tmp2,101)iort
c 101        FORMAT (I2)
c         ame=tmp1//tmp2//'.esr'
c      endif   
c        write(*,*)tmp1, tmp2,ame             
c        OPEN ( UNIT=1,FILE=ame)
c        write(1,*)''
c        write(1,*)''
c        tmp='Spectr:'
c        write(1,6)tmp
c  6     format(A7) 
c        DO 10 I=1,NTOCH
c        WRITE(1,*)SPEC1(I),SPEC2(I)
c 5     FORMAT (1X,4(G14.8,1X))
c   10    CONTINUE
c       CLOSE (UNIT=1)
c      return
c      end 
c######################################################################
