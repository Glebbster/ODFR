cAXV For ODFR
	 subroutine specread(FileErr,jSim,numSpec,keySp)
c       subroutine read the experimental spectrum
          IMPLICIT DOUBLE PRECISION (A-H),DOUBLE PRECISION (O-Z)
       include 'dimnls.inc'
       include 'name.inc'
      include 'expdat.inc'
      include 'mpi_state.inc'

	 Logical FileErr
	 character*30 nametemp,STR
       integer NptsExp
       integer keySp(MXSPC,2)
       integer jSim,numSpec


       double precision FieldExp,SpecExp,temp
	 common/iexp/NptsExp(MXSPC)
 	 common/exp/FieldExp(MXINP,MXSPC),SpecExp(MXINP,MXSPC)
c####################################################################################
c       
      DO j=1, numSpec
          nametemp=NamSpecInp(j,jSim)
	    irest=1
  20	    NptsExp(j)=1
	    str='No'
     	    IOCHECK=0
          OPEN ( UNIT=1,FILE=nametemp,err=30,status='old')
          DO WHILE ((LGT(STR,'Spectr:')).OR.(LLT(STR,'Spectr:')))
            READ(1,15)STR 
  15        FORMAT(A12)		                
          END DO
c  
	          READ(1,*,IOSTAT=iocheck) A,B
              if (iocheck.ne.0) goto 30
	          FieldExp(1,j)=A-ShiftExp(j)
          temp=A		 
          SPECEXP(1,j)=B
          NptsExp(j)=2
          DO WHILE ((IOCHECK).EQ.0)
            do 3 i=1,irest 
   	      READ(1,*,IOSTAT=iocheck) A,B
              if (iocheck.ne.0) goto 3
   	      if (NptsExp(j).lt.MXINP) then
               FieldExp(NptsExp(j),j)=A-ShiftExp(j)
               SPECEXP(NptsExp(j),j)=B
               if (FieldExp(NptsExp(j),j).EQ.temp) then
                  NptsExp(j)=NptsExp(j)-1  
                  SPECEXP(NptsExp(j),j)=(SPECEXP(NptsExp(j),j)+B)/2.d0
               endif
               temp=FieldExp(NptsExp(j),j)               
            endif    
   3	      continue
	      NptsExp(j)=NptsExp(j)+1
          END DO !WHILE     
	    NptsExp(j)=NptsExp(j)-1
          CLOSE (UNIT=1)
   	    if (NptsExp(j).gt.MXINP) then
	       irest=int(NptsExp(j)/MXINP+1)
	       if (mpi_rank.eq.0) write(*,*)
     #          ' Number points in spectrum are', NptsExp(j)
	       if (mpi_rank.eq.0) write (*,*)
     #          ' Every ',irest,'-d point of spectrum is taken'
	       goto 20
	    endif
	          DO WHILE ((NptsExp(j).gt.1).and.
     #              ((FieldExp(NptsExp(j),j)).EQ.0))
	            NptsExp(j)=NptsExp(j)-1
	          END DO 

          DO I=1,NptsExp(j)
	       SpecExp(I,j)=SpecExp(I,j)*keySp(j,1)
          END DO
c
          goto 31
c
 30       FileErr=.true.
          if (mpi_rank.eq.0) OPEN ( UNIT=2,FILE='history',
     #                              POSITION='APPEND')
          if (mpi_rank.eq.0) write(2,*)' FILE ',nametemp,
     #       ' NOT FOUND OR NOT PREPARED PROPERLY!'
          if (mpi_rank.eq.0) write(*,*)' FILE ',nametemp,
     #       ' NOT FOUND OR NOT PREPARED PROPERLY!'
          if (mpi_rank.eq.0) CLOSE ( UNIT=2)
c
 31	    continue
      END DO
       return
       end
