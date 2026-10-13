      module optimization_steps
c     Passive NL2SOL observer. Never writes solver/common-block state.
c     Two snapshots follow the solver's current R and saved R exactly.
      implicit none
      private
      public steps_begin,steps_end,steps_capture,steps_save,
     #       steps_restore,steps_accept
      logical, save :: enabled=.false.,valid(2)=.false.
      integer, save :: frame=0,last_iteration=-1
      integer, save :: evaluations(2)=0
      character(160), save :: directory=''
      double precision, allocatable, save :: fields(:,:),models(:,:),
     #       parameters(:,:),last_x(:)
      double precision, save :: losses(2,2)=0.d0
      contains

      subroutine steps_end()
      enabled=.false.
      valid=.false.
      if (allocated(fields)) deallocate(fields)
      if (allocated(models)) deallocate(models)
      if (allocated(parameters)) deallocate(parameters)
      if (allocated(last_x)) deallocate(last_x)
      end subroutine

      subroutine steps_fail(message)
      character(*), intent(in) :: message
      write(*,'(A)') 'Optimization steps disabled: '//message
      call steps_end()
      end subroutine

      subroutine steps_begin(n,p)
      include 'mpi_state.inc'
      integer, intent(in) :: n,p
      integer ios,cmdstat,exitstat,k
      logical exists
      character(32) option,os
      character(160) native
      call steps_end()
      if (mpi_rank.ne.0) return
      call get_environment_variable('ODFR_OPTIMIZATION_STEPS',option,
     #                              status=ios)
      if (ios.ne.0) return
      if (trim(option).ne.'1') return
      allocate(fields(n,2),models(n,2),parameters(p,2),last_x(p),
     #         stat=ios)
      if (ios.ne.0) then
         call steps_fail('cannot allocate snapshot buffers')
         return
      endif
      inquire(file='optimization_steps',exist=exists)
      if (.not.exists) then
         call execute_command_line('mkdir optimization_steps',
     #        exitstat=exitstat,cmdstat=cmdstat)
         if (cmdstat.ne.0.or.exitstat.ne.0) then
            call steps_fail('cannot create optimization_steps')
            return
         endif
      endif
      call get_environment_variable('OS',os,status=ios)
      do k=1,999999
         write(directory,'(A,I6.6)') 'optimization_steps/fit_',k
         inquire(file=trim(directory),exist=exists)
         if (exists) cycle
         native=directory
         if (trim(os).eq.'Windows_NT') native(19:19)=achar(92)
         call execute_command_line('mkdir "'//trim(native)//'"',
     #        exitstat=exitstat,cmdstat=cmdstat)
         if (cmdstat.eq.0.and.exitstat.eq.0) exit
c        Another process may have claimed this number in the meantime.
         inquire(file=trim(directory),exist=exists)
         if (exists) cycle
         call steps_fail('cannot create fit directory')
         return
      enddo
      if (k.gt.999999) then
         call steps_fail('fit directory numbers exhausted')
         return
      endif
      enabled=.true.
      frame=0
      last_iteration=-1
      write(*,'(A)') 'Optimization steps: '//trim(directory)
      call steps_experiment()
      end subroutine

      subroutine steps_experiment()
      include 'dimnls.inc'
      integer NptsExp,lu,ios,close_ios,isp,i
      double precision FieldExp,SpecExp
      common /iexp/ NptsExp(MXSPC)
      common /exp/ FieldExp(MXINP,MXSPC),SpecExp(MXINP,MXSPC)
      include 'expdat.inc'
      open(newunit=lu,file=trim(directory)//'/experiment.dat',
     #     status='new',action='write',iostat=ios)
      if (ios.ne.0) then
         call steps_fail('cannot open experiment.dat')
         return
      endif
      write(lu,'(A)',iostat=ios) '# spectrum field signal (original input)'
      do isp=1,nspc
         do i=1,NptsExp(isp)
            if (ios.ne.0) exit
            write(lu,'(I4,2(1X,ES24.16E3))',iostat=ios)
     #           isp,FieldExp(i,isp),SpecExp(i,isp)
         enddo
      enddo
      close(lu,iostat=close_ios)
      if (ios.ne.0.or.close_ios.ne.0)
     #   call steps_fail('cannot write experiment.dat')
      end subroutine

      subroutine steps_capture(p,x,nf)
      include 'dimnls.inc'
      include 'expdat.inc'
      include 'iterat.inc'
      integer, intent(in) :: p,nf
      double precision, intent(in) :: x(p)
      double precision Field
      common /Field/ Field(MXSPT,MXSPC)
      integer isp,i,k
      if (.not.enabled) return
      k=0
      do isp=1,nspc
         do i=1,npts(isp)
            k=k+1
            fields(k,1)=Field(i,isp)
            models(k,1)=ClSpecDat(i,isp,1)
         enddo
      enddo
      parameters(:,1)=x
      evaluations(1)=nf
      losses(:,1)=(/lossSseW,lossSseUn/)
      valid(1)=.true.
      end subroutine

      subroutine steps_save()
      if (.not.enabled) return
      fields(:,2)=fields(:,1)
      models(:,2)=models(:,1)
      parameters(:,2)=parameters(:,1)
      losses(:,2)=losses(:,1)
      evaluations(2)=evaluations(1)
      valid(2)=valid(1)
      end subroutine

      subroutine steps_restore()
      if (.not.enabled) return
      fields(:,1)=fields(:,2)
      models(:,1)=models(:,2)
      parameters(:,1)=parameters(:,2)
      losses(:,1)=losses(:,2)
      evaluations(1)=evaluations(2)
      valid(1)=valid(2)
      end subroutine

      subroutine steps_accept(iteration,p,x)
      include 'dimnls.inc'
      include 'expdat.inc'
      integer, intent(in) :: iteration,p
      double precision, intent(in) :: x(p)
      integer lu,ios,close_ios,isp,i,k
      character(32) name
      if (.not.enabled) return
      if (.not.valid(1)) then
         call steps_fail('accepted snapshot is unavailable')
         return
      endif
c     Fail closed rather than label a stale/Jacobian curve as accepted.
      if (any(parameters(:,1).ne.x)) then
         call steps_fail('accepted parameters do not match snapshot')
         return
      endif
      if (last_iteration.ge.0) then
         if (all(last_x.eq.x)) return
      endif
      write(name,'(A,I6.6,A)') 'step_',frame,'.dat'
      open(newunit=lu,file=trim(directory)//'/'//trim(name),
     #     status='new',action='write',iostat=ios)
      if (ios.ne.0) then
         call steps_fail('cannot open step table')
         return
      endif
      write(lu,'(A,I0)',iostat=ios) '# iteration ',iteration
      if (ios.eq.0) write(lu,'(A,I0)',iostat=ios)
     #    '# evaluation ',evaluations(1)
      if (ios.eq.0) write(lu,'(A,ES24.16E3)',iostat=ios)
     #    '# SSE_W ',losses(1,1)
      if (ios.eq.0) write(lu,'(A,ES24.16E3)',iostat=ios)
     #    '# SSE_UNW ',losses(2,1)
      if (ios.eq.0) write(lu,'(A,*(1X,ES24.16E3))',iostat=ios)
     #    '# x_optimizer',x
      if (ios.eq.0) write(lu,'(A)',iostat=ios)
     #    '# spectrum field model'
      k=0
      do isp=1,nspc
         do i=1,npts(isp)
            k=k+1
            if (ios.ne.0) exit
            write(lu,'(I4,2(1X,ES24.16E3))',iostat=ios)
     #           isp,fields(k,1),models(k,1)
         enddo
      enddo
c     A reader ignores an interrupted or incomplete table.
      if (ios.eq.0) write(lu,'(A)',iostat=ios) '# complete'
      close(lu,iostat=close_ios)
      if (ios.ne.0.or.close_ios.ne.0) then
         call steps_fail('cannot finish step table')
         return
      endif
      last_x=x
      last_iteration=iteration
      frame=frame+1
      end subroutine
      end module
