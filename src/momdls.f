cAXV based on MOMDLS.f from NLSL Version 1.5 beta 11/23/95
c    The orientation distribution explicitly presented as 
c     an expansion in a series of generalized spherical harmonics
c     is added.
c*********************************************************************
c
c                       ==================
c                        subroutine:MOMDLS
c                       ==================
c
c       Subroutine version of EPRLL family of programs.
c       This routine is intended for use with nonlinear least-squares
c       applications. The routine calculates a MOMD spectrum by
c       averaging spectra for a specified number of orientation using
c       the parameters passed in the fparm and iparm arrays. Spectra
c       are calculated using the Conjugate Gradients version of the Lanczos
c       algorithm.
c*********************************************************************
c
      subroutine momdls(fparm,iparm,icalc,spectr,nft,jsite,jsp,
     #                  numstor,tmpNams,DiscTrid,ierr,cache_mode,
     #                  reset_task,across_calls,first_site,last_site)
      use, intrinsic :: iso_fortran_env, only: int64
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none
      include 'mpif.h'
      include 'mpi_state.inc'
      include 'dimnls.inc'
      include 'eprprm.inc'
      include 'prmeqv.inc'
      include 'errmsg.inc'
      include 'pidef.inc'
      include 'rndoff.inc'
      include 'parODF.inc'      
c
      integer iparm(NIPRM),icalc,ierr,nft,jsite,jsp,ndone,	
     # numstor(MXSITE,MXSPC),zndone(MXSITE*MXSPC*MAXNORT)      
      integer ipnt,ioerr
      double precision fparm(NFPRM),spectr(MXSPT),
     #                 work(MXSPT),local_spectr(MXSPT)
      double complex store1(MXDIM,MAXStorMat),store2(MXDIM,MAXStorMat)
c   For preservation the values in course of
c   consecutive calling the subroutine lfunR.f:     	
      common /tempmomdls/store1,store2,zndone
c Euler angles SampleFrame -> Resonator frame:     
	 double precision alphaSF, betaSF, gammaSF
c
      double precision dlam,xi,cc,ss
      common /dfunc/ dlam,xi,cc,ss
c
      integer i,iort,id,j,nstp,nevl,jstore,iort_start,iort_step
      double precision cospsi,cosxi,delcos,dnorm,dwt,onorm,sinxi
      logical Full,DiscTrid,reset_task,across_calls,first_site,last_site
      integer cache_mode
c Each entry owns spectra and the completed CG prefixes for this rank.
c Legacy slots may be shared or overwritten on the next evaluation.
      type orientation_result
         double complex, allocatable :: coefficients(:,:)
         integer :: evaluation=0
      end type
      type orientation_cache_entry
         integer(int64) real_key(NFPRM),field_shift_key
         integer(int64) cg_real_key(NFPRM)
         integer(int64) :: last_access=0
         integer int_key(NIPRM),cg_int_key(NIPRM),fft_size
         integer first_angle,angle_stride
         integer :: bytes=0
         double precision, allocatable :: spectra(:,:)
         type(orientation_result) :: result(MAXNORT)
      end type
c Keep variants from neighboring evaluations, with a per-rank limit
c on their allocatable spectrum and coefficient payload.
      integer, parameter :: max_entries=2*MXSITE*MXSPC
      integer, parameter :: max_cache_bytes=112*1024*1024
      integer, parameter :: spectrum_bytes=storage_size(0.d0)/8
      integer, parameter :: cg_step_bytes=2*storage_size((0.d0,0.d0))/8
      integer, parameter :: CFONLY=0
      type(orientation_cache_entry), save :: cache(max_entries)
      integer(int64), save :: access_count=0
      integer, save :: cache_bytes=0,evaluation=0
      integer, save :: full_count,cf_count,hit_count,cross_count
      integer current_entry,coefficient_source,local_count
      logical use_cache
c
      double precision ZERO,ONE,SMALL,D180,acc
      parameter(ZERO=0.0D0,ONE=1.0D0,SMALL=1.0D-16,D180=180.0D0,
     #  acc=1.0d-4)
c
      double precision fu20,fu20phi
      external fu20,fu20phi
c Function from WignerL.for for calculation of ODF value:	 
      external ODFL_RoReFast_Integrated
      double precision ODFL_RoReFast_Integrated
      double precision psiFr 
c
       character*15  tmpNams(MXSITE*MXSPC*MAXNORT)
c
cFor debugging:
cdb       double precision TempField(MXINP), TempSpec(MXINP)
cdb       integer flagname                  
c#####################################################################
c
      Full = mod(icalc,1000)/100 .ne. 0
      use_cache=Full.and.(cache_mode.eq.1.or.cache_mode.eq.2)
      current_entry=0
      coefficient_source=0
      if (first_site) then
         if (reset_task) then
            evaluation=0
            access_count=0
         endif
         if (reset_task.or..not.across_calls) then
            do i=1,max_entries
               call clear_entry(cache(i))
            enddo
         endif
         evaluation=evaluation+1
         full_count=0
         cf_count=0
         hit_count=0
         cross_count=0
      endif
c     ----------------------------------------------
c     Load values from parameter array into /eprprm/
c     ----------------------------------------------
      do i=1,NFPRM
         fepr(i)=fparm(i)
      end do
c
      do i=1,NIPRM
         iepr(i)=iparm(i)
      end do
c
      if (Full) then
         call tchange(jsite)
      else         
c  *** Nevetheless w0 which is used in eprls ***
         w0=(wxx+wyy+wzz)/3.0D0
      endif         
c
      do j=1,nfld
         spectr(j)=ZERO
         local_spectr(j)=ZERO
      end do
c
      nstp=nstep	 !user defined max number of steps for CG tridiagonalization
      if (nstep.le.0) nstp=MXDIM
c
      jstore=numstor(jsite,jsp)
      if (jstore.lt.1) jstore=1   
c
      if (nort.le.1 .or. ipt.le.0) then
c       --------------------------
c       --- Single orientation ---
c       --------------------------
         if ((.not.mpi_enabled).or.(mpi_rank.eq.0)) then
            iort_start=1
            iort_step=1
            local_count=1
            if (use_cache.and.across_calls)
     #         call prepare_cache(.false.)
            call evaluate_orientation(1,local_spectr)
            if (ierr.ge.FATAL) then
               if (mpi_rank.eq.0) write(*,*) ' eprls gives FATAL error'
               if (mpi_enabled) call MPI_ABORT(MPI_COMM_WORLD,ierr,
     #              mpi_ierr)
               stop 1
            endif
         endif
         if (mpi_enabled) then
            call MPI_ALLREDUCE(local_spectr,spectr,nfld,
     #           MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpi_ierr)
         else
            do j=1,nfld
               spectr(j)=local_spectr(j)
            end do
         endif
c         
      else
c      --------------------------
c      ---- MOMD calculation ----
c      --------------------------
         onorm=ONE
         if (nort.gt.1) onorm=ONE/dfloat(nort-1)
c        ------------------------------------------------------
c        Partial MOMD
c        If there is ordering of the director, calculate 
c        normalization factor for the pseudopotential function
c        In the case of partial director ordering, the parameter
c        psi actually represents the angle between the director 
c        ordering axis and the field, whereas the calculation
c        subroutine always interprets psi as the director tilt.
c        Save the given psi value as the angle xi. 
c        ------------------------------------------------------
            dlam=ZERO
            dnorm=ONE
            xi=dmod(psi,D180)
         if (ODFtype(jsp).eq.1) then
            cosxi=dcos(RADIAN*xi)
            sinxi=dsin(RADIAN*xi)
            dlam=dc20
            call ccrint(ZERO,ONE,acc,SMALL,dnorm,nevl,fu20,id)
            if (dabs(dnorm).le.SMALL) dnorm=ONE
         end if
c        -----------------------
c        Loop over orientations
c        -----------------------
         cospsi=ZERO
         delcos=ZERO
         if (nort.gt.1) delcos=ONE/(nort-1)
c         
         iort_start=1
         iort_step=1
         if (mpi_enabled) then
            iort_start=1+mpi_rank
            iort_step=mpi_size
         endif
c
         if (use_cache.and.iort_start.le.nort) then
            local_count=1+(nort-iort_start)/iort_step
            call prepare_cache(.true.)
         endif

         do iort=iort_start,nort,iort_step
           if (mpi_enabled) then
              cospsi=min(dfloat(iort-1)*delcos,ONE)
           else
              cospsi=min(cospsi,ONE)
           endif
           psiFr=dacos(cospsi)
           fepr(IPSI)=psiFr/RADIAN
c
           call evaluate_orientation(iort,work)
c
            if (ierr.ge.FATAL) then 
               if (mpi_rank.eq.0) write(*,*) ' eprls gives FATAL error'
               if (mpi_enabled) call MPI_ABORT(MPI_COMM_WORLD,ierr,
     #              mpi_ierr)
               stop 1
            endif
c           --------------------------------------------------
c           If there is ordering of the director, calculate the
c           weighting coefficient for this tilt angle and
c           weight the spectrum
c           --------------------------------------------------
            if (ODFtype(jsp).eq.1) then
               if (xi.gt.RNDOFF) then
                  ss=sinxi*dsqrt(ONE-cospsi*cospsi)
                  cc=cosxi*cospsi
                  call ccrint(ZERO,PI,acc,SMALL,dwt,nevl,fu20phi,id)
                  dwt=dwt/(dnorm*PI)
               else
                  dwt=fu20(cospsi)/dnorm
               end if  
c
               do j=1,nfld
                  work(j)=work(j)*dwt
               end do
c			               
            end if   !ODFtype=1
c ODF
            if ((ODFtype(jsp).ge.2).and.(keyODF(jsite).ge.0)) then
c  
               alphaSF=0
               betaSF=RADIAN*xi
               gammaSF=0
c
               dwt=ODFL_RoReFast_Integrated
     *             (alphaSF,betaSF,gammaSF,psiFr)
               do j=1,nfld
                  work(j)=work(j)*dwt
               end do
cdb           write(*,*) jsp,iort,dwt              
cdb           call promout(jsp,iort,nfld,TempField,TempSpec)
             end if !	 ODFtype=2 or 3
c          --------------------------------------------------------
c          Trapezoidal rule: first and last pts. get factor of 1/2
c          --------------------------------------------------------
            if (iort.eq.1 .or. iort.eq.nort) then
               do j=1,nfld
                  work(j)=work(j)*0.5d0
               end do
            end if
c        ------------------------------------------------
c        Add spectrum for current orientation into total 
c        ------------------------------------------------
            do j=1,nfld
               local_spectr(j)=local_spectr(j)+work(j)
            end do
c
            if (.not.mpi_enabled) cospsi=cospsi+delcos
         end do !iort, loop over orientations 
c        -------------------
c        Normalize spectrum
c        -------------------
         do j=1,nfld
            local_spectr(j)=local_spectr(j)*onorm
         end do
         if (mpi_enabled) then
            call MPI_ALLREDUCE(local_spectr,spectr,nfld,
     #           MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpi_ierr)
         else
            do j=1,nfld
               spectr(j)=local_spectr(j)
            end do
         endif
c        -------------------------------------------
c        Return total number of steps actually taken
c        -------------------------------------------
      end if !(nort.le.1
c
      if (last_site.and.(cache_mode.eq.2.or.cache_mode.eq.3)) then
         write(*,'(7(A,I0))') 'ORIENTATION_CACHE rank=',mpi_rank,
     #        ' full=',full_count,' cfonly=',cf_count,
     #        ' spectra=',hit_count,' cross=',cross_count,
     #        ' eval=',evaluation,' bytes=',cache_bytes
      endif
      return

      contains

      subroutine prepare_cache(average_angles)
      logical average_angles
      double precision local_inputs(NFPRM)
      integer local_options(NIPRM),cg_options(NIPRM),entry
      integer(int64) real_key(NFPRM),cg_real_key(NFPRM)
      integer(int64) field_shift_key

c In MOMD the input PSI is used for ODF weights. EPRLS uses the local
c angle determined by NORT and the orientation index. A single-angle
c calculation uses PSI directly. ITYPESITE only aliases legacy slots.
      local_inputs=fparm
      if (average_angles) local_inputs(IPSI)=ZERO
      local_inputs(IDC20)=ZERO
      local_inputs(IAMPFAC)=ZERO
      if (.not.all(ieee_is_finite(local_inputs)).or.
     #    .not.ieee_is_finite(fldShftSite(jsite))) return
      local_options=iparm
      local_options(IItypeSite)=0
      real_key=transfer(local_inputs,0_int64,NFPRM)
      field_shift_key=transfer(fldShftSite(jsite),0_int64)

c Only CFS/GCONVL use these values. WXX/WYY/WZZ also set the
c anisotropic matrix broadening, so all three remain in the CG key.
      local_inputs(IPHASE:IGIBX)=ZERO
      local_inputs(ILBL:ILBG)=ZERO
      local_inputs(IRANGE:IDFLD)=ZERO
      local_inputs(ISPSHFT)=ZERO
      local_inputs(IAadd:NFPRM)=ZERO
      cg_real_key=transfer(local_inputs,0_int64,NFPRM)
      cg_options=local_options
      cg_options(INFLD:IIDERV)=0
      cg_options(IIFLGGIB)=0
      cg_options(iNAddN:NIPRM)=0

      access_count=access_count+1
      do entry=1,max_entries
         if (cache(entry)%last_access.eq.0) cycle
         if (cache(entry)%first_angle.ne.iort_start.or.
     #       cache(entry)%angle_stride.ne.iort_step) cycle
         if (all(cache(entry)%int_key.eq.local_options).and.
     #       all(cache(entry)%real_key.eq.real_key).and.
     #       cache(entry)%field_shift_key.eq.field_shift_key.and.
     #       cache(entry)%fft_size.eq.nft) then
            current_entry=entry
            coefficient_source=0
            cache(entry)%last_access=access_count
            return
         endif
         if (across_calls.and.coefficient_source.eq.0) then
            if (all(cache(entry)%cg_int_key.eq.cg_options).and.
     #          all(cache(entry)%cg_real_key.eq.cg_real_key).and.
     #          any(cache(entry)%result%evaluation.ne.0))
     #          coefficient_source=entry
         endif
      enddo

      current_entry=oldest_entry(.true.)
      call clear_entry(cache(current_entry))
      cache(current_entry)%real_key=real_key
      cache(current_entry)%int_key=local_options
      cache(current_entry)%cg_real_key=cg_real_key
      cache(current_entry)%cg_int_key=cg_options
      cache(current_entry)%field_shift_key=field_shift_key
      cache(current_entry)%fft_size=nft
      cache(current_entry)%first_angle=iort_start
      cache(current_entry)%angle_stride=iort_step
      cache(current_entry)%last_access=access_count
      end subroutine prepare_cache

      subroutine evaluate_orientation(orientation,values)
      integer orientation,slot
      double precision values(MXSPT)

      slot=jstore+orientation-1
      call calculate_orientation(orientation,store1(1,slot),
     #     store2(1,slot),zndone(slot),values)
      end subroutine evaluate_orientation

      subroutine calculate_orientation(orientation,al,be,steps,values)
      integer orientation,steps
      double complex al(MXDIM),be(MXDIM)
      double precision values(MXSPT)
      integer source_evaluation,local_index,method

      method=icalc
      if (current_entry.ne.0) then
         local_index=1+(orientation-iort_start)/iort_step
         source_evaluation=cache(current_entry)%result(orientation)%evaluation
         if (source_evaluation.ne.0) then
            values(1:nfld)=cache(current_entry)%spectra(:,local_index)
            call restore_coefficients(current_entry,orientation,al,be,steps)
            ierr=0
            hit_count=hit_count+1
            if (source_evaluation.lt.evaluation) cross_count=cross_count+1
            return
         endif
         if (coefficient_source.ne.0) then
            source_evaluation=cache(coefficient_source)%result(orientation)%evaluation
            if (source_evaluation.ne.0) then
               call restore_coefficients(coefficient_source,orientation,al,be,steps)
               cache(coefficient_source)%last_access=access_count
               method=CFONLY
               ierr=0
               if (source_evaluation.lt.evaluation) cross_count=cross_count+1
            endif
         endif
      endif
      call eprls(method,jsite,al,be,values,nft,nstp,steps,ierr)
      if (method.eq.CFONLY) then
         cf_count=cf_count+1
      else
         full_count=full_count+1
      endif
      if (current_entry.eq.0) return
      if (ierr.ne.0.or.steps.lt.1.or.steps.gt.MXDIM) return
      if (.not.all(ieee_is_finite(values(1:nfld)))) return
      if (method.ne.CFONLY) then
         if (.not.all(ieee_is_finite(dble(al(1:steps)))).or.
     #       .not.all(ieee_is_finite(dimag(al(1:steps)))).or.
     #       .not.all(ieee_is_finite(dble(be(1:steps)))).or.
     #       .not.all(ieee_is_finite(dimag(be(1:steps))))) return
      endif
      call save_orientation(cache(current_entry),orientation,local_index,
     #                      al,be,steps,values)
      end subroutine calculate_orientation

      subroutine restore_coefficients(entry,orientation,al,be,steps)
      integer entry,orientation,steps
      double complex al(MXDIM),be(MXDIM)
c Linked sites and later operations still see the completed prefix
c in the same work slot they would see after FULL.
      associate(coefficients=>cache(entry)%result(orientation)%coefficients)
         steps=size(coefficients,1)
         al(1:steps)=coefficients(:,1)
         be(1:steps)=coefficients(:,2)
      end associate
      end subroutine restore_coefficients

      subroutine save_orientation(entry,orientation,local_index,al,be,
     #                            steps,values)
      type(orientation_cache_entry) entry
      integer orientation,local_index,steps
      double complex al(*),be(*)
      double precision values(*)
      integer bytes,alloc_status

      bytes=steps*cg_step_bytes
      if (.not.allocated(entry%spectra))
     #    bytes=bytes+spectrum_bytes*nfld*local_count
      if (.not.reserve_bytes(bytes)) return

      if (.not.allocated(entry%spectra)) then
         allocate(entry%spectra(nfld,local_count),stat=alloc_status)
         if (alloc_status.ne.0) return
         entry%bytes=spectrum_bytes*nfld*local_count
         cache_bytes=cache_bytes+entry%bytes
      endif
      associate(result=>entry%result(orientation))
         allocate(result%coefficients(steps,2),stat=alloc_status)
         if (alloc_status.ne.0) return
         entry%spectra(:,local_index)=values(1:nfld)
         result%coefficients(:,1)=al(1:steps)
         result%coefficients(:,2)=be(1:steps)
         result%evaluation=evaluation
         entry%bytes=entry%bytes+steps*cg_step_bytes
         cache_bytes=cache_bytes+steps*cg_step_bytes
      end associate
      end subroutine save_orientation

      logical function reserve_bytes(bytes)
      integer bytes,entry
      reserve_bytes=.false.
      do while (cache_bytes+bytes.gt.max_cache_bytes)
         entry=oldest_entry(.false.)
         if (entry.eq.0) return
         call clear_entry(cache(entry))
      enddo
      reserve_bytes=.true.
      end function reserve_bytes

      integer function oldest_entry(allow_empty) result(oldest)
      logical allow_empty
      integer entry
      integer(int64) oldest_access
      oldest=0
      oldest_access=huge(access_count)
      do entry=1,max_entries
         if (entry.eq.current_entry.or.entry.eq.coefficient_source) cycle
         if (.not.allow_empty.and.cache(entry)%last_access.eq.0) cycle
         if (cache(entry)%last_access.ge.oldest_access) cycle
         oldest=entry
         oldest_access=cache(entry)%last_access
      enddo
      end function oldest_entry

      subroutine clear_entry(entry)
      type(orientation_cache_entry) entry
      integer orientation
      if (entry%last_access.eq.0) return
      if (allocated(entry%spectra)) deallocate(entry%spectra)
      do orientation=1,MAXNORT
         if (allocated(entry%result(orientation)%coefficients))
     #       deallocate(entry%result(orientation)%coefficients)
         entry%result(orientation)%evaluation=0
      enddo
      cache_bytes=cache_bytes-entry%bytes
      entry%bytes=0
      entry%last_access=0
      end subroutine clear_entry
      end

c                         function FU20
c                    =========================
c  Function defining Boltzmann population at a given orientation
c  resulting from an orienting pseudopotential using only the c20
c  term. Used for calculating weighting as a function of angle
c  between the average director and the spectrometer field, for
c  the "partial MOMD" model. The potential is intended to describe
c  partial macroscopic orientation of the directors of local domains.
c
c  This is required to be a single-argument function so that it 
c  may be called by ccrint to calculate the normalization factor
c  for orientational weighting coefficients. The second parameter, the 
c  director c20 coefficient, is passed in common /dfunc/
c----------------------------------------------------------------------
      function fu20(cost)
      implicit none
      double precision fu20,cost
c
      double precision dlam,xi,cc,ss
      common /dfunc/ dlam,xi,cc,ss
c
      fu20 = exp( dlam*(1.5d0*cost*cost-0.5d0 ) )
      return
      end
c
c  phi is expected in radians
      function fu20phi(phi)
      implicit none
      double precision costht,fu20phi,phi
c
      include 'pidef.inc'
      double precision dlam,xi,cc,ss
      common /dfunc/ dlam,xi,cc,ss
c
      costht = ss*cos(phi)+cc
      fu20phi = exp( dlam*(1.5d0*costht*costht-0.5d0 ) )
      return
      end
