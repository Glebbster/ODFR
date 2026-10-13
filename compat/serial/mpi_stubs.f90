subroutine MPI_INIT(ierr)
  implicit none
  integer, intent(out) :: ierr
  ierr = 0
end subroutine MPI_INIT

subroutine MPI_COMM_RANK(comm, rank, ierr)
  implicit none
  integer, intent(in) :: comm
  integer, intent(out) :: rank, ierr
  rank = 0
  ierr = 0
end subroutine MPI_COMM_RANK

subroutine MPI_COMM_SIZE(comm, size, ierr)
  implicit none
  integer, intent(in) :: comm
  integer, intent(out) :: size, ierr
  size = 1
  ierr = 0
end subroutine MPI_COMM_SIZE

subroutine MPI_BARRIER(comm, ierr)
  implicit none
  integer, intent(in) :: comm
  integer, intent(out) :: ierr
  ierr = 0
end subroutine MPI_BARRIER

subroutine MPI_ALLREDUCE(sendbuf, recvbuf, count, datatype, op, comm, ierr)
  implicit none
  integer, intent(in) :: count, datatype, op, comm
  integer, intent(out) :: ierr
  integer :: i
  double precision, intent(in) :: sendbuf(*)
  double precision, intent(out) :: recvbuf(*)

  do i = 1, count
    recvbuf(i) = sendbuf(i)
  end do
  ierr = 0
end subroutine MPI_ALLREDUCE

subroutine MPI_ABORT(comm, errorcode, ierr)
  implicit none
  integer, intent(in) :: comm, errorcode
  integer, intent(out) :: ierr
  ierr = errorcode
  stop 1
end subroutine MPI_ABORT

subroutine MPI_FINALIZE(ierr)
  implicit none
  integer, intent(out) :: ierr
  ierr = 0
end subroutine MPI_FINALIZE
