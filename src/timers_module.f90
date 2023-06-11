! Adapted from https://github.com/Onturenio 
! Timer class to provide tic()/toc() functionality like MATLAB
module timers_module

  implicit none
  private
  public :: Timer
  integer, parameter ::  KINT4 = selected_int_kind(6)

  type Timer
    private
    integer(KINT4) :: start, rate=-1
  contains
    procedure, public :: tic, toc
  end type Timer

contains

  subroutine tic(self)
    class (Timer), intent(inout) :: self
    integer(KINT4) :: start, rate

    call system_clock(count_rate=rate)
    call system_clock(start)
    self%start=start
    self%rate=rate
  end subroutine

  subroutine toc(self)
    class (Timer), intent(in) :: self
    integer(KINT4) :: finish

    if(self%rate<0) then
      print*, 'Call to ''toc'' subroutine must come after call to ''tic'''
      stop
    endif
    call system_clock(finish)
    print*, 'Elapsed time in seconds:', float(finish-self%start)/self%rate
  end subroutine

end module timers_module

! USAGE (from https://www.um.es/gmar/staff/gomez/computation/2016/06/28/timers.html):

! program testtimers
!   use Timers
!   implicit none

!   type(Timer) :: crono1

!   call crono1%Tic()

!   ! some heavy stuff to be benchmarked

!   call crono1%Toc()

! end program testtimers
