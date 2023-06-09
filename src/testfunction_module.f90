module testfunction_module
    implicit none
contains

function testfunction(a) result(ret)
    integer, intent(in) :: a
    integer :: ret

    ret = a * a
end function

end module testfunction_module