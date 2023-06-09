module utilities_module
    implicit none
contains

subroutine push_realloc_1d(the_array, the_item)
    ! VARIABLES
    real, allocatable, intent(inout) :: the_array(:)
    real, intent(in) :: the_item
    real, allocatable :: tmp_array(:)

    ! SUBROUTINE
    allocate(tmp_array(size(the_array) + 1))
    tmp_array(1:size(the_array)) = the_array
    tmp_array(size(tmp_array)) = the_item
    the_array = tmp_array ! req. fortran 2003?

end subroutine push_realloc_1d

subroutine setup_output_directory(dirname)
    ! VARIABLES
    character(len=*), intent(in) :: dirname
    character(len=60) :: cmd

    ! SUBROUTINE
    cmd = "bash ../src/out_dir.sh "//dirname
    call execute_command_line(cmd)

end subroutine setup_output_directory


end module utilities_module