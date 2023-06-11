module utilities_module
    implicit none
contains

subroutine setup_output_directory(dirname)
    ! VARIABLES
    character(len=*), intent(in) :: dirname
    character(len=60) :: cmd

    ! SUBROUTINE
    cmd = "bash ../src/out_dir.sh "//dirname
    call execute_command_line(cmd)

end subroutine setup_output_directory

end module utilities_module