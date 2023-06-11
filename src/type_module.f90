module type_module
    implicit none

    ! Basic particle type
    type particle
        integer :: id
        integer :: ast_id
        integer :: color
        real :: mass
        real :: radius
        real :: pos(3)
        real :: vel(3)
        real, allocatable :: grav_neighs(:) ! neighbor ids for gravity
    end type particle

    ! Basic collision type
    type collision_struct
        integer, dimension(2) :: particle_ids = (/ 0, 0 /)
        real :: collision_time = 10000000.0
    end type collision_struct

    ! For selecting force calculation scheme
    type force_calc_struct
        logical :: use_brute = .false.
        logical :: use_radial_prob = .false.
        real :: radius
        real :: prop_full ! proportion of the time to do the full force calculation
    end type force_calc_struct


contains

! ROUTINES FOR PARTICLES

subroutine print_particle(the_particle)
    type(particle), intent(in) :: the_particle

    print*, "PARTICLE START"
    print*, "id:"
    print*, the_particle%id
    print*, "mass:"
    print*, the_particle%mass
    print*, "position:"
    print*, the_particle%pos
    print*, "velocity:"
    print*, the_particle%vel
    print*, "PARTICLE END"
end subroutine print_particle

subroutine create_particle(theParticle)
    implicit none
    type(particle), intent(out) :: theParticle

    theParticle%id = 0
    theParticle%ast_id = 0
    theParticle%color = 0
    theParticle%mass = 1
    theParticle%radius = 1
    theParticle%pos = [0,0,0]
    theParticle%vel = [0,0,0]
    theParticle%grav_neighs = [0]
end subroutine

subroutine push_particle(the_array, the_item)
    ! VARIABLES
    type(particle), allocatable, intent(inout) :: the_array(:)
    type(particle), intent(in) :: the_item
    type(particle), allocatable :: tmp_array(:)

    ! SUBROUTINE
    allocate(tmp_array(1))
    tmp_array(1) = the_item
    if (allocated(the_array)) then
        the_array = [the_array(:), tmp_array(:)]
    else
        allocate(the_array(1))
        the_array = tmp_array
    endif
end subroutine push_particle

subroutine write_particle_list_to_file(particle_list, filename)
    use, intrinsic :: iso_fortran_env, only: error_unit
    ! Inputs
    type(particle), dimension(:), intent(in) :: particle_list
    character(len=*), intent(in) :: filename

    ! Internal variables
    integer file_unit, rc, i

    ! subroutine
    print*, "[write_particle_list_to_file] writing..."
    open(action='write', file=filename, iostat=rc, newunit=file_unit)

    if (rc /= 0) then
        write (error_unit, '(3a, i0)') 'Writing file "', filename, '" failed: ', rc
        stop
    end if
    
    write (file_unit, *, iostat=rc) size(particle_list)
    do i = 1, size(particle_list), 1
        write (file_unit, *, iostat=rc) particle_list(i)%id, particle_list(i)%ast_id, particle_list(i)%color
        if (rc /= 0) exit
        write (file_unit, *, iostat=rc) particle_list(i)%mass, particle_list(i)%radius
        write (file_unit, *, iostat=rc) particle_list(i)%pos
        write (file_unit, *, iostat=rc) particle_list(i)%vel
        !write (file_unit, *, iostat=rc) particle_list(i)%grav_neighs                         
    end do

    close (file_unit)
    print*, "[write_particle_list_to_file] DONE writing."

end subroutine write_particle_list_to_file

subroutine read_particle_list_from_file(particle_list, filename)
    use, intrinsic :: iso_fortran_env, only: error_unit
    ! Inputs
    type(particle), dimension(:), allocatable, intent(out) :: particle_list
    character(len=*), intent(in) :: filename

    ! Internal variables
    integer file_unit, rc, i
    integer num_particles

    ! subroutine
    print*, "[read_particle_list_from_file] reading..."
    open(action='read', file=filename, iostat=rc, newunit=file_unit)

    if (rc /= 0) then
        write (error_unit, '(3a, i0)') 'Reading file "', filename, '" failed: ', rc
        stop
    end if

    ! read in the number of particles and allocate
    read (file_unit, *, iostat=rc) num_particles
    allocate(particle_list(num_particles))
    
    do i = 1, size(particle_list), 1
        read (file_unit, *, iostat=rc) particle_list(i)%id, particle_list(i)%ast_id, particle_list(i)%color
        if (rc /= 0) exit
        read (file_unit, *, iostat=rc) particle_list(i)%mass, particle_list(i)%radius
        read (file_unit, *, iostat=rc) particle_list(i)%pos
        read (file_unit, *, iostat=rc) particle_list(i)%vel
        !write (file_unit, *, iostat=rc) particle_list(i)%grav_neighs                         
    end do

    close (file_unit)
    print*, "[read_particle_list_from_file] DONE reading."
end subroutine read_particle_list_from_file

character(len=20) function str(k)
!   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str

subroutine write_particle_list_for_paraview(particle_list, dirname, filenumber)
    use, intrinsic :: iso_fortran_env, only: error_unit
    ! Inputs
    type(particle), dimension(:), intent(in) :: particle_list
    character(len=*), intent(in) :: dirname
    integer, intent(in) :: filenumber

    ! Internal variables
    integer file_unit, rc, i
    character(len=:), allocatable :: filename
    character(len=10) num_str

    ! subroutine
    ! print*, "[write_particle_list_for_paraview] writing..."
    ! print*, dirname

    write(num_str, '(I10.10)') filenumber

    open(action='write', file='../OUT/'//dirname//'/ast.csv.'//num_str, iostat=rc, newunit=file_unit)

    if (rc /= 0) then
        write (error_unit, '(3a, i0)') 'Writing file "', filename, '" failed: ', rc
        stop
    end if

    write (file_unit, *, iostat=rc) "rx, ry, rz, col_ast"
    
    do i = 1, size(particle_list), 1
        write (file_unit, *, iostat=rc) particle_list(i)%pos(1), ', ', particle_list(i)%pos(2), ', ', &
            particle_list(i)%pos(3), ', ', particle_list(i)%color
        if (rc /= 0) exit                        
    end do

    close (file_unit)
    !print*, "[write_particle_list_for_paraview] DONE writing."

end subroutine write_particle_list_for_paraview


real function my_norm2(the_list) result(the_norm)
    real, dimension(:), intent(in) :: the_list
    integer :: i

    the_norm = 0
    do i = 1, 3, 1
        the_norm = the_norm + the_list(i)*the_list(i)
    end do
    the_norm = sqrt(the_norm)
    return
end function


end module type_module 
