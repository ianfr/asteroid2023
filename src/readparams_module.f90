module readparams_module
    implicit none

contains

! based on https://jblevins.org/log/control-file
subroutine read_params_from_file(filename, num_particles, num_timesteps, max_time, asteroid_masses, asteroid_radii, &
    asteroid_positions, asteroid_velocities, dt, particle_radius, out_dir, write_mod)
    ! VARIABLES
    ! For returning to main program
    integer, intent(out) :: num_particles, num_timesteps
    real, intent(out) :: max_time, dt, particle_radius
    real, dimension(:), allocatable, intent(out) :: asteroid_masses ! n masses
    real, dimension(:), allocatable, intent(out) :: asteroid_radii ! the radii of each asteroid
    real, dimension(:,:), allocatable, intent(out) :: asteroid_positions ! 3 rows, n columns
    real, dimension(:,:), allocatable, intent(out) :: asteroid_velocities ! 3 rows, n columns
    character(len=:), allocatable, intent(out) :: out_dir
    integer, intent(out) :: write_mod
    ! For the subroutine
    character(len=*), intent(in) :: filename
    character(len=100) :: buffer, label, out_dir_buffer
    integer :: pos
    integer, parameter :: fh = 15
    integer :: ios = 0
    integer :: line = 0
    character(len=5) ast_m_file, ast_rad_file, ast_pos_file, ast_vel_file
    real, dimension(:), allocatable :: temp3vector
    integer :: num_asteroids, n, i


    ! ROUTINE
    open(fh, file=filename)

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.

    do while (ios == 0)
        read(fh, '(A)', iostat=ios) buffer
        ! print*, "ios=", ios
        if (ios == 0) then
            line = line + 1

            ! Find the first instance of whitespace.  Split label and data.
            pos = scan(buffer, '    ')
            label = buffer(1:pos)
            ! print*, "label:", label
            buffer = buffer(pos+1:)
            ! print*, "buffer:", buffer

            select case (label)
            case ('num_particles')
                read(buffer, *, iostat=ios) num_particles
                print *, 'Read num_particles: ', num_particles
            case ('num_timesteps')
                read(buffer, *, iostat=ios) num_timesteps
                print *, 'Read num_timesteps: ', num_timesteps
            case ('max_time')
                read(buffer, *, iostat=ios) max_time
                print *, 'Read max_time: ', max_time
            case ('num_asteroids')
                read(buffer, *, iostat=ios) num_asteroids
                print *, 'Read num_asteroids: ', num_asteroids
            case ('ast_m_file')
                read(buffer, *, iostat=ios) ast_m_file
                print *, 'Read ast_m_file: ', ast_m_file
            case ('ast_rad_file')
                read(buffer, *, iostat=ios) ast_rad_file
                print *, 'Read ast_rad_file: ', ast_rad_file
            case ('ast_pos_file')
                read(buffer, *, iostat=ios) ast_pos_file
                print *, 'Read ast_pos_file: ', ast_pos_file
            case ('ast_vel_file')
                read(buffer, *, iostat=ios) ast_vel_file
                print *, 'Read ast_vel_file: ', ast_vel_file
            case ('dt')
                read(buffer, *, iostat=ios) dt
                print *, 'Read dt: ', dt
            case ('particle_radius')
                read(buffer, *, iostat=ios) particle_radius
                print *, 'Read particle_radius: ', particle_radius
            case ('out_dir')
                read(buffer, *, iostat=ios) out_dir_buffer
                out_dir = trim(out_dir_buffer)
                print *, 'Read out_dir: ', out_dir
            case ('write_mod')
                read(buffer, *, iostat=ios) write_mod
                print *, 'Read write_mod: ', write_mod
            case default
                print *, 'Skipping invalid label at line', line
            end select
        end if
    end do

    ! now, read in the matrices using the filenames stored in ast_[etc]
    ! used https://shocksolution.com/2010/05/21/reading-an-array-from-a-text-file-with-fortran-9095/

    print*, "Reading asteroid masses, radii, positions, and velocities..."
    n = num_asteroids

    open (unit=90, file="../config/"//ast_m_file, status='old', action='read')
    allocate(asteroid_masses(n))
    read(90,*) asteroid_masses
    close(90)

    open (unit=90, file="../config/"//ast_rad_file, status='old', action='read')
    allocate(asteroid_radii(n))
    read(90,*) asteroid_radii
    close(90)

    ! the matrices in the files get transposed since fortran is column-major in memory
    allocate(temp3vector(3))
    open (unit=91, file="../config/"//ast_pos_file, status='old', action='read')
    allocate(asteroid_positions(3,n))
    do i=1,n,1
        read(91,*) temp3vector
        asteroid_positions(:, i) = temp3vector
    end do
    close(91)

    open (unit=92, file="../config/"//ast_vel_file, status='old', action='read')
    allocate(asteroid_velocities(3,n))
    do i=1,n,1
        read(92,*) temp3vector
        asteroid_velocities(:, i) = temp3vector
    end do
    close(92)

    print*, "...done reading asteroid masses, radii, positions, and velocities."

end subroutine

end module readparams_module

