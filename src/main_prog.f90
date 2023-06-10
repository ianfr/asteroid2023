program main_prog

    ! IMPORTS
    use asteroid_module
    use type_module
    use utilities_module
    use gravity_module
    use collision_module

    ! VARIABLES
    implicit none
    integer :: i, file_counter, bbox_count
    real :: total_time, eff_dt
    type(particle), allocatable :: particle_list(:)

    ! Program parameters
    real, parameter :: DT  = 0.15 ! DT is an upper bound for timesteps
    character(len=*), parameter :: PARAM_FILE_NAME = "../config/control.txt"
    character(len=*), parameter :: OUT_DIR = "cmake-test"
    integer, parameter :: MAX_NUM_TSTEPS = 100
    integer, parameter :: WRITE_MOD = 10
    real, parameter :: MAX_TIME = 3
    real, parameter :: COLLISION_RANGE = 0.05 ! check with particle radii below to see if appropriate

    ! PROGRAM START
    print*, "SIMULATION OF RUBBLE-PILE ASTEROID COLLISION & FORMATION WITH:"
    print*, "    FORCES:"
    print*, "        * GRAVITATIONAL"
    print*, "    INTEGRATION SCHEME OPTIONS:"
    print*, "        * EULER"
    print*, "        * RK4 (RUNGE-KUTTA 4TH-ORDER)"
    print*, "    COLLISION SCHEME:"
    print*, "        IMPULSE-BASED W/ COEFF. OF RESTITUTION, AND"
    print*, "        'BOUNDING-SPHERE' COLLISION DETECTION, WITH"
    print*, "        VARIABLE TIMESTEP TO MINIMIZE MISSED COLLISIONS"
    print*, ""
    print*, "(c) Ian Friedrichs 2023"
    print*, ""

    ! create the asteroids
    print*, "[main_prog] adding asteroids..."

    call add_asteroid_ellipsoid(particle_list, 78e9, 280.0, 250.0, 250.0, 11.7, [0.0,0.0,0.0], [0.0,0.0,0.0], 10000, 0)
    call add_asteroid_ellipsoid(particle_list, 78e9, 280.0, 250.0, 250.0, 11.7, [700.0,0.0,0.0], [-10.0,0.0,0.0], 10000, 1)

    print*, "[main_prog] DONE adding asteroids."
    print*, "number of particles: ", size(particle_list)

    ! set up the output directory
    call setup_output_directory(OUT_DIR)

    total_time = 0.0
    call write_particle_list_for_paraview(particle_list, OUT_DIR, 0)
    i = 0 ! keep track of number of timesteps passed (counts contribs from colls)
    file_counter = 1

    eff_dt = calculate_next_dt(particle_list)
    if (eff_dt .gt. DT) then
      eff_dt = DT
    end if

    do while (total_time < MAX_NUM_TSTEPS * DT .and. total_time < MAX_TIME)
      call gravity_update_euler(particle_list, eff_dt)
      total_time = total_time + eff_dt
      i = i + 1
      if (mod(i, write_mod) .eq. 0) then
        call write_particle_list_for_paraview(particle_list, OUT_DIR, file_counter)
        file_counter = file_counter + 1
        print*, "Total sim time passed", total_time, "of", MAX_TIME
        print*, "Effective timestep:", eff_dt
      end if
      do bbox_count = 1,2,1
        call bbox_collisions(particle_list, COLLISION_RANGE)
      end do
      eff_dt = calculate_next_dt(particle_list)
      if (eff_dt .gt. DT) then
        eff_dt = DT
      end if
    end do
end program main_prog
