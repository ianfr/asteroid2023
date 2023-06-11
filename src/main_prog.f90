program main_prog

    ! IMPORTS
    use asteroid_module
    use type_module
    use utilities_module
    use gravity_module
    use collision_module
    use timers_module

    ! VARIABLES
    implicit none
    integer :: i, file_counter, bbox_count
    real :: total_time, eff_dt
    type(particle), allocatable :: particle_list(:)
    real, dimension(:,:), allocatable :: prev_for ! prev. forces for verlet update
    type(Timer) :: chrono

    ! Program parameters
    real, parameter :: DT  = 0.15 ! DT is an upper bound for timesteps
    real, parameter :: DT_LOWER_B = 1e-5 ! absolute lower bound for the adaptive timestep
    character(len=*), parameter :: OUT_DIR = "timing"
    integer, parameter :: MAX_NUM_TSTEPS = 3
    integer, parameter :: WRITE_MOD = 10
    real, parameter :: MAX_TIME = 100
    real, parameter :: COLLISION_RANGE = 0.05 ! check with particle radii below to see if appropriate
    integer, parameter :: BBOX_COLL_ITER = 2 ! iterations for collision checking/handling

    ! Pick an integration scheme
    logical, parameter :: USE_EULER = .false.
    logical, parameter :: USE_VERLET = .true.
    logical, parameter :: USE_RK4 = .false.

    ! Pick an adaptive timestepping scheme
    logical, parameter :: USE_ATS_MAX = .false. ! use maximum particle velocity, assume uniform radius
    logical, parameter :: USE_ATS_MEAN = .true. ! use mean particle velocity and radius

    ! PROGRAM START
    print*, "SIMULATION OF RUBBLE-PILE ASTEROID COLLISION AND FORMATION WITH:"
    print*, "    FORCES:"
    print*, "        * GRAVITATIONAL"
    print*, "    INTEGRATION SCHEME OPTIONS:"
    print*, "        * EULER"
    print*, "        * VERLET (SYMPLECTIC 2ND-ORDER)"
    print*, "        * RK4 (RUNGE-KUTTA 4TH-ORDER)"
    print*, "    COLLISION SCHEME:"
    print*, "        * IMPULSE-BASED W/ COEFF. OF RESTITUTION, AND"
    print*, "            'BOUNDING-SPHERE' COLLISION DETECTION"
    print*, "    ADAPTIVE TIME-STEPPING:"
    print*, "        * BASED ON MAXIMUM PARTICLE VELOCITY"
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

    ! Write the initial state
    call write_particle_list_for_paraview(particle_list, OUT_DIR, 0)

    total_time = 0.0 ! simulation time passed
    i = 0 ! keep track of number of timesteps passed (counts contribs from colls)
    file_counter = 1 ! for making filenames

    ! Adaptive timestep
    if (USE_ATS_MAX) then
      eff_dt = calculate_next_dt_max(particle_list)
    else ! USE_ATX_MEAN
      eff_dt = calculate_next_dt_mean(particle_list)
    end if
    if (eff_dt .gt. DT) then
      eff_dt = DT
    end if
    if (eff_dt .lt. DT_LOWER_B) then
      eff_dt = DT_LOWER_B
    end if

    do while (total_time < MAX_NUM_TSTEPS * DT .and. total_time < MAX_TIME)

      ! Select the appropriate integrator
      if (USE_EULER) then
        call gravity_update_euler(particle_list, eff_dt)
      else if (USE_VERLET) then
        call gravity_update_verlet_pos(particle_list, eff_dt, prev_for)
        call gravity_update_verlet_vel(particle_list, eff_dt, prev_for)
      else
        call gravity_update_rk4(particle_list, eff_dt)
      end if

      ! Update sim time passed, # of iters
      total_time = total_time + eff_dt
      i = i + 1

      ! Conditionally write to disk and print progress
      if (mod(i, write_mod) .eq. 0) then
        call write_particle_list_for_paraview(particle_list, OUT_DIR, file_counter)
        file_counter = file_counter + 1
        print*, "Total sim time passed", total_time, "of", MAX_TIME
        print*, "Effective timestep:", eff_dt
      end if

      ! Iterated collision handling for resolving closely packed situations
      do bbox_count = 1,BBOX_COLL_ITER,1
        call bbox_collisions(particle_list, COLLISION_RANGE)
      end do

      ! Adaptive timestep
      if (USE_ATS_MAX) then
        eff_dt = calculate_next_dt_max(particle_list)
      else ! USE_ATX_MEAN
        eff_dt = calculate_next_dt_mean(particle_list)
      end if
      if (eff_dt .gt. DT) then
        eff_dt = DT
      end if
      if (eff_dt .lt. DT_LOWER_B) then
        eff_dt = DT_LOWER_B
      end if
    end do
end program main_prog
