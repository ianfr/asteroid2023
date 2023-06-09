program main_prog

    ! IMPORTS
    use asteroid_module
    use type_module
    use readparams_module
    use utilities_module
    use gravity_module
    use collision_module

    ! VARIABLES
    implicit none
    integer :: i, fileCounter, bbox_count
    real :: total_time, accum_coll_time, eff_dt

    ! Variables for asteroid creation
    type(particle), allocatable :: particle_list(:)
    integer :: ast_ind
    integer :: num_asteroids
    type(logical), dimension(:,:), allocatable :: collision_matrix

    ! Varibles that are hard-coded parameters for now
    real :: PARTICLE_RADIUS
    real :: DT ! DT now acts as an upper bound for timesteps

    ! Variables for program parameters
    character(len=*), parameter :: PARAM_FILE_NAME = "../config/control.txt"
    character(len=:), allocatable :: OUT_DIR
    integer :: NUM_PARTICLES, NUM_TIMESTEPS, WRITE_MOD
    real :: MAX_TIME
    ! for now, we just use 2 asteroids
    real, dimension(:), allocatable :: ASTEROID_MASSES
    real, dimension(:), allocatable :: ASTEROID_RADII
    real, dimension(:,:), allocatable :: ASTEROID_POSITIONS ! every column is the initial asteroid position
    real, dimension(:,:), allocatable :: ASTEROID_VELOCITIES ! every column is the initial asteroid velocity

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
    print*, "(c) Ian Friedrichs 2022"
    print*, ""

    call read_params_from_file(PARAM_FILE_NAME, NUM_PARTICLES, NUM_TIMESTEPS, &
      MAX_TIME, ASTEROID_MASSES, ASTEROID_RADII, ASTEROID_POSITIONS, ASTEROID_VELOCITIES, &
      DT, PARTICLE_RADIUS, OUT_DIR, WRITE_MOD)

    ! create the asteroids
    print*, "[main_prog] adding asteroids..."
    num_asteroids = size(ASTEROID_MASSES)
    ! do ast_ind = 1, num_asteroids, 1
    !   call add_asteroid(particle_list, ASTEROID_MASSES(ast_ind), ASTEROID_RADII(ast_ind), &
    !     PARTICLE_RADIUS, ASTEROID_POSITIONS(:,ast_ind), ASTEROID_VELOCITIES(:,ast_ind), ast_ind) ! modifies particle_list
    ! end do
    ! based on Bennu
    call add_asteroid_ellipsoid(particle_list, 78e9, 280.0, 250.0, 250.0, 11.0, [0.0,0.0,0.0], [0.0,0.0,0.0], 10000, 0)
    call add_asteroid_ellipsoid(particle_list, 78e9, 280.0, 250.0, 250.0, 11.0, [700.0,0.0,0.0], [-10.0,0.0,0.0], 10000, 1)

    print*, "[main_prog] DONE adding asteroids."
    print*, "number of particles: ", size(particle_list)

    ! allocate the logical collision matrix to avoid memory realloc in bbox_collisions
    allocate(collision_matrix(size(particle_list), size(particle_list)))
    collision_matrix = .false.

    ! set up the output directory
    call setup_output_directory(OUT_DIR)

    total_time = 0.0
    accum_coll_time = 0.0
    call write_particle_list_for_paraview(particle_list, OUT_DIR, 0)
    i = 0 ! keep track of number of timesteps passed (counts contribs from colls)
    fileCounter = 1

    eff_dt = calculate_next_dt(particle_list)
    if (eff_dt .gt. DT) then
      eff_dt = DT
    end if

    do while (total_time < NUM_TIMESTEPS * DT .and. total_time < MAX_TIME)
      call gravity_update_euler(particle_list, eff_dt)
      total_time = total_time + eff_dt
      i = i + 1
      if (mod(i, write_mod) .eq. 0) then
        call write_particle_list_for_paraview(particle_list, OUT_DIR, fileCounter)
        fileCounter = fileCounter + 1
        print*, "Total sim time passed", total_time, "of", MAX_TIME
        print*, "effective dt:", eff_dt
      end if
      do bbox_count = 1,2,1
        call bbox_collisions(particle_list, collision_matrix)
      end do
      eff_dt = calculate_next_dt(particle_list)
      if (eff_dt .gt. DT) then
        eff_dt = DT
      end if
    end do
end program main_prog
