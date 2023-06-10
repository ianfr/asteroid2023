module gravity_module
    use type_module
    implicit none
    real, parameter :: BIG_G = 6.67430e-11
contains

! calculate a list of gravitational forces on each particle from each particle
subroutine calculate_forces(particle_list, force_list)
    ! Input/Output variables
    type(particle), dimension(:), intent(inout) :: particle_list
    real, dimension(:,:), allocatable, intent(out) :: force_list

    ! Internal variables for gravity
    integer :: i, j, n
    real :: Gmm
    real, dimension(3) :: r_vector, r_unit_vector, tmp_vector

    allocate( force_list(3,size(particle_list)) )
    force_list = 0

    n = size(particle_list)
    do i = 1, n, 1
        do j = 1, n, 1
            if (i .ne. j) then
                Gmm = BIG_G * particle_list(i)%mass  * particle_list(j)%mass
                r_vector = particle_list(i)%pos - particle_list(j)%pos
                r_unit_vector = r_vector / my_norm2(r_vector)
                tmp_vector = my_norm2(r_vector)**2
                
                force_list(:,i) = force_list(:,i) - 1.0 * (Gmm * r_unit_vector) / tmp_vector
            end if
        end do
    end do

end subroutine calculate_forces

! update the particles with motion uder gravity
! euler integration
! no COM speedups
subroutine gravity_update_euler(particle_list, timestep)
    ! Input variables
    type(particle), dimension(:), intent(inout) :: particle_list
    real, intent(in) :: timestep

    ! Internal variables
    integer :: i
    real, dimension(:,:), allocatable :: force_list

    ! subroutine

    ! force_list = calculate_forces(particle_list)
    call calculate_forces(particle_list, force_list)

    ! apply each calculated force to each particle for one euler step
    do i = 1, size(particle_list), 1
        particle_list(i)%vel = particle_list(i)%vel + (force_list(:,i) / particle_list(i)%mass) * timestep
        particle_list(i)%pos = particle_list(i)%pos + particle_list(i)%vel * timestep
    end do

end subroutine gravity_update_euler

subroutine gravity_update_rk4(particle_list, timestep)
    ! Input variables
    type(particle), dimension(:), intent(inout) :: particle_list
    real, intent(in) :: timestep

    ! Internal variables
    integer :: i
    real, dimension(:,:), allocatable :: force_list
    real, dimension(3) :: kv1, kv2, kv3, kv4, kr1, kr2, kr3, kr4, a, r, v
    real :: h

    ! Subroutine
    h = timestep
    call calculate_forces(particle_list, force_list)

    do i = 1, size(particle_list), 1
        r = particle_list(i)%pos
        v = particle_list(i)%vel
        a = force_list(:,i) / particle_list(i)%mass
        kr1 = v
        kv1 = a * r
        kr2 = v * kv1 * (h/2.0)
        kv2 = a * (r + kr1 * (h/2.0))
        kr3 = v * kv2 * (h/2.0)
        kv3 = a * (r + kr2 * (h/2.0))
        kr4 = v * kv3 * h
        kv4 = a * (r + kr3 * h)
        particle_list(i)%vel = v + (h/6.0) * (kv1 + 2*kv2 + 2*kv3 + kv4)
        particle_list(i)%pos = r + (h/6.0) * (kr1 + 2*kr2 + 2*kr3 + kr4)
    end do

end subroutine gravity_update_rk4


end module gravity_module