! based on https://gitlab.com/ifriedri/advanced-asteroid-sim/-/blob/master/aaparticle.h
module collision_module
    use type_module
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

contains

! adapted from https://physics.stackexchange.com/questions/681396/elastic-collision-3d-eqaution
subroutine collide(part1, part2)
    type(particle), intent(inout) :: part1, part2

    real :: m1, m2, m_eff, epsilon, v_imp, j
    real, dimension(3) :: pos1, pos2, v1, v2
    real, dimension(3) :: n, dv1, dv2

    m1 = part1%mass
    m2 = part2%mass
    pos1 = part1%pos
    pos2 = part2%pos
    v1 = part1%vel
    v2 = part2%vel

    ! coeff of restitution, set to 1 for perfectly elastic
    epsilon = 0.86 ! https://www.researchgate.net/publication/222983190_A_coefficient_of_restitution_of_rock_materials

    n = (pos2 - pos1) / norm2(pos2 - pos1)
    m_eff = 1.0 / ( (1.0/m1) + (1.0/m2) )
    v_imp = dot_product(n, v1 - v2)
    j = (1 + epsilon) * m_eff * v_imp

    dv1 = -(j/m1) * n
    dv2 =  (j/m2) * n
    
    part1%vel = part1%vel + dv1
    part2%vel = part2%vel + dv2

end subroutine

! fast-forward every particle without the influence of gravity until the next collision
    ! time, which is supplied as an argument
subroutine fast_forward(particle_list, DT)
    type(particle), dimension(:), intent(inout) :: particle_list
    real, intent(in) :: dt

    ! local variables
    integer :: i

    ! subroutine
    do i = 1, size(particle_list), 1
        particle_list(i)%pos = particle_list(i)%pos + particle_list(i)%vel * dt
    end do

end subroutine fast_forward

! if any particles are overlapping / very close, calculate the collision between them
! MODIFIES particle_list
subroutine bbox_collisions(particle_list, separation_buffer)
    type(particle), dimension(:), intent(inout) :: particle_list
    real, intent(in) :: separation_buffer
    integer :: i, j

    do i = 1, size(particle_list), 1
        do j = 1, size(particle_list), 1
            if (i .ne. j) then
                 if (norm2(particle_list(i)%pos - particle_list(j)%pos) <= &
                        particle_list(i)%radius + particle_list(j)%radius + separation_buffer) then
                    call collide(particle_list(i), particle_list(j))
                 end if
            end if
        end do
    end do

end subroutine bbox_collisions

! calculate a new timestep s.t.: (v_max + v_second_max) * new_dt < 2*particle_radius
! ASSUMES uniform radii, but it's also just a heuristic so not a huge deal, could also precompute average radius size
real function calculate_next_dt_max(particle_list) result(new_dt)
    type(particle), dimension(:), intent(in) :: particle_list

    integer :: i
    real, dimension(:), allocatable :: v_list 
    real :: v_max, v_second_max

    allocate(v_list(size(particle_list)))

    do i = 1, size(particle_list), 1
        v_list(i) = norm2(particle_list(i)%vel)
    end do

    v_max = maxval(v_list)
    v_second_max = maxval(v_list, mask = v_list .le. v_max)

    new_dt = (2.0 * particle_list(1)%radius) / ((v_max + v_second_max))

end function

! calculate a new timestep s.t.: v_mean * new_dt < mean_particle_radius
real function calculate_next_dt_mean(particle_list) result(new_dt)
    type(particle), dimension(:), intent(in) :: particle_list

    integer :: i
    real :: v_mean
    real :: radius_mean

    ! Avoiding sum()/size() to prevent overflows
    v_mean = 0
    radius_mean = 0
    do i = 1, size(particle_list), 1
        v_mean = v_mean + (norm2(particle_list(i)%vel) / size(particle_list))
        radius_mean = radius_mean + (particle_list(i)%radius / size(particle_list))
    end do

    ! subtract a small amount to keep the inequality true
    new_dt = (radius_mean / v_mean) - 1e-9

end function

end module collision_module