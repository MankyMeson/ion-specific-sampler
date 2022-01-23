module sampler
    implicit none
    integer, parameter :: dp = kind(1.d0)

contains

    real(dp) function rand_phi(uniform_rand)
        real(dp) :: uniform_rand
        real(dp), parameter :: pi = 4.d0*atan(1.d0)
        rand_phi = 2.d0 * pi * uniform_rand
    end function

    real(dp) function rand_cos_theta(uniform_rand)
        real(dp) :: uniform_rand
        rand_cos_theta = 2.d0 * uniform_rand - 1.d0
    end function

    integer function n_degen(n)
        integer :: n, l, d
        d = 0
        do l = 0,n-1
            d = d + 4*l + 2
        end do
        n_degen = d
    end function

    real(dp) function dist(r,Z,I)
        ! this function controls the distribution being sampled
        real(dp) :: r, I
        integer :: Z
        dist = r**2*(exp(-2.d0*Z*r) + exp(-2*(r**2)*sqrt(2*I)))
    end function

    real(dp) function norm_coeff(Z,I)
        ! Trapezoid integrator to find Z and I specific normalisation coefficient
        ! the integrand is asymptotic to 0 in the positive r direction so a suitable
        ! cutoff can be chosen

        real(dp) :: I, r1, r2, area, f1, f2, width
        real(dp), parameter :: cutoff_r = 100.d0
        integer :: Z, trpzd
        integer, parameter :: nr_trapezoids = 1000

        width = cutoff_r/nr_trapezoids

        area = 0.d0
        do trpzd = 1,nr_trapezoids
            r1 = (trpzd-1)*width
            r2 = trpzd*width
            f1 = dist(r1,Z,I)
            f2 = dist(r1,Z,I)
            area = area + width*(f1+f2)/2
        end do

        norm_coeff = 1/area
    end function

    real(dp) function rand_r(Z,I,norm)
        ! Samples r using a von Neumann acceptance-rejection method
        real(dp) :: I, cutoff_radius, r, norm
        real(dp), dimension(2) :: rands
        integer :: Z

        do
            call random_number(rands)
            cutoff_radius = 5.d0
            r = rands(1)*cutoff_radius
            if (cutoff_radius*rands(2) < dist(r,Z,I)*norm) then
                exit
            end if
        end do

        rand_r = r
    end function

end module sampler


program ionspecificsampler
    use sampler
    implicit none
    integer :: ele, n_ele, Z_ion
    real(dp) :: phi, cos_theta, I, C, r
    real(dp), dimension(:,:), allocatable :: ele_dist_cartesian, ele_dist_spherical, angular_rand_nums

    n_ele = 10000000
    Z_ion = 87
    I = 0.5d0

    C = norm_coeff(Z_ion,I)

    allocate (angular_rand_nums(2,n_ele))
    call random_number(angular_rand_nums)

    allocate (ele_dist_spherical(3,n_ele))

    do ele = 1, n_ele
        phi = rand_phi(angular_rand_nums(1,ele))
        cos_theta = rand_cos_theta(angular_rand_nums(2,ele))
        r = rand_r(Z_ion,I,C)

        print*, r, acos(cos_theta), phi

        ele_dist_spherical(1,ele) = r
        ele_dist_spherical(2,ele) = acos(cos_theta)
        ele_dist_spherical(3,ele) = phi

    end do

    deallocate (angular_rand_nums)

    deallocate (ele_dist_spherical)

end program ionspecificsampler
