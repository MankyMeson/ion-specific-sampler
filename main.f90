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

!       dist = r**2*(exp(-2.d0*Z*r) + exp(-2*(r**2)*sqrt(2*I)))
        dist = r**2*exp((-2.d0*Z*r*(1+sqrt(2.d0*I)*r))/(1.d0+Z*r)) ! specifically the one discussed during the meeting

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
            f2 = dist(r2,Z,I)
            area = area + width*(f1+f2)/2.d0
        end do

        norm_coeff = 1/area

    end function


    real(dp) function rand_r(Z,I,norm)
        ! Samples r using a von Neumann acceptance-rejection method

        real(dp) :: I, cutoff_radius, r, norm, max_r
        real(dp), dimension(2) :: rands
        integer :: Z, iter

        ! First need to find the maxvalue of dist(r,Z,I)
        max_r = 0.d0
        cutoff_radius = 5.d0
        do i = 1,1001
            r = dble((i-1))*cutoff_radius/1000.d0
            if (dist(r,Z,I) > dist(max_r,Z,I)) then
                max_r = r
            end if
        end do

        do
            call random_number(rands)
            r = rands(1)*cutoff_radius
!           if (cutoff_radius*rands(2) < dist(r,Z,I)*norm) then
            if (max_r*rands(2) < dist(r,Z,I)) then
                exit
            end if
        end do

        rand_r = r

    end function


    subroutine histogram(dist,dist_size,resolution,hist)

        integer, intent(in) :: dist_size
        integer, dimension(:), allocatable :: hist_count
        integer :: j, hist_elem, hist_size
        real(dp), dimension(dist_size), intent(in) :: dist
        real(dp), intent(in) :: resolution
        real(dp), dimension(:), allocatable, intent(out) :: hist

        hist_size = ceiling(maxval(dist)/resolution)
        allocate (hist_count(hist_size))
        hist_count(:) = 0

        do j = 1, dist_size
            hist_elem = ceiling(dist(j)/resolution)
            hist_count(hist_elem) = hist_count(hist_elem) + 1
        end do

        allocate (hist(hist_size))

        do j = 1, hist_size
            hist(j) = dble(hist_count(j)) / (dist_size*resolution)
        end do

        deallocate(hist_count)

    end subroutine


end module sampler



program ionspecificsampler

    use sampler

    implicit none
    integer :: ele, n_ele, Z_ion, j
    real(dp) :: phi, cos_theta, I, C, r, res
    real(dp), dimension(:,:), allocatable :: ele_dist_cartesian, ele_dist_spherical, angular_rand_nums
    real(dp), dimension(:), allocatable :: ele_dist_r, r_histogram

    n_ele = 100000
    Z_ion = 87
    I = 0.5d0

    C = norm_coeff(Z_ion,I)

    allocate (angular_rand_nums(2,n_ele))
    call random_number(angular_rand_nums)

    allocate (ele_dist_spherical(3,n_ele))
    allocate (ele_dist_r(n_ele))

    do ele = 1, n_ele
        phi = rand_phi(angular_rand_nums(1,ele))
        cos_theta = rand_cos_theta(angular_rand_nums(2,ele))
        r = rand_r(Z_ion,I,C)

!       print*, r, acos(cos_theta), phi

        ele_dist_r(ele) = r

        ele_dist_spherical(1,ele) = r
        ele_dist_spherical(2,ele) = acos(cos_theta)
        ele_dist_spherical(3,ele) = phi
    end do

    deallocate (angular_rand_nums)
    deallocate (ele_dist_spherical)

    res = 0.01d0
    call histogram(ele_dist_r(1), n_ele, res, r_histogram)

    deallocate (ele_dist_r)

    open (file="histogramdata", unit=8, status="replace")
        do j = 1,size(r_histogram)
            write (8,"(2f20.16)") (j - 0.5d0)*res, r_histogram(j)
        end do
    close (8)

    deallocate (r_histogram)

end program ionspecificsampler
