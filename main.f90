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

    function orbital_gen(ele)
        integer :: ele, i, n, l, nplusl
        integer, dimension(2) :: orbital_gen
        ! given the number of the electron used, this function
        ! returns the orbital in the form (n,l)
        i = ele
        n = 1
        l = 0
        nplusl = 1
        do nplusl = 1,10
            do n = ((nplusl/2)+1),nplusl
                l = nplusl - n
                i = i - (4*l + 2)
                if (i<=0) then
                    exit
                end if
            end do
        end do
!       do n = 1,10
!           i = i - n_degen(n)
!           if (i<=0) then
!               do l = n-1,0,-1
!                   i = i + 4*l + 2
!                   if (i>0) then
!                       exit
!                   end if
!               end do
!               exit
!           end if
!       end do
        orbital_gen = (/n,l/)
    end function

end module sampler


program ionspecificsampler
    use sampler
    implicit none
    integer :: ele, n_ele, Z_ion, l, max_l, n
    integer, dimension(2) :: orbital
    real(dp) :: phi, cos_theta, sin_theta, zeta
    real(dp), dimension(:,:), allocatable :: ele_dist_cartesian, ele_dist_spherical, angular_rand_nums

    n_ele = 100
    Z_ion = 2

    max_l = 3 ! highest l quantum number found in the system

    allocate (angular_rand_nums(2,n_ele))
    call random_number(angular_rand_nums)

    allocate (ele_dist_spherical(3,n_ele))

    do ele = 1, n_ele
        phi = rand_phi(angular_rand_nums(1,ele))
        cos_theta = rand_cos_theta(angular_rand_nums(2,ele))

        ! peaks at n/zeta, zeta is the effective Z_ion
        orbital = orbital_gen(ele)
        n = orbital(1)
        l = orbital(2)
        if (orbital(1)==1.and.orbital(2)==0) then
            zeta = Z_ion - 0.3d0 * ele
        end if
        print *, n, l

    end do

    deallocate (angular_rand_nums)

    deallocate (ele_dist_spherical)

end program ionspecificsampler
