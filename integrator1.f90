program itg_calculation

        implicit none

        real (kind=8), parameter :: pi=3.141592653589793
        complex (kind=8), parameter :: zi=(0.0,1.0)

        integer :: nvy, nvz
        integer :: iy, iz
        real (kind=8) :: lvy, lvz
        complex (kind=8) :: omega 
        real (kind=8) :: omn, omt, a_r
        real (kind=8) :: ky, kz        

        real (kind=8), dimension(:), allocatable :: vy_grid, vz_grid

        complex(kind=8) :: intgral

        nvy = 1000
        nvz = 2000
        lvy = 4.0
        lvz = 8.0
        omn = 0.5
        omt = 2.0
        a_r = 3.0

        allocate(vy_grid(nvy), vz_grid(nvz))

        do iy = 1, nvy
                vy_grid(iy) = (iy-1.0)*lvy/nvy
        end do

        do iz = 1, nvz
                vz_grid(iz) = -lvz/2.0 + (iz-1.0)*lvz/nvz
        end do
        
        ky = 0.3
        kz = 0.1

        omega = (0.5 + 0.1*zi)*omt
        intgral = 0.0

        do iy = 1, nvy
                do iz = 1, nvz
                        intgral = intgral + &
                        2.0*sqrt(1.0/pi)*vy_grid(iy)*exp(-vy_grid(iy)**2-&
                        vz_grid(iz)**2)*(omega-omn*ky-&
                        (-3.0/2.0+vy_grid(iy)**2+vz_grid(iz)**2)*omt*ky)/&
                        (omega-2.0*kz*vz_grid(iz))*lvy*lvz/nvy/nvz
                end do  
        end do

        print*, 2-intgral

        omega = (0.5 - 0.1*zi)*omt
        intgral = 0.0

        do iy = 1, nvy
                do iz = 1, nvz
                        intgral = intgral + &
                        2.0*sqrt(1.0/pi)*vy_grid(iy)*exp(-vy_grid(iy)**2-&
                        vz_grid(iz)**2)*(omega-omn*ky-&
                        (-3.0/2.0+vy_grid(iy)**2+vz_grid(iz)**2)*omt*ky)/&
                        (omega-2.0*kz*vz_grid(iz))*lvy*lvz/nvy/nvz
                end do  
        end do

        print*, 2-intgral


end program
