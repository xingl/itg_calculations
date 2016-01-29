program itg_calculation

        implicit none

        real (kind=8), parameter :: pi=3.141592653589793
        complex (kind=8), parameter :: zi=(0.0,1.0)

        integer :: nvy, nvz
        integer :: iy, iz
        real (kind=8) :: lvy, lvz
        complex (kind=8) :: omega 
        real (kind=8) :: omn, omt, a_r, theta
        real (kind=8) :: ky, kz        
        real (kind=8) :: vi
        integer :: nstep
        integer :: istep
        integer :: out_unit=101

        real (kind=8), dimension(:), allocatable :: vy_grid, vz_grid
        complex (kind=8), dimension(:), allocatable :: x, f

        complex(kind=8) :: intgral

        nstep = 100
        nvy = 1000
        nvz = 2000
        lvy = 4.0
        lvz = 8.0
        omn = 0.5
        omt = 2.0
        a_r = 3.0
        theta = 1.0
        vi  = 1.0

        allocate(vy_grid(nvy), vz_grid(nvz))

        do iy = 1, nvy
                vy_grid(iy) = (iy-1.0)*lvy/nvy
        end do

        do iz = 1, nvz
                vz_grid(iz) = -lvz/2.0 + (iz-1.0)*lvz/nvz
        end do
        
        ky = 0.3
        kz = 0.5

        allocate(x(nstep),f(nstep))
        x = 0.0
        f = 0.0

        x(1) = (0.5 + 0.1*zi)*omt
        call integrator(x(1),f(1))
        x(2) = (0.5 - 0.1*zi)*omt
        call integrator(x(2),f(2))
        
        do istep = 3, nstep
                x(istep) = x(istep-1) - f(istep-1) * &
                (x(istep-1)-x(istep-2))/(f(istep-1)-f(istep-2))
                call integrator(x(istep),f(istep))
                print*, x(istep), f(istep)
        end do

        open(unit=out_unit,file='root_trace.dat',action="write")
        do istep = 1, nstep
                write(out_unit,'(4e12.4)') x(istep), f(istep)
        end do
        close(unit=out_unit)
contains

  subroutine integrator (omega, intgral)

        implicit none

        complex (kind=8), intent(in) :: omega
        complex (kind=8), intent(out) :: intgral
        complex (kind=8) :: intgral_tmp
        
        intgral_tmp = - omega**2+1
        
        intgral = intgral_tmp

  end subroutine
end program
