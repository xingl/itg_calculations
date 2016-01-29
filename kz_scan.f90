program itg_calculation

        implicit none

        real (kind=8), parameter :: pi=3.141592653589793
        complex (kind=8), parameter :: zi=(0.0,1.0)

        integer :: nvy, nvz, nkz
        integer :: iy, iz, ikz
        real (kind=8) :: lvy, lvz
        complex (kind=8) :: omega 
        real (kind=8) :: omn, omt, a_r, theta
        real (kind=8), dimension(:), allocatable :: kz
        real (kind=8) :: ky        
        real (kind=8) :: vi, omd
        integer :: nstep
        integer :: istep
        integer :: out_unit=101

        real (kind=8), dimension(:), allocatable :: vy_grid, vz_grid
        complex (kind=8), dimension(:), allocatable :: x, f
        complex (kind=8) :: x_tmp, f_tmp, om_kz
        real (kind=8) :: kz_tmp

        complex(kind=8) :: intgral

        nstep = 30
        nvy = 1000
        nvz = 2000
        nkz = 10
        lvy = 4.0
        lvz = 8.0
        omn = 0.5
        omt = 2.0
        omd = 0.0
        a_r = 3.0
        vi  = 1.0
        theta = 1.0

        allocate(vy_grid(nvy), vz_grid(nvz))

        do iy = 1, nvy
                vy_grid(iy) = (iy-1.0)*lvy/nvy
        end do

        do iz = 1, nvz
                vz_grid(iz) = -lvz/2.0 + (iz-1.0)*lvz/nvz
        end do
        
        ky = 0.3

        allocate(kz(nkz))
        do ikz = 1, nkz
                kz(ikz) = ikz*0.1
        end do

        allocate(x(nstep),f(nstep))
        
        open(unit=out_unit,file='scan_kz.dat',action="write")
        do ikz = 5, 10
                x = 0.0
                f = 0.0
                kz_tmp = kz(ikz)
                print*, kz_tmp
                x(1) = (0.5 + 0.1*zi)*omt
                x_tmp = x(1)
                call integrator(x_tmp,kz_tmp,f_tmp)
                f(1) = f_tmp
                x(2) = (0.5 - 0.1*zi)*omt
                x_tmp = x(2)
                call integrator(x_tmp,kz_tmp,f_tmp)
                f(2) = f_tmp
                print*, x(1), f(1)
                print*, x(2), f(2)

                do istep = 3, nstep
                        x(istep) = x(istep-1) - f(istep-1) * &
                                (x(istep-1)-x(istep-2))/(f(istep-1)-f(istep-2))
                        x_tmp = x(istep)
                        call integrator(x_tmp,kz_tmp,f_tmp)
                        f(istep) = f_tmp
                        print*, x(istep), f(istep)
                        if (abs(f(istep))<0.000001) then
                                om_kz = x(istep)
                                exit
                        end if
                end do
                write(out_unit,'(3e12.4)') kz(ikz), om_kz 
        end do
        close(unit=out_unit)
contains

  subroutine integrator (omega, kz_tmp, intgral)

        implicit none

        complex (kind=8), intent(in) :: omega
        real (kind=8), intent(in) :: kz_tmp
        complex (kind=8), intent(out) :: intgral
        complex (kind=8) :: intgral_tmp
        
        intgral_tmp = 0.0

        do iy = 1, nvy
                do iz = 1, nvz
                        intgral_tmp = intgral_tmp + &
                        2.0*sqrt(1.0/pi)*vy_grid(iy)*exp(-vy_grid(iy)**2-&
                        (vz_grid(iz)**2-zi*2.0*vz_grid(iz)*vi-&
                        vi**2))*(omega-omn*ky-&
                        (-3.0/2.0+vy_grid(iy)**2+&
                        (vz_grid(iz)**2-zi*2.0*vz_grid(iz)*vi-vi**2))&
                        *omt*ky)/&
                        (omega-2.0*kz_tmp*(vz_grid(iz)-zi*vi)&
                        -(vy_grid(iy)**2+2*&
                        (vz_grid(iz)**2-zi*2.0*vz_grid(iz)*vi-vi**2))&
                        *omd*ky/a_r)*lvy*lvz/nvy/nvz
                end do  
        end do
        
        intgral = 1+theta-theta*intgral_tmp

  end subroutine
end program
