program itg_calculation

        implicit none

        real (kind=8), parameter :: pi=3.141592653589793
        complex (kind=8), parameter :: zi=(0.0,1.0)

        integer :: nvy, nvz
        integer :: iy, iz
        real (kind=8) :: lvy, lvz
        complex (kind=8) :: omega 
        real (kind=8) :: omn, omt, a_r, theta
        real (kind=8), dimension(:), allocatable :: theta_arr
        integer :: itheta
        real (kind=8) :: ky, kz        
        real (kind=8) :: vi
        integer :: nstep
        integer :: istep
        integer :: out_unit=101

        real (kind=8), dimension(:), allocatable :: vy_grid, vz_grid
        complex (kind=8), dimension(:), allocatable :: x, f
        complex (kind=8) :: x_tmp, f_tmp, om_kz

        complex(kind=8) :: intgral

        nstep = 100
        nvy = 1000
        nvz = 2000
        lvy = 4.0
        lvz = 8.0
        omn = 0.0
        omt = 0.0
        a_r = 3.0
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
        
        allocate(theta_arr(10))
        theta_arr(1:10) = (/0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0/)

        allocate(x(nstep),f(nstep))
        
        open(unit=out_unit,file='ionacc_theta.dat',action="write")
        do itheta = 1, 10
                print*, theta_arr(itheta)
                theta = theta_arr(itheta)
                x = 0.0
                f = 0.0
                if (itheta==1) then
                        x(1) = (0.1 - 0.8*zi)
                        x_tmp = x(1)
                        call integrator(x_tmp,f_tmp)
                        f(1) = f_tmp
                        x(2) = (0.1 - 0.1*zi)
                        x_tmp = x(2)
                        call integrator(x_tmp,f_tmp)
                        f(2) = f_tmp
                        print*, x(1), f(1)
                        print*, x(2), f(2)
                else
                        x(1) = om_kz
                        x_tmp = x(1)
                        call integrator(x_tmp,f_tmp)
                        f(1) = f_tmp
                        x(2) = om_kz - 0.2*zi
                        x_tmp = x(2)
                        call integrator(x_tmp,f_tmp)
                        f(2) = f_tmp
                        print*, x(1), f(1)
                        print*, x(2), f(2)
                end if

                do istep = 3, nstep
                        x(istep) = x(istep-1) - f(istep-1) * &
                                (x(istep-1)-x(istep-2))/(f(istep-1)-f(istep-2))
                        x_tmp = x(istep)
                        call integrator(x_tmp,f_tmp)
                        f(istep) = f_tmp
                        print*, x(istep), f(istep)
                        if (abs(f(istep))<0.000001) then
                                om_kz = x(istep)
                                exit
                        end if
                end do
                write(out_unit,'(3e12.4)') theta_arr(itheta)**(-1),-aimag(om_kz)/real(om_kz)
        end do
        close(unit=out_unit)
contains

  subroutine integrator (omega, intgral)

        implicit none

        complex (kind=8), intent(in) :: omega
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
                        (omega-2.0*kz*(vz_grid(iz)-zi*vi))&
                        *lvy*lvz/nvy/nvz
                end do  
        end do
        
        intgral = 1+theta-theta*intgral_tmp

  end subroutine
end program
