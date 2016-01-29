program itg_calculation

        use file_utils

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
        complex (kind=8) :: x_tmp, f_tmp, om_kz
        real (kind=8) :: om1_re, om1_im, om2_re, om2_im

        complex(kind=8) :: intgral

        call init_file_utils
        call read_input_file

        allocate(vy_grid(nvy), vz_grid(nvz))

        do iy = 1, nvy
                vy_grid(iy) = (iy-1.0)*lvy/nvy
        end do

        do iz = 1, nvz
                vz_grid(iz) = -lvz/2.0 + (iz-1.0)*lvz/nvz
        end do
        
        allocate(x(nstep),f(nstep))
        x = 0.0
        f = 0.0

        x(1) = om1_re + om1_im*zi
        call integrator(x(1),f(1))
        print*, x(1),f(1)
        x(2) = om2_re + om2_im*zi
        call integrator(x(2),f(2))
        print*, x(2),f(2)
        
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
        print*, om_kz
        
        call finish_file_utils
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

  subroutine read_input_file

        use file_utils, only: input_unit_exist, input_unit

        implicit none

        integer :: in_file
        logical :: exist

        namelist / parameters / nstep,nvy,nvz,lvy,lvz,omn,omt,a_r,theta,&
                                vi,ky,kz,om1_re,om1_im,om2_re,om2_im
        
        nstep = 100
        nvy = 1000
        nvz = 2000
        lvy = 4.0
        lvz = 8.0
        omn = 0.0
        omt = 0.0
        a_r = 3.0
        theta = 1.0
        vi  = 1.0
        ky = 0.3
        kz = 0.5
        om1_re = 1.0
        om1_im = 0.1
        om2_re = 1.0
        om2_im = 0.1

    in_file = input_unit_exist ("parameters", exist)
    if (exist) read (unit=input_unit("parameters"), nml=parameters)

  end subroutine read_input_file
end program
