program itg_calculation

        use file_utils

        implicit none

        real (kind=8), parameter :: pi=3.141592653589793
        complex (kind=8), parameter :: zi=(0.0,1.0)

        integer :: nkz, nky 
        integer :: nfprime, ntprime, nomd
        integer :: iy, iz, ikz, iky
        integer :: ifprime, itprime, iomd
        real (kind=8) :: vy_1, vy_2, vz_1,vz_2, dkz, dky
        real (kind=8) :: dfprime, dtprime, domd
        complex (kind=8) :: omega 
        real (kind=8) :: A, theta, fprime, tprime, omd, ky, kz
        real (kind=8) :: vi, fprime_start, ky_start,tol
        integer :: nstep
        integer :: istep
        integer :: out_unit=101

        real (kind=8), dimension(:), allocatable :: ky_grid, kz_grid
        real (kind=8), dimension(:), allocatable :: fprime_grid, tprime_grid, omd_grid

        complex (kind=8) :: root
        real (kind=8) :: om1_re, om1_im, om2_re, om2_im
        complex (kind=8) :: seed1, seed2
        
        character(len=2) :: ky_direction, fprime_direction
        
        call init_file_utils
        call read_input_file

        allocate(ky_grid(nky), kz_grid(nkz))
        do ikz = 1, nkz
                kz_grid(ikz) = (nkz+1-ikz)*dkz
        end do
        if (ky_direction=='bw') then
                do iky = 1, nky
                        ky_grid(iky) = (nky+1-iky)*dky+ky_start
                end do
        else
                do iky = 1, nky
                        ky_grid(iky) = iky*dky+ky_start
                end do
                
        end if

        allocate(fprime_grid(nfprime), tprime_grid(ntprime), omd_grid(nomd))
        if (fprime_direction=='bw') then
                do ifprime = 1, nfprime
                        fprime_grid(ifprime)=(nfprime+1-ifprime)*dfprime+fprime_start
                end do 
        else
                do ifprime = 1, nfprime
                        fprime_grid(ifprime) = ifprime*dfprime+fprime_start
                end do 
        end if
                
        do itprime = 1, ntprime
                tprime_grid(itprime) = (ntprime+1-itprime)*dtprime
        end do
        do iomd = 1, nomd
                omd_grid(iomd) = (nomd+1-iomd)*domd
        end do

        call open_output_file(out_unit,'.dat')
        write (out_unit,'(7a12)') "tprime","fprime","omd","ky","kz","frequency","growth rate"
        do iomd = 1, nomd
           omd = omd_grid(iomd)
           do ifprime = 1, nfprime
              fprime = fprime_grid(ifprime)
              do itprime = 1, ntprime
                 tprime = tprime_grid(itprime)
                 print*, "tprime, fprime, omd"
                 print*, tprime, fprime, omd
                 do iky = 1, nky
                    do ikz = 1, nkz
                       print*, "ky, kz"
                       print*, ky_grid(iky), kz_grid(ikz)
                       ky = ky_grid(iky)
                       kz = kz_grid(ikz)
                       if (ifprime ==1 .and. itprime==1 .and. iky==1 .and. ikz==1) then
                            print*, "Two initial omega values:"
                            seed1 = om1_re + om1_im*zi
                            seed2 = om2_re + om2_im*zi
                            print*, seed1, seed2
                            call rootfinder(seed1,seed2,root)
                            print*, "Root is found:"
                            print*, root    
                            write (out_unit, '(7e12.4)') tprime,fprime,omd/A,ky,kz,root
                       else
                            print*, "Two initial omega values:"
                            seed1 = root 
                            seed2 = root*0.9
                            print*, seed1, seed2
                            call rootfinder(seed1,seed2,root)
                            print*, "Root is found:"
                            print*, root
                            write (out_unit, '(7e12.4)') tprime,fprime,omd/A,ky,kz,root
                       end if
                    end do
                 end do
              end do
           end do
        end do
        call close_output_file (out_unit)
        call finish_file_utils

contains

 subroutine rootfinder(a,b,c)

        implicit none

        complex (kind=8), intent(in) :: a,b
        complex (kind=8), intent(out) :: c
        complex (kind=8), dimension(:), allocatable :: x, f
        complex (kind=8) :: x_tmp, f_tmp

        allocate(x(nstep),f(nstep))
        x = 0.0
        f = 0.0
        x(1) = a
        call simpson_adaptive(x(1),f(1),tol)
        x(2) = b
        call simpson_adaptive(x(2),f(2),tol)

        do istep = 3, nstep
                x(istep) = x(istep-1) - f(istep-1) * &
                           (x(istep-1)-x(istep-2))/(f(istep-1)-f(istep-2))
                x_tmp = x(istep)
                call simpson_adaptive(x_tmp,f_tmp,tol)
                f(istep) = f_tmp
                print*, istep
                print*, x(istep), f(istep)
                if (abs(f(istep))<1.0E-04) then
                        c = x(istep)
                        exit
                end if
        end do
        deallocate(x,f)

  end subroutine

subroutine simpson_adaptive(omega,integral,tol)

        implicit none

        real(kind=8), intent(in) :: tol
        complex(kind=8), intent(in) :: omega
        complex (kind=8), intent(out) :: integral
        real(kind=8) :: x_left, x_center, x_right, dx
        complex (kind=8) :: f_left, f_center, f_right
        complex (kind=8) :: i_trapezoid, i_simpson
        real(kind=8) :: a,b

        a = vy_1
        b = vy_2
        dx = 0.1
        x_left = a
        x_right = x_left + dx
        call vz_integral(x_left,f_left,tol,omega)
        integral = 0.0

        do while (x_right<b)
                x_center = 0.5*(x_left+x_right)
                call vz_integral(x_center,f_center,tol,omega)
                call vz_integral(x_right,f_right,tol,omega)
                i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                do while (abs(i_trapezoid-i_simpson)>tol)
                        dx=0.5*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                        if (dx>0.2) dx=0.2
                        x_right = x_left + dx
                        x_center = 0.5*(x_left+x_right)
                        call vz_integral(x_center,f_center,tol,omega)
                        call vz_integral(x_right,f_right,tol,omega)
                        i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                        i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                end do
                integral = integral + i_simpson
                x_left = x_right
                dx = 0.5*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                if (dx>0.2) dx=0.2
                x_right = x_left + dx
                if (x_right>b) x_right = b
                f_left = f_right

        end do

        integral = 1.0 + theta - theta*integral

  end subroutine

  subroutine vz_integral(vy,integral,tol,omega)

        implicit none

        real(kind=8), intent(in) :: tol, vy
        complex(kind=8), intent(in) :: omega
        real(kind=8) :: a, b
        complex (kind=8), intent(out) :: integral
        real(kind=8) :: x_left, x_center, x_right, dx
        complex (kind=8) :: f_left, f_center, f_right
        complex (kind=8) :: i_trapezoid, i_simpson

        a = vz_1
        b = vz_2
        dx = 0.1
        x_left = a
        x_right = x_left + dx
        f_left = integrand(x_left,vy,omega)
        integral = 0.0

        do while (x_right<b)
                x_center = 0.5*(x_left+x_right)
                f_center = integrand(x_center,vy,omega)
                f_right = integrand(x_right,vy,omega)
                i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                do while (abs(i_trapezoid-i_simpson)>tol)
                        dx=0.5*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                        if (dx>0.2) dx=0.2
                        x_right = x_left + dx
                        x_center = 0.5*(x_left+x_right)
                        f_center = integrand(x_center,vy,omega)
                        f_right = integrand(x_right,vy,omega)
                        i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                        i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                end do
                integral = integral + i_simpson
                x_left = x_right
                dx = 0.5*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                if (dx>0.2) dx=0.2
                x_right = x_left + dx
                if (x_right>b) x_right = b
                f_left = f_right

        end do

  end subroutine

  function integrand(vz,vy,omega)

        real(kind=8) :: vz,vy
        complex(kind=8) :: omega,integrand,int_tmp
        complex(kind=8) :: omega_star_n, omega_star_t, omega_star_d
        complex(kind=8) :: omega_parallel
        real(kind=8) :: rho,bj,dj,fj
        integer :: n

        omega_star_n = fprime*ky
        int_tmp = 0.0
        rho = ky*vy
        n = 0
        call bjndd (n, rho, bj, dj, fj)
        omega_star_t = 0.5*(-3.0+vy**2+&
                        (vz**2-zi*2.0*vz*vi-vi**2))*tprime*ky
        omega_star_d = (0.5*vy**2+&
                        vz**2-zi*2.0*vz*vi-vi**2)*omd*ky/A
        omega_parallel = kz*(vz-zi*vi)

        int_tmp = int_tmp + bj**2*&
                        sqrt(1.0/pi/2.0)*vy*exp(-0.5*vy**2-&
                        0.5*(vz**2-zi*2.0*vz*vi-&
                        vi**2))*(omega-omega_star_n-omega_star_t)/&
                        (omega-omega_parallel-omega_star_d)
        integrand = int_tmp

  end function

  subroutine read_input_file

        use file_utils, only: input_unit_exist, input_unit

        implicit none

        integer :: in_file
        logical :: exist

        namelist / parameters / nstep,vy_1,vy_2,vz_1,vz_2,A,theta,&
                                vi,om1_re,om1_im,om2_re,om2_im,nkz,dkz, &
                                nky,dky,nfprime,dfprime,ntprime,dtprime, &
                                nomd,domd,ky_direction,fprime_direction,&
                                ky_start, fprime_start, tol
        
        nstep = 100
        vy_1 = 0.0
        vy_2 = 0.0
        vz_1 = 8.0
        vz_2 = 8.0
        dfprime = 0.5
        nfprime = 1.0
        A = 3.0
        theta = 1.0
        vi  = 1.0
        om1_re = 1.0
        om1_im = 0.1
        om2_re = 1.0
        om2_im = 0.1
        nkz = 1
        dkz = 0.1
        nky = 1
        dky = 0.1
        ntprime = 1
        dtprime = 1.0
        nomd = 1
        domd = 1.0
        ky_direction = 'fw'
        fprime_direction = 'bw'
        ky_start = 0.0
        fprime_start = 0.0
        tol = 1.0E-05

    in_file = input_unit_exist ("parameters", exist)
    if (exist) read (unit=input_unit("parameters"), nml=parameters)

  end subroutine read_input_file
end program
