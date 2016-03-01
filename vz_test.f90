program vz_res_integral

        use file_utils

        implicit none

        real (kind=8), parameter :: pi=3.141592653589793
        complex (kind=8), parameter :: zi=(0.0,1.0)
        integer (kind=8) :: nkz, nky, nkx
        real (kind=8) :: dkz, dky, dkx
        integer (kind=8) :: ikz, iky, ikx
        real (kind=8) :: ky_start, kz_start, kx_start
        real (kind=8) :: vy_1, vz_1, vz_2
        complex (kind=8) :: omega,integral1,integral2,integral3,integral4
        real (kind=8) :: A, fprime, tprime, omd, ky, kz, kx, omd_kx
        real (kind=8) :: vi, tol, fr
        integer (kind=8) :: nstep, steps
        integer :: gamma_unit=101, Dmixing_unit=102, out_unit=103
        real (kind=8), dimension(:), allocatable :: ky_grid, kz_grid, kx_grid
        real (kind=8) :: om1_re, om1_im
        integer :: s

        call init_file_utils
        call read_input_file
        call init_grids

        call open_output_file(gamma_unit,'.dat')
        write (gamma_unit,'(13a12)') 'kx','ky','kz','omega','','int1','','int2','','int3','','int4',''
        call open_output_file(out_unit,'.res')
        write (out_unit,'(13a12)') 'kx','ky','kz','omega','','v1','',&
                                        'v2','','res1','','res2',''

        do iky = 1, nky
           ky = ky_grid(iky)
           do ikz = 1, nkz
              kz = kz_grid(ikz)
              do ikx = 1, nkx
                 kx = kx_grid(ikx)
                 omega = om1_re + zi*0.1
                 do while (aimag(omega)>-0.1)
                     print*, omega
                     call vz_integral(omega,vy_1,integral1,s=1)
                     call vz_integral(omega,vy_1,integral2,s=2)
                     call vz_integral(omega,vy_1,integral3,s=3)
                     call vz_integral(omega,vy_1,integral4,s=4)
                     
                     write(gamma_unit,'(13e12.4)') kx,ky,kz,omega,&
                                integral1,integral2,integral3,integral4
                     omega = omega - 0.001*zi
                 end do
              end do
           end do
        end do
        call close_output_file(gamma_unit)
        call close_output_file(out_unit)
        call finish_file_utils

contains

  subroutine vz_integral(omega,vy,integral,s)

        implicit none

        integer, intent(in) :: s
        real(kind=8), intent(in) :: vy
        complex(kind=8), intent(in) :: omega
        real(kind=8) :: int_a, int_b
        complex (kind=8), intent(out) :: integral
        real(kind=8) :: x_left, x_center, x_right, dx
        complex (kind=8) :: f_left, f_center, f_right
        complex (kind=8) :: i_trapezoid, i_simpson
        complex (kind=8) :: res_v1, res_v2
        real(kind=8) :: rho,bj,dj,fj
        integer :: n
        real(kind=8) :: kv,v0,coeff_a,coeff_b,coeff_d,dvcr
        complex(kind=8) :: dv,v1,v2
        complex(kind=8) :: coeff_c, coeff_e, coeff_f, coeff_g, tmp_f, tmp_g
        complex(kind=8) :: tmp_res1, tmp_res2

        n = 0
        rho = sqrt(ky**2+kx**2)*vy
        call bjndd (n, rho, bj, dj, fj)
        kv = (1.0-omd_kx)*ky + omd_kx*kx
        coeff_a = kv/A
        coeff_b = kz
        coeff_c = 0.5*kv/A*vy**2-omega
        coeff_d = 0.5*tprime*ky
        coeff_e = (-3./2.+vy**2/2.)*tprime*ky+fprime*ky-omega 
        v0 = -coeff_b/2./coeff_a
        dv = sqrt(coeff_b**2-4.*coeff_a*coeff_c)/coeff_a
        dvcr = coeff_b**2-4.0*coeff_a*(0.5*kv/A*vy**2-real(omega))
        print*, coeff_a,coeff_b,coeff_c,coeff_d,coeff_e,v0,dv,dvcr,vy
        if (coeff_a>0) then
                v1 = v0+0.5*dv
                v2 = v0-0.5*dv
        else
                v1 = v0-0.5*dv
                v2 = v0+0.5*dv
        end if
        coeff_f = (coeff_d/coeff_a*(coeff_b*v1+coeff_c)-coeff_e)/(v2-v1)
        coeff_g = (-coeff_d/coeff_a*(coeff_b*v2+coeff_c)+coeff_e)/(v2-v1)

        tmp_res1 = 2.*pi*zi/coeff_a*coeff_f
        tmp_res2 = 2.*pi*zi/coeff_a*coeff_g

        res_v1 =  sqrt(1.0/pi/2.0)*vy*exp(-0.5*vy**2-&
                0.5*v1**2)*tmp_res1
        res_v2 =  sqrt(1.0/pi/2.0)*vy*exp(-0.5*vy**2-&
                0.5*v2**2)*tmp_res2
        print*, coeff_f,coeff_g, res_v1,res_v2
        int_a = vz_1
        int_b = vz_2
        dx = 0.1
        x_left = int_a
        x_right = x_left + dx
        call integrand(x_left,vy,omega,v1,f_left,dvcr,s)
        integral = 0.0

        do while (x_left<int_b)
                x_center = 0.5*(x_left+x_right)
                call integrand(x_center,vy,omega,v1,f_center,dvcr,s)
                call integrand(x_right,vy,omega,v1,f_right,dvcr,s)
                i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                do while (abs(i_trapezoid-i_simpson)>tol)
                        dx=fr*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                        if (dx>0.2) dx=0.2
                        x_right = x_left + dx
                        if (x_right>int_b) then
                                x_right = int_b
                                dx = x_right - x_left
                        end if
                        x_center = 0.5*(x_left+x_right)
                        call integrand(x_center,vy,omega,v1,f_center,dvcr,s)
                        call integrand(x_right,vy,omega,v1,f_right,dvcr,s)
                        i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                        i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                end do
                integral = integral + i_simpson
                x_left = x_right
                dx = fr*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                if (dx>0.2) dx=0.2
                x_right = x_left + dx
                if (x_right>int_b) then 
                        x_right = int_b
                        dx = x_right - x_left
                end if
                f_left = f_right



        end do
        if (s==1) then
           if (real(dvcr)<0.) then 
              integral = bj**2*integral
           else 
              if (aimag(v1)>0.1) then
                 integral = bj**2*integral
              else if (aimag(v1)>-0.1) then
                 integral = bj**2*(integral - res_v2)
              else
                 integral = bj**2*(integral + res_v1 - res_v2)
              end if
           end if
        end if
        if (s==2) integral = bj**2*(integral + res_v1)
        if (s==3) integral = bj**2*(integral - res_v2)
        if (s==4) then 
           !if (aimag(v1)>0) integral = bj**2*integral
           !if (aimag(v1)<0) integral = bj**2*(integral + res_v1 - res_v2)
           integral = bj**2*integral
        end if

        write(out_unit,'(13e12.4)') kx,ky,kz,omega,v1,v2,res_v1,res_v2

  end subroutine

  subroutine integrand(vz,vy,omega,v1,int_tmp,cr,s)

        real(kind=8),intent(in) :: vz,vy,cr
        complex(kind=8),intent(in) :: omega,v1
        integer, intent(in) :: s
        complex(kind=8),intent(out) :: int_tmp
        complex(kind=8) :: omega_star_n, omega_star_t, omega_star_d
        complex(kind=8) :: omega_parallel
        complex(kind=8) :: vzprime

        if (s==1) then
           if (cr <0.) then
                vzprime = vz
           else 
                if (aimag(v1)>0.1) then
                        vzprime = vz
                else if (aimag(v1)>-0.1) then
                        vzprime = vz-zi*0.2
                else
                        vzprime = vz
                end if
           end if
                omega_star_n = fprime*ky
                omega_star_t = 0.5*(-3.0+vy**2+&
                        vzprime**2)*tprime*ky
                omega_star_d = (1.0-omd_kx)*(0.5*vy**2+&
                        vzprime**2)*omd*ky/A +&
                        omd_kx*(0.5*vy**2+vzprime**2)*&
                        omd*kx/A
                omega_parallel = kz*vzprime

                int_tmp = sqrt(1.0/pi/2.0)*vy*exp(-0.5*vy**2-&
                          0.5*vzprime**2)*&
                          (omega-omega_star_n-omega_star_t)/&
                          (omega-omega_parallel-omega_star_d)



        else if (s==2) then
                vzprime = vz+zi*vi
                omega_star_n = fprime*ky
                omega_star_t = 0.5*(-3.0+vy**2+vzprime**2)*tprime*ky
                omega_star_d = (1.0-omd_kx)*(0.5*vy**2+&
                        vzprime**2)*omd*ky/A +&
                        omd_kx*(0.5*vy**2+vzprime**2)*&
                        omd*kx/A
                omega_parallel = kz*vzprime
                int_tmp = sqrt(1.0/pi/2.0)*vy*exp(-0.5*vy**2-&
                  0.5*vzprime**2)*&
                  (omega-omega_star_n-omega_star_t)/&
                  (omega-omega_parallel-omega_star_d)

        else if (s==3) then
                vzprime = vz-zi*vi
                omega_star_n = fprime*ky
                omega_star_t = 0.5*(-3.0+vy**2+vzprime**2)*tprime*ky
                omega_star_d = (1.0-omd_kx)*(0.5*vy**2+&
                        vzprime**2)*omd*ky/A +&
                        omd_kx*(0.5*vy**2+vzprime**2)*&
                        omd*kx/A
                omega_parallel = kz*vzprime
                int_tmp = sqrt(1.0/pi/2.0)*vy*exp(-0.5*vy**2-&
                  0.5*vzprime**2)*&
                  (omega-omega_star_n-omega_star_t)/&
                  (omega-omega_parallel-omega_star_d)

        else if (s==4) then
                vzprime = vz
                omega_star_n = fprime*ky
                omega_star_t = 0.5*(-3.0+vy**2+vzprime**2)*tprime*ky
                omega_star_d = (1.0-omd_kx)*(0.5*vy**2+&
                        vzprime**2)*omd*ky/A +&
                        omd_kx*(0.5*vy**2+vzprime**2)*&
                        omd*kx/A
                omega_parallel = kz*vzprime
                int_tmp = sqrt(1.0/pi/2.0)*vy*exp(-0.5*vy**2-&
                  0.5*vzprime**2)*&
                  (omega-omega_star_n-omega_star_t)/&
                  (omega-omega_parallel-omega_star_d)
        end if
  end subroutine

  subroutine read_input_file

        use file_utils, only: input_unit_exist, input_unit

        implicit none

        integer :: in_file
        logical :: exist

        namelist / parameters / nstep,nkz,dkz,kz_start, &
                                nky,dky,ky_start, &
                                vy_1,vz_1,vz_2,om1_re,om1_im, &
                                A,vi,tol,fr, &
                                nkx, dkx, kx_start, &
                                tprime, fprime, omd, omd_kx
        
        nstep = 10

        nkz = 1
        dkz = 0.1
        kz_start = 0.0

        nky = 1
        dky = 0.1
        ky_start = 0.0

        vy_1 = 0.0
        vz_1 = 8.0
        vz_2 = 8.0

        om1_re = 1.0
        om1_im = 0.1

        A = 3.0
        vi  = 1.0
        tol = 1.0E-05
        fr = 0.5
        
        nkx = 1
        dkx = 0.5
        kx_start = 0.0

        tprime = 15.
        fprime = 7.5
        omd = 1.0
        omd_kx = 1.0
        
    in_file = input_unit_exist ("parameters", exist)
    if (exist) read (unit=input_unit("parameters"), nml=parameters)

  end subroutine read_input_file

subroutine init_grids

        implicit none

        allocate(ky_grid(nky), kz_grid(nkz), kx_grid(nkx))

        do ikz = 1, nkz
                kz_grid(ikz) = ikz*dkz+kz_start
        end do
        do iky = 1, nky
                ky_grid(iky) = iky*dky+ky_start
        end do
        do ikx = 1, nkx
                kx_grid(ikx) = ikx*dkx+kx_start
        end do

end subroutine

end program
