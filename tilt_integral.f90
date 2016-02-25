program imp_tor_itg

        use file_utils

        implicit none

        real (kind=8), parameter :: pi=3.141592653589793
        complex (kind=8), parameter :: zi=(0.0,1.0)

        integer (kind=8) :: nkz, nky, nkx
        real (kind=8) :: dkz, dky, dkx
        integer (kind=8) :: ikz, iky, ikx
        real (kind=8) :: ky_start, kz_start, kx_start

        integer (kind=8) :: nfprime, ntprime, nomd
        real (kind=8) :: dfprime, dtprime, domd
        integer (kind=8) :: ifprime, itprime, iomd
        real (kind=8) :: tprime_start, fprime_start, omd_start

        real (kind=8) :: vy_1, vy_2, vz_1, vz_2

        complex (kind=8) :: omega,integral_i 
        real (kind=8) :: A, fprime, tprime, omd, ky, kz, kx
        real (kind=8) :: vi, tol, fr, na_e, na_z, na_i, omd_kx
        real (kind=8) :: theta, Zeff, Z, mu_e, mu_z
        real (kind=8) :: dalpha
        integer (kind=8) :: nstep, steps
        integer :: gamma_unit=101, Dmixing_unit=102, out_unit=103

        real (kind=8), dimension(:), allocatable :: ky_grid, kz_grid, kx_grid
        real (kind=8), dimension(:), allocatable :: kx_lp, kx_lp_D
        real (kind=8), dimension(:), allocatable :: fprime_grid, tprime_grid
        real (kind=8), dimension(:), allocatable :: Dmixing_lp, omd_grid
        complex (kind=8), dimension(:), allocatable :: omega_lp
        integer (kind=8), dimension(:), allocatable :: ikx_lp_gamma,ikx_lp_Dmixing

        complex (kind=8) :: seed1, seed2
        real (kind=8) :: om1_re, om1_im, om2_re, om2_im
        integer :: nlp,s
        complex (kind=8) :: root, root_kz, root_ky_kz
        real(kind=8) :: kx_ref
        
        call init_file_utils
        call read_input_file
        call init_grids

        do iomd = 1, nomd
           omd = omd_grid(iomd)
           do ifprime = 1, nfprime
              fprime = fprime_grid(ifprime)
              do itprime = 1, ntprime
                 tprime = tprime_grid(itprime)
                 !print*, "tprime, fprime, omd"
                 !print*, tprime, fprime, omd
                 do iky = 1, nky
                       ky = ky_grid(iky)
                    do ikz = 1, nkz
                       kz = kz_grid(ikz)
                       do ikx = 1, nkx
                          kx = kx_grid(ikx)
!                         print*, "ky, kz"
!                         print*, ky, kz
                          omega = om1_re + om1_im*zi
                          call sa_integral(omega,integral_i,s=2)
                          print*, integral_i
                          call sa_integral(omega,integral_i,s=4)
                          print*, integral_i
                       end do
                    end do
                 end do
              end do
           end do
        end do
        call finish_file_utils

contains

subroutine sa_integral(omega,integral,s)

        implicit none

        integer, intent(in) :: s
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
        call vz_integral(x_left,f_left,s,omega)
        integral = 0.0

        do while (x_left<b)
                x_center = 0.5*(x_left+x_right)
                call vz_integral(x_center,f_center,s,omega)
                call vz_integral(x_right,f_right,s,omega)
                i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                do while (abs(i_trapezoid-i_simpson)>tol)
                        dx=fr*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                        if (dx>0.2) dx=0.2
                        x_right = x_left + dx
                        if (x_right>b) then
                                x_right = b
                                dx = x_right - x_left
                        end if
                        x_center = 0.5*(x_left+x_right)
                        call vz_integral(x_center,f_center,s,omega)
                        call vz_integral(x_right,f_right,s,omega)
                        i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                        i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                end do
                integral = integral + i_simpson
                x_left = x_right
                !print*, x_left, integral
                dx = fr*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                if (dx>0.2) dx=0.2
                x_right = x_left + dx
                if (x_right>b) then
                        x_right = b
                        dx = x_right - x_left
                end if
                f_left = f_right

        end do

  end subroutine

  subroutine vz_integral(vy,integral,s,omega)

        implicit none

        integer, intent(in) :: s
        real(kind=8), intent(in) :: vy
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
        f_left = integrand(x_left,vy,omega,s)
        integral = 0.0

        do while (x_left<b)
                x_center = 0.5*(x_left+x_right)
                f_center = integrand(x_center,vy,omega,s)
                f_right = integrand(x_right,vy,omega,s)
                i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                do while (abs(i_trapezoid-i_simpson)>tol)
                        dx=fr*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                        if (dx>0.2) dx=0.2
                        x_right = x_left + dx
                        if (x_right>b) then
                                x_right = b
                                dx = x_right - x_left
                        end if
                        x_center = 0.5*(x_left+x_right)
                        f_center = integrand(x_center,vy,omega,s)
                        f_right = integrand(x_right,vy,omega,s)
                        i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                        i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                end do
                integral = integral + i_simpson
                x_left = x_right
                !print*, x_left, integral
                dx = fr*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                if (dx>0.2) dx=0.2
                x_right = x_left + dx
                if (x_right>b) then 
                        x_right = b
                        dx = x_right - x_left
                end if
                f_left = f_right

        end do

  end subroutine

  function integrand(vz,vy,omega,s)

        real(kind=8) :: vz,vy
        complex(kind=8) :: omega,integrand,int_tmp
        complex(kind=8) :: omega_star_n, omega_star_t, omega_star_d
        complex(kind=8) :: omega_parallel
        real(kind=8) :: rho,bj,dj,fj
        integer :: s, n
        real(kind=8) :: v0,alpha,b,c
        complex(kind=8) :: dv,vzprime,d

        n = 0
        if (s==2) then
                b = (1.0-omd_kx)*ky/A+omd_kx*kx/A
                c = kz
                d = 0.5*((1.0-omd_kx)*ky/A+omd_kx*kx/A)*vy**2-omega
                v0 = -c/2./b
                dv = sqrt(c**2-4.*b*d)/b
                alpha = atan(aimag(dv)/real(dv))
                print*, 'alpha'
                print*, alpha*180./pi, alpha
                if (alpha>1e-2) then 
                        alpha =  0.0
                else
                        alpha = alpha+dalpha
                end if
                vzprime = (cos(alpha)+zi*sin(alpha))*(vz-v0)+v0
                !vzprime = vz+zi*alpha
                omega_star_n = fprime*ky
                omega_star_t = 0.5*(-3.0+vy**2+&
                        vzprime**2)*tprime*ky
                omega_star_d = (1.0-omd_kx)*(0.5*vy**2+&
                        vzprime**2)*omd*ky/A +&
                        omd_kx*(0.5*vy**2+vzprime**2)*&
                        omd*kx/A
                omega_parallel = kz*vzprime
                rho = sqrt(ky**2+kx**2)*vy
                call bjndd (n, rho, bj, dj, fj)
                int_tmp = bj**2*sqrt(1.0/pi/2.0)*vy*exp(-0.5*vy**2-&
                  0.5*vzprime**2)*&
                  (omega-omega_star_n-omega_star_t)/&
                  (omega-omega_parallel-omega_star_d)*&
                  (cos(alpha)+zi*sin(alpha))
        else if (s==4) then
                omega_star_n = fprime*ky
                omega_star_t = 0.5*(-3.0+vy**2+&
                        (vz**2-zi*2.0*vz*vi-vi**2))*tprime*ky
                omega_star_d = (1.0-omd_kx)*(0.5*vy**2+&
                        vz**2-zi*2.0*vz*vi-vi**2)*omd*ky/A +&
                        omd_kx*(0.5*vy**2+vz**2-zi*2.0*vz*vi-vi**2)*&
                        omd*kx/A
                omega_parallel = kz*(vz-zi*vi)
                rho = sqrt(ky**2+kx**2)*vy
                call bjndd (n, rho, bj, dj, fj)
                int_tmp = bj**2*sqrt(1.0/pi/2.0)*vy*exp(-0.5*vy**2-&
                  0.5*(vz**2-zi*2.0*vz*vi-vi**2))*&
                  (omega-omega_star_n-omega_star_t)/&
                  (omega-omega_parallel-omega_star_d)
        else if (s==1) then
                omega_star_n = -theta*fprime*ky
                omega_star_t = -theta*0.5*(-3.0+vy**2+&
                        (vz**2-zi*2.0*vz*vi-vi**2))*tprime*ky
                omega_star_d = -theta*((1.0-omd_kx)*(0.5*vy**2+&
                        vz**2-zi*2.0*vz*vi-vi**2)*omd*ky/A +&
                        omd_kx*(0.5*vy**2+vz**2-zi*2.0*vz*vi-vi**2)*&
                        omd*kx/A)
                omega_parallel = sqrt(theta/mu_e)*kz*(vz-zi*vi)
                rho = sqrt(theta*mu_e)*sqrt(ky**2+kx**2)*vy
                call bjndd (n, rho, bj, dj, fj)
        else if (s==3) then
                omega_star_n = 1.0/Z*fprime*ky
                omega_star_t = 1.0/Z*0.5*(-3.0+vy**2+&
                        (vz**2-zi*2.0*vz*vi-vi**2))*tprime*ky
                omega_star_d = 1.0/Z*((1.0-omd_kx)*(0.5*vy**2+&
                        vz**2-zi*2.0*vz*vi-vi**2)*omd*ky/A+&
                        omd_kx*(0.5*vy**2+vz**2-zi*2.0*vz*vi-vi**2)*&
                        omd*kx/A)
                omega_parallel = sqrt(1.0/mu_z)*kz*(vz-zi*vi)
                rho = sqrt(mu_z)/Z*sqrt(ky**2+kx**2)*vy
                call bjndd (n, rho, bj, dj, fj)
        end if
        integrand = int_tmp

  end function

  subroutine read_input_file

        use file_utils, only: input_unit_exist, input_unit

        implicit none

        integer :: in_file
        logical :: exist

        namelist / parameters / nstep,nkz,dkz,kz_start, &
                                nky,dky,ky_start, dfprime,nfprime, &
                                fprime_start,dtprime,ntprime, &
                                tprime_start, nomd,domd,omd_start, &
                                vy_1,vy_2,vz_1,vz_2,om1_re,om1_im, &
                                om2_re,om2_im,A,vi,tol,fr, &
                                theta,Zeff,Z,mu_e,mu_z,na_e,na_z,na_i, &
                                nkx, dkx, kx_start, omd_kx, dalpha
        
        nstep = 10

        nkz = 1
        dkz = 0.1
        kz_start = 0.0

        nky = 1
        dky = 0.1
        ky_start = 0.0

        dfprime = 0.5
        nfprime = 1.0
        fprime_start = 0.0

        ntprime = 1
        dtprime = 1.0
        tprime_start = 0.0

        nomd = 1
        domd = 1.0
        omd_start = 0.0

        vy_1 = 0.0
        vy_2 = 0.0
        vz_1 = 8.0
        vz_2 = 8.0

        om1_re = 1.0
        om1_im = 0.1
        om2_re = 1.0
        om2_im = 0.1

        A = 3.0
        vi  = 1.0
        tol = 1.0E-05
        fr = 0.5

        theta = 1.0
        Zeff = 1.65
        Z = 5
        mu_e = 5.45D-04
        mu_z = 10.8
        na_e = 0.0
        na_z = 1.0
        na_i = 1.0
        
        nkx = 1
        dkx = 0.5
        kx_start = 0.0
        omd_kx = 0.0
        
        dalpha = 0.1
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
                kx_grid(ikx) = -ikx*dkx+kx_start
        end do
        allocate(fprime_grid(nfprime), tprime_grid(ntprime), omd_grid(nomd))

        do ifprime = 1, nfprime
                fprime_grid(ifprime) = ifprime*dfprime+fprime_start
        end do 
        do itprime = 1, ntprime
                tprime_grid(itprime) = itprime*dtprime+tprime_start
        end do
        do iomd = 1, nomd
                omd_grid(iomd) = iomd*domd+omd_start
        end do

        allocate(kx_lp(4), omega_lp(4), kx_lp_D(4), Dmixing_lp(4))
        
        allocate(ikx_lp_gamma(1),ikx_lp_Dmixing(1))

end subroutine

end program
