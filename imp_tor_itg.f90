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

        complex (kind=8) :: omega 
        real (kind=8) :: A, fprime, tprime, omd, ky, kz, kx
        real (kind=8) :: vi, tol, fr, na_e, na_z, na_i, omd_kx
        real (kind=8) :: theta, Zeff, Z, mu_e, mu_z
        integer (kind=8) :: nstep, steps
        integer :: out_unit=101

        real (kind=8), dimension(:), allocatable :: ky_grid, kz_grid, kx_grid
        real (kind=8), dimension(:), allocatable :: fprime_grid, tprime_grid
        real (kind=8), dimension(:), allocatable :: omd_grid

        complex (kind=8) :: root, root_kz
        complex (kind=8) :: seed1, seed2
        real (kind=8) :: om1_re, om1_im, om2_re, om2_im
        
        call init_file_utils
        call read_input_file
        call init_grids

        call open_output_file(out_unit,'.dat')
        write (out_unit,'(8a12)') "tprime","fprime","omd","kx",&
                                "ky","kz","omega","gamma"
        do iomd = 1, nomd
           omd = omd_grid(iomd)
           do ifprime = 1, nfprime
              fprime = fprime_grid(ifprime)
              do itprime = 1, ntprime
                 tprime = tprime_grid(itprime)
                 print*, "tprime, fprime, omd"
                 print*, tprime, fprime, omd
                 do ikx = 1, nkx
                    do iky = 1, nky
                       do ikz = 1, nkz
                          ky = ky_grid(iky)
                          kz = kz_grid(ikz)
                          kx = kx_grid(ikx)
                          print*, "kx, ky, kz"
                          print*, kx,ky,kz
                          if (ifprime ==1 .and. itprime==1 .and. iomd==1 &
                            .and. iky==1 .and. ikz==1 .and. ikx==1) then
                            seed1 = om1_re + om1_im*zi
                            seed2 = om2_re + om2_im*zi
                            call rootfinder(seed1,seed2,root,steps)
                            if (steps==nstep) then
                                print*, 'No root is found.'
                            else 
                                print*, "Root is found:"
                                print*, root    
                                write (out_unit, '(8e12.4)') tprime,fprime,&
                                                        omd,kx,ky,kz,root
                            end if
                            root_kz = root
                          else if (ifprime==1 .and. itprime==1 .and. iomd==1&
                                .and. ikz==1) then
                            seed1 = root_kz
                            seed2 = root_kz*0.9
                            call rootfinder(seed1,seed2,root,steps)
                            if (steps==nstep) then
                                print*, 'No root is found.'
                            else 
                                print*, "Root is found:"
                                print*, root    
                                write (out_unit, '(8e12.4)') tprime,fprime,&
                                                        omd,kx,ky,kz,root
                            end if
                            root_kz = root
                          else
                            seed1 = root 
                            seed2 = root*0.9
                            call rootfinder(seed1,seed2,root,steps)
                            if (steps==nstep) then
                                print*, 'No root is found.'
                            else 
                                print*, "Root is found:"
                                print*, root    
                                write (out_unit, '(8e12.4)') tprime,fprime,&
                                                        omd,kx,ky,kz,root
                            end if
                          end if
                       end do
                    end do
                 end do
              end do
           end do
        end do
        call close_output_file (out_unit)
        call finish_file_utils

contains

 subroutine rootfinder(sd1,sd2,rt,istep)

        implicit none

        complex (kind=8), intent(in) :: sd1,sd2
        complex (kind=8), intent(out) :: rt
        integer (kind=8), intent(out) :: istep
        complex (kind=8) :: x_2,f_2,x_1,f_1
        complex (kind=8) :: x_tmp, f_tmp

        x_1 = sd1
        call dispersion_relation(x_1,f_1)
        x_2 = sd2
        call dispersion_relation(x_2,f_2)
        istep = 0
        do while (abs(f_1)>1.0E-04 .and. abs(f_2)>1.0E-04 .and. istep<nstep) 
                x_tmp = x_1 - f_1 * ((x_1-x_2)/(f_1-f_2))
                call dispersion_relation(x_tmp,f_tmp)
                f_2 = f_1
                f_1 = f_tmp
                x_2 = x_1
                x_1 = x_tmp
                istep = istep +1
                print*, istep
                print*, x_tmp, f_tmp
        end do
        rt = x_tmp

  end subroutine

subroutine dispersion_relation(omega, rhs)
        
        implicit none

        complex(kind=8), intent(in) :: omega
        complex (kind=8), intent(out) :: rhs
        complex (kind=8) :: integral_e, integral_i, integral_z
        integer :: s
        
        call sa_integral(omega,integral_e,s=1)
        call sa_integral(omega,integral_i,s=2)
        call sa_integral(omega,integral_z,s=3)
        rhs = 1.0 + theta*Zeff - na_e*integral_e - &
              na_i*theta*(Z-Zeff)/(z-1.0)*integral_i - &
              na_z*theta*(Zeff-1.0)/z/(z-1.0)*integral_z


end subroutine

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

        n = 0
        if (s==2) then
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
                omega_star_n = theta/Z*fprime*ky
                omega_star_t = theta/Z*0.5*(-3.0+vy**2+&
                        (vz**2-zi*2.0*vz*vi-vi**2))*tprime*ky
                omega_star_d = theta/Z*((1.0-omd_kx)*(0.5*vy**2+&
                        vz**2-zi*2.0*vz*vi-vi**2)*omd*ky/A+&
                        omd_kx*(0.5*vy**2+vz**2-zi*2.0*vz*vi-vi**2)*&
                        omd*kx/A)
                omega_parallel = sqrt(1.0/mu_z)*kz*(vz-zi*vi)
                rho = sqrt(mu_z)/Z*sqrt(ky**2+kx**2)*vy
                call bjndd (n, rho, bj, dj, fj)
        end if
                
        
        int_tmp = bj**2*sqrt(1.0/pi/2.0)*vy*exp(-0.5*vy**2-&
                  0.5*(vz**2-zi*2.0*vz*vi-vi**2))*&
                  (omega-omega_star_n-omega_star_t)/&
                  (omega-omega_parallel-omega_star_d)

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
                                nkx, dkx, kx_start, omd_kx
        
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

end subroutine

end program
