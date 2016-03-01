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
        integer :: gamma_unit=101, Dmixing_unit=102, out_unit=103

        real (kind=8), dimension(:), allocatable :: ky_grid, kz_grid, kx_grid
        real (kind=8), dimension(:), allocatable :: kx_lp, kx_lp_D
        real (kind=8), dimension(:), allocatable :: fprime_grid, tprime_grid
        real (kind=8), dimension(:), allocatable :: Dmixing_lp, omd_grid
        complex (kind=8), dimension(:), allocatable :: omega_lp
        integer (kind=8), dimension(:), allocatable :: ikx_lp_gamma,ikx_lp_Dmixing

        complex (kind=8) :: seed1, seed2
        real (kind=8) :: om1_re, om1_im, om2_re, om2_im
        integer :: nlp
        complex (kind=8) :: root, root_kz, root_ky_kz, root_kx0
        real(kind=8) :: kx_ref
        
        call init_file_utils
        call read_input_file
        call init_grids

        call open_output_file(gamma_unit,'.datgam')
        write (gamma_unit,'(8a12)') "tprime","fprime","omd","kx",&
                                "ky","kz","omega","gamma"
        call open_output_file(Dmixing_unit,'.datDmixing')
        write (Dmixing_unit,'(9a12)') "tprime","fprime","omd","kx",&
                                "ky","kz","omega","gamma","Dmixing"
        call open_output_file(out_unit,'.dat')
        write (out_unit,'(8a12)') "tprime","fprime","omd","kx",&
                                "ky","kz","omega","gamma"

        root = om1_re + om1_im*zi
        root_kz = om1_re + om1_im*zi
        root_ky_kz = om1_re + om1_im*zi
        root_kx0 = om1_re + om1_im*zi

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
                       print*, "ky, kz"
                       print*, ky, kz
                       kx_lp = 0.0; omega_lp = 0.0; Dmixing_lp = 0.0

                       do nlp = 1,4
                          call kx_scan(nlp,kx_lp(nlp),omega_lp(nlp),kx_lp_D(nlp),Dmixing_lp(nlp))
                       end do

                       ikx_lp_gamma = maxloc(aimag(omega_lp))
                       ikx_lp_Dmixing = maxloc(Dmixing_lp)

                       if (.not.maxval(aimag(omega_lp))==-9999.) write (gamma_unit,&
                           '(8e12.4)') tprime,fprime,omd, &
                           kx_lp(ikx_lp_gamma(1)),ky,kz,omega_lp(ikx_lp_gamma(1))
                       if (.not.maxval(Dmixing_lp)==-9999.) write (Dmixing_unit, &
                           '(9e12.4)') tprime,fprime,omd, &
                           kx_lp_D(ikx_lp_Dmixing(1)),ky,kz,omega_lp(ikx_lp_Dmixing(1)),&
                           Dmixing_lp(ikx_lp_Dmixing(1))
                       if (.not.maxval(aimag(omega_lp))==-9999.) then  
                           kx_ref = kx_lp(ikx_lp_gamma(1))
                           root_ky_kz = omega_lp(ikx_lp_gamma(1))
                           if (abs(kx_ref)<0.1) root_kx0 = omega_lp(ikx_lp_gamma(1))
                           if (ikz==1) root_kz = omega_lp(ikx_lp_gamma(1))
                       end if
                    end do
                    write(out_unit,*)
                 end do
                 write(out_unit,*)
              end do
           end do
        end do
        call close_output_file (gamma_unit)
        call close_output_file (Dmixing_unit)
        call close_output_file (out_unit)
        call finish_file_utils

contains

 subroutine kx_scan (iloop,kx_maxgam,omega,kx_maxDmix,Dmixing)

        implicit none

        integer, intent(in) :: iloop
        real (kind=8), intent(out) :: kx_maxgam, kx_maxDmix, Dmixing
        complex (kind=8), intent(out) :: omega
        complex (kind=8), dimension(:), allocatable :: omega_kx_scan
        integer (kind=8), dimension(:), allocatable :: ikx_maxgamma,ikx_maxDmixing
        real (kind=8), dimension(:), allocatable :: Dmixing_kx_scan
        real (kind=8) :: kx_init


        allocate(omega_kx_scan(nkx),Dmixing_kx_scan(nkx),ikx_maxgamma(1),ikx_maxDmixing(1))
        omega_kx_scan = -9999.0*zi
        Dmixing_kx_scan = -9999.0

        if (iloop==1) then
            kx_init = kx_ref
            if (iky==1.and.ikz==1) kx_init = kx_start
            do ikx = 1, nkx
                kx_grid(ikx) = (ikx-1)*dkx+kx_init
            end do
        else if (iloop==2) then
            kx_init = kx_ref
            if (iky==1.and.ikz==1) kx_init = kx_start
            do ikx = 1, nkx
                kx_grid(ikx) = -ikx*dkx+kx_init
            end do
        else if (iloop==3) then
            kx_init = 0.0
            do ikx = 1, nkx
                kx_grid(ikx) = (ikx-1)*dkx+kx_init
            end do
        else if (iloop==4) then
            kx_init = 0.0
            do ikx = 1, nkx
                kx_grid(ikx) = -ikx*dkx+kx_init
            end do

       end if 

       do ikx = 1, nkx
          kx = kx_grid(ikx)
          ! residue calculation doesn't allow kx==0.
          if (abs(kx)<0.01) kx=0.01
          print*, "kx, iloop"
          print*, kx,iloop
          if (ifprime ==1 .and. itprime==1 .and. iomd==1 &
              .and. iky==1 .and. ikz==1 .and. ikx==1) then
            seed1 = om1_re + om1_im*zi
            seed2 = om2_re + om2_im*zi
          else if (ikz==1 .and. ikx==1) then
            seed1 = root_kz
            seed2 = root_kz*0.9
          else if (ikx==1 .and. abs(kx)>0.1) then
            seed1 = root_ky_kz
            seed2 = root_ky_kz*0.9
          else if (abs(kx)<0.1 ) then
            seed1 = root_kx0
            seed2 = root_kx0*0.9
          else
            seed1 = root
            seed2 = root*0.9
          end if

          call rootfinder(seed1,seed2,root,steps)

          if (steps==nstep) then
              print*, 'No root is found.'
              root = seed1
          else if (abs(root)<5.0 .and. abs(aimag(root))>1.0E-06 ) then
              print*, "Root is found:"
              print*, root    
              write (out_unit, '(8e12.4)') tprime,fprime,&
                                        omd,kx,ky,kz,root
              omega_kx_scan(ikx) = root
              if (abs(kx)<0.1) root_kx0 = root
          else 
              print*, 'No root is found.'
              root = seed1
          end if
       end do

       ikx_maxgamma = maxloc(aimag(omega_kx_scan))
       Dmixing_kx_scan = aimag(omega_kx_scan)/(kx_grid**2+ky**2)
       ikx_maxDmixing = maxloc(Dmixing_kx_scan)

       kx_maxgam = kx_grid(ikx_maxgamma(1))
       omega = omega_kx_scan(ikx_maxgamma(1))
       kx_maxDmix = kx_grid(ikx_maxDmixing(1))
       Dmixing = Dmixing_kx_scan(ikx_maxDmixing(1))

       deallocate(omega_kx_scan,Dmixing_kx_scan,ikx_maxgamma,ikx_maxDmixing)

  end subroutine

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
                if (abs(f_1-f_2)<1.0E-05) then 
                        x_tmp  = x_1 - 1.0E-03
                else
                        x_tmp = x_1 - f_1 * ((x_1-x_2)/(f_1-f_2))
                end if
                if (aimag(x_tmp)<-1.0) then
                        istep = nstep
                        exit
                end if
                call dispersion_relation(x_tmp,f_tmp)
                f_2 = f_1
                f_1 = f_tmp
                x_2 = x_1
                x_1 = x_tmp
                istep = istep +1
        end do
        if ((.not.isnan(real(f_2))) .and. (.not.isnan(real(f_1))))  rt = x_tmp

  end subroutine

subroutine dispersion_relation(omega, rhs)
        
        implicit none

        complex(kind=8), intent(in) :: omega
        complex (kind=8), intent(out) :: rhs
        complex (kind=8) :: integral_e, integral_i, integral_z
        integer :: s
        
        call vy_integral(omega,integral_e,s=1)
        call vy_integral(omega,integral_i,s=2)
        call vy_integral(omega,integral_z,s=3)
        rhs = 1.0 + theta*Zeff - na_e*integral_e - &
              na_i*theta*(Z-Zeff)/(z-1.0)*integral_i - &
              na_z*theta*(Zeff-1.0)/z/(z-1.0)*integral_z


end subroutine

subroutine vy_integral(omega,integral,s)

        implicit none

        integer, intent(in) :: s
        complex(kind=8), intent(in) :: omega
        complex (kind=8), intent(out) :: integral
        real(kind=8) :: x_left, x_center, x_right, dx
        complex (kind=8) :: f_left, f_center, f_right
        complex (kind=8) :: i_trapezoid, i_simpson
        real(kind=8) :: int_a,int_b

        int_a = vy_1
        int_b = vy_2
        dx = 0.1
        x_left = int_a
        x_right = x_left + dx
        call vz_integral(omega,x_left,f_left,s)
        integral = 0.0

        do while (x_left<int_b)
                x_center = 0.5*(x_left+x_right)
                call vz_integral(omega,x_center,f_center,s)
                call vz_integral(omega,x_right,f_right,s)
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
                        call vz_integral(omega,x_center,f_center,s)
                        call vz_integral(omega,x_right,f_right,s)
                        i_trapezoid = 0.5*(f_left+2.0*f_center+f_right)*0.5*dx
                        i_simpson = 1.0/6.0*dx*(f_left+4.0*f_center+f_right)
                end do
                integral = integral + i_simpson
                x_left = x_right
                dx = fr*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                if (dx>0.2) dx=0.2
                x_right = x_left + dx
                dx = fr*dx*(abs(i_trapezoid-i_simpson)/tol)**(-1.0/3.0)
                if (dx>0.2) dx=0.2
                x_right = x_left + dx
                if (x_right>int_b) then
                        x_right = int_b
                        dx = x_right - x_left
                end if
                f_left = f_right

        end do

  end subroutine

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
        if (s==2) then
            rho = sqrt(ky**2+kx**2)*vy
            call bjndd (n, rho, bj, dj, fj)
            if (isnan(bj)) bj = 1.0
            kv = (1.0-omd_kx)*ky + omd_kx*kx
            coeff_a = kv/A
            coeff_b = kz
            coeff_c = 0.5*kv/A*vy**2-omega
            coeff_d = 0.5*tprime*ky
            coeff_e = (-3./2.+vy**2/2.)*tprime*ky+fprime*ky-omega
            v0 = -coeff_b/2./coeff_a
            dv = sqrt(coeff_b**2-4.*coeff_a*coeff_c)/coeff_a
            dvcr = coeff_b**2-4.0*coeff_a*(0.5*kv/A*vy**2-real(omega))
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
        else if (s==1) then
            rho = sqrt(theta*mu_e)*sqrt(ky**2+kx**2)*vy
            call bjndd (n, rho, bj, dj, fj)
            if (isnan(bj)) bj = 1.0
            kv = -theta*((1.0-omd_kx)*ky + omd_kx*kx)
            coeff_a = kv/A
            coeff_b = sqrt(theta/mu_e)*kz
            coeff_c = 0.5*kv/A*vy**2-omega
            coeff_d = -theta*0.5*tprime*ky
            coeff_e = -theta*((-3./2.+vy**2/2.)*tprime*ky+fprime*ky)-omega
            v0 = -coeff_b/2./coeff_a
            dv = sqrt(coeff_b**2-4.*coeff_a*coeff_c)/coeff_a
            dvcr = coeff_b**2-4.0*coeff_a*(0.5*kv/A*vy**2-real(omega))
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
        else if (s==3) then
            rho = sqrt(mu_z)/Z*sqrt(ky**2+kx**2)*vy
            call bjndd (n, rho, bj, dj, fj)
            if (isnan(bj)) bj = 1.0
            kv = 1./Z*((1.0-omd_kx)*ky + omd_kx*kx)
            coeff_a = kv/A
            coeff_b = sqrt(1./mu_z)*kz
            coeff_c = 0.5*kv/A*vy**2-omega
            coeff_d = 1./Z*0.5*tprime*ky
            coeff_e = 1./Z*((-3./2.+vy**2/2.)*tprime*ky+fprime*ky)-omega
            v0 = -coeff_b/2./coeff_a
            dv = sqrt(coeff_b**2-4.*coeff_a*coeff_c)/coeff_a
            dvcr = coeff_b**2-4.0*coeff_a*(0.5*kv/A*vy**2-real(omega))
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
        end if

        int_a = vz_1
        int_b = vz_2
        dx = 0.1
        x_left = int_a
        x_right = x_left + dx
        call integrand(x_left,vy,omega,v1,f_left,dvcr,coeff_a,s)
        integral = 0.0

        do while (x_left<int_b)
                x_center = 0.5*(x_left+x_right)
                call integrand(x_center,vy,omega,v1,f_center,dvcr,coeff_a,s)
                call integrand(x_right,vy,omega,v1,f_right,dvcr,coeff_a,s)
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
                        call integrand(x_center,vy,omega,v1,f_center,dvcr,coeff_a,s)
                        call integrand(x_right,vy,omega,v1,f_right,dvcr,coeff_a,s)
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

        !if (s==2) then
        if (dvcr<0.) then
           integral = bj**2*integral
        else if (coeff_a>0.) then
           if (aimag(v1)>0.1) then
              integral = bj**2*integral
           else if (aimag(v1)>-0.1) then
              integral = bj**2*(integral - res_v2)
           else
              integral = bj**2*(integral + res_v1 - res_v2)
           end if
         else
           if (aimag(v1)<-0.1) then
              integral = bj**2*integral
           else if (aimag(v1)<0.1) then
              integral = bj**2*(integral - res_v1)
           else
              integral = bj**2*(integral - res_v1 + res_v2)
           end if
        end if
        !end if

  end subroutine

  subroutine integrand(vz,vy,omega,v1,int_tmp,cr,coeff_a,s)

        real(kind=8),intent(in) :: vz,vy,cr,coeff_a
        complex(kind=8),intent(in) :: omega,v1
        integer, intent(in) :: s
        complex(kind=8),intent(out) :: int_tmp
        complex(kind=8) :: omega_star_n, omega_star_t, omega_star_d
        complex(kind=8) :: omega_parallel
        complex(kind=8) :: vzprime

        if (s==2) then
           if (cr <0.) then
                vzprime = vz
           else
             if (coeff_a>0) then   
                if (aimag(v1)>0.1) then
                        vzprime = vz
                else if (aimag(v1)>-0.1) then
                        vzprime = vz-zi*0.2
                else
                        vzprime = vz
                end if
             else 
                if (aimag(v1)<-0.1) then
                        vzprime = vz
                else if (aimag(v1)<0.1) then
                        vzprime = vz-zi*0.2
                else
                        vzprime = vz
                end if
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

        else if (s==1) then
           if (cr <0.) then
                vzprime = vz
           else
             if (coeff_a>0) then   
                if (aimag(v1)>0.1) then
                        vzprime = vz
                else if (aimag(v1)>-0.1) then
                        vzprime = vz-zi*0.2
                else
                        vzprime = vz
                end if
             else 
                if (aimag(v1)<-0.1) then
                        vzprime = vz
                else if (aimag(v1)<0.1) then
                        vzprime = vz-zi*0.2
                else
                        vzprime = vz
                end if
             end if
           end if
                omega_star_n = -theta*fprime*ky
                omega_star_t = -theta*0.5*(-3.0+vy**2+&
                        (vz**2-zi*2.0*vz*vi-vi**2))*tprime*ky
                omega_star_d = -theta*((1.0-omd_kx)*(0.5*vy**2+&
                        vz**2-zi*2.0*vz*vi-vi**2)*omd*ky/A +&
                        omd_kx*(0.5*vy**2+vz**2-zi*2.0*vz*vi-vi**2)*&
                        omd*kx/A)
                omega_parallel = sqrt(theta/mu_e)*kz*(vz-zi*vi)

        else if (s==3) then
           if (cr <0.) then
                vzprime = vz
           else
             if (coeff_a>0) then   
                if (aimag(v1)>0.1) then
                        vzprime = vz
                else if (aimag(v1)>-0.1) then
                        vzprime = vz-zi*0.2
                else
                        vzprime = vz
                end if
             else 
                if (aimag(v1)<-0.1) then
                        vzprime = vz
                else if (aimag(v1)<0.1) then
                        vzprime = vz-zi*0.2
                else
                        vzprime = vz
                end if
             end if
           end if
                omega_star_n = 1.0/Z*fprime*ky
                omega_star_t = 1.0/Z*0.5*(-3.0+vy**2+&
                        (vz**2-zi*2.0*vz*vi-vi**2))*tprime*ky
                omega_star_d = 1.0/Z*((1.0-omd_kx)*(0.5*vy**2+&
                        vz**2-zi*2.0*vz*vi-vi**2)*omd*ky/A+&
                        omd_kx*(0.5*vy**2+vz**2-zi*2.0*vz*vi-vi**2)*&
                        omd*kx/A)
                omega_parallel = sqrt(1.0/mu_z)*kz*(vz-zi*vi)
        end if
        int_tmp = sqrt(1.0/pi/2.0)*vy*exp(-0.5*vy**2-&
                  0.5*vzprime**2)*&
                  (omega-omega_star_n-omega_star_t)/&
                  (omega-omega_parallel-omega_star_d)
 end subroutine

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
