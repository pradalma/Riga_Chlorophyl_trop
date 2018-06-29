# Riga_Chlorophyl_trop
change param in tropics
-------------------------------------------
-------------------------------------------
Riga_Chlorophyl_trop/CM2Mc/src/mom4p1/ocean_param/sources/ocean_shortwave_gfdl.F90
-------------------------------------------
-------------------------------------------

New code change

  do i=isc,iec
            if (Grd%tmask(i,j,1) > 0.0) then

                k=1
                wrk1(i,j,k) =                                  &
                     0.5*exp(-Thickness%dzt(i,j,k)*(0.225+0.037*wrk3(i,j,k)**0.629))
                if (abs(Grd%yt(i,j))<30.0) then
                wrk2(i,j,k) =                                  &
                     0.5*exp(-Thickness%dzt(i,j,k)*(0.0232+0.0513*wrk3(i,j,k)**0.668+0.71*wrk5(i,j,k)**1.13))
                     !GKIM formulation
                else
                wrk2(i,j,k) =                                  &
                     0.5*exp(-Thickness%dzt(i,j,k)*(0.0232+0.071*wrk3(i,j,k)**0.674)   !Manizza outside Tropics
                endif

                sw_frac_zw(i,j,k) = f_vis(i,j)*(wrk1(i,j,k) + wrk2(i,j,k))
                sw_frac_zt(i,j,k) = 0.5*(f_vis(i,j) + sw_frac_zw(i,j,1))

                do k=2,nk-1
                   if (wrk4(i,j,k) < zmax_pen) then
                       wrk1(i,j,k) = wrk1(i,j,k-1) &
                            *exp(-Thickness%dzt(i,j,k)*(0.2250+0.037*wrk3(i,j,k)**0.629))
                if (abs(Grd%yt(i,j))<30.0) then
                       wrk2(i,j,k) = wrk2(i,j,k-1) &
                            *exp(-Thickness%dzt(i,j,k)*(0.0232+0.0513*wrk3(i,j,k)**0.668+0.71*wrk5(i,j,k)**1.13))
                            !GKIM formulation
                else
                wrk2(i,j,k) =                                  &
                     0.5*exp(-Thickness%dzt(i,j,k)*(0.0232+0.071*wrk3(i,j,k)**0.674)   !Manizza outside Tropics
                endif
                       sw_frac_zw(i,j,k)  = f_vis(i,j)*(wrk1(i,j,k) + wrk2(i,j,k))
                       sw_frac_zt(i,j,k)  = &
                            0.5*(sw_frac_zw(i,j,k-1) + sw_frac_zw(i,j,k))
                   endif
                enddo
                
                
    ! compute shortwave penetration using Manizza etal optics
!      do j=jsc,jec
!         do i=isc,iec
!            if (Grd%tmask(i,j,1) > 0.0) then

!                k=1
!                wrk1(i,j,k) =                                  &
!                     0.5*exp(-Thickness%dzt(i,j,k)*(0.225+0.037*wrk3(i,j,k)**0.629))
!                wrk2(i,j,k) =                                  &
!                     0.5*exp(-Thickness%dzt(i,j,k)*(0.0232+0.074*wrk3(i,j,k)**0.674))
!                sw_frac_zw(i,j,k) = f_vis(i,j)*(wrk1(i,j,k) + wrk2(i,j,k))
!                sw_frac_zt(i,j,k) = 0.5*(f_vis(i,j) + sw_frac_zw(i,j,1))
!
!                do k=2,nk-1
!                   if (wrk4(i,j,k) < zmax_pen) then
!                       wrk1(i,j,k) = wrk1(i,j,k-1) &
!                            *exp(-Thickness%dzt(i,j,k)*(0.2250+0.037*wrk3(i,j,k)**0.629))
!                       wrk2(i,j,k) = wrk2(i,j,k-1) &
!                            *exp(-Thickness%dzt(i,j,k)*(0.0232+0.074*wrk3(i,j,k)**0.674))
!                       sw_frac_zw(i,j,k)  = f_vis(i,j)*(wrk1(i,j,k) + wrk2(i,j,k))
!                       sw_frac_zt(i,j,k)  = &
!                            0.5*(sw_frac_zw(i,j,k-1) + sw_frac_zw(i,j,k))
!                   endif
!                enddo
!
!            endif
!         enddo
!      enddo
!
!  endif  ! endif for optics_manizza
                
 -------------------------------------------
 -------------------------------------------
 in the exec_hpcs
 
 edited  compile_CM2M_compile.csh to point to Riga_Chlorophyl_trop
 
 ------------------------------------------
 
make clean
--------------------------------------------
rm .x
-----------------------------------------------
rm Makefile*
-----------------------------------------------
execute csh
./compile_CM2M_compile.csh
