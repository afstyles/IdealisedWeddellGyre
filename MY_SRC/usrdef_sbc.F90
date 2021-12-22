MODULE usrdef_sbc
   !!======================================================================
   !!                       ***  MODULE  usrdef_sbc  ***
   !! 
   !!                      ===  CANAL configuration  ===
   !!
   !! User defined :   surface forcing of a user configuration
   !!======================================================================
   !! History :  4.0   ! 2017-11  (J.Chanut)  user defined interface
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_sbc    : user defined surface bounday conditions in OVERFLOW case
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE phycst          ! physical constants
   USE usrdef_nam      !, ONLY : rn_u10, rn_uofac, rn_windszy  &
                       !&  , ln_tau_acc, rn_tau_acc, rn_tau_wg, rn_tau_ext &
                       !&  , rn_0yratio, rn_dy, rn_sponge_ly
   USE usrdef_zgr
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_fortran     ! Fortran utilities (allows no signed zero when 'key_nosignedzero' defined) 
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usrdef_sbc_oce      ! routine called in sbcmod module
   PUBLIC   usrdef_sbc_ice_tau  ! routine called by icestp.F90 for ice dynamics
   PUBLIC   usrdef_sbc_ice_flx  ! routine called by icestp.F90 for ice thermo

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_sbc.F90 10074 2018-08-28 16:15:49Z nicolasmartin $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usrdef_sbc_oce( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE usr_def_sbc  ***
      !!              
      !! ** Purpose :   provide at each time-step the surface boundary
      !!              condition, i.e. the momentum, heat and freshwater fluxes.
      !!
      !! ** Method  :   all 0 fields, for CANAL case
      !!                CAUTION : never mask the surface stress field !
      !!
      !! ** Action  : - set to ZERO all the ocean surface boundary condition, i.e.   
      !!                   utau, vtau, taum, wndm, qns, qsr, emp, sfx
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      
      INTEGER  ::   ji, jj               ! dummy loop indices
      REAL(wp) :: zrhoair = 1.22     ! approximate air density [Kg/m3]
      REAL(wp) :: zcd = 1.13e-3      ! approximate drag coefficient
      REAL(wp) :: zrhocd             ! Rho * Cd
      REAL(wp), DIMENSION(jpi,jpj) :: zwndrel   ! relative wind
      REAL(wp) :: zminphi, zmaxphi, freeze_integral, melt_integral, sfx_integral            ! Minimum/maximum value of phit
      REAL(wp), DIMENSION(jpi,jpj) :: zH, z2d, shelf_frac, ice_mask, sfx_freeze, sfx_melt !Basin info 
      INTEGER, DIMENSION(jpi,jpj) :: pk_bot
      !!---------------------------------------------------------------------
      !
      zrhocd = zrhoair * zcd
      
      IF( kt == nit000 ) THEN
         !
         IF(lwp) WRITE(numout,*)' usr_sbc : EW_CANAL case: surface forcing'
         IF(lwp) WRITE(numout,*)' ~~~~~~~~~~~   vtau = taum = wndm = qns = qsr = emp = sfx = 0'
         !
         utau(:,:) = 0.
         IF(lwp) WRITE(numout,*)' utau(:,:) = 0.'
         
         zminphi = MINVAL(gphit)
         CALL mpp_min( 'usrdef_sbc', zminphi )
         !
         IF(lwp) WRITE(numout,*)' MPP_MIN AND MPP_MAX'
         !
         WHERE( gphit < -rn_tau_ext )
            utau(:,:) = - rn_tau_wg * sin( 0.5*rpi*(gphit - zminphi)/(rn_tau_ext + zminphi) )**2
         END WHERE
         !
         IF(ln_tau_acc) THEN
            WHERE( (gphit >= -rn_tau_ext) )
               utau(:,:) = - rn_tau_wg + (rn_tau_acc + rn_tau_wg)*sin( rpi*(gphit + rn_tau_ext)/(rn_sponge_ly + 2*rn_tau_ext))**2
            END WHERE
         ELSE
            WHERE( (gphit >= -rn_tau_ext).AND.(gphit<=0) )
               utau(:,:) = - rn_tau_wg*sin( 0.5*rpi*gphit/(rn_tau_ext))**2
            END WHERE
            !
            WHERE( gphit > 0 )
               utau(:,:) = 0.0
            END WHERE
            !
         END IF
         !
         vtau(:,:) = 0._wp
         taum(:,:) = 0._wp
         wndm(:,:) = 0._wp
         !
         CALL ideal_WG_bath(glamt, gphit, rn_h_flat, rn_h_ork, rn_domszz, rn_x1, rn_x2, rn_x3, &
              &             rn_x_ork, rn_y2, rn_d1, rn_d2, rn_d3, -zminphi, ln_orkney, ln_fh,  &
              &             rn_r0, rn_r1, ff_t, rn_H_mbump, ln_mbump, rn_d_mbump,   zH, z2d,   &
              &             shelf_frac )

         !WHERE( ((gphit <= 0.).AND.(glamt <= 0)).OR.((gphit >= rn_sponge_ly).AND.(glamt <= 0)) )
         !   z2d(:,:) = 0.
         !END WHERE
         
         !CALL lbc_lnk( 'usrdef_sbc', z2d, 'T', 1. )           ! set surrounding land to zero (here jperio=0 ==>> closed)

         !pk_bot = INT( z2d(:,:) )
              

         emp (:,:) = 0._wp
         sfx (:,:) = 0._wp
         qns (:,:) = 0._wp
         qsr (:,:) = 0._wp

  
        ice_mask(:,:) = 0.
        
        WHERE( (gphit <= -rn_y_ice).AND.(gphit >= -(rn_y_ice + rn_d_ice)) )
           !
           ice_mask(:,:) = SIN(0.5*rpi*(gphit + rn_y_ice)/rn_d_ice )**2
           !
        END WHERE 

        WHERE(gphit < -(rn_y_ice + rn_d_ice))
           !
           ice_mask(:,:) = 1.
           !
        END WHERE
         
        sfx_freeze(:,:) = 0.

        WHERE(shelf_frac < rn_H_freeze / rn_H_flat )
           !
           sfx_freeze(:,:) = - rn_sfx_max*ice_mask*COS(rpi*shelf_frac/(2*rn_H_freeze/rn_H_flat))**2
           !
        END WHERE

        !WHERE( pk_bot == 0 )
        !   !
        !   sfx_freeze(:,:) = 0.
        !   !
        !END WHERE

        !Sum all the salt release due to ice freezing
        !freeze_integral = SUM(sfx_freeze * tmask(:,:,1) )
        freeze_integral = glob_sum( 'closea', sfx_freeze * tmask(:,:,1) )
        IF(lwp) WRITE(numout,*) "freeze_integral = ", freeze_integral        

        !Define the normalised form for the melting zone in the domain
        sfx_melt(:,:) = ice_mask

 
        WHERE((shelf_frac >= rn_H_freeze / rn_H_flat).AND.(shelf_frac < 1.))
            !
            sfx_melt(:,:) = ice_mask*SIN(rpi*(shelf_frac-rn_H_freeze/rn_H_flat)/(2*(1-rn_H_freeze/rn_H_flat)))**2    
            !
        END WHERE

        WHERE(shelf_frac < rn_H_freeze / rn_H_flat)
            !
            sfx_melt(:,:) = 0.
            !
        END WHERE

        !WHERE( pk_bot == 0 )
        !   !
        !   sfx_melt(:,:) = 0.
        !   !
        !END WHERE 

        !Calculate the normalised integral of the melt contribution
        !melt_integral = SUM(sfx_melt * tmask(:,:,1))
        melt_integral = glob_sum('usrdef_sbc', sfx_melt * tmask(:,:,1) )
        IF(lwp) WRITE(numout,*) "normalised melt_integral = ", melt_integral


        !Scale the melt contribution so that freeze_integral = melt_integral
        sfx_melt(:,:) = ABS(freeze_integral/melt_integral)*sfx_melt

        sfx(:,:) = sfx_freeze + sfx_melt 

        !Test the integral of sfx is correct
        sfx_integral = SUM(sfx * tmask(:,:,1))
        CALL mpp_sum('usrdef_sbc', sfx_integral)
        IF(lwp) WRITE(numout,*) "sfx_integral = ", sfx_integral
         
         !qns = - shelf_frac * rn_q_ice * ice_mask
         !         
      ENDIF

      IF( rn_uofac /= 0. ) THEN
         
         WHERE( ABS(gphit) <= rn_windszy/2. )
            zwndrel(:,:) = rn_u10 - rn_uofac * un(:,:,1)
         ELSEWHERE
            zwndrel(:,:) =        - rn_uofac * un(:,:,1)
         END WHERE
         utau(:,:) = zrhocd * zwndrel(:,:) * zwndrel(:,:)

         zwndrel(:,:) = - rn_uofac * vn(:,:,1)
         vtau(:,:) = zrhocd * zwndrel(:,:) * zwndrel(:,:)
          
      ENDIF
      !
   END SUBROUTINE usrdef_sbc_oce

   SUBROUTINE usrdef_sbc_ice_tau( kt )
      INTEGER, INTENT(in) ::   kt   ! ocean time step
   END SUBROUTINE usrdef_sbc_ice_tau

   SUBROUTINE usrdef_sbc_ice_flx( kt )
      INTEGER, INTENT(in) ::   kt   ! ocean time step
   END SUBROUTINE usrdef_sbc_ice_flx

   !!======================================================================
END MODULE usrdef_sbc
