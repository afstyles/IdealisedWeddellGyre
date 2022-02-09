MODULE usrdef_istate
   !!======================================================================
   !!                     ***  MODULE usrdef_istate   ***
   !!
   !!                      ===  CANAL configuration  ===
   !!
   !! User defined : set the initial state of a user configuration
   !!======================================================================
   !! History :  NEMO ! 2017-11  (J. Chanut) Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  usr_def_istate : initial state in Temperature and salinity
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean space and time domain
   USE dom_oce        
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   !   
   USE usrdef_nam
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_istate   ! called by istate.F90

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_istate.F90 10425 2018-12-19 21:54:16Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
  
   SUBROUTINE usr_def_istate( pdept, ptmask, pts, pu, pv, pssh )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers
      !!                Here CANAL configuration 
      !!
      !! ** Method  :   Set a gaussian anomaly of pressure and associated
      !!                geostrophic velocities
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   pdept   ! depth of t-point               [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask             [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(  out) ::   pts     ! T & S fields      [Celsius ; g/kg]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pu      ! i-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pv      ! j-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    ! sea-surface height
      !
      INTEGER  :: ji, jj, jk, jl  ! dummy loop indices
      REAL(wp) :: zx, zy, zP0, zumax, zlambda, zr_lambda2, zn2, zf0, zH, zrho1, za, zf, zdzF
      REAL(wp) :: zpsurf, zdyPs, zdxPs
      REAL(wp) :: zdt, zdu, zdv
      REAL(wp) :: zjetx, zjety, zbeta
      REAL(wp), DIMENSION(jpi,jpj)  ::   zrandom
      !!----------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(jpi,jpj) :: pdept_jk
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate : CANAL configuration, analytical definition of initial state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   '
      !
      IF (ln_sshnoise) CALL RANDOM_NUMBER(zrandom)
      zjetx = ABS(rn_ujetszx)/2.
      zjety = ABS(rn_ujetszy)/2.
      !
      SELECT CASE(nn_initcase)
      CASE(0)    ! ACC satisfying thermal wind with a basin to the south at rest
         
         ! sea level:
         pssh(:,:) = 0.

         
         pu(:,:,:) = 0.
         pv(:,:,:) = 0.
         pts(:,:,:,jp_tem) = 0.
         pts(:,:,:,jp_sal) = rn_sponge_so

         IF( ln_sponge_uoconst ) THEN
            ! >>>>>>
            ! CASE 0a: ACC does not vary with space
            ! >>>>>>
   
            pu(:,:,:) = rn_sponge_uo
            pv(:,:,:) = rn_sponge_vo
            pts(:,:,:,jp_tem) = rn_sponge_to
   
         ELSE IF ( ln_sponge_uovar ) THEN
            ! >>>>>>
            ! CASE 0b : ACC varies sinusoidally with space
            ! >>>>>>

            !Calculate beta and f0 for the Coriolis parameter: f = zf0 + beta*(y-y0)
            zbeta = 1e3 * 2._wp * omega * COS( rad * rn_ppgphi0 ) / ra
            zf0   = 2._wp * omega * SIN( rad * rn_ppgphi0 )

            DO jk = 1,jpk
               
               pdept_jk = gdept_0(:,:,jk)          

               WHERE( ((gphiu > rn_d3 ).AND.(gphiu < rn_sponge_ly)) )
               ! Channel initial velocity >>>>>>>>>>>>>>
               ! In the space between the channel shelf and the northern boundary of the channel, 
               ! the velocity varies sinusoidally. The current is extended zonally across the
               ! whole domain
               ! _________________________________ y = rn_sponge_ly
               !     ->
               !     --->
               !     ----->
               !     ---> 
               !     ->                           ____ y = d3
               !%%%%%%%%% Channel shelf %%%%%%%%%%  \/  NO FLOW OVER
               !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  /\  THE CHANNEL SHELF
                  pu(:,:,jk) = rn_sponge_uomax * (SIN( rpi*(gphiu-rn_d3)/(rn_sponge_ly-rn_d3))**2)   &
                        &    * EXP(-pdept_jk / rn_depth_decay )

               END WHERE

               WHERE( gphiu >= rn_sponge_ly )
                  pu(:,:,jk) = 0.
               END WHERE    

               WHERE( (gphit > rn_d3).AND.(gphit < rn_sponge_ly) )
                  !
                  ! In the space between the channel shelf and the channel's northern boundary,
                  ! the target temperature satisfies the thermal wind relation of the target velocity
                  ! and has a stratification prescribed by rn_sponge_tomax and rn_depth_decay

                  ! _________ T_max_______________________ y = rn_sponge_ly
                  !           |
                  !          /
                  !         /
                  !        /
                  ! T_min |                          ____ y = d3
                  !%%%%%% |  % Channel shelf %%%%%%%%%%  \/ Temperature does not vary 
                  !%%%%%% |  %%%%%%%%%%%%%%%%%%%%%%%%%%  /\ horizontally over channel shelf
                  !
                  pts(:,:,jk,jp_tem) = + ( rn_sponge_uomax / (grav*rn_a0_user*rn_depth_decay*4*(rpi**2)) )            &
                                    & * EXP(-pdept_jk/rn_depth_decay) * (                                            &
                                    & 2*(rpi**2)*zf0*(rn_sponge_ly - gphit)                                          &
                                    & + (rpi**2)*zbeta*(rn_sponge_ly**2 - gphit**2)                                  &
                                    & + rpi*(zf0 + zbeta*gphit)*(rn_sponge_ly - rn_d3)*SIN(2*rpi*(gphit-rn_d3)/(rn_sponge_ly - rn_d3))   &
                                    & - zbeta*((rn_sponge_ly - rn_d3)**2)*(SIN(rpi*(gphit - rn_d3)/(rn_sponge_ly - rn_d3))**2) )         &
                                    & * 1e3  &  !Conversion for km --> m
                                    & + (rn_sponge_tomax-rn_sponge_tobot) * exp(-pdept_jk/rn_depth_decay)            &
                                    & + rn_sponge_tobot

               END WHERE

               WHERE( gphit >= rn_sponge_ly )
                  !
                  ! In the Northern sponge layer, the initial temperature does not vary horizontally
                  ! (equal to horizontal max of ACC temperature)
                  !
                  pts(:,:,jk,jp_tem) = (rn_sponge_tomax-rn_sponge_tobot) * exp(-pdept_jk/rn_depth_decay) + rn_sponge_tobot
                  !
               END WHERE
               !
               
               WHERE( (gphit <= rn_d3))
               !
               ! Above the channel shelf and in the basin, the temperature does not vary horizontally and is cooler
               ! than the temperature in the Northern sponge layer
               ! (equal to horizontal min of ACC temperature)
               !
                  pts(:,:,jk,jp_tem) = + ( rn_sponge_uomax / (grav*rn_a0_user*rn_depth_decay*4*(rpi**2)) )            &
                                       & * EXP(-pdept_jk/rn_depth_decay) * (                                         &
                                       & 2*(rpi**2)*zf0*(rn_sponge_ly - rn_d3)                                       & 
                                       & + (rpi**2)*zbeta*(rn_sponge_ly**2 - rn_d3**2) )                             & 
                                       & * 1e3  &  !Conversion for km --> m
                                       & + (rn_sponge_tomax-rn_sponge_tobot) * exp(-pdept_jk/rn_depth_decay)         &
                                       & + rn_sponge_tobot
               END WHERE
               !
            END DO
            !
         END IF
        !
      CASE(1)    ! Start from complete rest

         ! sea level:
         pssh(:,:) = 0.


         pu(:,:,:) = 0.
         pv(:,:,:) = 0.
         pts(:,:,:,jp_tem) = rn_sponge_to
         pts(:,:,:,jp_sal) = rn_sponge_so

      END SELECT
         
      IF (ln_sshnoise) THEN
         CALL RANDOM_NUMBER(zrandom)
         pssh(:,:) = pssh(:,:) + ( 0.1  * zrandom(:,:) - 0.05 )
      END IF
      CALL lbc_lnk( 'usrdef_istate', pssh, 'T',  1. )
      CALL lbc_lnk(  'usrdef_istate', pts, 'T',  1. )
      CALL lbc_lnk(   'usrdef_istate', pu, 'U', -1. )
      CALL lbc_lnk(   'usrdef_istate', pv, 'V', -1. )

   END SUBROUTINE usr_def_istate
     
   !!======================================================================
END MODULE usrdef_istate
