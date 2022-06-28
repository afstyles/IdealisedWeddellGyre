MODULE usrdef_zgr
   !!======================================================================
   !!                       ***  MODULE  usrdef_zgr  ***
   !!
   !!                      ===  CANAL configuration  ===
   !!
   !! User defined : vertical coordinate system of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2017-11  (J. Chanut)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_zgr   : user defined vertical coordinate system
   !!      zgr_z      : reference 1D z-coordinate 
   !!      zgr_top_bot: ocean top and bottom level indices
   !!      zgr_zco    : 3D verticl coordinate in pure z-coordinate case
   !!---------------------------------------------------------------------
   USE oce            ! ocean variables
   USE dom_oce        ! ocean domain
   USE phycst         ! physical constants
   USE usrdef_nam
   USE depth_e3       ! depth <=> e3
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp        ! distributed memory computing library
   USE iom   

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_zgr        ! called by domzgr.F90
   PUBLIC   ideal_WG_bath      ! called by usrdef_sbc.F90

  !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: usrdef_zgr.F90 10425 2018-12-19 21:54:16Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS             

   SUBROUTINE usr_def_zgr( ld_zco  , ld_zps  , ld_sco  , ld_isfcav,    &   ! type of vertical coordinate
      &                    pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d  ,    &   ! 1D reference vertical coordinate
      &                    pdept , pdepw ,                             &   ! 3D t & w-points depth
      &                    pe3t  , pe3u  , pe3v   , pe3f ,             &   ! vertical scale factors
      &                    pe3w  , pe3uw , pe3vw         ,             &   !     -      -      -
      &                    k_top  , k_bot,                             &   ! top & bottom ocean level
      &                    psponge_gamma_u, psponge_gamma_v, psponge_gamma_t , &
      &                    ptarget_uo, ptarget_vo, ptarget_to, ptarget_so         )    
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE usr_def_zgr  ***
      !!
      !! ** Purpose :   User defined the vertical coordinates
      !!
      !!----------------------------------------------------------------------
      LOGICAL                   , INTENT(out) ::   ld_zco, ld_zps, ld_sco      ! vertical coordinate flags
      LOGICAL                   , INTENT(out) ::   ld_isfcav                   ! under iceshelf cavity flag
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pdept_1d, pdepw_1d          ! 1D grid-point depth     [m]
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pe3t_1d , pe3w_1d           ! 1D grid-point depth     [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pdept, pdepw                ! grid-point depth        [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pe3t , pe3u , pe3v , pe3f   ! vertical scale factors  [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(out) ::   pe3w , pe3uw, pe3vw         ! i-scale factors 
      INTEGER , DIMENSION(:,:)  , INTENT(out) ::   k_top, k_bot                ! first & last ocean level
      !
      INTEGER  ::   inum   ! local logical unit
      REAL(WP) ::   z_zco, z_zps, z_sco, z_cav
      REAL(wp), DIMENSION(jpi,jpj) ::   z2d   ! 2D workspace
      !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(out) :: psponge_gamma_u, psponge_gamma_v, psponge_gamma_t
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(out) :: ptarget_uo, ptarget_vo, ptarget_to, ptarget_so
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_zgr : CANAL configuration (z-coordinate closed flat box ocean)'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      !
      !
      ! type of vertical coordinate
      ! ---------------------------
      ld_zco    = .FALSE.         ! CANAL case:  z-coordinate without ocean cavities
      !
      IF(nn_zgr_type .NE. 2) THEN
         ld_zps    = .TRUE.
         ld_sco    = .FALSE.
      ELSE
         ld_zps    = .FALSE.
         ld_sco    = .TRUE.
      END IF
      !
      ld_isfcav = .FALSE.
      !
      !
      ! Build the vertical coordinate system
      ! ------------------------------------
      IF(lwp) WRITE(numout,*) '    (1) Reference z coordinate system'
      CALL zgr_z( pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d )   ! Reference z-coordinate system
      !
      IF(lwp) WRITE(numout,*) '    (2) zgr_msk_top_bot   '
      CALL zgr_msk_top_bot( k_top , k_bot, z2d )                 ! masked top and bottom ocean t-level indices
      !
      IF(lwp) WRITE(numout,*) '    (3) zgr_zco           ' 
      !                                                     ! z-coordinate (3D arrays) from the 1D z-coord.
      CALL zgr_zco( pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d,   &   ! in  : 1D reference vertical coordinate
         &          pdept   , pdepw   ,                     &   ! out : 3D t & w-points depth
         &          pe3t    , pe3u    , pe3v   , pe3f   ,   &   !       vertical scale factors
         &          pe3w    , pe3uw   , pe3vw             )     !           -      -      -
      !
      !
      IF(lwp) WRITE(numout,*) '   (4) sponge_setup       '
      CALL sponge_setup( pdept, psponge_gamma_u, psponge_gamma_v, psponge_gamma_t,  &
               &   ptarget_uo, ptarget_vo, ptarget_to, ptarget_so)
      !
      !
   END SUBROUTINE usr_def_zgr


   SUBROUTINE zgr_z( pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d )   ! 1D reference vertical coordinate
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zgr_z  ***
      !!
      !! ** Purpose :   set the 1D depth of model levels and the resulting 
      !!              vertical scale factors.
      !!
      !! ** Method  :   1D z-coordinate system (use in all type of coordinate)
      !!       The depth of model levels is set from dep(k), an analytical function:
      !!                   w-level: depw_1d  = dep(k)
      !!                   t-level: dept_1d  = dep(k+0.5)
      !!       The scale factors are the discrete derivative of the depth:
      !!                   e3w_1d(jk) = dk[ dept_1d ] 
      !!                   e3t_1d(jk) = dk[ depw_1d ]
      !!           with at top and bottom :
      !!                   e3w_1d( 1 ) = 2 * ( dept_1d( 1 ) - depw_1d( 1 ) )
      !!                   e3t_1d(jpk) = 2 * ( dept_1d(jpk) - depw_1d(jpk) )
      !!       The depth are then re-computed from the sum of e3. This ensures 
      !!    that depths are identical when reading domain configuration file. 
      !!    Indeed, only e3. are saved in this file, depth are compute by a call
      !!    to the e3_to_depth subroutine.
      !!
      !!       Here the Madec & Imbard (1996) function is used.
      !!
      !! ** Action  : - pdept_1d, pdepw_1d : depth of T- and W-point (m)
      !!              - pe3t_1d , pe3w_1d  : scale factors at T- and W-levels (m)
      !!
      !! Reference : Marti, Madec & Delecluse, 1992, JGR, 97, No8, 12,763-12,766.
      !!             Madec and Imbard, 1996, Clim. Dyn.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pdept_1d, pdepw_1d   ! 1D grid-point depth        [m]
      REAL(wp), DIMENSION(:)    , INTENT(out) ::   pe3t_1d , pe3w_1d    ! 1D vertical scale factors  [m]
      !
      INTEGER  ::   jk       ! dummy loop indices
      REAL(wp) ::   zd, zsur, za0, za1, zkth, zacr, zt, zw       
      !!----------------------------------------------------------------------
      !
      zd = rn_domszz/FLOAT(jpkm1)
      !
      IF(lwp) THEN            ! Parameter print
         WRITE(numout,*)
         WRITE(numout,*) '    zgr_z   : Reference vertical z-coordinates '
         WRITE(numout,*) '    ~~~~~~~'
      ENDIF

      !
      ! 1D Reference z-coordinate    (using Madec and Imbard 1996 function)
      ! -------------------------
      !
      SELECT CASE(nn_zgr_type)
      
      CASE(0,2) ! Regular z grid
         IF(lwp) WRITE(numout,*)'        CANAL case : uniform vertical grid :'
         IF(lwp) WRITE(numout,*)'                     with thickness = ', zd
         pdepw_1d(1) = 0._wp
         pdept_1d(1) = 0.5_wp * zd
         !
         DO jk = 2, jpk          ! depth at T and W-points
            pdepw_1d(jk) = pdepw_1d(jk-1) + zd 
            pdept_1d(jk) = pdept_1d(jk-1) + zd 
         END DO
         !
      CASE(1) ! Madec and Imbard z grid
         !
         ! Set parameters of z(k) function
         ! -------------------------------
         zsur = -2033.194295283385_wp
         za0  =   155.8325369664153_wp
         za1  =   146.3615918601890_wp
         zkth =    17.28520372419791_wp
         zacr =     5.0_wp
         !
         IF(lwp) THEN
            WRITE(numout,*) '       CANAL case : MI96 function with the following coefficients :'
            WRITE(numout,*) '                 zsur = ', zsur
            WRITE(numout,*) '                 za0  = ', za0
            WRITE(numout,*) '                 za1  = ', za1
            WRITE(numout,*) '                 zkth = ', zkth
            WRITE(numout,*) '                 zacr = ', zacr
         END IF
         !
         DO jk = 1, jpk          ! depth at T and W-points
            zw = REAL( jk , wp )
            zt = REAL( jk , wp ) + 0.5_wp
            pdepw_1d(jk) = ( zsur + za0 * zw + za1 * zacr *  LOG( COSH( (zw-zkth) / zacr ) )  )
            pdept_1d(jk) = ( zsur + za0 * zt + za1 * zacr *  LOG( COSH( (zt-zkth) / zacr ) )  )
         END DO
         !
         !Scale the MI96 function to get as many interior vertical grid points as possible
         !
         pdept_1d(:) = pdept_1d(:) * rn_domszz/pdepw_1d(jpk)
         pdepw_1d(:) = pdepw_1d(:) * rn_domszz/pdepw_1d(jpk)
         !             
      END SELECT
      !
      !                       ! e3t and e3w from depth
      CALL depth_to_e3( pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d ) 
      !
      !                       ! recompute depths from SUM(e3)  <== needed
      CALL e3_to_depth( pe3t_1d, pe3w_1d, pdept_1d, pdepw_1d ) 
      !
      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*) '              Reference 1D z-coordinate depth and scale factors:'
         WRITE(numout, "(9x,' level  gdept_1d  gdepw_1d  e3t_1d   e3w_1d  ')" )
         WRITE(numout, "(10x, i4, 4f9.2)" ) ( jk, pdept_1d(jk), pdepw_1d(jk), pe3t_1d(jk), pe3w_1d(jk), jk = 1, jpk )
      ENDIF
      !
   END SUBROUTINE zgr_z


   SUBROUTINE zgr_msk_top_bot( k_top , k_bot, z2d )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE zgr_msk_top_bot  ***
      !!
      !! ** Purpose :   set the masked top and bottom ocean t-levels
      !!
      !! ** Method  :   CANAL case = closed flat box ocean without ocean cavities
      !!                   k_top = 1     except along north, south, east and west boundaries
      !!                   k_bot = jpk-1 except along north, south, east and west boundaries
      !!
      !! ** Action  : - k_top : first wet ocean level index
      !!              - k_bot : last  wet ocean level index
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(:,:), INTENT(out) ::   k_top , k_bot  ! first & last wet ocean level
      REAL, DIMENSION(:,:), INTENT(out)    ::   z2d            !
      REAL(wp), DIMENSION(jpi,jpj) ::   zH, znoise, z2d_tmp   ! 2D local workspace
      REAL(wp), DIMENSION(jpi,jpj) ::   shelf_frac ! depth of shelf divided by flat depth
      REAL(wp)                     ::   zmaxlam, zminlam, zscl, zminphi, iamp, iphase, ikx, iky, ilx, ily
      INTEGER                      ::   ip, np, noise_case, inum_noise
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '    zgr_top_bot : defines the top and bottom wet ocean levels.'
      IF(lwp) WRITE(numout,*) '    ~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*) '       CANAL case : closed flat box ocean without ocean cavities'
      !
      SELECT CASE(nn_botcase)
      CASE(2)
         zminphi = MINVAL(gphit)
         CALL mpp_min('usrdef_zgr', zminphi)
         !
         CALL ideal_WG_bath(glamt, gphit, rn_h_flat, rn_h_ork, rn_domszz, rn_hs,               &
              &             rn_x1, rn_x2, rn_x3,                                               &    
              &             rn_x_ork, rn_y2, rn_d1, rn_d2, rn_d3, rn_s1, rn_s2, -zminphi,      &
              &             ln_orkney, ln_fh,                                                  &
              &             rn_r0, rn_r1, ff_t, rn_H_mbump, ln_mbump, rn_d_mbump,   zH, z2d,   &
              &             shelf_frac )

      END SELECT
      !
      !
      IF(ln_bathy_noise) THEN
         !Add sinusoidal noise to the Weddell GYRE bathymetry
         
         CALL iom_open( trim(cn_noise_file), inum_noise)
         CALL iom_get( inum_noise, jpdom_autoglo_xy, cn_noise_varname, znoise, ldxios=lrxios)
         CALL iom_close( inum_noise )

         znoise = rn_noise_scale * znoise

      END IF
      !
      zmaxlam = MAXVAL(glamt)
      zminlam = MINVAL(glamt)
      CALL mpp_max( 'usrdef_zgr', zmaxlam )
      CALL mpp_min( 'usrdef_zgr', zminlam )
      !
      IF(lwp) WRITE(numout,*) 'zmaxlam =', zmaxlam
      IF(lwp) WRITE(numout,*) 'zminlam =', zminlam
      !
      !
      !
      IF(ln_bathy_noise) THEN
         !
         CALL H_to_ind(zH + znoise, z2d_tmp)
         !
         WHERE(( z2d > 1 ).AND.(z2d_tmp > 1))
            zH(:,:) = zH(:,:) + znoise(:,:)
         END WHERE
         !
         CALL H_to_ind(zH, z2d)
         !
      END IF

      WHERE( ((gphit <= 0.).AND.(glamt <= 0)).OR.((gphit >= rn_sponge_ly).AND.(glamt <= 0)) )
         z2d(:,:) = 0.
      END WHERE

      CALL lbc_lnk( 'usrdef_zgr', z2d, 'T', 1. )           ! set surrounding land to zero (here jperio=0 ==>> closed)
      !
      k_bot(:,:) = INT( z2d(:,:) )           ! =jpkm1 over the ocean point, =0 elsewhere
      !
      k_top(:,:) = MIN( 1 , k_bot(:,:) )     ! = 1    over the ocean point, =0 elsewhere
      !
   END SUBROUTINE zgr_msk_top_bot
   

   SUBROUTINE zgr_zco( pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d,   &   ! in : 1D reference vertical coordinate
      &                pdept   , pdepw   ,                     &   ! out: 3D t & w-points depth
      &                pe3t    , pe3u    , pe3v   , pe3f   ,   &   !      vertical scale factors
      &                pe3w    , pe3uw   , pe3vw             )     !          -      -      -
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zgr_zco  ***
      !!
      !! ** Purpose :   define the reference z-coordinate system
      !!
      !! ** Method  :   set 3D coord. arrays to reference 1D array 
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:)    , INTENT(in   ) ::   pdept_1d, pdepw_1d          ! 1D grid-point depth       [m]
      REAL(wp), DIMENSION(:)    , INTENT(in   ) ::   pe3t_1d , pe3w_1d           ! 1D vertical scale factors [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pdept, pdepw                ! grid-point depth          [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3t , pe3u , pe3v , pe3f   ! vertical scale factors    [m]
      REAL(wp), DIMENSION(:,:,:), INTENT(  out) ::   pe3w , pe3uw, pe3vw         !    -       -      -
      !
      INTEGER  ::   jk
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !Additional parameters due to partial cell addition (A. Styles)
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !
      INTEGER(wp) :: ji, jj, jj_glo, ji_glo
      REAL(wp), DIMENSION(jpi,jpj):: z2d, zH, shelf_frac
      INTEGER, DIMENSION(jpi, jpj) :: k_top, k_bot
      REAL(wp) :: pe3t_p, pe3u_p, pe3v_p, pe3f_p, pe3w_p, pe3vw_p, pe3uw_p
      REAL(wp) ::   zmaxlam, zminlam, zscl, zminphi
      !
      !-------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*) 'Adjusting depths of cells' 
      !
      !
      DO jk = 1, jpk
         pdept(:,:,jk) = pdept_1d(jk)
         pdepw(:,:,jk) = pdepw_1d(jk)
         pe3t (:,:,jk) = pe3t_1d (jk)
         pe3u (:,:,jk) = pe3t_1d (jk)
         pe3v (:,:,jk) = pe3t_1d (jk)
         pe3f (:,:,jk) = pe3t_1d (jk)
         pe3w (:,:,jk) = pe3w_1d (jk)
         pe3uw(:,:,jk) = pe3w_1d (jk)
         pe3vw(:,:,jk) = pe3w_1d (jk)
      END DO
      !
      !
      zminphi = MINVAL(gphit)
      CALL mpp_min('usrdef_zgr', zminphi)
      
      !Call analytic form of idealized bathymetry
      CALL ideal_WG_bath(glamt, gphit, rn_h_flat, rn_h_ork,  rn_domszz, rn_hs,   &                       
              &     rn_x1, rn_x2, rn_x3, &
              &     rn_x_ork, rn_y2, rn_d1, rn_d2, rn_d3, rn_s1, rn_s2,          &
              &     -zminphi, ln_orkney, ln_fh,  &
              &     rn_r0, rn_r1, ff_t, rn_H_mbump, ln_mbump, rn_d_mbump,   zH,  &
              &     z2d, shelf_frac )
 
      !Call bottom mask, top mask, and z2d
      CALL zgr_msk_top_bot( k_top , k_bot, z2d)
      !
      SELECT CASE( nn_zgr_type )
      !
      CASE(0,1)
      !
      !Addition by Andrew Styles >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !Adding a partial cell representation of the bathymetry
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !
      DO ji = 1,jpi
         DO jj = 1, jpj
            !
            jk = k_bot(ji,jj) !Load last wet index for column
            !
            !
            !Scale the last wet grid cell by fraction covered by bathymetry 
            pe3t_p = (z2d(ji,jj)-jk)*pe3t(ji,jj,jk)
            
            !Only proceed if thickness is positive
            IF ((pe3t_p > 0)) THEN
               !
               !If thickness is smaller than 20m this will be unstable
               !Therefore set minimum cell thickness to 20m
               IF ((pe3t_p < 20)) THEN
                  !
                  pe3t_p = 20.0
                  !
               END IF
               !
               !Translate centre point (pdept) and bottom point (pdepw) towards the top sea surface 
               !Centre (T point) upwards moves by half the change in thickness
               !
               pdept(ji,jj,jk) = pdept(ji,jj,jk) - (pe3t(ji,jj,jk) - pe3t_p)/2.0
               !
               !Save partial thicknesses as new cell thickness
               pe3t(ji,jj,jk) = pe3t_p
               !
               jj_glo = mjg(jj) !Determine global domain index
               ji_glo = mig(ji)
               !
               IF (jk .ne. 1) THEN
                   !
                   pe3w(ji,jj,jk) = (pe3t(ji,jj,jk) + pe3t(ji,jj,jk-1))/2
                   !
               ELSE
                   !
                   pe3w(ji,jj,jk) = pe3t(ji,jj,jk)
                   !
               END IF               

               IF ( ji .ne. jpi ) THEN
                  !
                  pe3u(ji,jj,jk) = (pe3t(ji,jj,jk) + pe3t(ji+1,jj,jk))/2
                  pe3uw(ji,jj,jk) = (pe3w(ji,jj,jk) + pe3w(ji+1,jj,jk))/2
                  !
               ELSE
                  !
                  pe3u(ji,jj,jk) = pe3t(ji,jj,jk)
                  pe3uw(ji,jj,jk) = pe3w(ji,jj,jk)
                  !
               END IF
               !
               IF (jj .ne. jpj) THEN
                  !
                  pe3v(ji,jj,jk) = (pe3t(ji,jj,jk) + pe3t(ji,jj+1,jk))/2
                  pe3vw(ji,jj,jk) = (pe3w(ji,jj,jk) + pe3w(ji,jj+1,jk))/2
                  !
               ELSE
                  !
                  pe3v(ji,jj,jk) = pe3t(ji,jj,jk)
                  pe3vw(ji,jj,jk) = pe3w(ji,jj,jk)
                  !
               END IF
               !
               IF ((ji .ne. jpi).AND.(jj .ne. jpj)) THEN
                  !
                  pe3f(ji,jj,jk) = (pe3t(ji,jj,jk) + pe3t(ji+1,jj,jk) &
                                & + pe3t(ji,jj+1,jk) + pe3t(ji+1,jj+1,jk) )/4
                  !
               ELSE IF ( (ji .ne. jpi).AND.(jj == jpj) ) THEN
                  !
                  pe3f(ji,jj,jk) = (pe3t(ji,jj,jk)+ pe3t(ji+1,jj,jk))/2
                  !
               ELSE IF ( (ji == jpi).AND.(jj .ne. jpj) ) THEN
                  !
                  pe3f(ji,jj,jk) = (pe3t(ji,jj,jk) + pe3t(ji,jj+1,jk) )/2
                  !
               ELSE
                  !
                  pe3f(ji,jj,jk) = pe3t(ji,jj,jk)
                  !
               END IF
               !
            END IF
         END DO
      END DO

      ! Work in Progress - Implementation of hybrid coordinates
      !CASE(2)
      !   !
      !   DO ji = 1, jpi
      !      DO jj= 1, jpj
      !         !
      !         !Scale all T cell thicknesses equally if an unmasked column
      !         IF (k_bot(ji,jj) .NE. 0) THEN
      !            pe3t(ji,jj,:) =  zH(ji,jj)/k_bot(ji,jj)
      !         ELSE
      !           pe3t(ji,jj,:)  =  0.
      !         END IF
      !         !
      !         pdept(ji,jj,1) = e3t(1)/2
      !         DO jk = 2, jpk
      !            pdept(ji,jj,jk) = pdept(ji,jj,jk-1) + 0.5*(e3t(jk-1) + e3t(jk))
      !         END DO 
      !         !
      !         DO jk = 1,jpk
      !            IF (jk .ne. 1) THEN
      !                !
      !                pe3w(ji,jj,jk) = (pe3t(ji,jj,jk) + pe3t(ji,jj,jk-1))/2
      !                !
      !            ELSE
      !                !
      !                pe3w(ji,jj,jk) = pe3t(ji,jj,jk)
      !                !
      !            END IF
      !         END DO
      !
      !         pdepw(ji,jj,1) = 0.
      !         DO jk = 2, jpk
      !            pdepw(ji,jj,jk) = pdepw(ji,jj,jk-1) + 0.5*(e3w(jk-1) + e3w(jk))
      !         END DO
      !
      !         IF ( ji .ne. jpi ) THEN
      !            !
      !            pe3u(ji,jj,:) = (pe3t(ji,jj,:) + pe3t(ji+1,jj,:))/2
      !            pe3uw(ji,jj,:) = (pe3w(ji,jj,:) + pe3w(ji+1,jj,:))/2
      !            !
      !         ELSE
      !            !
      !            pe3u(ji,jj,:) = pe3t(ji,jj,:)
      !            pe3uw(ji,jj,:) = pe3w(ji,jj,:)
      !            !
      !         END IF
      !         !
      !         IF (jj .ne. jpj) THEN
      !            !
      !            pe3v(ji,jj,:) = (pe3t(ji,jj,:) + pe3t(ji,jj+1,:))/2
      !            pe3vw(ji,jj,:) = (pe3w(ji,jj,:) + pe3w(ji,jj+1,:))/2
      !            !
      !         ELSE
      !            !
      !            pe3v(ji,jj,:) = pe3t(ji,jj,:)
      !            pe3vw(ji,jj,:) = pe3w(ji,jj,:)
      !            !
      !         END IF
      !         !
      !         IF ((ji .ne. jpi).AND.(jj .ne. jpj)) THEN
      !            !
      !            pe3f(ji,jj,:) = (pe3t(ji,jj,:) + pe3t(ji+1,jj,:) &
      !                          & + pe3t(ji,jj+1,:) + pe3t(ji+1,jj+1,:) )/4
      !            !
      !         ELSE IF ( (ji .ne. jpi).AND.(jj == jpj) ) THEN
      !            !
      !            pe3f(ji,jj,:) = (pe3t(ji,jj,:)+ pe3t(ji+1,jj,:))/2
      !            !
      !         ELSE IF ( (ji == jpi).AND.(jj .ne. jpj) ) THEN
      !            !
      !            pe3f(ji,jj,:) = (pe3t(ji,jj,:) + pe3t(ji,jj+1,:) )/2
      !            !
      !         ELSE
      !            !
      !            pe3f(ji,jj,:) = pe3t(ji,jj,:)
      !            !
      !         END IF
      !         ! 
      !      END DO
      !   END DO 
      END SELECT         
     
   END SUBROUTINE zgr_zco

   !!======================================================================

   SUBROUTINE ideal_WG_bath(x_grid, y_grid, H_flat, H_ork, H_max, H_s,                   &
              &             x1, x2, x3, x_ork, y2,                                       & 
              &             d1, d2, d3, s1, s2, yacc, ln_orkney, ln_fh, r0, r1,          & 
              &             ff_t, H_mbump,                                               &
              &             ln_mbump, d_mbump, H, z2d, shelf_frac)

            implicit none
            !
            real, intent(in) :: x_grid(:,:)
            real, intent(in) :: y_grid(SIZE(x_grid,1),SIZE(x_grid,2))
	         real, intent(in) :: ff_t(SIZE(x_grid,1),SIZE(x_grid,2))
            !
            real, intent(in) :: H_flat, H_ork, H_max, H_mbump, H_s
            real, intent(in) :: x1, x2, x3, x_ork, y2, d1, d2, d3, yacc, r0, r1, d_mbump
            real, intent(in) :: s1, s2
            logical, intent(in) :: ln_orkney, ln_fh, ln_mbump
            !
            real, intent(out) :: H(SIZE(x_grid,1), SIZE(x_grid,2))
            real, intent(out) :: z2d(SIZE(x_grid,1), SIZE(x_grid,2))
            real, intent(out) :: shelf_frac(SIZE(x_grid,1), SIZE(x_grid,2))
            !
            real :: r(SIZE(x_grid,1),SIZE(x_grid,2))
            real :: g_r(SIZE(x_grid,1),SIZE(x_grid,2))
            real :: sinx(SIZE(x_grid,1),SIZE(x_grid,2))
            real :: y_grid2(SIZE(x_grid,1),SIZE(x_grid,2))
            real :: x_cent, y_cent, f_cent, x_max
            integer :: ind_cent(2)

            real :: fn(SIZE(x_grid,1),SIZE(x_grid,2)) !Dummy functions for calculating North, South and West
            real :: fw(SIZE(x_grid,1),SIZE(x_grid,2)) !boundary values for subdomains
            real :: fs(SIZE(x_grid,1),SIZE(x_grid,2))
            !
            logical :: domain_tmp(SIZE(x_grid,1),SIZE(x_grid,2)) !Dummy mask for highlighting temporary 
                                                                 !domains of interest
            !
            real :: dx, dy  !Dummy x and y widths
            !
            H(:,:) = H_flat
            shelf_frac(:,:) = 0.

            x_max = MAXVAL(x_grid)
            CALL mpp_max( 'usrdef_zgr', x_max )

            y_grid2 = y_grid + yacc

            
            !Sinusoidal slope on Western boundary
            WHERE( (x_grid < d1).AND.(x_grid > s1).AND.(y_grid2 > x1).AND.(y_grid2 <= yacc))
                !
                H(:,:) = (H_flat-H_s)*SIN( 0.5*rpi*(x_grid-s1) / (d1-s1) )**2 + H_s
                !
            END WHERE

            WHERE( (x_grid <= s1).AND.(y_grid2 > x1).AND.(y_grid2 <= yacc))
            !
            H(:,:) = H_s
            !
            END WHERE
        
            !Circular curvature of basin on SW corner
            r = SQRT( (x_grid-x1)**2 + (y_grid2-x1)**2 )
            WHERE( ((r < x1-s1).AND.(r > (x1-d1) ).AND.(x_grid <= x1).AND.(y_grid2 <=x1)) )
                !
                H(:,:) = (H_flat-H_s)*SIN( 0.5*rpi*(r - x1 + s1)/(d1-s1))**2 + H_s
                !
            END WHERE

            WHERE( ((r <= x1).AND.(r >= (x1-s1) ).AND.(x_grid <= x1).AND.(y_grid2 <=x1)) )
            !
            H(:,:) = H_s
            !
            END WHERE

            WHERE( ((r > x1).AND.(x_grid < x1).AND.(y_grid2 <=x1)) )
            !
            H(:,:) = 0.
            !
            END WHERE
            !
	    !Transition from meridional shelf of thickness d1 to zonal shelf of thickness d2 >>>
            !
            sinx = SIN(0.5*rpi*(x_grid - x1)/(x2 - x1))**2

            WHERE( (x_grid < x2).AND.(x_grid >=x1).AND.(y_grid2 >= (y2+s2-s1)*sinx + s1).AND.(y_grid2 <= (y2+d2-d1)*sinx + d1) )
            !
            H(:,:) = (H_flat-H_s)*SIN( 0.5*rpi*(y_grid2 - (y2+s2-s1)*sinx - s1)/( (d2-s2-d1+s1)*sinx + d1-s1 ) )**2 + H_s
            !
            END WHERE

            WHERE( (x_grid < x2).AND.(x_grid >=x1).AND.(y_grid2 < (y2+s2-s1)*sinx + s1).AND.(y_grid2 >= y2*sinx) )
            !
            H(:,:) = H_s
            !
            END WHERE

            WHERE( (x_grid < x2).AND.(x_grid >=x1).AND.(y_grid2 < y2*sinx))
                !
                H(:,:) = 0.
                !
            END WHERE
            
            !Zonal shelf of thickness d2 >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

            WHERE( (x_grid >= x2).AND.(y_grid2 >= y2 + s2).AND.(y_grid2 < (y2+d2)) )
                !
                H(:,:) = (H_flat-H_s)*SIN( 0.5*rpi*(y_grid2-y2-s2) / (d2-s2) )**2 + H_s
                !
            END WHERE

            WHERE( (x_grid >= x2).AND.(y_grid2 < y2 + s2).AND.(y_grid2 >=y2) )
            !
            H(:,:) = H_s
            !
            END WHERE


            WHERE( (x_grid >= x2).AND.(y_grid2 < y2) )
                !
                H(:,:) = 0.
                !
            END WHERE

            !Circular curvature of basin on SE corner >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            r = SQRT( (x_grid - (x_max-x3)) **2 + (y_grid2-(y2+s2+x3))**2 )
            WHERE( ((r <= x3).AND.(r > x3+s2-d2 ).AND.(x_grid >= (x_max - x3) ).AND.(y_grid2 <= y2+ x3 + s2)) )
                !
                H(:,:) = (H_flat-H_s)*SIN( 0.5*rpi*(r-x3)/(d2-s2) ) **2 + H_s
                !
            END WHERE

            WHERE( ((r > x3).AND.(r < x3 + s2).AND.(x_grid > (x_max - x3) ).AND.(y_grid2 <= y2 + x3 + s2)) )
            !
            H(:,:) = H_s
            !
            END WHERE

            WHERE( ((r >= x3+s2).AND.(x_grid > (x_max - x3) ).AND.(y_grid2 <= y2 + x3 + s2)) )
            !
            H(:,:) = 0.
            !
            END WHERE

            !Sinusoidal slope on Eastern boundary >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            WHERE( (x_grid > x_max - d2 + s2).AND.(x_grid <= x_max).AND.(y_grid2 > y2 + x3 + s2).AND.(y_grid2 <= yacc))
                !
                H(:,:) = (H_flat-H_s)*SIN( 0.5*rpi*(x_grid - x_max) / (d2-s2) )**2 + H_s
                !
            END WHERE
            
            shelf_frac = ABS(H(:,:)/H_flat) !Fractional depth over basin shelves only

            !Optional f/H perturbation >>>>>>>>>>>>>>>>>>>>>>>>>>>>
            IF ( ln_fh ) THEN

                x_cent = x2 + r1
                y_cent = y2 + d2 + r1

                r = SQRT((x_grid-x_cent)**2 + (y_grid2-y_cent)**2)

                ind_cent = minloc( (x_grid-x_cent)**2 + (y_grid2-y_cent)**2 )
                f_cent = ff_t( ind_cent(1), ind_cent(2) )

                g_r = (f_cent/H_max)*(1-r/r0) + (f_cent/H_flat)*(r/r0)

                WHERE( (r <= r1).AND.(r > r0) )
                    !
                    H(:,:) = ff_t/(  (ff_t/H_flat)*SIN(0.5*rpi*(r-r0)/(r1-r0))**2 + g_r*COS(0.5*rpi*(r-r0)/(r1-r0))**2 )
                    !
                END WHERE

                WHERE( r <= r0)
                    !
                    H(:,:) = ff_t/g_r
                    !
                END WHERE

            END IF

            !Optional introduction of an Orkney passage >>>>>>>>>>>
            IF( ln_orkney) THEN

                WHERE( (x_grid < d1).AND.(x_grid > s1).AND.(y_grid2 > yacc - d3).AND.(y_grid2 <= yacc)  )
                !
                H(:,:) = ( (H_flat-H_s) * (SIN(0.5*rpi*(x_grid-s1)/(d1-s1))**2)  + H_s )*(SIN(0.5*rpi*(y_grid2-yacc)/d3)**2) &
                      &  +( H_ork * (SIN(0.5*rpi*(x_grid-s1)/(d1-s1))**2) )*(COS(0.5*rpi*(y_grid2-yacc)/d3)**2)
                !
                END WHERE

                WHERE( (x_grid <= s1).AND.(y_grid2 > yacc - d3).AND.(y_grid2 <= yacc)  )
                !
                H(:,:) = H_s
                !
                END WHERE
                
                WHERE( (x_grid < d1).AND.(x_grid > s1).AND.(y_grid2 > yacc).AND.(y_grid2 <= yacc+d3)  )
                !
                H(:,:) = (H_flat-H_s) * (SIN(0.5*rpi*(y_grid2-yacc)/d3)**2) &
                      &  + H_ork * (SIN(0.5*rpi*x_grid/d1)**2)*(COS(0.5*rpi*(y_grid2-yacc)/d3)**2)
                !
                END WHERE

                WHERE( (x_grid <= s1).AND.(y_grid2 > yacc).AND.(y_grid2 <= yacc+d3)  )
                !
                H(:,:) = (H_flat - H_s) * (SIN(0.5*rpi*(y_grid2-yacc)/d3)**2) + H_s 
                !
                END WHERE
                
                WHERE( (x_grid <= 0).AND.(y_grid2 >= yacc).AND.(y_grid2 <= yacc+d3)  )
                !
                H(:,:) = (H_flat-H_s) * (SIN(0.5*rpi*(y_grid2-yacc)/d3)**2) + H_s
                !
                END WHERE
                
                WHERE( ((x_grid >= d1).AND.(x_grid <= x_ork).AND.(y_grid2 >= yacc - d3).AND.(y_grid2 <= yacc + d3)) )   !North of Okney
                !
                    H(:,:) = H_flat *(SIN(0.5*rpi*(y_grid2-yacc)/d3)**2) &
                      &  + H_ork *(COS(0.5*rpi*(y_grid2-yacc)/d3)**2)
                !
                END WHERE
                
                r = SQRT( (x_grid-x_ork)**2 + (y_grid2 - yacc)**2 )
                
                WHERE( (x_grid > x_ork).AND.(x_grid <= x_ork + d3).AND.(r < d3) )   !North of Okney
                !
                    H(:,:) = H_flat *(SIN(0.5*rpi*r/d3)**2) &
                        &  + H_ork *(COS(0.5*rpi*r/d3)**2)
                !
                END WHERE 
                ! 
            ELSE

               !If an Orkney ridge isn't used, use an elliptical transition to the channel shelf of thickness d3
               
                r = SQRT( (x_grid/d1)**2 + ((y_grid2-yacc)/d3)**2 )

                WHERE( (r <= 1).AND.(x_grid >=0).AND.(y_grid2 >= yacc).AND.(y_grid2 <= yacc + d3))
                    !
                    H(:,:) = H_flat*(SIN(0.5*rpi*r)**2)
                    !
                END WHERE

                WHERE( (x_grid < 0).AND.(y_grid2 >= yacc).AND.(y_grid2 <= yacc + d3))
                    !
                    H(:,:) = H_flat*(SIN(0.5*rpi*(y_grid2-yacc)/d3)**2)
                    !
                END WHERE

            END IF

            !Optional meridional bump in the ACC channel
            IF( ln_mbump ) THEN
               !
               dx = d_mbump
               dy = d3
               !

               fn(:,:) = 1.
               fs(:,:) = 1.
               fw(:,:) = 1.

               domain_tmp = (x_grid < 0.).AND.(x_grid > -dx).AND.(y_grid > 0).AND.(y_grid < dy)
               ! Note that we are deliberately using y_grid rather than y_grid2 as y_grid == 0 at the
               ! bottom of the ACC channel

               WHERE( domain_tmp)
                   fn(:,:) = 1 - (1 - H_mbump/H_flat)*(SIN(rpi*(x_grid+dx)/dx)**2)
               END WHERE

               WHERE( domain_tmp )
                   fw(:,:) = SIN(0.5*rpi*y_grid/d3)**2
               END WHERE

               WHERE( domain_tmp )
                   H(:,:) = H_flat*( fw(:,:)*fn(:,:) )
               END WHERE
               ! 

               domain_tmp = (x_grid < 0.).AND.(x_grid > -dx).AND.(y_grid >=dy)

               WHERE( domain_tmp )
                   H(:,:) = H_flat * ( 1 - (1 - H_mbump/H_flat)*(SIN(rpi*(x_grid+dx)/dx)**2) )
               END WHERE
                   
           END IF

            !Elliptical curvature of basin on NE corner >>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            r = SQRT( ((x_grid-x_max)/(d2-s2))**2 + ((y_grid2 - yacc)/d3)**2 )
                !
            WHERE( (r < 1).AND.(x_grid <=x_max).AND.(y_grid2 >= yacc).AND.(y_grid2 <= yacc + d3))
            !    !
                H(:,:) = H_flat*(SIN( 0.5*rpi*r)**2)
            !    !
            END WHERE 

            CALL H_to_ind(H,z2d) 

      END SUBROUTINE ideal_WG_bath   

   SUBROUTINE sponge_setup(pdept, psponge_gamma_u, psponge_gamma_v, psponge_gamma_t,  &
                           ptarget_uo, ptarget_vo, ptarget_to, ptarget_so)

      REAL(wp), intent(in) , DIMENSION(jpi,jpj,jpk) :: pdept
      REAL(wp), intent(out), DIMENSION(jpi,jpj) ::   psponge_gamma_u, psponge_gamma_v, psponge_gamma_t
      REAL(wp), intent(out), DIMENSION(jpi,jpj,jpk) :: ptarget_uo, ptarget_vo, ptarget_to, ptarget_so
      REAL(wp), DIMENSION(jpi,jpj) :: pdept_jk, ptarget_uo_jk, ptarget_to_jk, dt_dy, pff_v
      REAL(wp) :: zf0, zbeta, gphimax, xs1, xs2 
      INTEGER :: jk, jj, ji
      REAL(wp) :: gamma_maxval

      !Relaxation parameter (gamma) configuration >>>>>>>>>>>>>>>>>>>>>>
      psponge_gamma_u(:,:) = 0._wp
      psponge_gamma_v(:,:) = 0._wp
      psponge_gamma_t(:,:) = 0._wp
      !
      IF(lwp) WRITE(numout,*) ' Set sponge values to 0. ' 
      !Left and right x coordinates of the channel sponge
      xs1 = -(rn_chan_lx + rn_sponge_lx)/2         !_____________________
      xs2 = -(rn_chan_lx - rn_sponge_lx)/2         !     | %%%%%% |       |
      !                                            !     | Sponge |       |
      !                                            !_____|________|______ |
      !                                                 xs1      xs2     x=0

      
      !Channel sponge relaxation parameters >>>>>>>>>>>>>>
      IF(ln_sponge_chan_mom) THEN
         ! The sponge is independent of y and varies sinusoidally wrt x over the length rn_sponge_lx
         IF(lwp) WRITE(numout,*) "Channel momentum sponge is on with max gamma = ", rn_sponge_gm
         WHERE( ((gphiu > 0 ).AND.(gphiu < rn_sponge_ly).AND.(glamu >= xs1).AND.(glamu <=xs2)) )
            psponge_gamma_u(:,:) = rn_sponge_gm*( SIN(rpi*(glamu-xs1)/rn_sponge_lx)**2 )
         END WHERE

         WHERE( ((gphiv > 0 ).AND.(gphiv < rn_sponge_ly).AND.(glamv >= xs1).AND.(glamv <= xs2)))
            psponge_gamma_v(:,:) = rn_sponge_gm*( SIN(rpi*(glamv-xs1)/rn_sponge_lx)**2 )
         END WHERE

      END IF

      IF(ln_sponge_chan_tra) THEN
         IF(lwp) WRITE(numout,*) "Channel tracer sponge is on with max gamma = ", rn_sponge_gm_t
         WHERE( ((gphit > 0 ).AND.(gphit < rn_sponge_ly).AND.(glamt >= xs1).AND.(glamt <= xs2)))
            psponge_gamma_t(:,:) = rn_sponge_gm_t*( SIN(rpi*(glamt-xs1)/rn_sponge_lx)**2 )
         END WHERE
 
      END IF
      
      !Northern sponge  >>>>>>>>>>>>>>>>>>>>>>>>
      ! The northern sponge lies above the ACC and only varies with y sinusoidally.
      !
      gphimax = MAXVAL(gphit)
      CALL mpp_max( 'usrdef_zgr', gphimax)

      IF( ln_sponge_nort_mom) THEN
         !
         IF(lwp) WRITE(numout,*) "Northern momentum sponge is on with max gamma = ", rn_sponge_gm2
         WHERE(gphiu >= rn_sponge_ly)
            psponge_gamma_u(:,:) = rn_sponge_gm2*( COS(0.5*rpi*(gphiu-gphimax)/(rn_sponge_ly-gphimax))**2 )
         END WHERE
         !
         WHERE( gphiv >= rn_sponge_ly )
            psponge_gamma_v(:,:) = rn_sponge_gm2*( COS(0.5*rpi*(gphiv-gphimax)/(rn_sponge_ly-gphimax))**2 )
         END WHERE
         !
      END IF

      IF( ln_sponge_nort_tra ) THEN
         !
         IF(lwp) WRITE(numout,*) "Northern tracer sponge is on with gamma = ", rn_sponge_gm_t2
         WHERE( gphit >= rn_sponge_ly )
            psponge_gamma_t(:,:) = rn_sponge_gm_t2*( COS(0.5*rpi*(gphit-gphimax)/(rn_sponge_ly-gphimax))**2 )
         END WHERE
         !
      END IF

      gamma_maxval = MAXVAL(ABS(sponge_gamma_u))
      CALL mpp_max('usrdef_zgr', gamma_maxval)
      IF(lwp) WRITE(numout,*) "u gamma maxval", gamma_maxval

      gamma_maxval = MAXVAL(ABS(sponge_gamma_v))
      CALL mpp_max('usrdef_zgr', gamma_maxval)
      IF(lwp) WRITE(numout,*) "u gamma maxval", gamma_maxval 

      gamma_maxval = MAXVAL(ABS(sponge_gamma_t))
      CALL mpp_max('usrdef_zgr', gamma_maxval)
      IF(lwp) WRITE(numout,*) "t gamma maxval", gamma_maxval 

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! Determine target/initial temperature and velocity fields
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !
      ptarget_uo(:,:,:) = 0._wp
      ptarget_vo(:,:,:) = 0._wp
      ptarget_to(:,:,:) = 0._wp
      ptarget_so(:,:,:) = rn_sponge_so
      !
      IF( ln_sponge_uoconst ) THEN
         ! >>>>>>
         ! CASE 1: ACC does not vary with space
         ! >>>>>>

         ptarget_uo(:,:,:) = rn_sponge_uo
         ptarget_vo(:,:,:) = rn_sponge_vo
         ptarget_to(:,:,:) = rn_sponge_to
         
      ELSE IF ( ln_sponge_uovar ) THEN
         ! >>>>>>
         ! CASE 2 : ACC varies sinusoidally with space
         ! >>>>>>
         ptarget_vo(:,:,:) = 0._wp
         ptarget_to(:,:,:) = -999.
         !
         !Calculate f0 and beta for the Coriolis parameter : f = zf0 + beta * (y-y0)
         zbeta = 1e3 * 2._wp * omega * COS( rad * rn_ppgphi0 ) / ra
         zf0   = 2._wp * omega * SIN( rad * rn_ppgphi0 )
         !

         !Loop over each model level
         DO jk = 1, jpk
            !
            pdept_jk = pdept(:,:,jk)
            ptarget_uo_jk(:,:) = 0._wp
            !
            WHERE( ((gphiu > rn_d3 ).AND.(gphiu < rn_sponge_ly).AND.(glamu >= xs1).AND.(glamu <= xs2)))
               !
               ! Channel target velocity >>>>>>>>>>>>>>
               ! In the space between the channel shelf and the northern boundary of the channel, 
               ! the velocity varies sinusoidally
               ! _________________________________ y = rn_sponge_ly
               !     ->
               !     --->
               !     ----->
               !     ---> 
               !     ->                           ____ y = d3
               !%%%%%%%%% Channel shelf %%%%%%%%%%  \/  NO FLOW OVER
               !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  /\  THE CHANNEL SHELF
               !
               ptarget_uo_jk(:,:) = rn_sponge_uomax * (SIN( rpi*(gphiu-rn_d3)/(rn_sponge_ly-rn_d3))**2)   &
                               &    * EXP(-pdept_jk / rn_depth_decay)
               !
            END WHERE

            WHERE( gphiu >= rn_sponge_ly )
               !
               ! Northern target velocity >>>>>>>>>>>>>>>>
               ! u is simply zero in this region
               !
               ptarget_uo_jk(:,:) = 0._wp
            END WHERE
            !
            !Save uo values to the 3d array
            ptarget_uo(:,:,jk) = ptarget_uo_jk(:,:)
            !

            !Calculate target temperatures >>>>>>>>>>>>>>
            WHERE( (gphit > rn_d3).AND.(gphit < rn_sponge_ly).AND.(glamt >= xs1).AND.(glamt <= xs2) )
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


                ptarget_to(:,:,jk) = + ( rn_sponge_uomax / (grav*rn_a0_user*rn_depth_decay*4*(rpi**2)) )            &
                                   & * EXP(-pdept_jk/rn_depth_decay) * (                                            &
                                   & 2*(rpi**2)*zf0*(rn_sponge_ly - gphit)                                          &
                                   & + (rpi**2)*zbeta*(rn_sponge_ly**2 - gphit**2)                                  &
                                   & + rpi*(zf0 + zbeta*gphit)*(rn_sponge_ly - rn_d3)*SIN(2*rpi*(gphit-rn_d3)/(rn_sponge_ly - rn_d3))   &
                                   & - zbeta*((rn_sponge_ly - rn_d3)**2)*(SIN(rpi*(gphit - rn_d3)/(rn_sponge_ly - rn_d3))**2) )         &
                                   & * 1e3  &  !Conversion for km --> m
                                   & + (rn_sponge_tomax-rn_sponge_tobot) * exp(-pdept_jk/rn_depth_decay)            &
                                   & + rn_sponge_tobot

            END WHERE

            WHERE( (gphit <= rn_d3).AND.(gphit >=0).AND.(glamt >= xs1).AND.(glamt <= xs2) )
            !
            ! Temperature over the channel shelf does not vary horizontally
            !
                 ptarget_to(:,:,jk) = + ( rn_sponge_uomax / (grav*rn_a0_user*rn_depth_decay*4*(rpi**2)) )           &
                                   & * EXP(-pdept_jk/rn_depth_decay) * (                                            &
                                   & 2*(rpi**2)*zf0*(rn_sponge_ly - rn_d3)                                          &
                                   & + (rpi**2)*zbeta*(rn_sponge_ly**2 - rn_d3**2) )                                &
                                   & * 1e3  &  !Conversion for km --> m
                                   & + (rn_sponge_tomax-rn_sponge_tobot) * exp(-pdept_jk/rn_depth_decay)            &
                                   & + rn_sponge_tobot
            END WHERE

            

 
            WHERE( gphit >= rn_sponge_ly )
                !
                ! Temperature over the northern sponge does not vary horizontally
                !
                ptarget_to(:,:,jk) = (rn_sponge_tomax - rn_sponge_tobot) * exp(-pdept_jk/rn_depth_decay) + rn_sponge_tobot
                !
            END WHERE
            !
         END DO     
         !
      END IF

   END SUBROUTINE sponge_setup

   SUBROUTINE H_to_ind(zH, z2d)
      !
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj) :: zH
      REAL(wp), INTENT(out), DIMENSION(jpi,jpj) :: z2d
      REAL(wp), DIMENSION(jpk) :: pdept_1d, pdepw_1d, pe3t_1d, pe3w_1d
      REAL(wp) :: zH_point
      INTEGER :: ind_min, ik1, ik2, ji, jj
      !
      z2d(:,:) = 0._wp
      !
      SELECT CASE(nn_zgr_type)
      !
      CASE(0) ! Regular z grid
         z2d(:,:) = zH(:,:) * jpkm1 / rn_domszz
      !
      CASE(1) ! Madec and Imbard z grid
         !
         CALL zgr_z( pdept_1d, pdepw_1d, pe3t_1d , pe3w_1d )                  
         !
         DO ji = 1, jpi
            DO jj =1, jpj
               !
               zH_point = zH(ji,jj)
               ind_min = MINLOC(ABS(pdept_1d - zH_point), DIM=1)
               !
               IF( zH_point - pdept_1d(ind_min) > 0) THEN
                  !
                  IF( ind_min .NE. jpk) THEN
                     ik1 = ind_min
                     ik2 = ind_min + 1 
                  ELSE
                     ik1 = ind_min - 1
                     ik2 = ind_min
                  END IF
                  !
                  z2d(ji,jj) = ik1 + ( zH_point - pdept_1d(ik1) )/( pdept_1d(ik2) - pdept_1d(ik1) )
                  !
               ELSE IF( zH_point - pdept_1d(ind_min) < 0 ) THEN
                  !
                  IF( ind_min .NE. 1 ) THEN
                    ik1 = ind_min - 1
                    ik2 = ind_min
                  ELSE
                    ik1 = ind_min
                    ik2 = ind_min + 1
                  END IF 
                  !
                  z2d(ji,jj) = ik1 + ( zH_point - pdept_1d(ik1) )/( pdept_1d(ik2) - pdept_1d(ik1) )
                  !
               ELSE !zH_point lies exactly on an existing coordinate in pdept_1d
                  !
                  z2d(ji,jj) = ind_min
                  !
               END IF           
               !
            END DO
         END DO      
         !
         WHERE( z2d > jpkm1 )
            z2d = jpkm1
         END WHERE     
         !
         WHERE(z2d < 0 )
            z2d = 0._wp
         END WHERE              
         !
      !Work in progress - Hybrid (terrain following) coordinates
      !CASE(2) ! Terrain (sigma) following coordinate
      !   !
      !   DO ji = 1, jpi
      !      DO jj = 1, jpj
      !         !
      !         zH_point = zH(ji,jj)
      !         !
      !         z2d(ji,jj) = jpkm1
      !         !
      !         IF( zH_point < 20. ) THEN
      !            z2d(ji,jj) = 0.
      !         ELSE
      !            DO WHILE( zH_point / z2d(ji,jj) < 20.)
      !               z2d(ji,jj) = z2d(ji,jj) - 1
      !            END DO
      !         END IF
      !      END DO
      !   END DO
      !    !
      !   WHERE(z2d < 0 )
      !      z2d = 0._wp
      !   END WHERE

         
      END SELECT      
      !
   END SUBROUTINE H_to_ind   


END MODULE usrdef_zgr
