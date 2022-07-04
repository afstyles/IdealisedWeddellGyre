MODULE traldf
   !!======================================================================
   !!                       ***  MODULE  traldf  ***
   !! Ocean Active tracers : lateral diffusive trends 
   !!=====================================================================
   !! History :  9.0  ! 2005-11  (G. Madec)  Original code
   !!  NEMO      3.0  ! 2008-01  (C. Ethe, G. Madec)  merge TRC-TRA 
   !!            3.7  ! 2013-12  (G. Madec) remove the optional computation from T & S anomaly profiles and traldf_bilapg
   !!             -   ! 2013-12  (F. Lemarie, G. Madec)  triad operator (Griffies) + Method of Stabilizing Correction
   !!             -   ! 2014-01  (G. Madec, S. Masson)  restructuration/simplification of lateral diffusive operators
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_ldf       : update the tracer trend with the lateral diffusion trend
   !!   tra_ldf_init  : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE ldftra         ! lateral diffusion: eddy diffusivity & EIV coeff.
   USE ldfslp         ! lateral diffusion: iso-neutral slope
   USE traldf_lap_blp ! lateral diffusion: laplacian iso-level            operator  (tra_ldf_lap/_blp   routines)
   USE traldf_iso     ! lateral diffusion: laplacian iso-neutral standard operator  (tra_ldf_iso        routine )
   USE traldf_triad   ! lateral diffusion: laplacian iso-neutral triad    operator  (tra_ldf_triad      routine )
   USE trd_oce        ! trends: ocean variables
   USE trdtra         ! ocean active tracers trends
   USE usrdef_nam, ONLY: ln_sponge_chan_tra, ln_sponge_nort_tra   
   !
   USE prtctl         ! Print control
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE timing         ! Timing
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_ldf        ! called by step.F90 
   PUBLIC   tra_ldf_init   ! called by nemogcm.F90 
   
   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: traldf.F90 10068 2018-08-28 14:09:04Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_ldf( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf  ***
      !! 
      !! ** Purpose :   compute the lateral ocean tracer physics.
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   ztrdt, ztrds
      INTEGER :: ji, jj, jk
      REAL(wp) :: zvol, ztn
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: ztn_weighted
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('tra_ldf')
      !
      IF( l_trdtra )   THEN                    !* Save ta and sa trends
         ALLOCATE( ztrdt(jpi,jpj,jpk) , ztrds(jpi,jpj,jpk) ) 
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem) 
         ztrds(:,:,:) = tsa(:,:,:,jp_sal)
      ENDIF
      !
      SELECT CASE ( nldf_tra )                 !* compute lateral mixing trend and add it to the general trend
      CASE ( np_lap   )                                  ! laplacian: iso-level operator
         CALL tra_ldf_lap  ( kt, nit000,'TRA', ahtu, ahtv, gtsu, gtsv, gtui, gtvi, tsb,      tsa, jpts,  1   )
      CASE ( np_lap_i )                                  ! laplacian: standard iso-neutral operator (Madec)
         CALL tra_ldf_iso  ( kt, nit000,'TRA', ahtu, ahtv, gtsu, gtsv, gtui, gtvi, tsb, tsb, tsa, jpts,  1   )
      CASE ( np_lap_it )                                 ! laplacian: triad iso-neutral operator (griffies)
         CALL tra_ldf_triad( kt, nit000,'TRA', ahtu, ahtv, gtsu, gtsv, gtui, gtvi, tsb, tsb, tsa, jpts,  1   )
      CASE ( np_blp , np_blp_i , np_blp_it )             ! bilaplacian: iso-level & iso-neutral operators
         CALL tra_ldf_blp  ( kt, nit000,'TRA', ahtu, ahtv, gtsu, gtsv, gtui, gtvi, tsb      , tsa, jpts, nldf_tra )
      END SELECT
      !
      ! Add temperature sponge relaxation to the ldf trend.
      !
      IF( ln_sponge_chan_tra.OR.ln_sponge_nort_tra ) THEN 
         !
         DO jk = 1, jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  !
                  tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem) - tmask(ji,jj,jk) * sponge_gamma_t(ji,jj) *   &
                                     &   ( tsb(ji,jj,jk,jp_tem) - target_to(ji,jj,jk) )
                  !
                  tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal) - tmask(ji,jj,jk) * sponge_gamma_t(ji,jj) *   &
                                     &   ( tsb(ji,jj,jk,jp_sal) - target_so(ji,jj,jk) )
               END DO
            END DO
         END DO
         !
      END IF 

      !
      ! Output sponge parameters for diagnostic purposes
      !

      IF( iom_use("sponge_gamma_t") )  CALL iom_put( "sponge_gamma_t" , sponge_gamma_t(:,:) )
      IF( iom_use("target_to") )       CALL iom_put( "target_to", target_to(:,:,:) )

      IF( iom_use("mean_temperature").OR.iom_use("mean_salinity") ) THEN
         ALLOCATE( ztn_weighted(jpi,jpj,jpk) )
         ztn_weighted(:,:,:) = 0.
        
         !Calculate T cell volume of the domain
         DO jk = 1,jpk
               ztn_weighted(:,:,jk) = e1t(:,:) * e2t(:,:) * e3t_n(:,:,jk) * tmask(:,:,jk) 
         END DO
         
         zvol =  glob_sum('traldf', ztn_weighted) ! T cell volume of the domain
         
         !Calculate the xyz mean temperature
         IF( iom_use("mean_temperature") ) THEN
            DO jk = 1,jpk
               ztn_weighted(:,:,jk) = tsb(:,:,jk,jp_tem) * e1t(:,:) * e2t(:,:) * e3t_n(:,:,jk) * tmask(:,:,jk)
            END DO
            ztn = glob_sum('traldf', ztn_weighted)/zvol
            CALL iom_put( "mean_temperature", ztn )
         END IF 

         !Calculate the xyz mean salinity
         IF( iom_use("mean_salinity") ) THEN
            DO jk = 1,jpk
               ztn_weighted(:,:,jk) = tsb(:,:,jk,jp_sal) * e1t(:,:) * e2t(:,:) * e3t_n(:,:,jk) * tmask(:,:,jk) 
            END DO
            ztn = glob_sum('traldf', ztn_weighted)/zvol 
            CALL iom_put( "mean_salinity", ztn ) 
         END IF

         DEALLOCATE( ztn_weighted )

      END IF

      IF( l_trdtra )   THEN                    !* save the horizontal diffusive trends for further diagnostics
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem) - ztrdt(:,:,:)
         ztrds(:,:,:) = tsa(:,:,:,jp_sal) - ztrds(:,:,:)
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_ldf, ztrdt )
         CALL trd_tra( kt, 'TRA', jp_sal, jptra_ldf, ztrds )
         DEALLOCATE( ztrdt, ztrds ) 
      ENDIF
      !                                        !* print mean trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' ldf  - Ta: ', mask1=tmask,               &
         &                       tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      !
      IF( ln_timing )   CALL timing_stop('tra_ldf')
      !
   END SUBROUTINE tra_ldf


   SUBROUTINE tra_ldf_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_init  ***
      !! 
      !! ** Purpose :   Choice of the operator for the lateral tracer diffusion
      !!
      !! ** Method  :   set nldf_tra from the namtra_ldf logicals
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio, ierr   ! temporary integers 
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                     !==  Namelist print  ==!
         WRITE(numout,*)
         WRITE(numout,*) 'tra_ldf_init : lateral tracer diffusive operator'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtra_ldf: already read in ldftra module'
         WRITE(numout,*) '      see ldf_tra_init report for lateral mixing parameters'
         WRITE(numout,*)
         !
         SELECT CASE( nldf_tra )             ! print the choice of operator
         CASE( np_no_ldf )   ;   WRITE(numout,*) '   ==>>>   NO lateral diffusion'
         CASE( np_lap    )   ;   WRITE(numout,*) '   ==>>>   laplacian iso-level operator'
         CASE( np_lap_i  )   ;   WRITE(numout,*) '   ==>>>   Rotated laplacian operator (standard)'
         CASE( np_lap_it )   ;   WRITE(numout,*) '   ==>>>   Rotated laplacian operator (triad)'
         CASE( np_blp    )   ;   WRITE(numout,*) '   ==>>>   bilaplacian iso-level operator'
         CASE( np_blp_i  )   ;   WRITE(numout,*) '   ==>>>   Rotated bilaplacian operator (standard)'
         CASE( np_blp_it )   ;   WRITE(numout,*) '   ==>>>   Rotated bilaplacian operator (triad)'
         END SELECT
      ENDIF
      !
   END SUBROUTINE tra_ldf_init

   !!======================================================================
END MODULE traldf
