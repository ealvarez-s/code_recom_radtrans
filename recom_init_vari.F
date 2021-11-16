C $Header: /csys/software/cvsroot/REcoM/MITgcm/recom/recom_init_vari.F,v 1.11 2008/04/17 14:25:52 mlosch Exp $
C $Name:  $

#include "CPP_OPTIONS.h"
C#include "PTRACERS_OPTIONS.h"
#include "GCHEM_OPTIONS.h"
#include "RECOM_OPTIONS.h"
#ifdef ALLOW_EXF
# include "EXF_OPTIONS.h"
#endif

CBOP
C !ROUTINE: RECOM_INIT_VARI

C !INTERFACE: ============================================================
      SUBROUTINE RECOM_INIT_VARI( myThid )

C !DESCRIPTION:
C     Initialize RECOM data structures that could potentially be subject
C     to adjoint computations (to be called from S/R GCHEM_INIT_VARI)

C !USES: ================================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#ifdef ALLOW_MNC
#include "MNC_PARAMS.h"
#endif
#ifdef ALLOW_PTRACERS
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "RECOM.h"
#endif /* ALLOW_PTRACERS */
#ifdef ALLOW_EXF
# include "EXF_PARAM.h"
# include "RECOM_EXF.h"
#endif /* ALLOW_EXF */
CTW increased sinking speed
#ifdef ALLOW_SINK_INCREASE
#include "RECOM_PARAMS.h"
#endif /* ALLOW_SINK_INCREASE */

#ifdef RECOM_WAVEBANDS
#include "SPECTRAL_SIZE.h" /* tlam */
#include "SPECTRAL.h"      /* oasim_ed, oasim_es */
#include "WAVEBANDS_PARAMS.h"  /* QYmax  QYmax_d */  
#endif      

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none

#ifdef ALLOW_RECOM

C !LOCAL VARIABLES: ====================================================
C     i,j,k,bi,bj - Loop counters

      INTEGER i, j, k, bi, bj
      LOGICAL useVariableK
#ifdef ALLOW_EXF
      INTEGER interp_method
      PARAMETER (interp_method=1)
      _RL loc0, loc1
#endif /* ALLOW_EXF */
#ifdef RECOM_WAVEBANDS
      INTEGER ilam
      INTEGER nl
      _RL cu_area
      _RL cu_area_phy
      _RL cu_area_dia      
      _RL mQY
      _RL mQY_dia
#endif      
CEOP

C
C     initialise specific forcing fields
C      
      K = 1
      DO bj=myByLo(myThid),myByHi(myThid)
       DO bi=myBxLo(myThid),myBxHi(myThid)
        DO J=1-Oly,sNy+Oly
         DO I=1-Olx,sNx+Olx
          dicFlux  (i,j,bi,bj)      = 0. _d 0
          windSpeed(i,j,bi,bj)      = 0. _d 0
          atmospCO2(i,j,bi,bj)      = 0. _d 0
          atmospCO2(i,j,bi,bj)      = 380. _d 0 * maskC(I,J,K,bi,bj)
          feDustDeposit(i,j,bi,bj)  = 0. _d 0
C     initialise hPlus with some sensible values here
          hPlus(i,j,bi,bj)          = 1. _d -08 * maskC(I,J,K,bi,bj)
C     
          benthicLayerC(i,j,bi,bj)  = 0. _d 0
          benthicLayerN(i,j,bi,bj)  = 0. _d 0
CCV#ifdef ALLOW_RECOM_SILICATE
          benthicLayerSi(i,j,bi,bj) = 0. _d 0
CCV#endif /* ALLOW_RECOM_SILICATE */
          benthicLayerCaCO3(i,j,bi,bj) = 0. _d 0
         ENDDO
        ENDDO
C
C
#if defined(RECOM_WAVEBANDS) && defined(OASIM)
CEA initialize OASIM(x,y) forcing fields
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
             tmp01ed(i,j,bi,bj) = 0. _d 0
             tmp02ed(i,j,bi,bj) = 0. _d 0
             tmp03ed(i,j,bi,bj) = 0. _d 0
             tmp04ed(i,j,bi,bj) = 0. _d 0
             tmp05ed(i,j,bi,bj) = 0. _d 0
             tmp06ed(i,j,bi,bj) = 0. _d 0
             tmp07ed(i,j,bi,bj) = 0. _d 0
             tmp08ed(i,j,bi,bj) = 0. _d 0
             tmp09ed(i,j,bi,bj) = 0. _d 0
             tmp10ed(i,j,bi,bj) = 0. _d 0
             tmp11ed(i,j,bi,bj) = 0. _d 0
             tmp12ed(i,j,bi,bj) = 0. _d 0
             tmp13ed(i,j,bi,bj) = 0. _d 0
             tmp01es(i,j,bi,bj) = 0. _d 0
             tmp02es(i,j,bi,bj) = 0. _d 0
             tmp03es(i,j,bi,bj) = 0. _d 0
             tmp04es(i,j,bi,bj) = 0. _d 0
             tmp05es(i,j,bi,bj) = 0. _d 0
             tmp06es(i,j,bi,bj) = 0. _d 0
             tmp07es(i,j,bi,bj) = 0. _d 0
             tmp08es(i,j,bi,bj) = 0. _d 0
             tmp09es(i,j,bi,bj) = 0. _d 0
             tmp10es(i,j,bi,bj) = 0. _d 0
             tmp11es(i,j,bi,bj) = 0. _d 0
             tmp12es(i,j,bi,bj) = 0. _d 0
             tmp13es(i,j,bi,bj) = 0. _d 0
            ENDDO !SL
           ENDDO  !SL
C#ifndef ALLOW_EXF
           DO ilam=1,tlam
            DO j=1-Oly,sNy+Oly  !SL
             DO i=1-Olx,sNx+Olx  !SL 
               oasim_ed(i,j,ilam,bi,bj) = 0. _d 0
C               oasim_ed1(i,j,ilam,bi,bj) = 0. _d 0
               oasim_es(i,j,ilam,bi,bj) = 0. _d 0
C               oasim_es1(i,j,ilam,bi,bj) = 0. _d 0               
             ENDDO
            ENDDO
           ENDDO
C#endif
#endif /* OASIM */
C
C     initialise total matter fields (diagnostic)
C
        DO K=1, Nr
         DO J=1,sNy
          DO I=1,sNx
           recomC(i,j,k,bi,bj)  = 0. _d 0
           recomN(i,j,k,bi,bj)  = 0. _d 0
           recomSi(i,j,k,bi,bj) = 0. _d 0
          ENDDO
         ENDDO
        ENDDO
       ENDDO  !bi
      ENDDO   !bj
C
CTW  calculate increased sinking speed
C
#ifdef ALLOW_SINK_INCREASE
      DO k=1,Nr
       Vdetfast(k)= Vdet_a*abs(rF(k))+Vdet_b
       print *, 'layer & vdetfast',k,Vdetfast(k)
      ENDDO
#endif /* ALLOW_SINK_INCREASE */

      
#ifdef ALLOW_EXF
      loc0 = 0. _d 0
      loc1 = 1. _d 0
      IF ( recom_windFile .NE. ' '  ) THEN
       call exf_init_gen (
     &     recom_windfile, wspeedperiod, loc1, wspeedmask,
     &     wspeedconst, windSpeed, windSpeed0, windSpeed1, 
#ifdef USE_EXF_INTERPOLATION
     &     wspeed_lon0, wspeed_lon_inc, 
     &     wspeed_lat0, wspeed_lat_inc,
     &     wspeed_nlon, wspeed_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF
      IF ( recom_pCO2File .NE. ' ' ) THEN
       call exf_init_gen (
     &     recom_pCO2file, pCO2period, loc1, sfluxmask,
     &     loc0, atmospCO2, atmosp0, atmosp1, 
#ifdef USE_EXF_INTERPOLATION
     &     pCO2_lon0, pCO2_lon_inc, 
     &     pCO2_lat0, pCO2_lat_inc,
     &     pCO2_nlon, pCO2_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF
CCV#ifdef ALLOW_RECOM_FEINPUT
      IF ( recom_ironFile .NE. ' ' ) THEN
       call exf_init_gen (
     &     recom_ironFile, fedustperiod, loc1, sfluxmask,
     &     loc0, fedustDeposit, fedustDeposit0, fedustDeposit1, 
#ifdef USE_EXF_INTERPOLATION
     &     fedust_lon0, fedust_lon_inc, 
     &     fedust_lat0, fedust_lat_inc,
     &     fedust_nlon, fedust_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF
CCV#endif /* ALLOW_RECOM_FEINPUT */

#if defined(RECOM_WAVEBANDS) && defined(OASIM)
CEA ----------- Direct ----------------------
CEA Ed01 
      IF ( darwin_oasim_edFile01 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_edFile01, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp01ed, tmp01ed0, tmp01ed1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Ed02 
      IF ( darwin_oasim_edFile02 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_edFile02, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp02ed, tmp02ed0, tmp02ed1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Ed03
      IF ( darwin_oasim_edFile03 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_edFile03, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp03ed, tmp03ed0, tmp03ed1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Ed04
      IF ( darwin_oasim_edFile04 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_edFile04, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp04ed, tmp04ed0, tmp04ed1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Ed05
      IF ( darwin_oasim_edFile05 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_edFile05, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp05ed, tmp05ed0, tmp05ed1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Ed06
      IF ( darwin_oasim_edFile06 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_edFile06, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp06ed, tmp06ed0, tmp06ed1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Ed07
      IF ( darwin_oasim_edFile07 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_edFile07, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp07ed, tmp07ed0, tmp07ed1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Ed08
      IF ( darwin_oasim_edFile08 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_edFile08, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp08ed, tmp08ed0, tmp08ed1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Ed09
      IF ( darwin_oasim_edFile09 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_edFile09, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp09ed, tmp09ed0, tmp09ed1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Ed10
      IF ( darwin_oasim_edFile10 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_edFile10, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp10ed, tmp10ed0, tmp10ed1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Ed11
      IF ( darwin_oasim_edFile11 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_edFile11, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp11ed, tmp11ed0, tmp11ed1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Ed12
      IF ( darwin_oasim_edFile12 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_edFile12, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp12ed, tmp12ed0, tmp12ed1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Ed13
      IF ( darwin_oasim_edFile13 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_edFile13, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp13ed, tmp13ed0, tmp13ed1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA ----------  end direct -------------------
C
CEA ----------- Diffuse ----------------------
CEA Es01 
      IF ( darwin_oasim_esFile01 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_esFile01, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp01es, tmp01es0, tmp01es1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Es02 
      IF ( darwin_oasim_esFile02 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_esFile02, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp02es, tmp02es0, tmp02es1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Es03
      IF ( darwin_oasim_esFile03 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_esFile03, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp03es, tmp03es0, tmp03es1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Es04
      IF ( darwin_oasim_esFile04 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_esFile04, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp04es, tmp04es0, tmp04es1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Es05
      IF ( darwin_oasim_esFile05 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_esFile05, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp05es, tmp05es0, tmp05es1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Es06
      IF ( darwin_oasim_esFile06 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_esFile06, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp06es, tmp06es0, tmp06es1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Es07
      IF ( darwin_oasim_esFile07 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_esFile07, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp07es, tmp07es0, tmp07es1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Es08
      IF ( darwin_oasim_esFile08 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_esFile08, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp08es, tmp08es0, tmp08es1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Es09
      IF ( darwin_oasim_esFile09 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_esFile09, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp09es, tmp09es0, tmp09es1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Ed10
      IF ( darwin_oasim_esFile10 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_esFile10, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp10es, tmp10es0, tmp10es1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Es11
      IF ( darwin_oasim_esFile11 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_esFile11, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp11es, tmp11es0, tmp11es1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Es12
      IF ( darwin_oasim_esFile12 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_esFile12, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp12es, tmp12es0, tmp12es1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA Es13
      IF ( darwin_oasim_esFile13 .NE. ' '  ) THEN
       call exf_init_gen ( 
     &     darwin_oasim_esFile13, OASIMperiod, loc1, sfluxmask,
     &     loc0, tmp13es, tmp13es0, tmp13es1,
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc,
     &     OASIM_lat0, OASIM_lat_inc,
     &     OASIM_nlon, OASIM_nlat, xC, yC, interp_method,
#endif
     &     mythid )
      ENDIF   
CEA ----------  end diffuse -------------------
#endif /* OASIM */     

#else /* ALLOW_EXF */      
      _BEGIN_MASTER(myThid)
      IF ( recom_WindFile .NE. ' '  ) THEN
       CALL READ_FLD_XY_RS( recom_WindFile, ' ',windSpeed, 
     &      0, myThid )
      ENDIF
      IF ( recom_pCO2File .NE. ' '  ) THEN
       CALL READ_FLD_XY_RS( recom_pCO2File, ' ',atmospCO2,
     &      0, myThid )
      ENDIF
      IF ( recom_ironFile .NE. ' '  ) THEN
       CALL READ_FLD_XY_RS( recom_ironFile, ' ',feDustDeposit,
     &      0, myThid )
      ENDIF
#if defined(RECOM_WAVEBANDS) && defined(OASIM) 
      IF ( darwin_oasim_edFile01 .NE. ' '  ) THEN
       CALL READ_FLD_XY_RS( darwin_oasim_edFile01, ' ', tmp01ed,
     &      0, myThid )
      ENDIF
CC needs repetition or change the function to 
#endif /* OASIM */       
      _END_MASTER(myThid)
#endif /* ALLOW_EXF */
C
C
      _EXCH_XY_RS(windSpeed    , myThid )
      _EXCH_XY_RS(atmospCO2    , myThid )
      _EXCH_XY_RS(feDustDeposit, myThid )
      _EXCH_XY_RS(hPlus        , myThid )
#if defined(RECOM_WAVEBANDS) && defined(OASIM) 
      _EXCH_XY_RL(tmp01ed, myThid )
      _EXCH_XY_RL(tmp02ed, myThid )
      _EXCH_XY_RL(tmp03ed, myThid )
      _EXCH_XY_RL(tmp04ed, myThid )
      _EXCH_XY_RL(tmp05ed, myThid )
      _EXCH_XY_RL(tmp06ed, myThid )
      _EXCH_XY_RL(tmp07ed, myThid )
      _EXCH_XY_RL(tmp08ed, myThid )
      _EXCH_XY_RL(tmp09ed, myThid )
      _EXCH_XY_RL(tmp10ed, myThid )
      _EXCH_XY_RL(tmp11ed, myThid )
      _EXCH_XY_RL(tmp12ed, myThid )
      _EXCH_XY_RL(tmp13ed, myThid )
      _EXCH_XY_RL(tmp01es, myThid )
      _EXCH_XY_RL(tmp02es, myThid )
      _EXCH_XY_RL(tmp03es, myThid )
      _EXCH_XY_RL(tmp04es, myThid )
      _EXCH_XY_RL(tmp05es, myThid )
      _EXCH_XY_RL(tmp06es, myThid )
      _EXCH_XY_RL(tmp07es, myThid )
      _EXCH_XY_RL(tmp08es, myThid )
      _EXCH_XY_RL(tmp09es, myThid )
      _EXCH_XY_RL(tmp10es, myThid )
      _EXCH_XY_RL(tmp11es, myThid )
      _EXCH_XY_RL(tmp12es, myThid )
      _EXCH_XY_RL(tmp13es, myThid )
#endif /* OASIM */
C      
#ifdef RECOM_IRON_HYDROTHERMAL
C     if a file for hydrothermal iron input is defined, read it in; else
C     initialize it to zeros
      _BEGIN_MASTER(myThid)
      IF ( recom_hydroFile .NE. ' '  ) THEN
         CALL READ_FLD_XYZ_RL(recom_hydroFile,' ',
     &        feHydrothermal(1-OLx,1-OLy,1,1,1),0,myThid)
      ENDIF
      _END_MASTER(myThid)

C shifted the exchange routine outside of master thread!
      _EXCH_XYZ_RL(feHydrothermal , myThid)

C     Apply mask
      DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
            DO k=1,Nr
               DO j=1-OLy,sNy+OLy
                  DO i=1-OLx,sNx+OLx
                     IF (maskC(i,j,k,bi,bj).EQ.0.)
     &           feHydrothermal(i,j,k,bi,bj)=0. _d 0
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
#endif /* RECOM_IRON_HYDROTHERMAL */					 

C     This part can go to recom_sms to take the value of Qyield,
C     for now, we gave constant values to mQY and mQY_dia
#ifdef RECOM_WAVEBANDS
C      mQY = QYmax / 4.596        ! mmolC uE-1
C      mQY_dia = QYmax_d / 4.596  ! mmolC uE-1     
      mQY = 0.
      mQY_dia = 0.      
           do nl = 1,tlam          
              mQY = mQY + ((QYmax/tlam) / WtouEins(nl))
              mQY_dia = mQY_dia + ((QYmax_d/tlam) / WtouEins(nl))
           end do
      mQY = max(tiny, mQY)
      mQY_dia = max(tiny,mQY_dia)
C     
      call WAVEBANDS_INIT_VARI(MyThid)
          do nl = 1,tlam
             alphachl_nl(nl) =  mQY * aphy_chl_ps(nl)
             alphachl_nl_dia(nl) = mQY_dia * aphy_chl_ps_dia(nl)      
          end do
c find mean 
          cu_area = 0.d0
          do nl = 1,tlam
             cu_area_phy = cu_area + wb_width(nl) * alphachl_nl(nl)
             cu_area_dia = cu_area + wb_width(nl) * alphachl_nl_dia(nl)
          end do
          alpha_mean = cu_area_phy / wb_totalWidth
          alpha_mean_dia = cu_area_dia / wb_totalWidth
          
C          chl2cmin=chl2cmax/
C     &            (1+(chl2cmax(np)* alpha_mean(np) *2000. _d 0)/
C     &               (2*pcmax(np)))      
#endif
            
      
C Read benthic layer variables from a pickup file if needed
C (analogous to ptracers_read_checkpoint)
      IF (nIter0.GT.PTRACERS_Iter0) THEN
       CALL RECOM_READ_PICKUP( nIter0,myThid )
      ENDIF
#endif /* ALLOW_RECOM */

      RETURN
      END
