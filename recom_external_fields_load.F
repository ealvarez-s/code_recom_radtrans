C $Header: /csys/software/cvsroot/REcoM/MITgcm/recom/recom_external_fields_load.F,v 1.14 2008/04/17 14:25:52 mlosch Exp $
C $Name:  $

#include "GCHEM_OPTIONS.h"
#include "RECOM_OPTIONS.h"
#ifdef ALLOW_EXF
# include "EXF_OPTIONS.h"
#endif /* ALLOW_EXF */

CBOP
C     !ROUTINE: RECOM_EXTERNAL_FIELDS_LOAD
C     !INTERFACE:
      SUBROUTINE RECOM_EXTERNAL_FIELDS_LOAD( 
     ImyTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==================================================================
C     | SUBROUTINE RECOM_EXTERNAL_FIELDS_LOAD
C     | o Control reading of fields from external source 
C     |   for a bio-geo-chemical model
C     *==================================================================
C     | External source field loading routine.                    
C     | This routine is called every time we want to              
C     | load a a set of external fields. The routine decides      
C     | which fields to load and then reads them in.              
C     | This routine needs to be customised for particular        
C     | experiments.                                              
C     |
C     | RECOM_EXTERNAL_FIELDS_LOAD is called from GCHEM_FORCING, before 
C     | calling the dynamics of the biogeochemical model.
C     |
C     *==================================================================
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "GCHEM.h"
#ifdef ALLOW_RECOM
# include "RECOM.h"
#endif /* ALLOW_RECOM */
#ifdef ALLOW_EXF
# include "EXF_PARAM.h"
# include "EXF_FIELDS.h"
# include "EXF_CONSTANTS.h"
# include "RECOM_EXF.h"
#endif /* ALLOW_EXF */
#ifdef RECOM_ATMOSPCO2_HISTORY
# include "RECOM_PARAMS.h"
#endif /*RECOM_ATMOSPCO2_HISTORY*/
#if defined(RECOM_WAVEBANDS) && defined(OASIM)
#include "SPECTRAL_SIZE.h"       /* tlam */
#include "SPECTRAL.h"            /* oasim_ed, oasim_es*/
C#include "WAVEBANDS_PARAMS.h"        
#endif /* WB & OASIM */      

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     note: under the multi-threaded model myiter and 
C           mytime are local variables passed around as routine 
C           arguments. Although this is fiddly it saves the need to 
C           impose additional synchronisation points when they are 
C           updated.
C     myiter - iteration counter for this thread
C     mytime - time counter for this thread
C     mythid - thread number for this instance of the routine.
      integer mythid 
      integer myiter
      _RL     myTime

#ifdef ALLOW_RECOM

C     !LOCAL VARIABLES:
C     == Local variables ==
C     aWght, bWght :: Interpolation weights
      INTEGER bi,bj,i,j,intime0,intime1
      _RL aWght,bWght,rdt
      _RL tmp1Wght, tmp2Wght
      _RL tmpRepeatPeriod
#ifdef RECOM_ATMOSPCO2_HISTORY
      INTEGER ntim, ioUnit
      _RL tmp
#endif
      INTEGER nForcingPeriods,Imytm,Ifprd,Ifcyc,Iftm
      CHARACTER*(MAX_LEN_MBUF) msgBuf

#if defined(RECOM_WAVEBANDS) && defined(OASIM)
C      COMMON /SPECTRAL_INPUT_AUX/ tmp0, tmp1, tmp2
C#ifndef ALLOW_EXF
C     &  , oasim_ed0, oasim_ed1 , oasim_es0,  oasim_es1, 
C      _RS oasim_ed0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,tlam,nSx, nSy)
C      _RS oasim_ed1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,tlam,nSx, nSy)
C      _RS oasim_es0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,tlam,nSx, nSy)
C      _RS oasim_es1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,tlam,nSx, nSy)
C#endif
C      _RS tmp0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
C      _RS tmp1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
C      _RS tmp2(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      INTEGER ilam, fp, nj0, nj1
#endif

      
#ifdef ALLOW_EXF
      INTEGER interp_method
      PARAMETER(interp_method=1)
      _RL loc0, loc1    
#endif /* ALLOW_EXF */

CEOP

      
#if defined(RECOM_WAVEBANDS) && defined(OASIM)
CEA initialize OASIM (x,y,tlam) forcing fields
C        DO bj = myByLo(myThid), myByHi(myThid)
C         DO bi = myBxLo(myThid), myBxHi(myThid)
C           DO j=1-Oly,sNy+Oly
C            DO i=1-Olx,sNx+Olx
C             tmp0(i,j,bi,bj) = 0. _d 0
C             tmp1(i,j,bi,bj) = 0. _d 0
C             tmp2(i,j,bi,bj) = 0. _d 0
C            ENDDO !SL
C           ENDDO  !SL
C#ifndef ALLOW_EXF
C           DO ilam=1,tlam
C            DO j=1-Oly,sNy+Oly  !SL
C             DO i=1-Olx,sNx+Olx  !SL
C               oasim_ed0(i,j,ilam,bi,bj) = 0. _d 0
C               oasim_ed1(i,j,ilam,bi,bj) = 0. _d 0
C               oasim_es0(i,j,ilam,bi,bj) = 0. _d 0
C               oasim_es1(i,j,ilam,bi,bj) = 0. _d 0
C             ENDDO
C            ENDDO
C           ENDDO
C#endif
C          ENDDO
C         ENDDO
#endif /* OASIM */

         
#ifdef ALLOW_EXF
      loc0 = 0. _d 0
      loc1 = 1. _d 0

C-    Wind speed
C      CALL EXF_SET_FLD(
C     I     'wspeed', wspeedfile, wspeedmask,
C     I     wspeedStartTime, wspeedperiod, wspeedRepCycle,
C     I     exf_inscal_wspeed,
C     I     wspeed_exfremo_intercept, wspeed_exfremo_slope,
C     U     wspeed, wspeed0, wspeed1,
C#ifdef USE_EXF_INTERPOLATION
C     I     wspeed_lon0, wspeed_lon_inc,
C     I     wspeed_lat0, wspeed_lat_inc,
C     I     wspeed_nlon, wspeed_nlat, xC, yC, wspeed_interpMethod,
C#endif
C     I     myTime, myIter, myThid )

C     10-m winds
      IF (recom_windFile .NE. ' ') CALL EXF_SET_GEN  ( 
     &     recom_windFile, wspeedstartdate, wspeedperiod, 
     &     loc1, loc0, loc0, 
     &     windSpeed, windSpeed0, windSpeed1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     wspeed_lon0, wspeed_lon_inc, 
     &     wspeed_lat0, wspeed_lat_inc,
     &     wspeed_nlon, wspeed_nlat, xC, yC, interp_method,
#endif
     &     myTime, myIter, myThid )
C      ENDIF
#ifndef RECOM_ATMOSPCO2_HISTORY
C     atmospheric CO2
      IF (recom_pCO2File .NE. ' ') THEN
C     this is a hack to allow non-periodic atmospheric pCO2
C     with periodic (climatological) forcing: save repeatperiod
C     and restore it to original value after the call of exf_set_gen
       tmpRepeatPeriod = repeatPeriod
       repeatPeriod    = pCO2RepeatPeriod
       CALL EXF_SET_GEN  ( 
     &     recom_pCO2File, pCO2startdate, pCO2period, 
     &     loc1, loc0, loc0, 
     &     atmospCO2, atmosp0, atmosp1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     pCO2_lon0, pCO2_lon_inc, pCO2_lat0,
     &     pCO2_lat_inc, pCO2_nlon, pCO2_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
C     restore repeatPeriod to original value
       repeatPeriod    = tmpRepeatPeriod
      ENDIF
#endif
CCV#ifdef ALLOW_RECOM_FEINPUT
C     iron input (through dust)
      IF (recom_ironFile .NE. ' ') THEN
CMLC     this is a hack to allow non-periodic dust input
CMLC     with periodic (climatological) forcing: save repeatperiod
CMLC     and restore it to original value after the call of exf_set_gen
CML       tmpRepeatPeriod = repeatPeriod
CML       repeatPeriod    = feDustRepeatPeriod
       CALL EXF_SET_GEN  ( 
     &     recom_ironFile, feduststartdate, fedustperiod, 
     &     loc1, loc0, loc0, 
     &     feDustDeposit, feDustDeposit0, feDustDeposit1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     fedust_lon0, fedust_lon_inc, fedust_lat0,
     &     fedust_lat_inc, fedust_nlon, fedust_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
CMLC     restore repeatPeriod to original value
CML       repeatPeriod    = tmpRepeatPeriod
      ENDIF
CCV#endif /* ALLOW_RECOM_FEINPUT  */


#if defined(RECOM_WAVEBANDS) && defined(OASIM)  
CEA ---------------------- DIRECT --------------------------
CEA 1      
      IF ( darwin_oasim_edFile01 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_edFile01, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp01ed, tmp01ed0, tmp01ed1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )       
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
             DO i=1-Olx,sNx+Olx
       oasim_ed(i,j,1,bi,bj) = tmp01ed(i,j,bi,bj) 
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 2      
      IF ( darwin_oasim_edFile02 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_edFile02, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp02ed, tmp02ed0, tmp02ed1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_ed(i,j,2,bi,bj) = tmp02ed(i,j,bi,bj) 
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 3      
      IF ( darwin_oasim_edFile03 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_edFile03, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp03ed, tmp03ed0, tmp03ed1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_ed(i,j,3,bi,bj) = tmp03ed(i,j,bi,bj) 
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 4      
      IF ( darwin_oasim_edFile04 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_edFile04, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp04ed, tmp04ed0, tmp04ed1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_ed(i,j,4,bi,bj) = tmp04ed(i,j,bi,bj) 
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 5      
      IF ( darwin_oasim_edFile05 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_edFile05, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp05ed, tmp05ed0, tmp05ed1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_ed(i,j,5,bi,bj) = tmp05ed(i,j,bi,bj) 
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 6      
      IF ( darwin_oasim_edFile06 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_edFile06, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp06ed, tmp06ed0, tmp06ed1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_ed(i,j,6,bi,bj) = tmp06ed(i,j,bi,bj) 
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 7      
      IF ( darwin_oasim_edFile07 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_edFile07, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp07ed, tmp07ed0, tmp07ed1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_ed(i,j,7,bi,bj) = tmp07ed(i,j,bi,bj) 
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 8      
      IF ( darwin_oasim_edFile08 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_edFile08, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp08ed, tmp08ed0, tmp08ed1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_ed(i,j,8,bi,bj) = tmp08ed(i,j,bi,bj) 
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 9      
      IF ( darwin_oasim_edFile09 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_edFile09, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp09ed, tmp09ed0, tmp09ed1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_ed(i,j,9,bi,bj) = tmp09ed(i,j,bi,bj) 
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 10      
      IF ( darwin_oasim_edFile10 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_edFile10, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp10ed, tmp10ed0, tmp10ed1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_ed(i,j,10,bi,bj) = tmp10ed(i,j,bi,bj) 
             ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 11      
      IF ( darwin_oasim_edFile11 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_edFile11, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp11ed, tmp11ed0, tmp11ed1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_ed(i,j,11,bi,bj) = tmp11ed(i,j,bi,bj) 
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 12      
      IF ( darwin_oasim_edFile12 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_edFile12, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp12ed, tmp12ed0, tmp12ed1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_ed(i,j,12,bi,bj) = tmp12ed(i,j,bi,bj) 
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 13      
      IF ( darwin_oasim_edFile13 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_edFile13, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp13ed, tmp13ed0, tmp13ed1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_ed(i,j,13,bi,bj) = tmp13ed(i,j,bi,bj) 
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF

C       CALL EXCH_3D_RL (oasim_ed, tlam, myThid)

C --------------------- end direct ------------------------------      
CEA -------------------- DIFFUSE --------------------------------      
CEA 1
      IF ( darwin_oasim_esFile01 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_esFile01, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp01es, tmp01es0, tmp01es1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_es(i,j,1,bi,bj) = tmp01es(i,j,bi,bj)   
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 2
      IF ( darwin_oasim_esFile02 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_esFile02, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp02es, tmp02es0, tmp02es1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_es(i,j,2,bi,bj) = tmp02es(i,j,bi,bj)   
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 3
      IF ( darwin_oasim_esFile03 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_esFile03, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp03es, tmp03es0, tmp03es1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_es(i,j,3,bi,bj) = tmp03es(i,j,bi,bj)   
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 4
      IF ( darwin_oasim_esFile04 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_esFile04, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp04es, tmp04es0, tmp04es1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_es(i,j,4,bi,bj) = tmp04es(i,j,bi,bj)   
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 5
      IF ( darwin_oasim_esFile05 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_esFile05, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp05es, tmp05es0, tmp05es1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_es(i,j,5,bi,bj) = tmp05es(i,j,bi,bj)   
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 6
      IF ( darwin_oasim_esFile06 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_esFile06, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp06es, tmp06es0, tmp06es1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_es(i,j,6,bi,bj) = tmp06es(i,j,bi,bj)   
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 7
      IF ( darwin_oasim_esFile07 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_esFile07, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp07es, tmp07es0, tmp07es1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_es(i,j,7,bi,bj) = tmp07es(i,j,bi,bj)   
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 8
      IF ( darwin_oasim_esFile08 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_esFile08, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp08es, tmp08es0, tmp08es1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_es(i,j,8,bi,bj) = tmp08es(i,j,bi,bj)   
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 9
      IF ( darwin_oasim_esFile09 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_esFile09, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp09es, tmp09es0, tmp09es1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_es(i,j,9,bi,bj) = tmp09es(i,j,bi,bj) 
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 10
      IF ( darwin_oasim_esFile10 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_esFile10, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp10es, tmp10es0, tmp10es1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_es(i,j,10,bi,bj) = tmp10es(i,j,bi,bj)   
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 11
      IF ( darwin_oasim_esFile11 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_esFile11, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp11es, tmp11es0, tmp11es1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_es(i,j,11,bi,bj) = tmp11es(i,j,bi,bj)   
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 12
      IF ( darwin_oasim_esFile12 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_esFile12, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp12es, tmp12es0, tmp12es1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_es(i,j,12,bi,bj) = tmp12es(i,j,bi,bj)   
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 13
      IF ( darwin_oasim_esFile13 .NE. ' '  ) THEN
       CALL EXF_SET_GEN  ( 
     &     darwin_oasim_esFile13, OASIMstartdate, OASIMperiod, 
     &     loc1, loc0, loc0, tmp13es, tmp13es0, tmp13es1, sfluxmask, 
#ifdef USE_EXF_INTERPOLATION
     &     OASIM_lon0, OASIM_lon_inc, OASIM_lat0,
     &     OASIM_lat_inc, OASIM_nlon, OASIM_nlat, xC, yC,
     &     interp_method,
#endif
     &     myTime, myIter, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
       oasim_es(i,j,13,bi,bj) = tmp13es(i,j,bi,bj)   
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF

C       CALL EXCH_3D_RL (oasim_es, tlam, myThid)

C ---------------------- end diffuse -----------------------------      
#endif /* RECOM_WAVEBANDS && OASIM */


      
#else /* not ALLOW_EXF */

      IF ( periodicExternalForcing ) THEN

C First call requires that we initialize everything to zero for safety
cQQQ need to check timing
       IF ( myIter .EQ. nIter0 ) THEN
         CALL LEF_ZERO( windSpeed0,myThid )
         CALL LEF_ZERO( windSpeed1,myThid )
#ifndef RECOM_ATMOSPCO2_HISTORY
         CALL LEF_ZERO( atmosp0,myThid )
         CALL LEF_ZERO( atmosp1,myThid )
#endif
CCV#ifdef ALLOW_RECOM_FEINPUT 
         CALL LEF_ZERO( feDustDeposit0,myThid )
         CALL LEF_ZERO( feDustDeposit1,myThid )
CCV#endif /* ALLOW_RECOM_FEINPUT */
#ifdef OASIM
         CALL LEF_ZERO( oasim_ed0 ,myThid )
         CALL LEF_ZERO( oasim_ed1 ,myThid )
         CALL LEF_ZERO( oasim_es0 ,myThid )
         CALL LEF_ZERO( oasim_es1 ,myThid )
         CALL LEF_ZERO( tmp1 ,myThid )
         CALL LEF_ZERO( tmp2 ,myThid )
#endif /* OASIM */
       ENDIF

       rdt = 1. _d 0 / deltaTclock
       nForcingPeriods = NINT(externForcingCycle/externForcingPeriod)
       Imytm = NINT(myTime*rdt)
       Ifprd = NINT(externForcingPeriod*rdt)
       Ifcyc = NINT(externForcingCycle*rdt)
       Iftm  = MOD( Imytm+Ifcyc-Ifprd/2, Ifcyc)
C Now calculate whether it is time to update the forcing arrays
       intime0 = 1 + INT(Iftm/Ifprd)
       intime1 = 1 + MOD(intime0,nForcingPeriods)
C-jmc: with some option of g77, FLOAT results in real*4 evaluation
C      of aWght; using DFLOAT always force real*8 computation:
C-ph: however, TAF doesnt recognize DFLOAT, 
C     so I put it back to FLOAT for now 
c      aWght =  FLOAT( Iftm-Ifprd*(intime0 - 1) ) /  FLOAT( Ifprd )
c      aWght = DFLOAT( Iftm-Ifprd*(intime0 - 1) ) / DFLOAT( Ifprd )
CC-jmc: so let's try this:
       tmp1Wght = FLOAT( Iftm-Ifprd*(intime0 - 1) )
       tmp2Wght = FLOAT( Ifprd )
       aWght =  tmp1Wght / tmp2Wght
       bWght = 1. _d 0 - aWght
       
       IF (
     &   Iftm-Ifprd*(intime0-1).EQ. 0
     &   .OR. myIter .EQ. nIter0
     &    ) THEN

        _BEGIN_MASTER(myThid)

C      If the above condition is met then we need to read in
C      data for the period ahead and the period behind myTime.
        WRITE(msgBuf,'(A,E12.6,I11)') 
     &       'S/R RECOM_EXTERNAL_FIELDS_LOAD: Reading new REcoM data ',
     &       myTime,myIter
        CALL PRINT_MESSAGE( msgBuf, standardMessageUnit, 
     &       SQUEEZE_RIGHT , 1)

       IF ( recom_WindFile .NE. ' '  ) THEN
C         CALL MDSREADFIELD( recom_WindFile, readBinaryPrec, 
C     &       'RS', 1, windSpeed0, intime0, myThid )
          CALL MDS_READ_FIELD( recom_WindFile, readBinaryPrec, 
      &       .FALSE., 'RS', 1, 1, 1, windSpeed0, intime0, myThid )
C         CALL MDSREADFIELD( recom_WindFile, readBinaryPrec, 
C     &       'RS', 1, windSpeed1, intime1, myThid )
          CALL MDS_READ_FIELD( recom_WindFile, readBinaryPrec, 
      &       .FALSE., 'RS', 1, 1, 1, windSpeed1, intime1, myThid )
       ENDIF
#ifndef RECOM_ATMOSPCO2_HISTORY
       IF ( recom_pCO2File .NE. ' '  ) THEN
C         CALL MDSREADFIELD( recom_pCO2File, readBinaryPrec, 
C     &       'RS', 1, atmosp0, intime0, myThid )
         CALL MDS_READ_FIELD( recom_pCO2File, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1, atmosp0, intime0, myThid )
C         CALL MDSREADFIELD( recom_pCO2File, readBinaryPrec, 
C     &       'RS', 1, atmosp1, intime1, myThid )
         CALL MDS_READ_FIELD( recom_pCO2File, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1, atmosp1, intime1, myThid )
       ENDIF 
#endif
CCV#ifdef ALLOW_RECOM_FEINPUT 
       IF ( recom_IronFile .NE. ' '  ) THEN
C         CALL MDSREADFIELD( recom_IronFile, readBinaryPrec, 
C     &       'RS', 1, feDustDeposit0, intime0, myThid )
         CALL MDS_READ_FIELD( recom_IronFile, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1, feDustDeposit0, intime0, myThid )
C         CALL MDSREADFIELD( recom_IronFile, readBinaryPrec, 
C     &       'RS', 1, feDustDeposit1, intime1, myThid )
         CALL MDS_READ_FIELD( recom_IronFile, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1, feDustDeposit1, intime1, myThid )
       ENDIF 
CCV#endif /* ALLOW_RECOM_FEINPUT */

#if defined(RECOM_WAVEBANDS) && defined(OASIM)
CEA   Theta(x,y,z,time) is read with READ_FLD_XYZ_RL,
CEA   Is it possible to use it here?       
CEA ----------------------- DIRECT --------------------------
CEA 1
      IF ( darwin_oasim_edFile01 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_edFile01, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp1, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_edFile01, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp2, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_ed0(i,j,1,bi,bj) = tmp1(i,j,bi,bj)
          oasim_ed1(i,j,1,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 2
      IF ( darwin_oasim_edFile02 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_edFile02, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp1, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_edFile02, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp2, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_ed0(i,j,2,bi,bj) = tmp1(i,j,bi,bj)
          oasim_ed1(i,j,2,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 3
      IF ( darwin_oasim_edFile03 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_edFile03, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp1, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_edFile03, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp2, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_ed0(i,j,3,bi,bj) = tmp1(i,j,bi,bj)
          oasim_ed1(i,j,3,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 4
      IF ( darwin_oasim_edFile04 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_edFile04, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp1, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_edFile04, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp2, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_ed0(i,j,4,bi,bj) = tmp1(i,j,bi,bj)
          oasim_ed1(i,j,4,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 5
      IF ( darwin_oasim_edFile05 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_edFile05, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp1, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_edFile05, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp2, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_ed0(i,j,5,bi,bj) = tmp1(i,j,bi,bj)
          oasim_ed1(i,j,5,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 6
      IF ( darwin_oasim_edFile06 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_edFile06, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp1, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_edFile06, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp2, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_ed0(i,j,6,bi,bj) = tmp1(i,j,bi,bj)
          oasim_ed1(i,j,6,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 7
      IF ( darwin_oasim_edFile07 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_edFile07, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp1, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_edFile07, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp2, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_ed0(i,j,7,bi,bj) = tmp1(i,j,bi,bj)
          oasim_ed1(i,j,7,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 8
      IF ( darwin_oasim_edFile08 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_edFile08, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp1, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_edFile08, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp2, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_ed0(i,j,8,bi,bj) = tmp1(i,j,bi,bj)
          oasim_ed1(i,j,8,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 9
      IF ( darwin_oasim_edFile09 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_edFile09, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp1, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_edFile09, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp2, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_ed0(i,j,9,bi,bj) = tmp1(i,j,bi,bj)
          oasim_ed1(i,j,9,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 10
      IF ( darwin_oasim_edFile10 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_edFile10, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp1, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_edFile10, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp2, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_ed0(i,j,10,bi,bj) = tmp1(i,j,bi,bj)
          oasim_ed1(i,j,10,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 11
      IF ( darwin_oasim_edFile11 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_edFile11, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp1, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_edFile11, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp2, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_ed0(i,j,11,bi,bj) = tmp1(i,j,bi,bj)
          oasim_ed1(i,j,11,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 12
      IF ( darwin_oasim_edFile12 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_edFile12, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp1, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_edFile12, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp2, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_ed0(i,j,12,bi,bj) = tmp1(i,j,bi,bj)
          oasim_ed1(i,j,12,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 13
      IF ( darwin_oasim_edFile13 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_edFile13, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp1, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_edFile13, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  tmp2, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_ed0(i,j,13,bi,bj) = tmp1(i,j,bi,bj)
          oasim_ed1(i,j,13,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF     
CEA --------------------- end direct ------------------------
CEA ---------------------- DIFFUSE ---------------------------      
CEA 1
      IF ( darwin_oasim_esFile01 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_esFile01, readBinaryPrec, 
     &       .FALSE., 'RS', tlam, 1, tlam,  oasim_es0, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_esFile01, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  oasim_es1, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_es0(i,j,1,bi,bj) = tmp1(i,j,bi,bj)
          oasim_es1(i,j,1,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 2
      IF ( darwin_oasim_esFile02 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_esFile02, readBinaryPrec, 
     &       .FALSE., 'RS', tlam, 1, tlam,  oasim_es0, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_esFile02, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  oasim_es1, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_es0(i,j,2,bi,bj) = tmp1(i,j,bi,bj)
          oasim_es1(i,j,2,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 3
      IF ( darwin_oasim_esFile03 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_esFile03, readBinaryPrec, 
     &       .FALSE., 'RS', tlam, 1, tlam,  oasim_es0, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_esFile03, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  oasim_es1, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_es0(i,j,3,bi,bj) = tmp1(i,j,bi,bj)
          oasim_es1(i,j,3,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 4
      IF ( darwin_oasim_esFile04 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_esFile04, readBinaryPrec, 
     &       .FALSE., 'RS', tlam, 1, tlam,  oasim_es0, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_esFile04, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  oasim_es1, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_es0(i,j,4,bi,bj) = tmp1(i,j,bi,bj)
          oasim_es1(i,j,4,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 5
      IF ( darwin_oasim_esFile05 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_esFile05, readBinaryPrec, 
     &       .FALSE., 'RS', tlam, 1, tlam,  oasim_es0, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_esFile05, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  oasim_es1, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_es0(i,j,5,bi,bj) = tmp1(i,j,bi,bj)
          oasim_es1(i,j,5,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 6
      IF ( darwin_oasim_esFile06 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_esFile06, readBinaryPrec, 
     &       .FALSE., 'RS', tlam, 1, tlam,  oasim_es0, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_esFile06, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  oasim_es1, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_es0(i,j,6,bi,bj) = tmp1(i,j,bi,bj)
         oasim_es1(i,j,6,bi,bj) = tmp2(i,j,bi,bj)
      ENDIF
CEA 7
      IF ( darwin_oasim_esFile07 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_esFile07, readBinaryPrec, 
     &       .FALSE., 'RS', tlam, 1, tlam,  oasim_es0, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_esFile07, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  oasim_es1, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_es0(i,j,7,bi,bj) = tmp1(i,j,bi,bj)
          oasim_es1(i,j,7,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 8
      IF ( darwin_oasim_esFile08 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_esFile08, readBinaryPrec, 
     &       .FALSE., 'RS', tlam, 1, tlam,  oasim_es0, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_esFile08, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  oasim_es1, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_es0(i,j,8,bi,bj) = tmp1(i,j,bi,bj)
          oasim_es1(i,j,8,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 9
      IF ( darwin_oasim_esFile09 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_esFile09, readBinaryPrec, 
     &       .FALSE., 'RS', tlam, 1, tlam,  oasim_es0, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_esFile09, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  oasim_es1, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_es0(i,j,9,bi,bj) = tmp1(i,j,bi,bj)
          oasim_es1(i,j,9,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 10
      IF ( darwin_oasim_esFile10 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_esFile10, readBinaryPrec, 
     &       .FALSE., 'RS', tlam, 1, tlam,  oasim_es0, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_esFile10, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  oasim_es1, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_es0(i,j,10,bi,bj) = tmp1(i,j,bi,bj)
          oasim_es1(i,j,10,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 11
      IF ( darwin_oasim_esFile11 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_esFile11, readBinaryPrec, 
     &       .FALSE., 'RS', tlam, 1, tlam,  oasim_es0, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_esFile11, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  oasim_es1, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_es0(i,j,11,bi,bj) = tmp1(i,j,bi,bj)
          oasim_es1(i,j,11,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 12
      IF ( darwin_oasim_esFile12 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_esFile12, readBinaryPrec, 
     &       .FALSE., 'RS', tlam, 1, tlam,  oasim_es0, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_esFile12, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  oasim_es1, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_es0(i,j,12,bi,bj) = tmp1(i,j,bi,bj)
          oasim_es1(i,j,12,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF
CEA 13
      IF ( darwin_oasim_esFile13 .NE. ' '  ) THEN
         nj0= (intime0-1)*tlam
         nj1= (intime1-1)*tlam
          CALL MDS_READ_FIELD( darwin_oasim_esFile13, readBinaryPrec, 
     &       .FALSE., 'RS', tlam, 1, tlam,  oasim_es0, intime0, myThid )
          CALL MDS_READ_FIELD( darwin_oasim_esFile13, readBinaryPrec, 
     &       .FALSE., 'RS', 1, 1, 1,  oasim_es1, intime1, myThid )
       DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
           DO j=1-Oly,sNy+Oly
            DO i=1-Olx,sNx+Olx
          oasim_es0(i,j,13,bi,bj) = tmp1(i,j,bi,bj)
          oasim_es1(i,j,13,bi,bj) = tmp2(i,j,bi,bj)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
      ENDIF

CEA -------------------- end diffuse -------------------------
#endif /* RECOM_WAVEBANDS && OASIM */

       _END_MASTER(myThid)
C
       _EXCH_XY_RS(windSpeed0, myThid )
       _EXCH_XY_RS(windSpeed1, myThid )
#ifndef RECOM_ATMOSPCO2_HISTORY
       _EXCH_XY_RS(atmosp0, myThid )
       _EXCH_XY_RS(atmosp1, myThid )
#endif
CCV#ifdef ALLOW_RECOM_FEINPUT 
       _EXCH_XY_RS(feDustDeposit0, myThid )
       _EXCH_XY_RS(feDustDeposit1, myThid )
CCV#endif /* ALLOW_RECOM_FEINPUT */

#if defined(RECOM_WAVEBANDS) && defined(OASIM)
       CALL EXCH_3D_RS   (oasim_ed0, tlam, myThid)
       CALL EXCH_3D_RS   (oasim_ed1, tlam, myThid)
       CALL EXCH_3D_RS   (oasim_es0, tlam, myThid)
       CALL EXCH_3D_RS   (oasim_es1, tlam, myThid)
#endif /* OASIM */      
C
       ENDIF
      
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO j=1-Oly,sNy+Oly
          DO i=1-Olx,sNx+Olx
           IF ( recom_WindFile .NE. ' '  ) THEN
             windSpeed(i,j,bi,bj) = bWght*windSpeed0(i,j,bi,bj)
     &                             +aWght*windSpeed1(i,j,bi,bj)
           ELSE
             windSpeed(i,j,bi,bj) = 1. _d 0*maskC(i,j,1,bi,bj)
           ENDIF
#ifndef RECOM_ATMOSPCO2_HISTORY
           IF ( recom_pCO2File .NE. ' '  ) THEN
             atmospCO2(i,j,bi,bj) = Bwght*atmosp0(i,j,bi,bj)
     &                             +aWght*atmosp1(i,j,bi,bj)
           ELSE
             atmospCO2(i,j,bi,bj) = 370. _d 0*maskC(i,j,1,bi,bj)
           ENDIF
#endif /* RECOM_ATMOSPCO2_HISTORY */
           
CCV#ifdef ALLOW_RECOM_FEINPUT 
           IF ( recom_IronFile .NE. ' ' ) THEN
            feDustDeposit(i,j,bi,bj) = bWght*feDustDeposit0(i,j,bi,bj)
     &                                +aWght*feDustDeposit1(i,j,bi,bj)
           ELSE
            feDustDeposit(i,j,bi,bj) = 0. _d 0*maskC(i,j,1,bi,bj)
           ENDIF
CCV#else
CCV           feDustDeposit(i,j,bi,bj)  = 0. _d 0*maskC(i,j,1,bi,bj)
CCV#endif /* ALLOW_RECOM_FEINPUT */
          ENDDO
         ENDDO
        ENDDO
       ENDDO

#if defined(RECOM_WAVEBANDS) && defined(OASIM)
    IF ( darwin_oasim_edFile01 .NE. ' '  ) THEN
c oasim data (load as W/m2 per band)
       DO ilam=1,tlam
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          DO j=1-Oly,sNy+Oly
           DO i=1-Olx,sNx+Olx       
               oasim_ed(i,j,ilam,bi,bj) = 
     &                          bWght*DBLE(oasim_ed0(i,j,ilam,bi,bj))
     &                         +aWght*DBLE(oasim_ed1(i,j,ilam,bi,bj))

          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
    ENDIF
    IF ( darwin_oasim_esFile01 .NE. ' '  ) THEN
       DO ilam=1,tlam
        DO bj = myByLo(myThid), myByHi(myThid)
         DO bi = myBxLo(myThid), myBxHi(myThid)
          DO j=1-Oly,sNy+Oly
           DO i=1-Olx,sNx+Olx            
               oasim_es(i,j,ilam,bi,bj) =
     &                          bWght*DBLE(oasim_es0(i,j,ilam,bi,bj))
     &                         +aWght*DBLE(oasim_es1(i,j,ilam,bi,bj))
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
    ENDIF
       
#endif /* RECOM_WAVEBANDS && OASIM */
     
C     ELSE
C not  periodicForcing
C
C endif for periodicForcing
      ENDIF
#endif /* ALLOW_EXF */


      
      
#ifdef RECOM_ATMOSPCO2_HISTORY
C read from a file (note:
C                   recom_pco2_int2=number entries to read
C                   recom_pco2_int3=start timestep,
C                   recom_pco2_int4=timestep between file entries)
C linearly interpolate between file entries
      ntim=int((myIter-recom_pco2_int3)/recom_pco2_int4)+1
      aWght = FLOAT(myIter-recom_pco2_int3)
      bWght = FLOAT(recom_pco2_int4)

CJH   taken out 0.5 in the following line as this would assume 
CJH   that for January in the first year the mean of the last 
CJH   year would have been read in -> the interpolation for the 
CJH   first time step would always be the mean of the first and 
CJH   second CO2 entry which is wrong for both cases: (a) mean 
CJH   centered in January is being read in, (b) mean of the 
CJH   actual calendar year (centered in July) is read in
CJH    aWght = 0.5 _d 0 +  aWght/bWght - FLOAT(ntim-1)
CJH   new now: uses start time step (recom_pco2_int3) as given
CJH   in data.recom for first CO2 entry. Set recom_pco2_int3 to 
CJH   two weeks in time stepping units to use yearly means centered
CJH   in January for correct interpolation

      aWght =  aWght/bWght - FLOAT(ntim-1)
      IF (aWght.GT.1. _d 0) THEN
         ntim=ntim+1
         aWght=aWght-1. _d 0
      ENDIF
      bWght = 1. _d 0 - aWght
      tmp=co2atmos(ntim,2)*bWght+co2atmos(ntim+1,2)*aWght

      ioUnit = standardMessageUnit
      WRITE(ioUnit,*) 'weights',ntim, aWght, bWght, tmp

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j=1-Oly,sNy+Oly
         DO i=1-Olx,sNx+Olx
          atmospCO2(i,j,bi,bj) = tmp
         ENDDO
        ENDDO
       ENDDO
      ENDDO
#endif /* RECOM_ATMOSPCO2_HISTORY */

     
#endif /* ALLOW_RECOM */

      RETURN
      END

