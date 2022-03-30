C $Header: /csys/software/cvsroot/REcoM/MITgcm/recom/recom_init_fixed.F,v 1.4 2008/04/17 14:25:52 mlosch Exp $
C $Name:  $

#include "GCHEM_OPTIONS.h"
#include "RECOM_OPTIONS.h"

CBOP
C !ROUTINE: RECOM_INIT_FIXED

C !INTERFACE: ============================================================
      SUBROUTINE RECOM_INIT_FIXED( myThid )

C !DESCRIPTION:
C     Initialize fixed (not subject to adjoint) RECOM data structures
C     (to be called from S/R GCHEM_INIT_FIXED)

C !USES: ================================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#ifdef ALLOW_MNC
#include "MNC_PARAMS.h"
#endif
#ifdef ALLOW_EXF
#include "RECOM.h"
#include "RECOM_EXF.h"
#include "cal.h"
#include "EXF_PARAM.h"
#include "EXF_CONSTANTS.h"
#endif 
#ifdef RECOM_ATMOSPCO2_HISTORY
#include "RECOM_PARAMS.h"
#endif
#ifdef RECOM_WAVEBANDS
#include "SPECTRAL_SIZE.h"
#include "SPECTRAL.h"
#include "WAVEBANDS_PARAMS.h"
#endif  

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none

#ifdef ALLOW_RECOM

C !LOCAL VARIABLES: ====================================================
      LOGICAL useVariableK
#ifdef ALLOW_EXF
      INTEGER date_array(4), difftime(4), yearStartDate(4)
#endif
#ifdef RECOM_ATMOSPCO2_HISTORY
      INTEGER k, iUnit
#endif
CEOP

#ifdef RECOM_WAVEBANDS
      INTEGER ilam
C      _RL planck, c, hc, oavo, hcoavo, rlamm
C#endif

C#ifdef RECOM_WAVEBANDS
C#ifdef OASIM
Cc     Quanta conversion
C      planck = 6.6256 _d -34         !Plancks constant J sec
C      c = 2.998 _d 8                 !speed of light m/sec
C      hc = 1.0/(planck*c)
C      oavo = 1.0/6.023 _d 23         ! 1/Avogadros number
C      hcoavo = hc*oavo
C      do ilam = 1,tlam
C       rlamm = darwin_waves(ilam)*1 _d -9       !lambda in m
C       WtouEins(ilam) = 1 _d 6*rlamm*hcoavo     !Watts to uEin/s conversion
C      enddo
c#endif
c     read the files
        call wavebands_init_fixed(myThid)
#endif

        
#ifdef RECOM_O2
C define Schmidt no. coefficients for O2
C based on Keeling et al [GBC, 12, 141, (1998)]
      sox1 = 1638.0 _d 0
      sox2 = -81.83 _d 0
      sox3 =    1.483 _d 0
      sox4 =   -0.008004 _d 0

C coefficients for determining saturation O2
      oA0=  2.00907 _d 0
      oA1=  3.22014 _d 0
      oA2=  4.05010 _d 0
      oA3=  4.94457 _d 0
      oA4= -2.56847 _d -1
      oA5=  3.88767 _d 0
      oB0= -6.24523 _d -3
      oB1= -7.37614 _d -3
      oB2= -1.03410 _d -2
      oB3= -8.17083 _d -3
      oC0= -4.88682 _d -7
#endif
C
C     initialise additional output variables
C
#ifdef ALLOW_MNC     
      IF (useMNC) THEN
       CALL MNC_CW_ADD_VNAME('WCtave', 'Cen_xy_Hn__C__t', 4,5,myThid)
       CALL MNC_CW_ADD_VATTR_TEXT('WCtave','units','mmolC/m^2/s', 
     &      myThid)
       CALL MNC_CW_ADD_VNAME('WNtave', 'Cen_xy_Hn__C__t', 4,5, myThid)
       CALL MNC_CW_ADD_VATTR_TEXT('WNtave','units','mmolN/m^2/s',  
     &      myThid)
CCV#ifdef ALLOW_RECOM_SILICATE
       CALL MNC_CW_ADD_VNAME('WSitave', 'Cen_xy_Hn__C__t', 4,5, myThid)
       CALL MNC_CW_ADD_VATTR_TEXT('WSitave','units','mmolSi/m^2/s',
     &      myThid)
CCV#endif /* ALLOW_RECOM_SILICATE */
       useVariableK = useKPP .OR. usePP81 .OR. useMY82 .OR. useGGL90
     &      .OR. useGMredi .OR. ivdc_kappa.NE.0.
       IF ( useVariableK ) THEN
        CALL MNC_CW_ADD_VNAME('Cdiftave', 'Cen_xy_Hn__C__t',4,5,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('Cdiftave','units','mmolC/m^2/s',
     &       myThid)
        CALL MNC_CW_ADD_VNAME('Ndiftave', 'Cen_xy_Hn__C__t',4,5,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('Ndiftave','units','mmolN/m^2/s',
     &       myThid)
CCV#ifdef ALLOW_RECOM_SILICATE
        CALL MNC_CW_ADD_VNAME('Sidiftave', 'Cen_xy_Hn__C__t',4,5,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('Sidiftave','units','mmolSi/m^2/s',
     &       myThid)
CCV#endif /* ALLOW_RECOM_SILICATE */

        CALL MNC_CW_ADD_VNAME('BenC','Cen_xy_Hn__-__t',3,4,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('BenC','units','mmol/m^3', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('BenC','description',
     &       'benthic carbon', myThid)
        CALL MNC_CW_ADD_VNAME('BenN','Cen_xy_Hn__-__t',3,4,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('BenN','units','mmol/m^3', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('BenN','description',
     &       'benthic nitrogen', myThid)
CCV#ifdef ALLOW_RECOM_SILICATE
        CALL MNC_CW_ADD_VNAME('BenSi','Cen_xy_Hn__-__t',3,4,myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('BenSi','units','mmol/m^3', myThid)
        CALL MNC_CW_ADD_VATTR_TEXT('BenSi','description',
     &       'benthic silicate', myThid)
CCV#endif /* ALLOW_RECOM_SILICATE */

       ENDIF
      ENDIF
#endif /* ALLOW_MNC */

#ifdef ALLOW_CAL
      IF ( recom_pCO2file .NE. ' ' .and. pCO2period .gt. 0. ) THEN
       call cal_FullDate  ( pCO2startdate1,   pCO2startdate2,
     &      date_array, mythid )
       call cal_TimePassed(modelstartdate,date_array,difftime,mythid)
       call cal_ToSeconds ( difftime, pCO2startdate,   mythid )
       pCO2startdate    = modelstart  + pCO2startdate
      ENDIF
      IF ( recom_ironfile .NE. ' ' .and. fedustperiod .gt. 0. ) THEN
       call cal_FullDate  ( feduststartdate1, feduststartdate2,
     &      date_array, mythid )
       call cal_TimePassed(modelstartdate,date_array,difftime,mythid)
       call cal_ToSeconds ( difftime, feduststartdate, mythid )
       feduststartdate    = modelstart  + feduststartdate
      ENDIF
#endif /* ALLOW_CAL */

#ifdef ALLOW_EXF
      IF ( useEXF ) THEN
      _BEGIN_MASTER( myThid )
      IF ( recom_pCO2file .NE. ' ' .AND. pCO2period .GT. 0. ) THEN
       CALL cal_FullDate  ( pCO2startdate1, pCO2startdate2,
     &      date_array, myThid )
       IF ( useExfYearlyFields ) THEN
        yearStartDate(1) = int(date_array(1)/10000.) * 10000 + 101
        yearStartDate(2) = 0
        yearStartDate(3) = date_array(3)
        yearStartDate(4) = date_array(4)
        CALL cal_TimePassed( yearStartDate,date_array,difftime,myThid)
        CALL cal_ToSeconds ( difftime,   pCO2startdate      ,myThid)
       ELSE
        CALL cal_TimePassed(modelstartdate,date_array,difftime,myThid)
        CALL cal_ToSeconds ( difftime,   pCO2startdate      ,myThid)
        pCO2startdate  =  modelstart + pCO2startdate
       ENDIF
      ENDIF
      IF ( recom_ironfile .NE. ' ' .AND. fedustperiod .GT. 0. ) THEN
       CALL cal_FullDate  ( feduststartdate1, feduststartdate2,
     &      date_array, myThid )
       IF ( useExfYearlyFields ) THEN
        yearStartDate(1) = int(date_array(1)/10000.) * 10000 + 101
        yearStartDate(2) = 0
        yearStartDate(3) = date_array(3)
        yearStartDate(4) = date_array(4)
        CALL cal_TimePassed( yearStartDate,date_array,difftime,myThid)
        CALL cal_ToSeconds ( difftime,   feduststartdate      ,myThid)
       ELSE
        CALL cal_TimePassed(modelstartdate,date_array,difftime,myThid)
        CALL cal_ToSeconds ( difftime,   feduststartdate      ,myThid)
        feduststartdate  =  modelstart + feduststartdate
       ENDIF
      ENDIF
#ifdef OASIM
      IF (darwin_oasim_edFile01.NE.' ' .AND.OASIMperiod .GT. 0.) THEN
       CALL cal_FullDate  ( OASIMstartdate1, OASIMstartdate2,
     &      date_array, myThid )
       IF ( useExfYearlyFields ) THEN
        yearStartDate(1) = int(date_array(1)/10000.) * 10000 + 101
        yearStartDate(2) = 0
        yearStartDate(3) = date_array(3)
        yearStartDate(4) = date_array(4)
        CALL cal_TimePassed( yearStartDate,date_array,difftime,myThid)
        CALL cal_ToSeconds ( difftime,   OASIMstartdate      ,myThid)
       ELSE
        CALL cal_TimePassed(modelstartdate,date_array,difftime,myThid)
        CALL cal_ToSeconds ( difftime,   OASIMstartdate      ,myThid)
        OASIMstartdate  =  modelstart + OASIMstartdate
       ENDIF
      ENDIF
#endif     
      _END_MASTER( myThid )
      _BARRIER
      ENDIF
#endif /* ALLOW_EXF */

#ifdef ALLOW_DIAGNOSTICS
C     Define diagnostics Names :
      IF ( useDiagnostics ) THEN
        CALL RECOM_DIAGNOSTICS_INIT( myThid )
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#ifdef RECOM_ATMOSPCO2_HISTORY
      _BEGIN_MASTER( myThid )
C     Read in a history of atmopheric pCO2
      IF ( recom_pco2_int1.EQ.2 ) THEN
        CALL MDSFINDUNIT( iUnit, mythid )
        OPEN(UNIT=iUnit,FILE='co2atmos.dat',STATUS='old')
        DO k=1,recom_pco2_int2
          READ(iUnit,*) co2atmos(k,1),co2atmos(k,2)
          WRITE(standardMessageUnit,*) 'co2atmos',
     &         co2atmos(k,1),co2atmos(k,2)
        ENDDO
        CLOSE(iUnit)
      ENDIF
      _END_MASTER( myThid )
      _BARRIER
#endif /* RECOM_ATMOSPCO2_HISTORY */

#endif /* ALLOW_RECOM */

      RETURN
      END
