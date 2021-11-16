C $Header: /csys/software/cvsroot/REcoM/MITgcm/recom/recom_insolation.F,v 1.9 2008/02/06 13:24:24 mlosch Exp $
C $Name:  $

#include "GCHEM_OPTIONS.h"
#include "RECOM_OPTIONS.h"

CBOP
C     !ROUTINE: RECOM_INSOLATION
C     !INTERFACE:
      SUBROUTINE RECOM_INSOLATION(
     O     cosAngleOfIncidence, lengthOfDay, PARadiation,
#ifdef RECOM_RADTRANS
     O     solz,
#endif      
     I     iMin, iMax, jMin, jMax, bi, bj, myTime, myIter, myThid )
C     !DESCRIPTION: \bv
C     *==================================================================
C     | SUBROUTINE RECOM_INSOLATION
C     | o diurnal cycle for solar radiation
C     *==================================================================
C     |
C     | RECOM_INSOLATION computes the diurnal cycle of incoming
C     | solar radiation (light), in fact, the photosynthetically 
C     | available radiation (PAR), as a Gaussian bell curve with a 
C     | halfwidth determined by the local length of the day times 0.28.
C     | Adopted from REcoM related code by Markus Schartau.
C     | The incoming radiation is assumed to be a daily mean, if it is
C     | read from the short wave radiation, and is rescaled to be the
C     | insolation at noon.
C     | Alternatively, if the run time flag brockReedRadiation=.true.,
C     | the routines computes the insolation at noon from atronomical 
C     | parameters according to Brock (1981) and Reed (1977).
C     |
C     | NOTE: So far, this routine computes everything at each time
C     |       step. It may be sufficient to compute things like
C     |       lengthOfDay and InclinationAngle only once per day, but
C     |       I am not sure whether this reduces the computational
C     |       sufficiently.
C     |
C     | RECOM_INSOLATION is called from RECOM_FORCING.
C     |
C     | created by M. Losch, July, 2003
C     | modified Oct 21, 2003 by Martin Losch
C     | - REcoM needs the photosynthetically available light, so the
C     |   parameter parFrac is used to compute PAR as a fraction of
C     |   insolation; this causes renaming a few variables
C     | - REcoM needs either the mean insolation (e.g., in the case of 
C     |   Evans+Parslow, 1985) or a diurnal cycle. This is now possible 
C     |   for both brookReedInsolation=.true. and .false.
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
#ifdef ALLOW_RECOM
#include "RECOM.h"
#endif /* ALLOW_RECOM */
#ifdef SHORTWAVE_HEATING
# include "FFIELDS.h"
# ifdef ALLOW_LONGSTEP
#  include "LONGSTEP.h"
# endif
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     note: under the multi-threaded model myiter and 
C           mytime are local variables passed around as routine 
C           arguments. Although this is fiddly it saves the need to 
C           impose additional synchronisation points when they are 
C           updated.
C     iMin, iMax, jMin, jMax - loop-ranges
C     bi, bj - tile indices
C     myiter - iteration counter for this thread
C     mytime - time counter for this thread
C     mythid - thread number for this instance of the routine.
      INTEGER iMin, iMax, jMin, jMax
      integer bi, bj
      integer mythid 
      integer myiter
      _RL     myTime
C     cosAngleOfIncidence - of incoming solar irradiation in water (after
C                           refraction according to Snell's Law)
C     lengthOfDay  - length of the day (function of latitude and time of year)
C     PARadiation  - photosynthetically available radiation (PAR)
      _RS     cosAngleOfIncidence(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RS     lengthOfDay(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RS     PARadiation(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#ifdef RECOM_RADTRANS
C     solz  - solar zenith angle (radians)      
      _RS     solz(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#endif
      
#ifdef ALLOW_RECOM
C     !LOCAL VARIABLES:
C     == Local variables ==
      INTEGER i,j,k
C     dayfrac, yDay :: fraction of the current year, in radians
C     dayOfYear     :: dayfrac in days
C     locLat        :: local latitude, backed out from coriolis parameter
C     declination   :: of sun
      _RS     locLat       (1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RS     cosAngleNoon (1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RS     dayfrac, yDay, dayOfYear, declination, arg
      _RS     radiation, zenith, cos_zenith
      _RS     recipRadiusVector, facReduction
      _RS     integratedCycle
C     insolation    :: effective solar radiation at myTime
      _RS     insolation   (1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL     timeOfDay
      _RL     timeSinceNoon(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
C     widthFac  :: factor that sets the width of the Gaussian bell curve for 
C                  the diurnal cycle, 0.28 corresponds to 5% of the mid-day 
C                  (noon) insolation at sunrise, 0.3 corresponds to 6%
C     halfWidth :: widthFac*lengthOfDay, half width of Gaussian bell curve
      _RS     widthFac, halfWidth
      _RS     r_halfWidth  (1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      PARAMETER  ( widthFac = 0.28 )
C     newDay :: flag, only true during the first call to this routine
C               after midnight
      LOGICAL newDay
C     msgBuf    - Informational/error meesage buffer
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

#ifndef SHORTWAVE_HEATING
      IF ( .NOT. brockReedInsolation ) THEN
       WRITE(msgBuf,'(A)') 'S/R RECOM_INSOLATION'
       CALL PRINT_ERROR( msgBuf, 1 ) 
       WRITE(msgBuf,'(A)')
     &      'The insolation is meant to be computed from '//
     &      'the short wave heat,'
       CALL PRINT_ERROR( msgBuf, 1 ) 
       WRITE(msgBuf,'(A)')
     &      'but "SHORTWAVE_HEATING" is undefined.'
       CALL PRINT_ERROR( msgBuf, 1 ) 
       STOP 'ABNORMAL END: S/R RECOM_INSOLATION'
      ENDIF
#endif /* SHORTWAVE_HEATING */
C     some initialisation
      DO j=1-Oly,sNy+Oly
       DO i=1-Olx,sNx+Olx
        insolation(I,J) = 0. _d 0
       ENDDO
      ENDDO
C
C     is this a new day?
C
      IF ( int(myTime/recom_secondsPerDay) .gt. 
     &     int( (myTime-deltaTClock)/recom_secondsPerDay ) 
     &     .or. myIter .lt. 2 ) then
       newDay = .true.
      ELSE
       newDay = .false.
      ENDIF
C
C     find day (****NOTE for year starting in winter*****)
C     fraction of the year
      dayfrac = mod(myTime,daysPerYear*recom_secondsPerDay)
     &                           /(daysPerYear*recom_secondsPerDay)
      dayOfYear = dayfrac*daysPerYear
C     convert to radians
      yDay = 2.*pi*dayfrac
C     cosine zenith angle (paltridge and platt)
      declination = 0.006918 
     &     -     0.399912*cos(yDay)
     &     +     0.070257*sin(yDay)
     &     -     0.006758*cos(2.0*yDay)
     &     +     0.000907*sin(2.0*yDay)
     &     -     0.002697*cos(3.0*yDay)
     &     +     0.001480*sin(3.0*yDay) 
c$$$C     declination as a function of Julian day COOPER(1969)
c$$$      declination=23.45*sin(2.*pi*(equinox+dayOfYear)/daysPerYear)
c$$$     &     *deg2rad
      IF ( diurnalInsolation ) THEN
C
C     compute the time of the day, that is how many seconds since
C     midnight
C     
       timeOfDay   = mod(myTime,recom_secondsPerDay)
       IF ( usingCartesianGrid ) THEN
        DO j=1-Oly,sNy+Oly
         DO i=1-Olx,sNx+Olx
          timeSinceNoon(I,J) = timeOfDay-noonTime
         ENDDO
        ENDDO
       ELSE 
C     For computations on the sphere we can make the local time 
C     dependent on longitude
        DO j=1-Oly,sNy+Oly
         DO i=1-Olx,sNx+Olx
          timeSinceNoon(I,J) = abs(timeOfDay  - 0.5*recom_secondsPerDay
     &         + XC(I,J,bi,bj)*recom_secondsPerDay/360.)
          IF ( timeSinceNoon(I,J) .GT. 0.5*recom_secondsPerDay ) THEN
           timeSinceNoon(I,J) = recom_secondsPerDay-timeSinceNoon(I,J)
          ENDIF
         ENDDO
        ENDDO
       ENDIF
      ENDIF

CCV useSphereF=.true. has been replaced by selectCoriMap=2 as of version
CCV 62o      
CCV      IF ( .NOT. useSphereF ) THEN
      IF ( selectCoriMap .LT. 2 ) THEN
       DO j=jMin,jMax
        DO i=iMin,iMax
C     latitute in radians, backed out from coriolis parameter
C     (makes latitude independent of grid)
         locLat(i,j) = asin( fCori(i,j,bi,bj)/(2. _d 0 *omega) )
        ENDDO
       ENDDO
      ELSE
       DO j=jMin,jMax
        DO i=iMin,iMax
C     this is more accurate for computations on sphere
         locLat(i,j) = YC(I,J,bi,bj)*deg2rad
        ENDDO
       ENDDO
      ENDIF	
      DO j=jMin,jMax
       DO i=iMin,iMax
C     compute length of day
C     as a function of the hour-angle(earth=15 deg/hour) 
C     in seconds
        arg = -1.*( tan(locLat(i,j))*tan(declination) )
        IF ( (arg-1.) .gt. 0. ) THEN
         lengthOfDay(i,j) = 0. _d 0
        ELSEIF ( (arg+1.) .lt. 0. ) THEN
         lengthOfDay(i,j) = recom_secondsPerDay
        ELSE
         lengthOfDay(i,j) = ( acos(arg)/deg2rad/15.
     &        )*2.*recom_secondsPerHour
        ENDIF
C     compute COSINE(angle of incidence of solar radiation)
        cosAngleNoon(i,j) = SIN(locLat(i,j))*SIN(declination) 
     &       +         COS(locLat(i,j))*COS(declination)
C        cosAngleOfIncidence(i,j) = SQRT( 1.
C     &       - ( (1. - cosAngleNoon(i,j)**2)/nWater )**2 )
        cosAngleOfIncidence(i,j) = SQRT( 1.
     &       - (1. - cosAngleNoon(i,j)**2)/nWater**2 )
#ifdef RECOM_RADTRANS         
C        cosine of zenith angle at noon !! MILANKOVITCH(1930)
         cos_zenith=cosAngleNoon(i,j)       
         zenith=acos(cos_zenith)
         solz(i,j)=zenith
#endif   
       ENDDO
      ENDDO
      DO j=jMin,jMax
       DO i=iMin,iMax
C     compute half width of gaussian bell curve and its reciprocal
C     in seconds
C     factor that sets the with of the Gaussian bell curve for 
C     the diurnal cycle
        halfWidth = lengthOfDay(i,j)*widthFac
        IF ( halfWidth .LE. 0. _d 0 ) THEN
         r_halfWidth(i,j) = 0. _d 0 
        ELSE
         r_halfWidth(i,j) = 1. _d 0 / halfWidth
        ENDIF
       ENDDO
      ENDDO
C
C
C      
      IF ( brockReedInsolation ) THEN
       DO j=jMin,jMax
        DO i=iMin,iMax
C     --- calculation of radiation at the top of the ocean at noon,
C     --- according to BROCK(1981)+REED(1977) 
C     
C     radius vector NICHOLLS and CHILD(1979)
C     radiusVector=1./(sqrt(1.+(0.033*cos((2.*pi*time)/daysPerYear))))
         recipRadiusVector=
     &        sqrt( 1.  + (0.033*cos(2.*pi*dayOfYear/daysPerYear)) )
         
C     cosine of zenith angle at noon !! MILANKOVITCH(1930)
         cos_zenith=cosAngleNoon(i,j)       
CML                    sin(locLat(i,j))*sin(declination)
CML     &             +cos(locLat(i,j))*cos(declination)
         zenith=acos(cos_zenith)
C 
C     radiation at top of atmosphere
C
         radiation=solarConstant*cos_zenith*recipRadiusVector**2
C
C     ------------- factor of reduction of insolation ----
C     estimation of the insolation over the ocean REED(1977)
         facReduction=1.
     &        - 0.62*cloudCover
     &        + 0.0019*(.5*pi-zenith)

C     radiation at the top of the ocean surface at noon
         insolation(I,J) = facReduction*radiation

         IF ( .NOT. diurnalInsolation ) THEN
C     compute (approximate) normalization constant (integral 
C     of a Gaussian bell curve)
          integratedCycle = 0. _d 0
          IF ( r_halfWidth(i,j) .GT. 0. _d 0 ) THEN
           DO K = 1, 24
            integratedCycle = integratedCycle
     &           + EXP( -( real( k-12 )*recom_secondsPerHour 
     &           * r_halfWidth(i,j) )**2 )
           ENDDO	
          ENDIF
C     normalize assuming a Gaussian bell curve as the diurnal cycle
C     to obtain the mean irradiation from the irradiation at noon
          insolation(I,J) = insolation(I,J)*(integratedCycle)/24.
         ENDIF
        ENDDO
       ENDDO
C
C       
C
      ELSE 
       DO j=jMin,jMax
        DO i=iMin,iMax
C     not brookReedInsolation:
C     daily averaged insolation is the short wave heat flux Qsw
#ifdef SHORTWAVE_HEATING
# ifdef ALLOW_LONGSTEP
        insolation(I,J) = - LS_Qsw(i,j,bi,bj)
# else
        insolation(I,J) = - Qsw(i,j,bi,bj)
# endif /* ALLOW_LONGSTEP */
#else
        insolation(I,J) =  0.
#endif /* SHORTWAVE_HEATING */
         IF ( diurnalInsolation ) THEN
C     compute (approximate) normalization constant (integral 
C     of a Gaussian bell curve)
          integratedCycle = 0. _d 0
          IF ( r_halfWidth(i,j) .GT. 0. _d 0 ) THEN
           DO K = 1, 24
            integratedCycle = integratedCycle
     &           + EXP( -( real( k-12 )*recom_secondsPerHour
     &           * r_halfWidth(I,J) )**2 )
           ENDDO	
          ENDIF
C     normalize assuming a Gaussian bell curve as the diurnal cycle
C     to obtain irradiation at noon from the mean irradiation
          IF ( integratedCycle .GT. 0. _d 0 ) THEN
           insolation(I,J) = 24.*insolation(I,J)/(integratedCycle)
          ELSE
           insolation(I,J) = 0. _d 0
          ENDIF
         ENDIF
        ENDDO
       ENDDO
      ENDIF

      DO j=jMin,jMax
       DO i=iMin,iMax
C     compute the diurnal cycle at time myTime, that is timeOfDay
C     the diurnal cylce is assumed to be a Gaussian bell curve with peak
C     at noon, background radiation of 1 W/m^2
        IF ( diurnalInsolation ) THEN
         IF ( r_halfWidth(I,J) .le. 0. _d 0 ) THEN
          insolation(I,J) = 1. _d 0
         ELSE
          insolation(I,J) = 1. _d 0 + insolation(I,J)
     &        *EXP( -( ( timeSinceNoon(I,J) ) * r_halfWidth(I,J) )**2 )
         ENDIF
        ELSE
         insolation(I,J) = 1. _d 0 + insolation(I,J)
        ENDIF
C     compute the photosynthetically available radiation as a fraction
C     of the insolation (short wave heat)
        PARadiation(i,j) = parFrac*insolation(i,j)
C
C     -- end j/i-loops
       ENDDO
      ENDDO
C
#endif /* ALLOW_RECOM */
      RETURN
      END



