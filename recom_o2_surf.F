C     Adapted from the MITgcm DIC package routine o2_surfforcing.F
#include "RECOM_OPTIONS.h"
CBOP
C !ROUTINE: RECOM_O2_SURF

C !INTERFACE: ==========================================================
      SUBROUTINE RECOM_O2_SURF( PTR_O2,windSpeedLoc, FluxO2Loc,
     I           bi,bj,iMin,iMax,jMin,jMax,
     I           myIter, myTime, myThid )

C !DESCRIPTION:
C     Calculate the oxygen air-sea flux terms
C
C     input:
C     PTR_O2
C     Windspeedloc
C     
C     output:

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "FFIELDS.h"
C#include "DIC_VARS.h"
#include "RECOM_PARAMS.h"

c  !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
C  myIter               :: current timestep
C  myTime               :: current time
C  PTR_O2               :: oxygen tracer field
      _RL myTime
      _RL PTR_O2(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL windSpeedLoc  (1-olx:snx+olx,1-oly:sny+oly)
      INTEGER iMin,iMax,jMin,jMax, bi, bj
      INTEGER myIter, myThid

c  !OUTPUT PARAMETERS: ===================================================
C  SGO2                  :: air-sea exchange of oxygen
      _RL  FluxO2Loc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

#ifdef ALLOW_PTRACERS
#ifdef RECOM_O2
C#ifdef ALLOW_O2

C !LOCAL VARIABLES: ===================================================
C I, J, K - Loop counters
      INTEGER I,J,K
      EXTERNAL piston_velocity
      _RL piston_velocity
C Solubility relation coefficients
      _RL SchmidtNoO2(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL O2sat(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL Kwexch(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
C      _RL FluxO2Loc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  aTT
      _RL  aTK
      _RL  aTS
      _RL  aTS2
      _RL  aTS3
      _RL  aTS4
      _RL  aTS5
      _RL  o2s
      _RL  ttemp
      _RL  stemp
      _RL  oCnew
      _RL piston_vel
CEOP

C changes only at uppermost model box
      K=1

C calculate SCHMIDT NO. for O2
        DO j=jmin,jmax
          DO i=imin,imax
            IF (maskC(i,j,k,bi,bj).NE.0.) THEN
              ttemp = theta(i,j,k,bi,bj)
              stemp = salt(i,j,k,bi,bj)

              SchmidtNoO2(i,j) =
     &            sox1
     &          + sox2 * ttemp
     &          + sox3 * ttemp*ttemp
     &          + sox4 * ttemp*ttemp*ttemp

C Determine surface flux of O2
C exchange coeff accounting for Schmidt no Scaling for ice cover outside.
C Kwexch_Pre= pisvel*(1-fice): previously computed in dic_surfforcing.F
C pisten_vel [m/d]
              piston_vel =  piston_velocity(ttemp,windSpeedLoc(i,j))
     &                    / sqrt(SchmidtNoO2(i,j)/660.0 _d 0)*0.24

C determine saturation O2
C using Garcia and Gordon (1992), L&O (mistake in original ?)
              aTT  = 298.15 _d 0 -ttemp
              aTK  = 273.15 _d 0 +ttemp
              aTS  = log(aTT/aTK)
              aTS2 = aTS*aTS
              aTS3 = aTS2*aTS
              aTS4 = aTS3*aTS
              aTS5 = aTS4*aTS

              oCnew  = oA0 + oA1*aTS + oA2*aTS2 + oA3*aTS3 +
     &            oA4*aTS4 + oA5*aTS5
     &          + stemp*(oB0 + oB1*aTS + oB2*aTS2 + oB3*aTS3)
     &          + oC0*(stemp*stemp)

              o2s = EXP(oCnew)

C Convert from ml/l to mmol/m^3
              O2sat(i,j) = o2s/22391.6 _d 0 * 1. _d 6

C Determine flux, inc. correction for local atmos surface pressure
              FluxO2Loc(i,j) = piston_vel*
     &                     ( O2sat(i,j)
     &                      - PTR_O2(i,j))
            ELSE
              FluxO2Loc(i,j) = 0. _d 0
            ENDIF


          END DO
        END DO

C update surface tendencies
C        DO j=jmin,jmax
C          DO i=imin,imax
C           SGO2(i,j)= FluxO2Loc(i,j)
C     &         *recip_drF(K) * recip_hFacC(i,j,K,bi,bj)
C          ENDDO
C         ENDDO
#endif
#endif


        RETURN
        END

