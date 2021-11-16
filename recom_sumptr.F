C $Header: /csys/software/cvsroot/REcoM/MITgcm/recom/recom_sumptr.F,v 1.4 2008/04/17 14:46:15 mlosch Exp $
C $Name:  $

#include "GCHEM_OPTIONS.h"
#include "RECOM_OPTIONS.h"

CBOP
C !ROUTINE: RECOM_SUMPTR

C !INTERFACE: ============================================================
      SUBROUTINE RECOM_SUMPTR( 
     I     pTracer, nTracers, bi, bj, myThid )

C !DESCRIPTION:
C     copy one or the sum of a few passive tracer field onto recomPtr

C !USES: ================================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#ifdef ALLOW_RECOM
#include "RECOM.h"
#endif /* ALLOW_RECOM */

C !INPUT PARAMETERS: ===================================================
C     pTracer - input Fields
C     nTracers - number of tracers
C     myThid   - Thread number for this instance of the routine.
C     bi, bj   - tile numbering
      INTEGER myThid
      INTEGER nTracers, bi, bj
      _RL pTracer(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,nTracers)

C !OUTPUT PARAMETERS: ==================================================

#ifdef ALLOW_RECOM
C !LOCAL VARIABLES ====================================================
C     i,j,k - Loop counters
      INTEGER i, j, k
CEOP

      DO K=1, Nr
       DO J=1,sNy
        DO I=1,sNx
C     this is the total POM = particulate organic matter 
C     = the sum of detritus, phytoplankton, zooplankton, 
C     and maybe a fraction of D/EOM
         recomC(i,j,k,bi,bj) = pTracer(i,j,k,bi,bj,iphyc)
     &                       + pTracer(i,j,k,bi,bj,idetc)           
     &                       + pTracer(i,j,k,bi,bj,ihetc)
#ifdef NEXT_STEP
     &                       + pTracer(i,j,k,bi,bj,idiac) 
#endif          
     &                       + fracTEP*pTracer(i,j,k,bi,bj,idoc)
         recomN(i,j,k,bi,bj) = pTracer(i,j,k,bi,bj,iphyn)
     &                       + pTracer(i,j,k,bi,bj,idetn)           
     &                       + pTracer(i,j,k,bi,bj,ihetn)
#ifdef NEXT_STEP           
     &                       + pTracer(i,j,k,bi,bj,idian)  
#endif         
     &                       + fracTEP*pTracer(i,j,k,bi,bj,idon)
#ifdef NEXT_STEP
CCV#ifdef ALLOW_RECOM_SILICATE
         recomSi(i,j,k,bi,bj)= pTracer(i,j,k,bi,bj,idiasi)
     &                       + pTracer(i,j,k,bi,bj,idetsi)           
CCV#endif /* ALLOW_RECOM_SILICATE */
#endif
        ENDDO
       ENDDO
      ENDDO

#endif /* ALLOW_RECOM */
      RETURN
      END

CBOP
C !ROUTINE: RECOM_CUMUL_DIF_1T

C !INTERFACE: ==========================================================
      SUBROUTINE RECOM_CUMUL_DIF_1T( 
     I     kappaR, dT, bi,bj,myThid )

C !DESCRIPTION:
C     Calls the the time averaging cumulation routine for selected 
C     tracers

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_FIELDS.h"
#include "RECOM.h"
C !INPUT PARAMETERS: ===================================================
C  bi,bj   :: tile indices
C  KappaR  :: vertical diffusion coefficient
C  myThid  :: thread number
      INTEGER bi,bj
      _RL KappaR(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL dT
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none

#if ( (defined ALLOW_TIMEAVE) && (defined ALLOW_PTRACERS) )
#ifdef ALLOW_RECOM
C !LOCAL VARIABLES: ====================================================
CEOP

      IF (implicitDiffusion) THEN
         CALL RECOM_SUMPTR(    gPtr, PTRACERS_num, bi, bj, myThid )
      ELSE
         CALL RECOM_SUMPTR( pTracer, PTRACERS_num, bi, bj, myThid )
      ENDIF

      CALL TIMEAVE_CUMUL_DIF_1T(CDiffRtave, recomC, kappaR,
     I     Nr, 3, dT, bi, bj, myThid)
      CALL TIMEAVE_CUMUL_DIF_1T(NDiffRtave, recomN, kappaR,
     I     Nr, 3, dT, bi, bj, myThid)
#ifdef NEXT_STEP
CCV#ifdef ALLOW_RECOM_SILICATE
      CALL TIMEAVE_CUMUL_DIF_1T(SiDiffRtave, recomSi, kappaR,
     I     Nr, 3, dT, bi, bj, myThid)
CCV#endif /* ALLOW_RECOM_SILICATE */
#endif

#endif /* ALLOW_RECOM */
#endif /* ALLOW_PTRACERS && ALLOW_TIMEAVE */

      RETURN
      END

