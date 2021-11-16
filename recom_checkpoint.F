C$Header: /csys/software/cvsroot/REcoM/MITgcm/recom/recom_checkpoint.F,v 1.4 2008/04/17 14:25:52 mlosch Exp $
C$Name:  $

#include "CPP_OPTIONS.h"
#include "GCHEM_OPTIONS.h"
#ifdef ALLOW_RECOM
#include "RECOM_OPTIONS.h"
#endif /* ALLOW_RECOM */

CBOP
C !ROUTINE: RECOM_WRITE_PICKUP

C !INTERFACE: ==========================================================
      SUBROUTINE RECOM_WRITE_PICKUP( permPickup,
     I                               suff,myTime,myIter,myThid )

C !DESCRIPTION:
C     Writes current state of REcoM specific tracers to a pickup file

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#ifdef ALLOW_RECOM
#include "RECOM.h"
#endif /* ALLOW_RECOM */

C !INPUT PARAMETERS: ===================================================
C     permPickup :: write a permanent pickup (not used)
C     suff       :: suffix for pickup file (eg. ckptA or 0000000010)
C     myTime     :: Current time in simulation
C     myIter     :: Current iteration number in simulation
C     myThid     :: My Thread Id number
      LOGICAL permPickup
      CHARACTER*(*) suff
      INTEGER myIter
      _RL myTime
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none

#ifdef ALLOW_RECOM
C !LOCAL VARIABLES: ====================================================
C  iTracer		:: loop indices
C  iRec    		:: record number
C  fn   		:: character buffer for creating filename
C  prec 		:: precision of pickup files
C  lgf  		:: flag to write "global" files
      INTEGER prec,iRec,iChar,lChar
      CHARACTER*(MAX_LEN_FNAM) fn
      LOGICAL lgf
CEOP

      lChar = 0
      DO iChar = 1,len(suff)
       IF ( suff(iChar:iChar) .NE. ' ') lChar=iChar
      ENDDO
      WRITE(fn,'(A,A)') 'pickup_recom.',suff(1:lChar)
      prec = precFloat64
      lgf = globalFiles

C Write fields as consecutive records
C      CALL MDSWRITEFIELD(fn,prec,lgf,'RL',1,hPlus,1,myIter,myThid)
      CALL WRITE_REC_3D_RL(fn,prec,1,hPlus,1,myIter,myThid)
C      CALL MDSWRITEFIELD(fn,prec,lgf,'RL',1,
C     &     benthicLayerN,2,myIter,myThid)
      CALL WRITE_REC_3D_RL(fn,prec,1,
     &     benthicLayerN,2,myIter,myThid)
C      CALL MDSWRITEFIELD(fn,prec,lgf,'RL',1,
C     &      benthicLayerC,3,myIter,myThid)
      CALL WRITE_REC_3D_RL(fn,prec,1,
     &      benthicLayerC,3,myIter,myThid)
CCV#ifdef ALLOW_RECOM_SILICATE
C      CALL MDSWRITEFIELD(fn,prec,lgf,'RL',1,
C     &      benthicLayerSi,4,myIter,myThid)
      CALL WRITE_REC_3D_RL(fn,prec,1,
     &      benthicLayerSi,4,myIter,myThid)
CCV#endif /* ALLOW_RECOM_SILICATE */
      CALL WRITE_REC_3D_RL(fn,prec,1,
     &      benthicLayerCaCO3,5,myIter,myThid)

#endif /* ALLOW_RECOM */
      RETURN
      END

CBOP
C !ROUTINE: RECOM_READ_PICKUP

C !INTERFACE: ==========================================================
      SUBROUTINE RECOM_READ_PICKUP( myIter,myThid )

C !DESCRIPTION:
C     Reads current state of recom specific tracers from a pickup file

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#ifdef ALLOW_RECOM
#include "RECOM.h"
#endif /* ALLOW_RECOM */

C !INPUT PARAMETERS: ===================================================
C  myIter               :: time-step number
C  suff                 :: filename suffix
C  myThid               :: thread number
      INTEGER myIter
      CHARACTER*(10) suff
      INTEGER myThid

C !OUTPUT PARAMETERS: ==================================================
C  none

#ifdef ALLOW_RECOM
C !LOCAL VARIABLES: ====================================================
C  iTracer		:: loop indices
C  iRec    		:: record number
C  fn   		:: character buffer for creating filename
C  prec 		:: precision of pickup files
C  fn                   :: suffix for pickup file (eg. 0000000010)
      INTEGER prec,iRec
      CHARACTER*(MAX_LEN_MBUF) fn
CEOP

      WRITE(fn,'(A,I10.10)') 'pickup_recom.',myIter
      prec = precFloat64

C Read fields as consecutive records
C      CALL MDSREADFIELD(fn,prec,'RL',1,hPlus,1,myThid)
       CALL READ_REC_3D_RL(fn,prec,1,hPlus,1,myIter,myThid)
C      CALL MDSREADFIELD(fn,prec,'RL',1,benthicLayerN,2,myThid)
       CALL READ_REC_3D_RL(fn,prec,1,benthicLayerN,2,myIter,myThid)
C      CALL MDSREADFIELD(fn,prec,'RL',1,benthicLayerC,3,myThid)
       CALL READ_REC_3D_RL(fn,prec,1,benthicLayerC,3,myIter,myThid)
CCV#ifdef ALLOW_RECOM_SILICATE
C      CALL MDSREADFIELD(fn,prec,'RL',1,benthicLayerSi,4,myThid)
       CALL READ_REC_3D_RL(fn,prec,1,benthicLayerSi,4,myIter,myThid)
CCV#endif /* ALLOW_RECOM_SILICATE */
C      CALL MDSREADFIELD(fn,prec,'RL',1,benthicLayerCaCO3,5,myThid)
       CALL READ_REC_3D_RL(fn,prec,1,benthicLayerCaCO3,5,myIter,myThid)

#endif /* ALLOW_RECOM */
      RETURN
      END
