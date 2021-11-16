C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/monod/monod_acdom.F,v 1.3 2013/04/16 20:21:57 jahn Exp $
C $Name:  $

#include "RECOM_OPTIONS.h"
      
CBOP
C !ROUTINE: MONOD_ACDOM

C !INTERFACE: ==========================================================
      SUBROUTINE RECOM_APHYTO_PHY(
     I                       D1,QYm,Drel,
     I                       aphy_PSpigm,aphy_ALLpigm,
     O                       aphytolocal,
     I                       myThid )

C !DESCRIPTION:
C     computes aphyto from D1 state variable

C !USES: ===============================================================
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "SPECTRAL_SIZE.h"
#ifdef RECOM_WAVEBANDS
#include "WAVEBANDS_PARAMS.h"
C#include "RECOM_PARAMS.h"
#endif
C !INPUT PARAMETERS: ===================================================
C     D1         :: Fraction of PSII active for photochemistry
C     QYm        :: Maximum quantum yield of photochemistry
C     Drel       :: Minimum D1 that keeps QY=QYm
C     aphy_chl_ps:: absorption spectrum photosynthetic pigments (ap_ps, m2 mgChla-1)
C     aphy_chl   :: absorption spectrum all pigments (ap, m2 mgChla-1)  
C     myTime     :: Current time in simulation
C     myIter     :: Current iteration number in simulation
C     myThid     :: My Thread Id number
C !OUTPUT PARAMETERS: ==================================================
C     aphytolocal :: absortpion spectra for phytoplankton group per level
    
      _RL D1(Nr)
      _RL QYm
      _RL Drel
      _RL aphy_PSpigm(tlam)
      _RL aphy_ALLpigm(tlam)
      _RL aphytolocal(Nr,tlam)
C     _RL     myTime
C     INTEGER myIter
      INTEGER myThid
CEOP

#ifdef RECOM_CALC_APHYT
C     !LOCAL VARIABLES:
C     == Local variables ==
      INTEGER k, ilam
      _RL QY, NPQ, astar
      _RL ALPHAmar, PPC
      _RL slope1, slope2
      
C      Drel = 0.75d0
      do k = 1,Nr
C Compute alpha from D1
         if (D1(k) .gt. 1.0d0) then
            QY = QYm
         else   
            QY = min(((QYm/Drel) * D1(k)) ,  QYm)
         endif
         NPQ = 1-(QY/QYm)
         astar = 0.d0
       do ilam = 1,tlam
           astar = astar + (wb_width(ilam) * aphy_PSpigm(ilam))
       end do
         astar = astar / wb_totalWidth   
         ALPHAmar = astar * (1-NPQ) * QY * 86400

C Use alpha as proxy for non-photosynthetic pigments ratio
         PPC = 1-(ALPHAmar / (astar * QYm * 86400) )

C Compute slopes from PPC
C Slope 488-532nm, (Eisner et al 2003 L&O 48) linear model with NOMAD data
         slope1 = (-1.93d-3 * PPC) - 1.011e-4
 
C Modify absorption spectra with slopes
C TO DO: make tlam indexes flexible, find wbs closest to 488 and 532        
      do ilam = 1,tlam
         aphytolocal(k,ilam) = aphy_PSpigm(ilam)
      enddo !l
      aphytolocal(k,5) = aphytolocal(k,6) + ((-slope1)*(525-500))
      aphytolocal(k,4) = aphytolocal(k,6) + ((-slope1)*(525-475))
      aphytolocal(k,3) = aphytolocal(k,3)
     &                    + (aphytolocal(k,4) - aphy_PSpigm(4))
      aphytolocal(k,2) = aphytolocal(k,2)
     &                    + (aphytolocal(k,4) - aphy_PSpigm(4))    
      aphytolocal(k,1) = aphytolocal(k,1)
     &                    + (aphytolocal(k,4) - aphy_PSpigm(4))     
      enddo  !k
#endif /* RECOM_CALC_APHYT */
      RETURN
      END

