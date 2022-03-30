C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/darwin/wavebands_init_vari.F,v 1.2 2013/12/04 21:19:12 jahn Exp $
C $Name:  $

c wavebands_init_vari.F assigns actual values to phyto types (two separate)
c wavebands_init_vari.F creates output files

C#include "DARWIN_OPTIONS.h"
#include "RECOM_OPTIONS.h"      

CBOP
C     !ROUTINE: WAVEBANDS_INIT_VARI
C     !INTERFACE:
       SUBROUTINE WAVEBANDS_INIT_VARI(myThid)

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE WAVEBANDS_INIT_VARI
C     | o assigns actual values to choices made by coin-flips
C     |   in darwin_generate_phyto.F; creates output files
C     *==========================================================*
C     \ev

C     !USES:
       implicit none
C     == Global variables ===
#include "SIZE.h"
#include "SPECTRAL_SIZE.h"
#ifdef RECOM_WAVEBANDS
#include "WAVEBANDS_PARAMS.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     myThid     :: my Thread Id number
       INTEGER myThid
CEOP

#ifdef RECOM_WAVEBANDS

C     !FUNCTIONS:
C      LOGICAL MASTER_CPU_THREAD
C      EXTERNAL MASTER_CPU_THREAD

C     !LOCAL VARIABLES:
C     == Local variables ==
C       INTEGER npmax
C       _RL outfile(2,tlam)
C       _RL outfile_ps(2,tlam)

       INTEGER np,nl,i,ilam, nap
C       INTEGER IniUnit1, IniUnit2, IniUnit3
C       INTEGER IniUnit4, IniUnit5, IniUnit6
C       INTEGER IniUnit7

CEA Here groups=types and datafile has 1=small phyto, 2=diatoms (keep the order)       
          do nap=1, tnabp
              ap_type(nap) = nap
          enddo
       
C       npmax = 2
C          do nap=1, tnabp
C          if (ap_type(np).eq.1) then
           do i = 1,tlam
              aphy_chl(i) = ap(1,i)
              aphy_chl_ps(i) = ap_ps(1,i)
#ifdef RECOM_RADTRANS
              bphy_chl(i) = bp(1,i)
#ifdef DAR_NONSPECTRAL_BACKSCATTERING_RATIO
              bbphy_chl(i) = bp(1,i)*darwin_bbphy(1)
#else
              bbphy_chl(i) = bbp(1,i)
#endif
#endif
C           enddo
C          endif
        
C          if (ap_type(np).eq.2) then
C           do i = 1,tlam
              aphy_chl_dia(i) = ap(2,i)
              aphy_chl_ps_dia(i) = ap_ps(2,i)
#ifdef RECOM_RADTRANS
              bphy_chl_dia(i) = bp(2,i)
#ifdef DAR_NONSPECTRAL_BACKSCATTERING_RATIO
              bbphy_chl_dia(i) = bp(2,i) * darwin_bbphy(2)
#else
              bbphy_chl_dia(i) = bbp(2,i)
#endif
#endif
           enddo
C          endif  
C         enddo
       
          if(ap_type(1).eq.0) then
           do i=1,tlam
           aphy_chl(i) = 9.9d2
           aphy_chl_ps(i) = 9.9d2
#ifdef RECOM_RADTRANS
           bphy_chl(i) = 9.9d2
           bbphy_chl(i) = 9.9d2
#endif
           enddo
          endif
          if(ap_type(2).eq.0) then
           do i=1,tlam
           aphy_chl_dia(i) = 8.9d2
           aphy_chl_ps_dia(i) = 8.9d2
#ifdef RECOM_RADTRANS
           bphy_chl_dia(i) = 8.9d2
           bbphy_chl_dia(i) = 8.9d2
#endif
           enddo
          endif
#endif /* RECOM_WAVEBANDS */
         return
         end
