C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/monod/monod_acdom.F,v 1.3 2013/04/16 20:21:57 jahn Exp $
C $Name:  $

C#include "DARWIN_OPTIONS.h"
#include "RECOM_OPTIONS.h"
      
CBOP
C !ROUTINE: RTWB_ACDOM

c Only used if RECOM_CALC_ACDOM (if not, uses the file)
C !INTERFACE: ==========================================================
      SUBROUTINE MONOD_ACDOM(
#ifdef RECOM_CDOM
     I                       cdomlocal,
#else
     I                       P_chl, aclocal, awlocal,
#endif
     O                       acdomlocal,
     I                       myThid )

C !DESCRIPTION:
C     computes CDOM from Chlorophyll
C
c AH changed this a bit from WG code
c ------------- this is WG code:
c      do k = 1,km
c       actot450 = 0.0
c       atot450 = 0.0
c       do n = 1,nchl
c        actot450 = actot450  + P(k,n)*ac(n,nl450)
c       enddo
c       atot450 = aw(nl450) + actot450
c       do nl = 1,nlt
c        acdom(k,nl) = 0.2*atot450*excdom(nl)
c       enddo
c      enddo
c      do nl = 1,nlt
c       write(6,*)'nl,lam,aw,ac,acdom = ',nl,lam(nl),
c    *aw(nl),ac(1,nl),acdom(1,nl)
c      enddo
c this implies P biomass is being scaled by the absorption of phytoplankton
c at 450nm. ie. the waveband on the main chl-a peak.
c so (i think) this is effectively getting the absorption by the phytoplankton
c In our version, since we have phyto chl we need to scale the 
c chlorophyll normalised absorption spectra (the inputs)
c by phyt chl rather than phyto biomass 

C !USES: ===============================================================
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "SPECTRAL_SIZE.h"
#ifdef RECOM_WAVEBANDS
#include "WAVEBANDS_PARAMS.h"      
#endif

C !INPUT PARAMETERS: ===================================================
C     P_chl      :: Chlorophyll per species and level
c     ac         :: absorption spectra for each phtyo (chl normalised)
c     awlocal    :: absoprtion spectra for water (m-1)
C     myTime     :: Current time in simulation
C     myIter     :: Current iteration number in simulation
C     myThid     :: My Thread Id number

C !OUTPUT PARAMETERS: ==================================================
C     acdom      :: absortpion spectra for CDOM per level

#ifdef RECOM_CDOM
      _RL cdomlocal(NR)
#else
      _RL P_chl(2,Nr)
      _RL aclocal(2,tlam)
      _RL awlocal(tlam)
#endif   
      _RL acdomlocal(Nr,tlam)
C      _RL     myTime
C      INTEGER myIter
      INTEGER myThid
CEOP

C#ifdef RECOM_WAVEBANDS
#ifdef RECOM_CALC_ACDOM

C     !LOCAL VARIABLES:
C     == Local variables ==
      INTEGER npmax
      INTEGER k, np, ilam
      _RL actot450, atot450
c params in common block
c     nl450, excdom, 

#ifdef RECOM_CDOM
c use cdom-like tracer
      do k = 1,Nr
       do ilam = 1,tlam
        acdomlocal(k,ilam) = cdomcoeff*cdomlocal(k)*excdom(ilam)
       enddo
      enddo
#else
      npmax = 2      
      do k = 1,Nr
       actot450 = 0.0 _d 0
       atot450 = 0.0 _d 0
       do np = 1,npmax
        actot450 = actot450  + P_chl(np,k) * aclocal(np,nlaCDOM) !nb. n,k swapped from WG
       enddo
       atot450 = awlocal(nlaCDOM) + actot450
       do ilam = 1,tlam
        acdomlocal(k,ilam) = darwin_aCDOM_fac * atot450 * excdom(ilam)
       enddo
      enddo
c      do nl = 1,tlam
c       write(6,*)'nl,lam,aw,ac,acdom = ',nl,lam(nl),
c    *aw(nl),ac(1,nl),acdom(1,nl)
c      enddo
#endif /* RECOM_CDOM */
#endif /* RECOM_CALC_ACDOM */
C#endif /* RECOM_WAVEBANDS */
      RETURN
      END

