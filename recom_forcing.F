C $Header: /csys/software/cvsroot/REcoM/MITgcm/recom/recom_forcing.F,v 1.49 2008/04/17 14:25:52 mlosch Exp $
C $Name:  $

#include "GCHEM_OPTIONS.h"
#include "RECOM_OPTIONS.h"
#ifdef ALLOW_EXF
# include "EXF_OPTIONS.h"
#endif /* ALLOW_EXF */

CBOP
C     !ROUTINE: RECOM_FORCING
C     !INTERFACE:
      SUBROUTINE RECOM_FORCING( 
     I     iMin, iMax, jMin, jMax, bi, bj, myTime, myIter, myThid )
C     !DESCRIPTION: \bv
C     *==================================================================
C     | SUBROUTINE RECOM_FORCING
C     | o Calling interface for REcoM
C     *==================================================================
C     |
C     | RECOM_FORCING is the interface for calling a the bio-geo-chemical
C     | REcoM model (Markus Schartau)
C     | that provides source and sink terms for the passive tracers.
C     | the ptracers-package must be enabled for this interface to work
C     |
C     | RECOM_FORCING is called from GCHEM_FORCING, analogously to
C     | Stephanie Dutkiewicz's dic_forcing.
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
#include "DYNVARS.h"
#ifdef ALLOW_PTRACERS
# include "PTRACERS_SIZE.h"
# include "PTRACERS_PARAMS.h"
# include "PTRACERS_FIELDS.h"
# ifdef ALLOW_GCHEM
#  include "GCHEM.h"
# endif
# ifdef ALLOW_RECOM
#  include "RECOM.h"
#  include "RECOM_PARAMS.h"      
#  ifdef RECOM_WAVEBANDS
#  include "SPECTRAL_SIZE.h"
#  include "SPECTRAL.h"
#  include "WAVEBANDS_PARAMS.h"
#  endif      
# endif /* ALLOW_RECOM */
#endif /* ALLOW_PTRACERS */
#ifdef ALLOW_SEAICE
# include "SEAICE_SIZE.h"
# include "SEAICE_PARAMS.h"
# include "SEAICE.h"
#endif /* ALLOW_SEAICE */
#ifdef ALLOW_THSICE
# include "THSICE_VARS.h"
#endif /* ALLOW_THSICE */
#ifdef ALLOW_EXF
# include "EXF_FIELDS.h"
#endif /* ALLOW_EXF */
#ifdef ALLOW_LONGSTEP
# include "LONGSTEP.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     note: under the multi-threaded model myiter and 
C           mytime are local variables passed around as routine 
C           arguments. Although this is fiddly it saves the need to 
C           impose additional synchronisation points when they are 
C           updated.
C     iMin, iMax, jMin, jMax - loop ranges
C     bi, bj - tile indices
C     myiter - iteration counter for this thread
C     mytime - time counter for this thread
C     mythid - thread number for this instance of the routine.
      INTEGER iMin, iMax, jMin, jMax
      integer bi, bj
      integer mythid 
      integer myiter
      _RL     myTime

#ifdef ALLOW_RECOM
C     !LOCAL VARIABLES :
C     == Local variables ==
C     i,j,k,iTracer - loop indices
C     kSurface      - surface index
C     krLoc         - local vertical loop boundary klowC
C     recom_dic     - dissolved inorganic carbon (intermediate variable)
C     recom_alk     - total alkalinity (intermediate variable)

C     benthicLayerOld/New - auxillary varialbes
C     angleOfInclination - of incoming solar irradiation
C     daylen        - length of the day (function of latitude and time of year)
C     radForce      - short wave radiation (needs to be replaced by a proper
C                     field !
      INTEGER i,j,k,iTracer,krLoc, kSurface
CCV      INTEGER iDIC, iALK, iFe

      _RL     recom_dic   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     recom_alk   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     dicFluxLoc  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     hPlusLoc    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     windSpeedLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     atmospCO2Loc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     saltLoc     (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     thetaLoc    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     co2_ppm     (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     co2_mol     (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#ifdef RECOM_O2
      _RL     recom_oxy   (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     FluxO2Loc  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#endif
      _RL     benthicLayerloc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,1:4)
C     these need to get the proper dimension ... etc.
      _RS     cosAngleOfIncidence(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RS     lengthOfDay(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RS     PARadiation(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#ifdef RECOM_RADTRANS      
      _RS     solz(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#endif      
      _RS     feDustFluxLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy) 
      _RS     surfpco2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     resuspension variables
      _RL     bottomStress
      _RL     erosionCoeff
      _RL     erosionFlux(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
CMS   implementing LithPartFlux
#ifdef RECOM_IRON_LITHOGENIC_PARTICLES
      _RS     ironLithParticleFluxLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
#endif  /* RECOM_IRON_LITHOGENIC_PARTICLES */
C     diagnostics auxilliary fields
      INTEGER idiag, ndiags3d_used
      _RL     recom_diag_integral(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &                          nDiagsREcoM3D)
      _RL     recom_diag_3Dfields(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &                          Nr,nDiagsREcoM3D)
      _RL     recom_diag_2Dfields(1-OLx:sNx+OLx,1-OLy:sNy+OLy,
     &                          nDiagsREcoM2D)
C    Local gptr :: work-space for time-stepping
      _RL     gPtrLoc(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,
     &               PTRACERS_num)
CEOP

#ifdef RECOM_WAVEBANDS
      integer ilam
       _RL PARwup(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam)
       _RL PARwdn(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam)       
       _RL PARw_k(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam, Nr)       
       _RL PARl(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL PURl(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)       
       _RL PARwup_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam)
       _RL PARwup_total(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
       _RL Edwsf(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam)
       _RL Eswsf(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam)
#ifdef RECOM_CALC_REFLEC
       INTEGER index
       _RL PARw_kwb(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
#endif       
       _RL C_phot_nl(tlam)
       _RL C_phot_nl_dia(tlam)
       _RL Ek_nl(tlam)
       _RL Ek_nl_dia(tlam)
       integer iday,iyr,imon,isec,lp,wd,mydate(4)
       _RL a_k(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr,tlam)
       _RL acdom_k(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr, tlam)
       _RL aphy_chl_k(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr, tlam)
       _RL aphy_chl_dia_k(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr, tlam)       
       _RL apart_k(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr, tlam)    
       _RL actot(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr, tlam)  
       _RL aclocal(2,tlam)    
       _RL atten 
       _RL discEs, discEu
       _RL dz_k(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)      
       INTEGER idiscEs,jdiscEs,kdiscEs,ldiscEs
       INTEGER idiscEu,jdiscEu,kdiscEu,ldiscEu
CCEA    BIOMASS
       _RL part_k(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL Phy_k(1-OLx:sNx+OLx,1-OLy:sNy+OLy, 2, Nr)
       _RL phychl_k(1-OLx:sNx+OLx,1-OLy:sNy+OLy, 2, Nr)             
#ifdef RECOM_CDOM
       _RL cdom_k(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
#endif
#if defined(RECOM_CALC_APHYT) && defined(RECOM_MARSHALL)
       _RL phyD1_k(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL diaD1_k(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)       
#endif       
CCEA for diagnostics
       _RL a_kave(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL acdom_kave(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL apart_kave(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL actot_ave(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
c  iops surface
       _RL a_ksur(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam)
       _RL acdom_ksur(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam)
       _RL apart_ksur(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam) 
       _RL actot_sur(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam)     
#ifdef RECOM_CALC_REFLEC
       _RL a_kwb(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL acdom_kwb(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL apart_kwb(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL actot_wb(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)   
#endif
#ifdef RECOM_RADTRANS   
       _RL bt_k(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr, tlam)
       _RL bb_k(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr, tlam)       
       _RL bpart_k(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr, tlam)
       _RL bbpart_k(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr,tlam) 
       _RL bctot(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr, tlam)
       _RL bbctot(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr, tlam)          
       _RL Edz(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam, Nr)
       _RL Esz(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam, Nr)
       _RL Euz(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam, Nr)
       _RL Estop(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam, Nr)
       _RL Eutop(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam, Nr)
       _RL tirrq(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL tirrwq(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam, Nr)
       _RL amp1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam, Nr)
       _RL amp2(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam, Nr)
       _RL rmud(1-OLx:sNx+OLx,1-OLy:sNy+OLy)      
       _RL sinszaw, szaw, rmudl, rn
       data rn /1.341 _d 0/     !refractive index of seawater
c  only for diagnostics
       _RL bt_kave(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL bpart_kave(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL bctot_ave(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr) 
       _RL bb_kave(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL bbpart_kave(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL bbctot_ave(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr) 
       _RL Eupwel(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam)       
       _RL Reflec(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam)
c  iops surface
       _RL bt_ksur(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam)
       _RL bpart_ksur(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam)
       _RL bctot_sur(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam) 
       _RL bb_ksur(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam)
       _RL bbpart_ksur(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam)
       _RL bbctot_sur(1-OLx:sNx+OLx,1-OLy:sNy+OLy, tlam)      
#ifdef RECOM_CALC_REFLEC
       _RL bt_kwb(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL bpart_kwb(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL bctot_wb(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr) 
       _RL bb_kwb(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL bbpart_kwb(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL bbctot_wb(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL Edz_wb(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL Esz_wb(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL Euz_wb(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL Estop_wb(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL Eutop_wb(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL amp1_wb(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
       _RL amp2_wb(1-OLx:sNx+OLx,1-OLy:sNy+OLy, Nr)
#endif    
#endif /* RADTRANS */         
#endif /* RECOM_WAVEBANDS */

CCV: tracer indices are now set in RECOM.h
CCV      iDIC = 2
CCV      iALK = 3
CCV      iFe  = PTRACERS_numInUse
C     surface layer index
      ksurface = 1
C     resuspension
      erosionCoeff = recom_FeErosionRate*(1.-recom_porosity)
     &     * recom_sizeFraction / recomCritShearStress

C     initialise some local fields
      DO j=jMin,jMax
       DO i=iMin,iMax
        feDustFluxLoc(i,j) = 0. _d 0
        erosionFlux  (i,j) = 0. _d 0
CCMS   initialise LithParticleFlux too
#ifdef RECOM_IRON_LITHOGENIC_PARTICLES
        ironLithParticleFluxLoc(i,j) = 0. _d 0
#endif /* RECOM_IRON_LITHOGENIC_PARTICLES */
       ENDDO
      ENDDO
      DO idiag = 1,nDiagsREcoM3D
       DO k=1,Nr
        DO j=jMin,jMax
         DO i=iMin,iMax
          recom_diag_3Dfields(i,j,k,idiag) = 0. _d 0
         ENDDO 
        ENDDO
       ENDDO
       DO j=jMin,jMax
        DO i=iMin,iMax
         recom_diag_integral(i,j,idiag) = 0. _d 0
        ENDDO 
       ENDDO
      ENDDO

C     solar radiation
      CALL RECOM_INSOLATION(
     O     cosAngleOfIncidence, lengthOfDay, PARadiation,
#ifdef RECOM_RADTRANS
     O     solz,
#endif      
     I     iMin, iMax, jMin, jMax, bi, bj, myTime, myIter, myThid )

      
#ifdef RECOM_WAVEBANDS
       idiscEs = 0
       jdiscEs = 0
       kdiscEs = 0
       ldiscEs = 0
       idiscEu = 0
       jdiscEu = 0
       kdiscEu = 0
       ldiscEu = 0
       discEs = 0.
       discEu = 0.
   
C ========================== k,j,i loop =================================

       DO k=1,Nr 
        DO j=jMin,jMax
         DO i=iMin,iMax
             
C ------------ GET BIOMASSES  ----------------
        dz_k(i,j,k) = drF(k)*hFacC(i,j,k,bi,bj) 
        part_k(i,j,k) = MAX(pTracer(i,j,k,bi,bj,idetc), recom_tiny)     !mmolC           
#ifdef RECOM_CDOM
        cdom_k(i,j,k) = MAX(pTracer(i,j,k,bi,bj,icdom), recom_tiny)     !mmolC
#endif
        Phy_k(i,j,1,k) = MAX(pTracer(i,j,k,bi,bj,iphyc), recom_tiny)    !mmolC
        phychl_k(i,j,1,k) = MAX(pTracer(i,j,k,bi,bj,ipchl), recom_tiny) !mg Chla
        Phy_k(i,j,2,k) = MAX(pTracer(i,j,k,bi,bj,idiac), recom_tiny)    !mmolC
        phychl_k(i,j,2,k) = MAX(pTracer(i,j,k,bi,bj,idchl), recom_tiny) !mg Chla
#if defined(RECOM_CALC_APHYT) && defined(RECOM_MARSHALL)
        phyD1_k(i,j,k) = MAX(pTracer(i,j,k,bi,bj,id1), recom_tiny)      !rel
        diaD1_k(i,j,k) = MAX(pTracer(i,j,k,bi,bj,id1d), recom_tiny)     !rel
#endif

C ------ GET constant acdom_k -------
#ifndef RECOM_CALC_ACDOM
       do ilam = 1,tlam
        acdom_k(i,j,k,ilam) = acdom(ilam)
       enddo         
#endif
C ------ GET constant aphy_chl_k & aphy_chl_dia_k -------
#ifndef RECOM_CALC_APHYT
       do ilam = 1,tlam
        aphy_chl_k(i,j,k,ilam) = aphy_chl(ilam)
        aphy_chl_dia_k(i,j,k,ilam) = aphy_chl_dia(ilam)
       enddo 
#endif
#ifndef RECOM_2GROUPS
       do ilam = 1,tlam
        aphy_chl_dia_k(i,j,k,ilam) = aphy_chl_dia(ilam)
       enddo 
#endif
         ENDDO !i
       ENDDO  !j
      ENDDO  !k

     
C ========================== j,i loops =================================
#ifdef RECOM_CALC_ACDOM     
        DO j=jMin,jMax
         DO i=iMin,iMax    
C ------------ COMPUTE ACDOM_k  ----------------------------------
#ifdef RECOM_CDOM
        call MONOD_ACDOM(cdom_k(i,j,1:Nr),
     O                   acdom_k(i,j, 1:Nr, 1:tlam),
     I                   myThid)        
#else
        aclocal(1,1:tlam)=aphy_chl
        aclocal(2,1:tlam)=aphy_chl_dia        
         call MONOD_ACDOM(phychl_k(i,j, 2, 1:Nr), aclocal, aw,
     O                    acdom_k(i,j, 1:Nr, 1:tlam),
     I                    myThid)
#endif /* RECOM_CDOM */
        ENDDO
       ENDDO
#endif /* RECOM_CALC_ACDOM */

C ------------ COMPUTE aphy_chl_k & aphy_chl_dia_k ----------------
C Slope 488-532nm, (Eisner et al 2003 L&O 48) linear model 
C with NOMAD data <50% fDiatoms and >50% fDiatoms
       
#if defined(RECOM_CALC_APHYT) && defined(RECOM_MARSHALL)
        DO j=jMin,jMax
         DO i=iMin,iMax
C        call RECOM_APHYTO_SLOPES(phyD1_k(i,j,1:Nr),QYmax,Drel,
C     I                    aphy_chl_ps(1:tlam),aphy_chl(1:tlam),
C     I                    -1.9288E-3, -1.011E-4,
C     O                    aphy_chl_k(i,j, 1:Nr, 1:tlam),
C     I                    myThid)

        call RECOM_APHYTO_PHY(phyD1_k(i,j,1:Nr),QYmax,Drel,
     I                    aphy_chl_ps(1:tlam),aphy_chl(1:tlam),
     O                    aphy_chl_k(i,j, 1:Nr, 1:tlam),
     I                    myThid)             

#ifdef RECOM_2GROUPS             
C        call RECOM_APHYTO_SLOPES(diaD1_k(i,j,1:Nr),QYmax_d,Drel,
C     I                    aphy_chl_ps_dia(1:tlam),aphy_chl_dia(1:tlam),
C     I                    -1.418E-3, -5.8E-5,
C     O                    aphy_chl_dia_k(i,j, 1:Nr, 1:tlam),
C     I                    myThid)

        call RECOM_APHYTO_DIA(diaD1_k(i,j,1:Nr),QYmax_d,Drel,
     I                    aphy_chl_ps_dia(1:tlam),aphy_chl_dia(1:tlam),
     O                    aphy_chl_dia_k(i,j, 1:Nr, 1:tlam),
     I                    myThid)

#endif /* RECOM_2GROUPS */
        ENDDO
       ENDDO
#endif /* RECOM_CALC_APHYT */
       
C========================== l,k,j,i loop =================================
C
C ------------ GET PART_k FOR WAVEBANDS_3D and RADTRANS ----------------
C     In Darwin-MONOD particulate matter is calculated in P units.
C     We use either detC and a C per particle factor (Stramski 2001) or
C     detC with a biomass-specific spectrum (Gallegos 2011).     
 
      DO ilam = 1,tlam    
       DO k=1,Nr
        DO j=jMin,jMax
         DO i=iMin,iMax
             
#ifdef RECOM_CALC_APART
      apart_k(i,j,k,ilam) = part_k(i,j,k) * aparcoeff * exapar(ilam)
#ifdef RECOM_RADTRANS
      bpart_k(i,j,k,ilam) = part_k(i,j,k) * bparcoeff * exbpar(ilam)
      bbpart_k(i,j,k,ilam) = part_k(i,j,k) * bparcoeff
     &                     * exbpar(ilam) * bb_to_b      
#endif
#else
      apart_k(i,j,k,ilam) = part_k(i,j,k)*apart_P(ilam)
#ifdef RECOM_RADTRANS
      bpart_k(i,j,k,ilam) = part_k(i,j,k)*bpart_P(ilam)
      bbpart_k(i,j,k,ilam) = part_k(i,j,k)*bbpart_P(ilam)      
#endif
#endif /* RECOM_CALC_APART*/
          ENDDO !i
         ENDDO  !j      
        ENDDO   !k
       ENDDO    !l
       
C ------------- GET SPECTRAL LIGHT -----------------
       do ilam = 1,tlam
        DO j=jMin,jMax
         DO i=iMin,iMax      
#ifdef OASIM
c add direct and diffuse, convert to uEin/m2/s/nm             
          Edwsf(i,j,ilam) = oasim_ed(i,j,ilam,bi,bj) 
          Eswsf(i,j,ilam) = oasim_es(i,j,ilam,bi,bj)         
          PARwup(i,j,ilam) = WtouEins(ilam) * (Edwsf(i,j,ilam)
     &                          + Eswsf(i,j,ilam))
#else
c sf is per nm; convert to per waveband
          Edwsf(i,j,ilam) = wb_width(ilam) * sf(ilam)
     &                    * PARadiation(i,j)
          Eswsf(i,j,ilam) = tiny
          PARwup(i,j,ilam) = Edwsf(i,j,ilam) * WtouEins(ilam)
#endif
          PARwup_diag(i,j,ilam) = PARwup(i,j,ilam)          
          ENDDO !i
         ENDDO  !j    
        enddo   !l

       DO j=jMin,jMax
        DO i=iMin,iMax
          PARwup_total(i,j) = 0.
         do ilam = 1,tlam
          PARwup_total(i,j) = PARwup_total(i,j) + PARwup(i,j,ilam)
         enddo !l
        ENDDO  !i
       ENDDO   !j

#ifdef RECOM_RADTRANS         
c     Compute 1/cos(zenith) for direct light below surface given solar zenith
c     angle (in radians) at surface (solz=zenith_deg from RECOM_INSOLATION.F)       

       DO j=jMin,jMax
        DO i=iMin,iMax  
         sinszaw = sin(solz(i,j))/rn
         szaw = asin(sinszaw)
         rmudl = 1.0/cos(szaw)    !avg cosine direct (1 over)
         rmud(i,j) = min(rmudl,1.5)
         rmud(i,j) = max(rmud(i,j),0.0)         
       ENDDO !i
      ENDDO  !j
#endif   
   
       
#ifndef RECOM_RADTRANS
c ------------ WAVEBANDS W/O RADTRANS ----------------------------------

         do k=1,Nr              
           do ilam = 1,tlam
             DO j=jMin,jMax
               DO i=iMin,iMax             
               if (hFacC(i,j,k,bi,bj).gt.0. _d 0) then
c get total attenuation (absorption) by phyto at each wavelength
              actot(i,j,k,ilam) = 0.
              actot(i,j,k,ilam) = actot(i,j,k,ilam)
     &                      + (phychl_k(i,j,1,k)*aphy_chl_k(i,j,k,ilam))
     &                  + (phychl_k(i,j,2,k)*aphy_chl_dia_k(i,j,k,ilam))    
              a_k(i,j,k,ilam) = aw(ilam)
     &                      + actot(i,j,k,ilam)
     &                      + acdom_k(i,j,k,ilam)
     &                      + apart_k(i,j,k,ilam)
              atten = a_k(i,j,k,ilam) * drF(k)    ! why not dz_k ?
              PARwdn(i,j,ilam) = PARwup(i,j,ilam)*exp(-atten)
              endif
             ENDDO !i
            ENDDO !j             
           enddo !l
           
c find for the midpoint of the gridcell (gridcell mean)
           do ilam = 1,tlam
             DO j=jMin,jMax
                DO i=iMin,iMax
                   if (hFacC(i,j,k,bi,bj).gt.0. _d 0) then                   
                   PARw_k(i,j,ilam,k)=sqrt(PARwup(i,j,ilam)
     &                     * PARwdn(i,j,ilam))
                   else
                   PARw_k(i,j,ilam,k) = 0. _d 0                         
                   endif    
              ENDDO !i
             ENDDO  !j            
            enddo   !l

c cycle
           do ilam=1,tlam
              DO j=jMin,jMax
               DO i=iMin,iMax
                if (hFacC(i,j,k,bi,bj).gt.0. _d 0) then
                   PARwup(i,j,ilam) = PARwdn(i,j,ilam)
                endif    
               ENDDO !i
              ENDDO  !j      
             enddo   !j

c sum wavebands for total PAR at the mid point of the gridcell (PARl)
CCEA sum wavebands for total PAR absorbed by phytoplankton (PURl)
        DO j=jMin,jMax
          DO i=iMin,iMax
             PARl(i,j,k) = 0.
             PURl(i,j,k) = 0.
            do ilam = 1,tlam            
               PARl(i,j,k) = PARl(i,j,k) + PARw_k(i,j,ilam,k)
               PURl(i,j,k) = PURl(i,j,k) + (PARw_k(i,j,ilam,k)
     &                                    * actot(i,j,k,ilam))  
            enddo !l
         ENDDO    !i
       ENDDO      !j          


CCEA compute averages and wb=index for exporting
        DO j=jMin,jMax
          DO i=iMin,iMax
             a_kave(i,j,k) = 0.d0
             acdom_kave(i,j,k) = 0.d0
             apart_kave(i,j,k) = 0.d0
             actot_ave(i,j,k) = 0.d0 
              do ilam = 1,tlam  
                 a_kave(i,j,k) = a_kave(i,j,k)
     &                         + wb_width(ilam) * a_k(i,j,k,ilam)
                 acdom_kave(i,j,k) = acdom_kave(i,j,k)
     &                     + wb_width(ilam) * acdom_k(i,j,k,ilam)
                 apart_kave(i,j,k) = apart_kave(i,j,k)
     &                     + wb_width(ilam) * apart_k(i,j,k,ilam)
                 actot_ave(i,j,k) = actot_ave(i,j,k)
     &                       + wb_width(ilam) * actot(i,j,k,ilam)
              enddo   !l
             a_kave(i,j,k) = a_kave(i,j,k) / wb_totalWidth
             acdom_kave(i,j,k) = acdom_kave(i,j,k) / wb_totalWidth
             apart_kave(i,j,k) = apart_kave(i,j,k) / wb_totalWidth
             actot_ave(i,j,k) = actot_ave(i,j,k) / wb_totalWidth
#ifdef RECOM_CALC_REFLEC
             index = darwin_diag_acdom_ilam
             PARw_kwb(i,j,k) = PARw_k(i,j,index,k)  
             a_kwb(i,j,k) = a_k(i,j,k,index) 
             acdom_kwb(i,j,k) = acdom_k(i,j,k,index)
             apart_kwb(i,j,k) = apart_k(i,j,k,index) 
             actot_wb(i,j,k) = actot(i,j,k,index)
#endif            
          ENDDO    !i
       ENDDO      !j
      enddo       !k

C iops surface for diagnostics   
        DO ilam = 1,tlam
          DO j=jMin,jMax
            DO i=iMin,iMax
              IF( dz_k(i,j,kSurface) .GT. 0.0 )THEN
                IF(Eswsf(i,j,ilam).GE.tiny .OR.
     &             Edwsf(i,j,ilam).GE.tiny ) THEN
             a_ksur(i,j,ilam) = a_k(i,j,kSurface,ilam) 
             acdom_ksur(i,j,ilam) = acdom_k(i,j,kSurface,ilam)
             apart_ksur(i,j,ilam) = apart_k(i,j,kSurface,ilam) 
             actot_sur(i,j,ilam) = actot(i,j,kSurface,ilam)           
                ENDIF !light
              ENDIF !depth         
            ENDDO !i
           ENDDO !j
        ENDDO  !l
        
#else /* RECOM_RADTRANS */

c ------------ FULL RADIATIVE TRANSFER CODE ----------------------------
CEA Direct and difusse irradiance: both with OASIM, only Ed without
C           
C Compute total absorption/scattering coefficients           
      DO k=1,Nr
       DO ilam = 1,tlam
        DO j=jMin,jMax
         DO i=iMin,iMax                    
c   absorption by phyto: if RADTRANS actot is computed twice
            actot(i,j,k,ilam) = 0.0
            bctot(i,j,k,ilam) = 0.0
            bbctot(i,j,k,ilam) = 0.0
C            DO np = 1,npmax
            actot(i,j,k,ilam)  = actot(i,j,k,ilam)
     &                      + phychl_k(i,j,1,k) * aphy_chl_k(i,j,k,ilam)
     &                  + phychl_k(i,j,2,k) * aphy_chl_dia_k(i,j,k,ilam)
#ifdef RECOM_BMASS
            bctot(i,j,k,ilam)  = bctot(i,j,k,ilam)
     &                        + Phy_k(i,j,1,k) * bphy_chl(ilam)
     &                    + Phy_k(i,j,2,k) * bphy_chl_dia(ilam)       
            bbctot(i,j,k,ilam) = bbctot(i,j,k,ilam)
     &                       + Phy_k(i,j,1,k) * bbphy_chl(ilam)
     &                   + Phy_k(i,j,2,k) * bbphy_chl_dia(ilam)
#else
            bctot(i,j,k,ilam)  = bctot(i,j,k,ilam)
     &                      + phychl_k(i,j,1,k) * bphy_chl(ilam)
     &                  + phychl_k(i,j,2,k) * bphy_chl_dia(ilam)       
            bbctot(i,j,k,ilam) = bbctot(i,j,k,ilam)
     &                      + phychl_k(i,j,1,k) * bbphy_chl(ilam)
     &                  + phychl_k(i,j,2,k) * bbphy_chl_dia(ilam)
#endif       
C            ENDDO
c   total: water, CDOM, phyto, particles
            a_k(i,j,k,ilam) = aw(ilam) + acdom_k(i,j,k,ilam)
     &                   + actot(i,j,k,ilam) + apart_k(i,j,k,ilam)
            bt_k(i,j,k,ilam) = bw(ilam)
     &                   + bctot(i,j,k,ilam) + bpart_k(i,j,k,ilam)
            bb_k(i,j,k,ilam) = darwin_bbw * bw(ilam)
     &                 + bbctot(i,j,k,ilam) + bbpart_k(i,j,k,ilam)
            bb_k(i,j,k,ilam) = MAX(darwin_bbmin, bb_k(i,j,k,ilam))
c   initialize output variables
            Edz(i,j,ilam,k) = 0.0
            Esz(i,j,ilam,k) = 0.0
            Euz(i,j,ilam,k) = 0.0
            Estop(i,j,ilam,k) = 0.0
            Eutop(i,j,ilam,k) = 0.0
            amp1(i,j,ilam,k) = 0.0
            amp2(i,j,ilam,k) = 0.0          
           ENDDO !i
         ENDDO   !j           
        ENDDO    !l
       ENDDO     !k

       
C ------ Propagate three-beam light in the water column -------  
        DO j=jMin,jMax
         DO i=iMin,iMax
CEA Some of the routines use drF and others dz_k, why?                         
         IF (darwin_radtrans_niter.GE.0) THEN
           call MONOD_RADTRANS_ITER(
     I               dz_k(i,j,1:Nr),rmud(i,j),
     I               Edwsf(i,j,1:tlam),
     I               Eswsf(i,j,1:tlam),
     I               a_k(i,j,1:Nr,1:tlam),
     I               bt_k(i,j,1:Nr,1:tlam),
     I               bb_k(i,j,1:Nr,1:tlam),          
     I               darwin_radtrans_kmax,darwin_radtrans_niter,
     O               Edz(i,j,1:tlam,1:Nr),
     O               Esz(i,j,1:tlam,1:Nr),
     O               Euz(i,j,1:tlam,1:Nr),
     O               Eutop(i,j,1:tlam,1:Nr),
     O               tirrq(i,j,1:Nr),
     O               tirrwq(i,j,1:tlam,1:Nr),
     O               amp1(i,j,1:tlam,1:Nr),amp2(i,j,1:tlam,1:Nr),         
     I               myThid)
         ELSEIF (darwin_radtrans_niter.EQ.-1) THEN
           call MONOD_RADTRANS(
     I               drF(1:Nr),rmud(i,j),
     I               Edwsf(i,j,1:tlam),Eswsf(i,j,1:tlam),
     I               a_k(i,j,1:Nr,1:tlam),
     I               bt_k(i,j,1:Nr,1:tlam),
     I               bb_k(i,j,1:Nr,1:tlam),           
     O               Edz(i,j,1:tlam,1:Nr),Esz(i,j,1:tlam,1:Nr),
     O               Euz(i,j,1:tlam,1:Nr),Eutop(i,j,1:tlam,1:Nr),
     O               tirrq(i,j,1:Nr),
     O               tirrwq(i,j,1:tlam,1:Nr),          
     I               myThid)
         ELSE
            call MONOD_RADTRANS_DIRECT(
     I               dz_k(i,j,1:Nr),rmud(i,j),
     I               Edwsf(i,j,1:tlam),Eswsf(i,j,1:tlam),
     I               a_k(i,j,1:Nr,1:tlam),
     I               bt_k(i,j,1:Nr,1:tlam),
     I               bb_k(i,j,1:Nr,1:tlam),
     I               darwin_radtrans_kmax,
     O               Edz(i,j,1:tlam,1:Nr),Esz(i,j,1:tlam,1:Nr),
     O               Euz(i,j,1:tlam,1:Nr),
     O               Estop(i,j,1:tlam,1:Nr),Eutop(i,j,1:tlam,1:Nr),
     O               tirrq(i,j,1:Nr),
     O               tirrwq(i,j,1:tlam,1:Nr),
     O               amp1(i,j,1:tlam,1:Nr),amp2(i,j,1:tlam,1:Nr),          
     I               myThid)
          
         ENDIF

      ENDDO  !i
      ENDDO  !j

      
c     Uses chl from prev timestep (as wavebands does) keep like this in case
c     need to consider upwelling irradiance as affecting the grid box above
c     Pass to sms: PARw_k only, but will be for this timestep for RADTRANST
c     and previous timestep for WAVEBANDS.
c
c     Now copy

       DO k=1,Nr
               DO j=jMin,jMax
                   DO i=iMin,iMax
                       PARl(i,j,k) = tirrq(i,j,k)
                   ENDDO
               ENDDO
           
           DO ilam = 1,tlam
               DO j=jMin,jMax
                   DO i=iMin,iMax
                       PARw_k(i,j,ilam,k) = tirrwq(i,j,ilam,k)
                   ENDDO
               ENDDO                      
           ENDDO  !l

        DO j=jMin,jMax
          DO i=iMin,iMax
             PURl(i,j,k) = 0.
            do ilam = 1,tlam
               PURl(i,j,k) = PURl(i,j,k) + (tirrwq(i,j,ilam,k)
     &                                    * actot(i,j,k,ilam))  
            enddo !l
         ENDDO    !i
       ENDDO      !j          

CCEA Compute averages and wb=index for exporting
      
        DO j=jMin,jMax
         DO i=iMin,iMax
             a_kave(i,j,k) = 0.d0
             acdom_kave(i,j,k) = 0.d0
             apart_kave(i,j,k) = 0.d0
             actot_ave(i,j,k) = 0.d0 

             bt_kave(i,j,k) = 0.d0
             bpart_kave(i,j,k) = 0.d0
             bctot_ave(i,j,k) = 0.d0

             bb_kave(i,j,k) = 0.d0
             bbpart_kave(i,j,k) = 0.d0
             bbctot_ave(i,j,k) = 0.d0
             
              do ilam = 1,tlam  
                 a_kave(i,j,k) = a_kave(i,j,k)
     &                           + wb_width(ilam) * a_k(i,j,k,ilam)
                 acdom_kave(i,j,k) = acdom_kave(i,j,k)
     &                       + wb_width(ilam) * acdom_k(i,j,k,ilam)
                 apart_kave(i,j,k) = apart_kave(i,j,k)
     &                       + wb_width(ilam) * apart_k(i,j,k,ilam)
                 actot_ave(i,j,k) = actot_ave(i,j,k)
     &                         + wb_width(ilam) * actot(i,j,k,ilam)

                 bt_kave(i,j,k) = bt_kave(i,j,k)
     &                           + wb_width(ilam) * bt_k(i,j,k,ilam)
                 bpart_kave(i,j,k) = bpart_kave(i,j,k)
     &                        + wb_width(ilam) * bpart_k(i,j,k,ilam)
                 bctot_ave(i,j,k) = bctot_ave(i,j,k)
     &                          + wb_width(ilam) * bctot(i,j,k,ilam)

                 bb_kave(i,j,k) = bb_kave(i,j,k)
     &                          + wb_width(ilam) * bb_k(i,j,k,ilam)
                 bbpart_kave(i,j,k) = bbpart_kave(i,j,k)
     &                      + wb_width(ilam) * bbpart_k(i,j,k,ilam)
                 bbctot_ave(i,j,k) = bbctot_ave(i,j,k)
     &                        + wb_width(ilam) * bbctot(i,j,k,ilam)
              
              enddo   !l
           
             a_kave(i,j,k) = a_kave(i,j,k) / wb_totalWidth
             acdom_kave(i,j,k) = acdom_kave(i,j,k) / wb_totalWidth
             apart_kave(i,j,k) = apart_kave(i,j,k) / wb_totalWidth
             actot_ave(i,j,k) = actot_ave(i,j,k) / wb_totalWidth

             bt_kave(i,j,k) = bt_kave(i,j,k) / wb_totalWidth
             bpart_kave(i,j,k) = bpart_kave(i,j,k) / wb_totalWidth
             bctot_ave(i,j,k) = bctot_ave(i,j,k) / wb_totalWidth

             bb_kave(i,j,k) = bb_kave(i,j,k) / wb_totalWidth
             bbpart_kave(i,j,k) = bbpart_kave(i,j,k) / wb_totalWidth
             bbctot_ave(i,j,k) = bbctot_ave(i,j,k) / wb_totalWidth
             
#ifdef RECOM_CALC_REFLEC
             index = darwin_diag_acdom_ilam
             PARw_kwb(i,j,k) = PARw_k(i,j,index,k)  
             a_kwb(i,j,k) = a_k(i,j,k,index) 
             acdom_kwb(i,j,k) = acdom_k(i,j,k,index)
             apart_kwb(i,j,k) = apart_k(i,j,k,index) 
             actot_wb(i,j,k) = actot(i,j,k,index)

             bt_kwb(i,j,k) = bt_k(i,j,k,index) 
             bpart_kwb(i,j,k) = bpart_k(i,j,k,index) 
             bctot_wb(i,j,k) = bctot(i,j,k,index)

             bb_kwb(i,j,k) = bb_k(i,j,k,index) 
             bbpart_kwb(i,j,k) = bbpart_k(i,j,k,index) 
             bbctot_wb(i,j,k) = bbctot(i,j,k,index)

             Edz_wb(i,j,k) = Edz(i,j,index,k)
             Esz_wb(i,j,k) = Esz(i,j,index,k)
             Euz_wb(i,j,k) = Euz(i,j,index,k)
             Estop_wb(i,j,k) = Estop(i,j,index,k)
             Eutop_wb(i,j,k) = Eutop(i,j,index,k)
             amp1_wb(i,j,k) = amp1(i,j,index,k)
             amp2_wb(i,j,k) = amp2(i,j,index,k)
#endif            
          ENDDO    !i
        ENDDO      !j 
      ENDDO        !k

c     PARw and PARwup from non-spectral RECOM are from previous timestep
c     (attenuation done in recom_sms) but PARw and PARwup from WAVEBANDS
c     and RADTRANS are for the current timestep.


c ----------------- Reflectance -------------------------         

        DO ilam = 1,tlam
          DO j=jMin,jMax
            DO i=iMin,iMax
              IF( dz_k(i,j,kSurface) .GT. 0.0 )THEN
                IF(Eswsf(i,j,ilam).GE.darwin_radmodThresh .OR.
     &             Edwsf(i,j,ilam).GE.darwin_radmodThresh ) THEN
                Eupwel(i,j,ilam) = Eutop(i,j,ilam,kSurface)
                Reflec(i,j,ilam) = Eupwel(i,j,ilam) /
     &                  (Eswsf(i,j,ilam) + Edwsf(i,j,ilam))
C     iops surface
             a_ksur(i,j,ilam) = a_k(i,j,kSurface,ilam) 
             acdom_ksur(i,j,ilam) = acdom_k(i,j,kSurface,ilam)
             apart_ksur(i,j,ilam) = apart_k(i,j,kSurface,ilam) 
             actot_sur(i,j,ilam) = actot(i,j,kSurface,ilam)
             bt_ksur(i,j,ilam) = bt_k(i,j,kSurface,ilam) 
             bpart_ksur(i,j,ilam) = bpart_k(i,j,kSurface,ilam) 
             bctot_sur(i,j,ilam) = bctot(i,j,kSurface,ilam)
             bb_ksur(i,j,ilam) = bb_k(i,j,kSurface,ilam)
             bbpart_ksur(i,j,ilam) = bbpart_k(i,j,kSurface,ilam)
             bbctot_sur(i,j,ilam) = bbctot(i,j,kSurface,ilam)            
                ENDIF !light
              ENDIF !depth         
            ENDDO !i
           ENDDO  !j
          ENDDO   !l   
#endif /* RECOM_RADTRANS */ 
#endif /* RECOM_WAVEBANDS */

       
#ifdef ALLOW_EXF
      IF ( recom_windFile .EQ. ' '  ) THEN
C     if no wind has been read copy wind speed of exf package (which
C     may or may not have been read from a file)
       DO j=jMin,jMax
        DO i=iMin,iMax
         windSpeed(i,j,bi,bj) = wspeed(i,j,bi,bj)
        ENDDO
       ENDDO  
      ENDIF 
#endif /* ALLOW_EXF */
C     Compute DIC Flux
      DO j=1-OLy, sNy+OLy
       DO i=1-OLx, sNx+OLx
C     initialise local arrays
C     recom_co2flux should not have to deal with dic/alk=0, so make
C     sure that recom_tiny is larger than zero.
        recom_dic(I,J) = MAX(recom_tiny,
     &       pTracer(I,J,kSurface,bi,bj,iDIC))
        recom_alk(I,J) = MAX(recom_tiny,
     &       pTracer(I,J,kSurface,bi,bj,iALK))
        dicFluxLoc(I,J)   = 0. _d 0
        co2_ppm(I,J)      = 0. _d 0
        co2_mol(I,J)      = 0. _d 0
        hPlusLoc(I,J)     = hPlus(I,J,bi,bj)
        windSpeedLoc(I,J) = windSpeed(I,J,bi,bj)
        atmospCO2Loc(I,J) = atmospCO2(I,J,bi,bj)
C        print*,'atmCO2', atmospCO2Loc
#ifdef RECOM_O2
        recom_oxy(I,J) = MAX(recom_tiny,
     &       pTracer(I,J,kSurface,bi,bj,ioxy))
        FluxO2Loc(I,J)   = 0. _d 0
#endif
#ifdef ALLOW_LONGSTEP
        thetaLoc(I,J)     = LS_theta(I,J,kSurface,bi,bj)
C     this is a terrible hack and we should never have negative
C     salinities
        saltLoc(I,J)      = MAX(0.,LS_salt(I,J,kSurface,bi,bj) )
#else
        thetaLoc(I,J)     = theta(I,J,kSurface,bi,bj)
C     this is a terrible hack and we should never have negative
C     salinities
        saltLoc(I,J)      = MAX(0.,salt(I,J,kSurface,bi,bj) )
#endif /*  ALLOW_LONGSTEP */
CML        saltLoc(I,J)      = salt (I,J,kSurface,bi,bj)
CML        print '(A,2I4,3E10.2)','ml-dic',i,j,hPlus(I,J,bi,bj),
CML     &       pTracer(I,J,kSurface,bi,bj,iDIC),
CML     &       pTracer(I,J,kSurface,bi,bj,iALK)
       ENDDO
      ENDDO
      CALL RECOM_CO2FLUX(
     I     atmospCO2loc, recom_dic, recom_alk,
     I     thetaLoc, saltLoc, windSpeedLoc,
     O     dicFluxLoc, co2_ppm, co2_mol,
     U     hPlusLoc,
     I     iMin, iMax, jMin, jMax )

#ifdef RECOM_O2
      CALL RECOM_O2_SURF(
     &     recom_oxy, windSpeedLoc,FluxO2Loc,
     &     bi,bj,iMin,iMax,jMin,jMax,
     I     myIter, myTime, myThid )
#endif
CML      DO j=jMin,jMax
CML       DO i=iMin,iMax
CML        IF ( hFacC(i,j,ksurface,bi,bj) .gt. 0. _d 0 ) THEN
CML         CALL RECOM_CO2FLUX(
CML     I        atmospCO2loc(I,J), recom_dic(I,J), recom_alk(I,J),
CML     I        thetaLoc(I,J), saltLoc(I,J), windSpeedLoc(I,J),
CML     O        dicFluxLoc(I,J), co2_ppm(I,J), co2_mol(I,J),
CML     U        hPlusLoc(I,J),
CML     I        iMin, iMax, jMin, jMax )
CMLCMLC co2_ppm_dummy is micro atmosheres indeed (compare recom_co2_flux)
CMLCML         print '(A,8E10.2)','ml-dic',
CMLCML     &        dicFluxLoc,atmospCO2(I,J,bi,bj),
CMLCML     &        recom_dic,recom_alk, 
CMLCML     &        theta(I,J,kSurface,bi,bj),salt(I,J,kSurface,bi,bj),
CMLCML     &        windSpeed(I,J,bi,bj), hPlus(I,J,bi,bj)    
CML        ENDIF
CML       ENDDO
CML      ENDDO
C     now copy back some variables
      DO j=jMin,jMax
       DO i=iMin,iMax
        IF ( hFacC(i,j,ksurface,bi,bj) .gt. 0. _d 0 ) THEN
         surfpco2(I,J,bi,bj) = co2_ppm(I,J)   *maskC(I,J,kSurface,bi,bj)
         dicFlux (I,J,bi,bj) = dicFluxLoc(I,J)*maskC(I,J,kSurface,bi,bj)
         hPlus   (I,J,bi,bj) = hPlusLoc(I,J)  *maskC(I,J,kSurface,bi,bj)
CCV#ifdef ALLOW_RECOM_FEINPUT
         feDustFluxLoc(i,j) = feDustDeposit(i,j,bi,bj)
CCV#endif /* ALLOW_RECOM_FEINPUT */
CMS RECOM_IRON_LITHOGENIC_PARTICLES
#ifdef RECOM_O2
         FluxO2(I,J,bi,bj) = FluxO2Loc(I,J)*maskC(I,J,kSurface,bi,bj)
#endif /* RECOM_O2 */
#ifdef RECOM_IRON_LITHOGENIC_PARTICLES
CMS calculate ironLithogenicParticleFluxLoc
CMS feDustDeposit :: iron dust deposition at the surface, in mumol/m^2/s
CMS vgl.: RECOM.h l.150
CMS so get it in mug/m^2/s by multipliing with M = 56 g/mol
CMS fe is 3.5% of dust so reverse it by multipliing with 100/3.5
CMS divide by 1000 to get it in mg/m^2/s
	 ironLithParticleFluxLoc(i,j) = feDustDeposit(i,j,bi,bj)*56*0.1/3.5
#endif  /* RECOM_IRON_LITHOGENIC_PARTICLES */
        ENDIF
       ENDDO
      ENDDO
C     Apply sea-ice mask according to the two different sea ice models
C     if necessary
#ifdef ALLOW_SEAICE
      IF ( useSeaice .AND. .NOT.useThSice ) THEN
       DO j=jMin,jMax
        DO i=iMin,iMax
         dicFlux(i,j,bi,bj) = dicFlux(i,j,bi,bj)
     &        *(1.0 _d 0 - AREA(i,j,bi,bj))
CCV#ifdef ALLOW_RECOM_FEINPUT
         feDustFluxLoc(i,j) = feDustFluxLoc(i,j)
     &        *(1.0 _d 0 - AREA(i,j,bi,bj))
CCV#endif /* ALLOW_RECOM_FEINPUT */
CMS RECOM_IRON_LITHOGENIC_PARTICLES
#ifdef RECOM_O2
         FluxO2(i,j,bi,bj) = FluxO2(i,j,bi,bj)
     &        *(1.0 _d 0 - AREA(i,j,bi,bj))
#endif /* RECOM_O2 */
#ifdef RECOM_IRON_LITHOGENIC_PARTICLES
	 ironLithParticleFluxLoc(i,j) = ironLithParticleFluxLoc(i,j)
     &        *(1.0 _d 0 - AREA(i,j,bi,bj))
#endif /* RECOM_IRON_LITHOGENIC_PARTICLES */
        ENDDO
       ENDDO
      ENDIF
#endif /* ALLOW_SEAICE */
#ifdef ALLOW_THSICE
      IF ( useThSice ) THEN
       DO j=jMin,jMax
        DO i=iMin,iMax
         dicFlux(i,j,bi,bj) = dicFlux(i,j,bi,bj)
     &        *(1.0 _d 0 - iceMask(i,j,bi,bj))
CCV#ifdef ALLOW_RECOM_FEINPUT
         feDustFluxLoc(i,j) = feDustFluxLoc(i,j)
     &        *(1.0 _d 0 - iceMask(i,j,bi,bj))
CCV#endif /* ALLOW_RECOM_FEINPUT */
CMS RECOM_IRON_LITHOGENIC_PARTICLES
#ifdef RECOM_O2
         FluxO2(i,j,bi,bj) = FluxO2(i,j,bi,bj)
     &        *(1.0 _d 0 - iceMask(i,j,bi,bj))
#endif /* RECOM_O2 */
#ifdef RECOM_IRON_LITHOGENIC_PARTICLES
	 ironLithParticleFluxLoc(i,j) = ironLithParticleFluxLoc(i,j)
     &        *(1.0 _d 0 - iceMask(i,j,bi,bj))
#endif /* RECOM_IRON_LITHOGENIC_PARTICLES */
        ENDDO
       ENDDO
      ENDIF	
#endif /* ALLOW_THSICE */
      DO j=jMin,jMax
       DO i=iMin,iMax
        IF ( hFacC(i,j,ksurface,bi,bj) .gt. 0. _d 0 ) THEN
C     determine loop ranges
         krLoc = klowC(i,j,bi,bj)
C     only if wet point:
C     surface forcing for individual tracers:
C     DIC
C     use surfaceForcingPtr to store forcing for the next time step to be
C     used with kpp (and pray that it is correct); IMPORTANT: convert to
C     units per second (dicFlux is in DIC/m^2/day)
         surfaceForcingPtr(i,j,bi,bj,iDIC) = 
     &        - dicFlux(i,j,bi,bj)/recom_secondsPerDay
CCV#ifdef ALLOW_RECOM_FEINPUT
C     Iron:
C     feDustFlux is iron, needs to be converted into dissolved iron, also
C     obey MITgcm sign conventions (positive = upward flux)
         surfaceForcingPtr(i,j,bi,bj,iFe) = 
     &        - maskC(i,j,kSurface,bi,bj)
     &       * constantIronSolubility*feDustFluxLoc(i,j)
         pTracer(i,j,ksurface,bi,bj,iFe) =
     &        pTracer(i,j,ksurface,bi,bj,iFe)
     &        - surfaceForcingPtr(i,j,bi,bj,iFe)
     &        * recip_drF(kSurface)*recip_hFacC(i,j,kSurface,bi,bj)
     &        * dTtracerLev(kSurface)
         
C     resuspension of sediment leads to iron flux into the ocean at 
C     the bottom, which is assumed to be proportional to an *imaginary* 
C     sediment erosion flux according to Ariathurai and Arulanandan (1978).
C     Stress units are m^2/s^2 which is a bit unusual for MITgcm, but
C     we stick with them for now. Also bottomStress has been scaled by
C     the bottomDragLinear/Quadratic already
#ifdef ALLOW_LONGSTEP 
         bottomStress = 0.25 _d 0 *
     &        ((LS_uVel(I,J,krLoc,bi,bj)+LS_uVel(I+1,J,krLoc,bi,bj))**2 
     &        +(LS_vVel(I,J,krLoc,bi,bj)+LS_vVel(I,J+1,krLoc,bi,bj))**2)
#else
         bottomStress = 0.25 _d 0 *
     &        ( (uVel(I,J,krLoc,bi,bj)+uVel(I+1,J,krLoc,bi,bj))**2 
     &        + (vVel(I,J,krLoc,bi,bj)+vVel(I,J+1,krLoc,bi,bj))**2 )
#endif /* ALLOW_LONGSTEP */
C     limit the stress and thus the erosion flux (may be useful in
C     very shallow regions with high velocities)
         bottomStress = MIN( bottomStress, recomShearStressMax )
         erosionFlux(i,j) = erosionCoeff
     &        * MAX( 0. _d 0, bottomStress-recomCritShearStress )
         pTracer(i,j,krLoc,bi,bj,iFe) =
     &        pTracer(i,j,krLoc,bi,bj,iFe)
     &        + erosionFlux(i,j)*dTtracerLev(krLoc)
     &        * recip_drF(krLoc) * _recip_hFacC(i,j,krLoc,bi,bj)
C     MS Implementing RECOM_IRON_LITHOGENIC_PARTICLES
#ifdef RECOM_O2
         surfaceForcingPtr(i,j,bi,bj,ioxy) = 
     &        - FluxO2(i,j,bi,bj)/recom_secondsPerDay
#endif /* RECOM_O2 */
#ifdef RECOM_IRON_LITHOGENIC_PARTICLES
         surfaceForcingPtr(i,j,bi,bj,iDust) = 
     &        - maskC(i,j,kSurface,bi,bj)
     &       * ironLithParticleFluxLoc(i,j)
c         surfaceForcingPtr(i,j,bi,bj,iDust) = 
c     &        - maskC(i,j,kSurface,bi,bj) * 1.0
         pTracer(i,j,ksurface,bi,bj,iDust) =
     &        pTracer(i,j,ksurface,bi,bj,iDust)
     &        - surfaceForcingPtr(i,j,bi,bj,iDust)
     &        * recip_drF(kSurface)*recip_hFacC(i,j,kSurface,bi,bj)
     &        * dTtracerLev(kSurface)
#endif /* RECOM_IRON_LITHOGENIC_PARTICLES */
CCV#endif /* ALLOW_RECOM_FEINPUT */
C     copy global fields to local arrays
         BenthicLayerloc(i,j,1) = benthicLayerN(i,j,bi,bj)
         benthicLayerloc(i,j,2) = benthicLayerC(i,j,bi,bj)
CCV#ifdef ALLOW_RECOM_SILICATE
         benthicLayerloc(i,j,3) = benthicLayerSi(i,j,bi,bj)
CCV#else /* not ALLOW_RECOM_SILICATE */
CCV         benthicLayer(i,j,3) = 0.
CCV#endif /* ALLOW_RECOM_SILICATE */
         benthicLayerloc(i,j,4) = benthicLayerCaCO3(i,j,bi,bj)

C     end if wet point
        ENDIF
       ENDDO
      ENDDO
      
#ifdef RECOM_IRON_HYDROTHERMAL
C add hydrothermal flux to iron field
      DO k=1,Nr
       DO j=jMin,jMax
        DO i=iMin,iMax
         pTracer(i,j,k,bi,bj,iFe)=
     &         pTracer(i,j,k,bi,bj,iFe)
     &         + feHydrothermal(i,j,k,bi,bj) 
     &         * maskC(i,j,k,bi,bj)
     &         * dTtracerLev(kSurface)
#ifdef RECOM_IRON_VARIABLE_LIGANDS
         pTracer(i,j,k,bi,bj,iligw)=
     &         pTracer(i,j,k,bi,bj,iligw)
     &         + hydro_lig2fe * feHydrothermal(i,j,k,bi,bj) 
     &         * maskC(i,j,k,bi,bj)
     &         * dTtracerLev(kSurface)
#endif /* RECOM_IRON_VARIABLE_LIGANDS */
        ENDDO
       ENDDO
      ENDDO
#endif /* RECOM_IRON_HYDROTHERMAL */

C     call the bio-geo-chemical model
      CALL RECOM_SMS(iMin,iMax,jMin,jMax,bi,bj,
     I     Nr, klowC, PTRACERS_num, nDiagsREcoM3D, nDiagsREcoM2D,
     I     dTtracerLev(kSurface), 
     I     maskC, hFacC, recip_hFacC, drF, recip_drF, recip_drC,
#ifdef ALLOW_LONGSTEP
     I     LS_theta, 
#else
     I     theta, 
#endif /* ALLOW_LONGSTEP */
     I     dicFlux, PARadiation, cosAngleOfIncidence, 
     I     pTracer,
     U     benthicLayerloc, 
     O     gPtrLoc, recom_diag_3Dfields, recom_diag_2Dfields
#ifdef RECOM_WAVEBANDS
     I     ,PARl      
     I     ,PARw_k      
     O     ,C_phot_nl 
     O     ,C_phot_nl_dia                   
     O     ,Ek_nl
     O     ,Ek_nl_dia
#endif
     &     ) 
C     copy new local tendency term to global array and perform
C     explicit (Eulerian) time step
      DO iTracer=1,PTRACERS_numInUse
       DO k=1,Nr
        DO j=jMin,jMax
         DO i=iMin,iMax
          pTracer(i,j,k,bi,bj,iTracer)=max(recom_tiny,
     &         pTracer(i,j,k,bi,bj,iTracer)
     &         +gPtrLoc(i,j,k,bi,bj,iTracer))*maskC(i,j,k,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
C     store benthic layer concentrations to global array
      DO j=jMin,jMax
       DO i=iMin,iMax
        benthicLayerN(i,j,bi,bj) = benthicLayerloc(i,j,1)
        benthicLayerC(i,j,bi,bj) = benthicLayerloc(i,j,2)
CCV#ifdef ALLOW_RECOM_SILICATE
        benthicLayerSi(i,j,bi,bj) = benthicLayerloc(i,j,3)
CCV#endif /* ALLOW_RECOM_SILICATE */
        benthicLayerCaCO3(i,j,bi,bj) = benthicLayerloc(i,j,4)
       ENDDO
      ENDDO
#ifdef ALLOW_DIAGNOSTICS
      IF (useDiagnostics) THEN
C     vertical integral of 3D diagnostics
       DO idiag = 1,nDiagsRecom3D
        DO k=1,Nr
         DO j=jMin,jMax
          DO i=iMin,iMax
           recom_diag_integral(i,j,idiag) = 
     &          recom_diag_integral(i,j,idiag) 
     &          + recom_diag_3Dfields(i,j,k,idiag)
     &          * drF(k) * _hFacC(I,J,k,bi,bj)
          Enddo
         ENDDO
        ENDDO
       ENDDO
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

      
#ifdef ALLOW_DIAGNOSTICS
      IF (useDiagnostics) THEN
       IF ( nDiagsREcoM3D .GE. 9 ) THEN
        CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,1),
     &                       'net_pps ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_integral(1-OLx,1-OLy,1),
     &                       'NETPPVIS', 0,  1, 2, bi, bj, myThid)

        CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,2),
     &                       'net_ppd ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_integral(1-OLx,1-OLy,2),
     &                       'NETPPVID', 0,  1, 2, bi, bj, myThid)

        CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,3),
     &                       'gr_pps  ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_integral(1-OLx,1-OLy,3),
     &                       'GRPPVIS ', 0,  1, 2, bi, bj, myThid)

        CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,4),
     &                       'gr_ppd  ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_integral(1-OLx,1-OLy,4),
     &                       'GRPPVID ', 0,  1, 2, bi, bj, myThid)

        CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,5),
     &                       'net_nass', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_integral(1-OLx,1-OLy,5),
     &                       'NETNAVIS', 0,  1, 2, bi, bj, myThid)

        CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,6),
     &                       'net_nasd', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_integral(1-OLx,1-OLy,6),
     &                       'NETNAVID', 0,  1, 2, bi, bj, myThid)

        CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,7),
     &                       'n_assims', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_integral(1-OLx,1-OLy,7),
     &                       'GRNAVIS ', 0,  1, 2, bi, bj, myThid)

        CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,8),
     &                       'n_assimd', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_integral(1-OLx,1-OLy,8),
     &                       'GRNAVID ', 0,  1, 2, bi, bj, myThid)

        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,9),
     &            'free_Fe ', 0, Nr, 2, bi, bj, myThid)
       ENDIF

#ifdef RECOM_MAREMIP
       ndiags3d_used = 9
       if ( nDiagsREcom3D . GE. (ndiags3d_used+28) ) then
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+1),
     &            'reminer ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+2),
     &            'graztot ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+3),
     &            'resptot ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+4),
     &            'calcprod', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+5),
     &            'calcdiss', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+6),
     &            'rdoczoo ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+7),
     &            'rpoczoo ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+8),
     &            'cdiapoc ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+9),
     &            'cphypoc ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+10),
     &            'cdiadoc ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+11),
     &            'cphydoc ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+12),
     &            'cresdia ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+13),
     &            'cresphy ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+14),
     &            'creszoo ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+15),
     &            'grazdia ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+16),
     &            'grazphy ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+17),
     &            'par3d   ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+18),
     &            'nlimdia ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+19),
     &            'nlimphy ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+20),
     &            'llimdia ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+21),
     &            'llimphy ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+22),
     &            'felimdia', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+23),
     &            'felimphy', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+24),
     &            'silimdia', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+25),
     &            'sink_pon', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+26),
     &            'sink_poc', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+27),
     &            'sink_bsi', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+28),
     &            'sinkcalc', 0, Nr, 2, bi, bj, myThid)
       endif
#endif
#ifdef RECOM_IRON_VARIABLE_LIGANDS
       ndiags3d_used = 9
#ifdef RECOM_MAREMIP
     &    + 28
#endif
#ifdef RECOM_TWO_LIGANDS
       if ( nDiagsREcom3D . GE. (ndiags3d_used+11) ) then
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+1),
     &            'ligw_rem', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+2),
     &            'ligw_pho', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+3),
     &            'ligw_deg', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+4),
     &            'ligw_phd', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+5),
     &            'ligw_upt', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+6),
     &            'ligw_coa', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+7),
     &            'ligs_doc', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+8),
     &            'ligs_deg', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+9),
     &            'ligs_phd', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+10),
     &            'ligs_upt', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+11),
     &            'ligs_coa', 0, Nr, 2, bi, bj, myThid)
       endif       
#else
       if ( nDiagsREcom3D . GE. (ndiags3d_used+6) ) then
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+1),
     &            'ligsrrem', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+2),
     &            'ligsrdoc', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+3),
     &            'ligremin', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+4),
     &            'ligphot ', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+5),
     &            'liguptak', 0, Nr, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(
     &            recom_diag_3Dfields(1-OLx,1-OLy,1,ndiags3d_used+6),
     &            'lig_coag', 0, Nr, 2, bi, bj, myThid)
       endif
#endif /* RECOM_TWO_LIGANDS */
#endif /* RECOM_VARIABLE_LIGANDS */ 

       
       IF ( nDiagsREcoM2D .GE. 4) THEN
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,1),
     &                       'SEDFN   ', 0,  1, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,2),
     &                       'SEDFC   ', 0,  1, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,3),
     &                       'SEDFSI  ', 0,  1, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,4),
     &                       'SEDFCALC', 0,  1, 2, bi, bj, myThid)
       ENDIF

       IF ( nDiagsREcoM2D .GE. 9) THEN
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,5),
     &                       'DIAFN   ', 0,  1, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,6),
     &                       'DIAFC   ', 0,  1, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,7),
     &                       'DIAFALK ', 0,  1, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,8),
     &                       'DIAFSI  ', 0,  1, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,9),
     &                       'DIAFFE  ', 0,  1, 2, bi, bj, myThid)
       ENDIF

#ifdef RECOM_EXPORT_DIAGNOSTICS
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,10),
     &                       'EXPORTN ', 0,  1, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,11),
     &                       'EXPORTC ', 0,  1, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,12),
     &                       'EXPCALC ', 0,  1, 2, bi, bj, myThid)
        CALL DIAGNOSTICS_FILL(recom_diag_2Dfields(1-OLx,1-OLy,13),
     &                       'EXPORTSI', 0,  1, 2, bi, bj, myThid)
#endif

       CALL DIAGNOSTICS_FILL(PARadiation,
     &                      'PARSURF ', 0,  1, 2, bi, bj, myThid)      
       CALL DIAGNOSTICS_FILL(surfaceForcingPtr(1-Olx,1-Oly,1,1,iDIC),
     &                      'CO2Flx  ', 0,  1, 1, bi, bj, myThid) 
       CALL DIAGNOSTICS_FILL(surfpco2,
     &                      'pCO2surf ', 0,  1, 1, bi, bj, myThid)
       CALL DIAGNOSTICS_FILL(hPlus,
     &                      'HPlus   ', 0,  1, 1, bi, bj, myThid)
       CALL DIAGNOSTICS_FILL(surfaceForcingPtr(1-Olx,1-Oly,1,1,iFe),
     &                      'FeSrfFlx', 0,  1, 1, bi, bj, myThid)
       CALL DIAGNOSTICS_FILL(erosionFlux,
     &                      'FeBtmFlx', 0,  1, 2, bi, bj, myThid)
       CALL DIAGNOSTICS_FILL(benthicLayerN,
     &                      'NBENTHOS', 0,  1, 1, bi, bj, myThid)
       CALL DIAGNOSTICS_FILL(benthicLayerC,
     &                      'CBENTHOS', 0,  1, 1, bi, bj, myThid)
CCV#ifdef ALLOW_RECOM_SILICATE
       CALL DIAGNOSTICS_FILL(benthicLayerSi,
     &                      'SBENTHOS', 0,  1, 1, bi, bj, myThid)
CCV#endif /* ALLOW_RECOM_SILICATE */
       CALL DIAGNOSTICS_FILL(benthicLayerCaCO3,
     &                      'CALCBENT', 0,  1, 1, bi, bj, myThid)
      ENDIF


    
#ifdef RECOM_PHOTODAMAGE
C Diagnostics photodamage      
       ndiags3d_used = 9
#ifdef RECOM_MAREMIP
     &    + 28
#endif
#ifdef RECOM_IRON_VARIABLE_LIGANDS
#ifdef RECOM_IRON_TWO_LIGANDS
     &    + 11
#else        
     &    + 6
#endif
#endif
       if ( nDiagsREcom3D . GE. (ndiags3d_used+8) ) then       
      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                        ndiags3d_used+1),
     &                        'pphotphy', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                        ndiags3d_used+2),
     &                        'pphotdia', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                        ndiags3d_used+3),
     &                        'ppmaxphy', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                         ndiags3d_used+4),
     &                        'ppmaxdia', 0, Nr, 2, bi, bj, myThid)      
      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                        ndiags3d_used+5),
     &                        'chlsyphy', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                         ndiags3d_used+6),
     &                        'chlsydia', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                        ndiags3d_used+7),
     &                        'kochlphy', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                         ndiags3d_used+8),
     &                        'kochldia', 0, Nr, 2, bi, bj, myThid)
      endif
#endif /* RECOM_PHOTODAMAGE */




#ifdef RECOM_MARSHALL
C Diagnostics Non-photosynthetic pigments      
       ndiags3d_used = 9
#ifdef RECOM_MAREMIP
     &    + 28
#endif
#ifdef RECOM_IRON_VARIABLE_LIGANDS
#ifdef RECOM_IRON_TWO_LIGANDS
     &    + 11
#else        
     &    + 6
#endif     
#endif        
#ifdef RECOM_PHOTODAMAGE
     &    + 8
#endif
      IF ( nDiagsREcom3D . GE. (ndiags3d_used+4) ) then         
      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                        ndiags3d_used+1),
     &                       'alphaphy', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                         ndiags3d_used+2),
     &                       'alphadia', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                        ndiags3d_used+3),
     &                       'ppcphy  ', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                         ndiags3d_used+4),
     &                       'ppcdia  ', 0, Nr, 2, bi, bj, myThid)
      ENDIF
#endif     

      
C     Diagnostics spectral light      
#ifdef RECOM_WAVEBANDS
       ndiags3d_used = 9
#ifdef RECOM_MAREMIP
     &    + 28
#endif
#ifdef RECOM_IRON_VARIABLE_LIGANDS
#ifdef RECOM_IRON_TWO_LIGANDS
     &    + 11
#else        
     &    + 6
#endif     
#endif        
#ifdef RECOM_PHOTODAMAGE
     &    + 8
#endif
#ifdef RECOM_MARSHALL
     &    + 4
#endif       
CEA 3D (x,y,z) from sms      
       IF ( nDiagsREcom3D . GE. (ndiags3d_used+3) ) then        
      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                        ndiags3d_used+1),
     &                       'par3dw  ', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                        ndiags3d_used+4),
     &                       'Ek_phy  ', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL(recom_diag_3Dfields(1-OLx,1-OLy,1,
     &                        ndiags3d_used+5),
     &                       'Ek_dia  ', 0, Nr, 2, bi, bj, myThid)
      ENDIF

CEA 3D (x,y,z) from forcing
      CALL DIAGNOSTICS_FILL( a_kave ,
     &                      'abtotave', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( acdom_kave,
     &                      'acdomave', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( apart_kave,
     &                      'abparave', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( actot_ave,
     &                      'abphyave', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( PURl,
     &                      'pur3dw  ', 0, Nr, 2, bi, bj, myThid)
#ifdef RECOM_CALC_REFLEC
      CALL DIAGNOSTICS_FILL( PARw_kwb,
     &                      'par3d_wb', 0, Nr, 2, bi, bj, myThid)      
      CALL DIAGNOSTICS_FILL( a_kwb,
     &                      'abtot_wb', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( acdom_kwb,
     &                      'acdom_wb', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( apart_kwb,
     &                      'abpar_wb', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( actot_wb,
     &                      'abphy_wb', 0, Nr, 2, bi, bj, myThid)
#endif
#ifdef RECOM_RADTRANS
      CALL DIAGNOSTICS_FILL( bt_kave,
     &                      'bttotave', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( bpart_kave,
     &                      'btparave', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( bctot_ave,
     &                      'btphyave', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( bb_kave,
     &                      'bbtotave', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( bbpart_kave,
     &                      'bbparave', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( bbctot_ave,
     &                      'bbphyave', 0, Nr, 2, bi, bj, myThid)

#ifdef RECOM_CALC_REFLEC
      CALL DIAGNOSTICS_FILL( bt_kwb,
     &                      'bttot_wb', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( bpart_kwb,
     &                      'btpar_wb', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( bctot_wb,
     &                      'btphy_wb', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( bb_kwb,
     &                      'bbtot_wb', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( bbpart_kwb,
     &                      'bbpar_wb', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( bbctot_wb,
     &                      'bbphy_wb', 0, Nr, 2, bi, bj, myThid)     
      CALL DIAGNOSTICS_FILL( Edz_wb,
     &                      'Edz_wb  ', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( Esz_wb,
     &                      'Esz_wb  ', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( Euz_wb,
     &                      'Euz_wb  ', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( Estop_wb,
     &                      'Estop_wb', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( Eutop_wb,
     &                      'Eutop_wb', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( amp1_wb,
     &                      'amp1_wb ', 0, Nr, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( amp2_wb,
     &                      'amp2_wb ', 0, Nr, 2, bi, bj, myThid)
#endif
#endif /* RADTRANS */


CEA 3D (x,y,l)     
      CALL DIAGNOSTICS_FILL(PARwup_diag,
     &                     'par3doa ', 0, tlam, 2, bi, bj, myThid)   
      CALL DIAGNOSTICS_FILL(Edwsf ,
     &                     'edirect ', 0, tlam, 2, bi, bj, myThid)      
      CALL DIAGNOSTICS_FILL(Eswsf ,
     &                     'ediffuse', 0, tlam, 2, bi, bj, myThid)
C iops surface
      CALL DIAGNOSTICS_FILL( a_ksur ,
     &                      'abtotsur', 0, tlam, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( acdom_ksur,
     &                      'acdomsur', 0, tlam, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( apart_ksur,
     &                      'abparsur', 0, tlam, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( actot_sur,
     &                      'abphysur', 0, tlam, 2, bi, bj, myThid)
#ifdef RECOM_RADTRANS
      CALL DIAGNOSTICS_FILL(Eupwel,
     &                     'euptop  ', 0, tlam, 2, bi, bj, myThid)      
      CALL DIAGNOSTICS_FILL(Reflec,
     &                     'reflecta', 0, tlam, 2, bi, bj, myThid)       
C iops surface
      CALL DIAGNOSTICS_FILL( bt_ksur,
     &                      'bttotsur', 0, tlam, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( bpart_ksur,
     &                      'btparsur', 0, tlam, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( bctot_sur,
     &                      'btphysur', 0, tlam, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( bb_ksur,
     &                      'bbtotsur', 0, tlam, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( bbpart_ksur,
     &                      'bbparsur', 0, tlam, 2, bi, bj, myThid)
      CALL DIAGNOSTICS_FILL( bbctot_sur,
     &                      'bbphysur', 0, tlam, 2, bi, bj, myThid)
CEA 2D (x,y)       
      CALL DIAGNOSTICS_FILL(rmud,
     &                     'rmudave ', 0,  1, 2, bi, bj, myThid)
#endif /* RADTRANS */
      CALL DIAGNOSTICS_FILL(PARwup_total,
     &                     'par2doa ', 0,  1, 2, bi, bj, myThid)      
#endif /* WAVEBANDS */        

       
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_RECOM */

      RETURN
      END

#ifdef ALLOW_RECOM
      SUBROUTINE RECOM_FORCING_DUMMY( 
     I     iMin, iMax, jMin, jMax, bi, bj, myTime, myIter, myThid )
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#include "RECOM.h"
      INTEGER iMin, iMax, jMin, jMax
      integer bi, bj
      integer mythid 
      integer myiter
      _RL     myTime
      INTEGER i,j,k,iTracer,krLoc, kSurface
      ksurface = 1
      DO j=jMin,jMax
       DO i=iMin,iMax
        IF ( hFacC(i,j,ksurface,bi,bj) .gt. 0. _d 0 ) THEN
C--------------------------------------
CCV The meaning of the following part is completely unclear to me. Delete?
CCV#ifdef ALLOW_RECOM_FEINPUT
C     Iron:
         iTracer = 4
C     feDustFlux is iron, needs to be converted into dissolved iron, also
C     obey MITgcm sign conventions (positive = upward flux)
         surfaceForcingPtr(i,j,bi,bj,iTracer) = 
     &        - maskC(i,j,kSurface,bi,bj)
     &       * constantIronSolubility*feDustDeposit(i,j,bi,bj)
         pTracer(i,j,ksurface,bi,bj,iTracer) =
     &        pTracer(i,j,ksurface,bi,bj,iTracer)
     &        - surfaceForcingPtr(i,j,bi,bj,iTracer)
     &        * recip_drF(kSurface)*recip_hFacC(i,j,kSurface,bi,bj)
     &        * dTtracerLev(kSurface)
         pTracer(i,j,ksurface,bi,bj,5) =
     &        pTracer(i,j,ksurface,bi,bj,5)
     &        - surfaceForcingPtr(i,j,bi,bj,iTracer)
     &        * recip_drF(kSurface)*recip_hFacC(i,j,kSurface,bi,bj)
     &        * dTtracerLev(kSurface)
CCV#endif /* ALLOW_RECOM_FEINPUT */
C--------------------------------------
        ENDIF
       ENDDO
      ENDDO
      RETURN
      END
#endif /* ALLOW_RECOM */
#ifdef ALLOW_RECOM
CBOP
C     !ROUTINE: BGC_DUMMY_MODEL
C     !INTERFACE:
      SUBROUTINE BGC_DUMMY_MODEL( 
     INr, kr, nTracer, rC, drF, bgc_temp, bgc_ptracer, 
     Obgc_gPtr, 
     ImyTime, myIter, myThid 
     &)

C     !DESCRIPTION: \bv
C     *==================================================================
C     | SUBROUTINE BGC_DUMMY_MODEL
C     | o Dummy routine for a bio-geo-chemical model
C     *==================================================================
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global variables ==

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     note: under the multi-threaded model myiter and 
C           mytime are local variables passed around as routine 
C           arguments. Although this is fiddly it saves the need to 
C           impose additional synchronisation points when they are 
C           updated.
C     Nr - number of vertical points
C     kr - index of lowest wet tracer grid cell (k>kr means below bottom)
C     nTracer - number biological tracers
C     rC - depth of grid cell centers
C     drF - thickness of vertical grid cells
C     bgc_temp - local vertical array of temperature (potential!)
C     bgc_pTracer - local vertical array of biological tracers
C     bgc_gPtr - local vertical array of biological tracer tendencies
C     myiter - iteration counter for this thread
C     mytime - time counter for this thread
C     mythid - thread number for this instance of the routine.
      INTEGER Nr, kr, nTracer
      _RS     rC(1:Nr),drF(1:Nr)
      _RL     bgc_temp(1:Nr)
      _RL     bgc_pTracer(1:Nr,1:nTracer)
      _RL     bgc_gPtr(1:Nr,1:nTracer)
      integer mythid 
      integer myiter
      _RL     myTime

C     !LOCAL VARIABLES :
C     == Local variables ==
      INTEGER k, iTracer
      CHARACTER*(3) NVAR
      CHARACTER*(128) FMT
CEOP

C--   the bio-geo-chemical model goes here
C     this is just a place-holder
      print *, 'ml-bgc: myIter = ', myIter, ', myTime = ', myTime
      print *, 'ml-bgc:   k, bgc_pTracer'
      write(NVAR,'(I3)') nTracer
      FMT = '(A,I3,2F10.4,'//NVAR//'E10.2)'
      print *, 'ml-bgc:   k     depth  bgc_temp  bgc_pTracer'
      DO k = 1,Nr
       write(*,FMT) ' ml-bgc: ', k, rC(k), bgc_temp(k),
     &      (bgc_pTracer(k,iTracer),iTracer=1,nTracer)
      ENDDO 
      print *, 'ml-bgc:   k     depth  bgc_temp  bgc_gpTr'
      DO k = 1,Nr
       DO iTracer = 1,nTracer
        bgc_gPtr(k,iTracer) = float(k+iTracer)
       ENDDO
       write(*,FMT) ' ml-bgc: ', k, rC(k), bgc_temp(k),
     &     (bgc_gPtr(k,iTracer),iTracer=1,nTracer)
      ENDDO
C--

      RETURN
      END
#endif /* ALLOW_RECOM */

