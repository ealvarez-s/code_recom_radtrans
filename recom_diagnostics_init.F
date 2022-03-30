C $Header: /csys/software/cvsroot/REcoM/MITgcm/recom/recom_diagnostics_init.F,v 1.5 2007/02/23 14:41:11 mlosch Exp $
C $Name:  $

#include "RECOM_OPTIONS.h"
    
CBOP
C     !ROUTINE: RECOM_DIAGNOSTICS_INIT
C     !INTERFACE:
      SUBROUTINE RECOM_DIAGNOSTICS_INIT( myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE RECOM_DIAGNOSTICS_INIT
C     | o Routine to initialize list of all available diagnostics
C     |   for GM/Redi package
C     *==========================================================*
C     \ev
C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "EEPARAMS.h"
c #include "SIZE.h"
c #include "PARAMS.h"
c #include "RECOM.h"
#ifdef RECOM_WAVEBANDS
#include "SPECTRAL_SIZE.h" /* tlam */   
#endif
      
C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

#ifdef ALLOW_DIAGNOSTICS
C     !LOCAL VARIABLES:
C     === Local variables ===
C     diagNum   :: diagnostics number in the (long) list of available diag.
C     diagName  :: local short name (8c) of a diagnostics
C     diagCode  :: local parser field with characteristics of the diagnostics
C              cf head of S/R DIAGNOSTICS_INIT_EARLY or DIAGNOSTICS_MAIN_INIT
C     diagUnits :: local string (16c): physical units of a diagnostic field
C     diagTitle :: local string (80c): description of field in diagnostic
      INTEGER        diagNum
      CHARACTER*8    diagName
      CHARACTER*16   diagCode
      CHARACTER*16   diagUnits
      CHARACTER*(80) diagTitle
      CHARACTER*4    diagSufx

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C----------------------------------------
C primary production
C----------------------------------------
      diagName  = 'net_pps '
      diagTitle = 'small Phy net primary production'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0,myThid )

      diagName  = 'net_ppd '
      diagTitle = 'diatom net primary production'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'gr_pps  '
      diagTitle = 'small Phy gross primary production'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'gr_ppd  '
      diagTitle = 'diatom gross primary production'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      
#ifdef RECOM_MAREMIP
C----------------------------------------
C Other carbon fluxes between ecosystem compartments
C----------------------------------------
      diagName  = 'reminer '
      diagTitle = 'remineralisation of DOC'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      
      diagName  = 'graztot '
      diagTitle = 'total grazing by zooplankton'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      
      diagName  = 'resptot '
      diagTitle = 'total respiration by ecosystem'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      
      diagName  = 'calcprod'
      diagTitle = 'calcite production'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'calcdiss'
      diagTitle = 'calcite dissolution'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'rdoczoo '
      diagTitle = 'DOC production by zooplankton'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'rpoczoo '
      diagTitle = 'POC production by zooplankton'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'cdiapoc '
      diagTitle = 'POC production by diatoms'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'cphypoc '
      diagTitle = 'POC production by small phytoplankton'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'cdiadoc '
      diagTitle = 'DOC production by diatoms'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'cphydoc '
      diagTitle = 'DOC production by small phytoplankton'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'cresdia '
      diagTitle = 'respiration by diatoms'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'cresphy '
      diagTitle = 'respiration by small phytoplankton'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'creszoo '
      diagTitle = 'respiration by zooplankton'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'grazdia '
      diagTitle = 'zooplankton grazing on diatoms'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'grazphy '
      diagTitle = 'zooplankton grazing on small phytoplankton'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C----------------------------------------
C Other 3D diagnostics (PAR, limitations)
C----------------------------------------

      diagName  = 'par3d   '
      diagTitle = 'photosynthetically available radiation'
      diagUnits = 'W/m^2           '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'nlimdia '
      diagTitle = 'diatom nitrogen limitation'
      diagUnits = '                '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      
      diagName  = 'nlimphy '
      diagTitle = 'small phytoplankton nitrogen limitation'
      diagUnits = '                '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      
      diagName  = 'llimdia '
      diagTitle = 'diatom light limitation'
      diagUnits = '                '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      
      diagName  = 'llimphy '
      diagTitle = 'small phytoplankton light limitation'
      diagUnits = '                '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      
      diagName  = 'felimdia'
      diagTitle = 'diatom iron limitation'
      diagUnits = '                '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      
      diagName  = 'felimphy'
      diagTitle = 'small phytoplankton iron limitation'
      diagUnits = '                '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'silimdia'
      diagTitle = 'diatom silicon limitation'
      diagUnits = '                '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )   
      
C----------------------------------------
C sinking flux diagnostics
C----------------------------------------

      diagName  = 'sink_pon'
      diagTitle = 'sinking flux of particulate organic nitrogen'
      diagUnits = 'mmol m^-2 d^-1  '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      
      diagName  = 'sink_poc'
      diagTitle = 'sinking flux of particulate organic carbon'
      diagUnits = 'mmol m^-2 d^-1  '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      
      diagName  = 'sink_bsi'
      diagTitle = 'sinking flux of biogenic silica'
      diagUnits = 'mmol m^-2 d^-1  '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      
      diagName  = 'sinkcalc'
      diagTitle = 'sinking flux of caco3'
      diagUnits = 'mmol m^-2 d^-1  '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      
#endif


#ifdef RECOM_PHOTODAMAGE
CC----------------------------------------
CC Synthesis diagnostics with Photodamage
CC----------------------------------------
      
      diagName  = 'pphotphy'
      diagTitle = 'photosynthesis rate Phyto'
      diagUnits = 'd-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'pphotdia'
      diagTitle = 'photosynthesis rate Diatoms'
      diagUnits = 'd-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'ppmaxphy'
      diagTitle = 'maximum production rate Phyto'
      diagUnits = 'd-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'ppmaxdia'
      diagTitle = 'maximum production rate Diatoms'
      diagUnits = 'd-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'chlsyphy'
      diagTitle = 'chlorophyll synthesis rate Phyto'
      diagUnits = 'gchl/molC/d     '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'chlsydia'
      diagTitle = 'chlorophyll synthesis rate Diatoms'
      diagUnits = 'gchl/molC/d     '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      
      diagName  = 'kochlphy'
      diagTitle = 'chlorophyll degradation rate Phyto'
      diagUnits = 'gchl/molC/d     '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'kochldia'
      diagTitle = 'chlorophyll degradation rate Diatoms'
      diagUnits = 'gchl/molC/d     '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )      

#endif      


#ifdef RECOM_MARSHALL
      diagName  = 'alphaphy'
      diagTitle = 'initial slope PE Phyto Marshall'
      diagUnits = 'm2 molC gChl-1 J-1 d-1'
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'alphadia'
      diagTitle = 'initial slope PE Diatoms Marshall'
      diagUnits = 'm2 molC gChl-1 J-1 d-1'
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'ppcphy  '
      diagTitle = 'ratio PPC:Chla small phyto'
      diagUnits = '                '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'ppcdia  '
      diagTitle = 'ratio PPC:Chla diatoms'
      diagUnits = '                '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif      


      
#ifdef RECOM_WAVEBANDS
CC----------------------------------------
CC WAVEBANDS 3D (x,y,z) diagnostics
CC----------------------------------------

      diagName  = 'par3dw  '
      diagTitle = 'average spectral PAR  '
      diagUnits = 'uE/m^2/s-1      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'Ek_phy  '
      diagTitle = 'phytoplankton spectral light saturation'
      diagUnits = '                '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'Ek_dia  '
      diagTitle = 'diatom spectral light saturation'
      diagUnits = '                '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'pur3dw  '
      diagTitle = 'average photon flux through phytoplankton'
      diagUnits = 'uE/m^3/s-1      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

CCC Absorption averages
      diagName  = 'abtotave'
      diagTitle = 'absorption total average'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'acdomave'
      diagTitle = 'absorption by cdom average'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'abparave'
      diagTitle = 'absorption by detritus average'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'abphyave'
      diagTitle = 'absorption by phytoplankton average'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

#ifdef RECOM_CALC_REFLEC
CCC Light at index waveband
      diagName  = 'par3d_wb'
      diagTitle = 'Light at selected wb'
      diagUnits = 'uE/m^2/s-1      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
CCC   Absorption at index waveband
      diagName  = 'abtot_wb'
      diagTitle = 'absorption total at selected wb'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'acdom_wb'
      diagTitle = 'absorption by cdom at selected wb'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'abpar_wb'
      diagTitle = 'absorption by detritus at selected wb'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'abphy_wb'
      diagTitle = 'absorption by phytoplankton at selected wb'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif /* RECOM_CALC_REFLEC */     

#ifdef RECOM_RADTRANS      
CCC Scatter and Backs. averages
      diagName  = 'bttotave'
      diagTitle = 'scatter total average'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'btparave'
      diagTitle = 'scatter by detritus average'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'btphyave'
      diagTitle = 'scatter by phytoplankton average'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )      
      diagName  = 'bbtotave'
      diagTitle = 'backscatter total average'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'bbparave'
      diagTitle = 'backscatter by detritus average'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'bbphyave'
      diagTitle = 'backscatter by phytoplankton average'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )      

#ifdef RECOM_CALC_REFLEC
CCC Scatter and Backs. at index waveband
      diagName  = 'bttot_wb'
      diagTitle = 'scatter total at selected wb'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'btpar_wb'
      diagTitle = 'scatter by detritus at selected wb'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'btphy_wb'
      diagTitle = 'scatter by phytoplankton at selected wb'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )      
      diagName  = 'bbtot_wb'
      diagTitle = 'backscatter total at selected wb'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'bbpar_wb'
      diagTitle = 'backscatter by detritus at selected wb'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'bbphy_wb'
      diagTitle = 'backscatter by phytoplankton at selected wb'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'Edz_wb  '
      diagTitle = 'Ed at the bottom of the layer at selected wb'
      diagUnits = 'W m-2           '
      diagCode  = 'SM      LR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'Esz_wb  '
      diagTitle = 'Es at the bottom of the layer at selected wb'
      diagUnits = 'W m-2           '
      diagCode  = 'SM      LR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'Euz_wb  '
      diagTitle = 'Eu at the bottom of the layer at selected wb'
      diagUnits = 'W m-2           '
      diagCode  = 'SM      LR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'Estop_wb'
      diagTitle = 'Es at the top of the layer at selected wb'
      diagUnits = 'W m-2           '
      diagCode  = 'SM      UR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'Eutop_wb'
      diagTitle = 'Eu at the top of the layer at selected wb'
      diagUnits = 'W m-2           '
      diagCode  = 'SM      UR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'amp1_wb '
      diagTitle = 'amplitude of downward increasing mode at selected wb'
      diagUnits = '                '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'amp2_wb '
      diagTitle = 'amplitude of downward decreasing mode at selected wb'
      diagUnits = '                '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif /* RECOM_CALC_REFLEC */ 

#endif /* RECOM_RADTRANS */
     
#endif /* RECOM_WAVEBANDS */



#ifdef RECOM_WAVEBANDS
CC----------------------------------------
CC WAVEBANDS 3D (x,y,tlam) diagnostics
CC----------------------------------------
      diagName  = 'par3doa '
      diagTitle = 'spectral PAR below surface      '
      diagUnits = 'uE/m-2/s-1      '
      diagCode  = 'SM      UX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      CALL DIAGNOSTICS_SETKLEV( diagName, tlam, myThid )
      diagName  = 'edirect '
      diagTitle = 'direct irradiance below surface  '
      diagUnits = 'W/m-2           '
      diagCode  = 'SM      UX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )         
      CALL DIAGNOSTICS_SETKLEV( diagName, tlam, myThid )
      diagName  = 'ediffuse'
      diagTitle = 'diffuse irradiance below surface '
      diagUnits = 'W/m-2           '
      diagCode  = 'SM      UX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )  
      CALL DIAGNOSTICS_SETKLEV( diagName, tlam, myThid )
CCC Absorption surface
      diagName  = 'abtotsur'
      diagTitle = 'absorption total in surface'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )  
      CALL DIAGNOSTICS_SETKLEV( diagName, tlam, myThid )  
      diagName  = 'acdomsur'
      diagTitle = 'absorption by cdom in surface'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )  
      CALL DIAGNOSTICS_SETKLEV( diagName, tlam, myThid )  
      diagName  = 'abparsur'
      diagTitle = 'absorption by detritus in surface'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )  
      CALL DIAGNOSTICS_SETKLEV( diagName, tlam, myThid )  
      diagName  = 'abphysur'
      diagTitle = 'absorption by phytoplankton in surface'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )  
      CALL DIAGNOSTICS_SETKLEV( diagName, tlam, myThid )

#ifdef RECOM_RADTRANS    
      diagName  = 'euptop  '
      diagTitle = 'water-leaving irradiance  '
      diagUnits = 'W/m-2           '
      diagCode  = 'SM      UX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      CALL DIAGNOSTICS_SETKLEV( diagName, tlam, myThid )
      diagName  = 'reflecta'
      diagTitle = 'Reflectance     '
      diagUnits = '                '
      diagCode  = 'SM      UX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      CALL DIAGNOSTICS_SETKLEV( diagName, tlam, myThid )

CCC Scatter and Backs. surface
      diagName  = 'bttotsur'
      diagTitle = 'scatter total in surface'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      CALL DIAGNOSTICS_SETKLEV( diagName, tlam, myThid )
      diagName  = 'btparsur'
      diagTitle = 'scatter by detritus in surface'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      CALL DIAGNOSTICS_SETKLEV( diagName, tlam, myThid )
      diagName  = 'btphysur'
      diagTitle = 'scatter by phytoplankton in surface'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      CALL DIAGNOSTICS_SETKLEV( diagName, tlam, myThid )   
      diagName  = 'bbtotsur'
      diagTitle = 'backscatter total in surface'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      CALL DIAGNOSTICS_SETKLEV( diagName, tlam, myThid )
      diagName  = 'bbparsur'
      diagTitle = 'backscatter by detritus in surface'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      CALL DIAGNOSTICS_SETKLEV( diagName, tlam, myThid )
      diagName  = 'bbphysur'
      diagTitle = 'backscatter by phytoplankton in surface'
      diagUnits = 'm-1             '
      diagCode  = 'SM      MX      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      CALL DIAGNOSTICS_SETKLEV( diagName, tlam, myThid )

#endif /* RECOM_RADTRANS */   
#endif /* RECOM_WAVEBANDS */      
      
C----------------------------------------
C free iron
C----------------------------------------
      diagName  = 'free_Fe '
      diagTitle = 'Fe available for scavenging'
      diagUnits = 'micromol/m^3    '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )      
C----------------------------------------
C Nitrogen assimilation
C----------------------------------------
      diagName  = 'net_nass'
      diagTitle = 'small Phy net nitrogen assimilation'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'net_nasd'
      diagTitle = 'diatom net nitrogen assimilation'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'n_assims'
      diagTitle = 'small Phy nitrogen assimilation'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'n_assimd'
      diagTitle = 'diatom nitrogen assimilation'
      diagUnits = 'mmol/m^3/d      '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C-----------------------
C     vertical integrals of the above
C-----------------------
      diagName  = 'NETPPVIS'
      diagTitle = 'Net PP small Phy, vertically integrated'
      diagUnits = 'mmol/m^2/d      '
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'NETPPVID'
      diagTitle = 'Net PP diatoms, vertically integrated'
      diagUnits = 'mmol/m^2/d      '
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'GRPPVIS '
      diagTitle = 'Gross PP small Phy, vertically integrated'
      diagUnits = 'mmol/m^2/d      '
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'GRPPVID '
      diagTitle = 'Gross PP diatoms, vertically integrated'
      diagUnits = 'mmol/m^2/d      '
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'NETNAVIS'
      diagTitle = 'Nitrogen assimilation minus excretion, small Phy '//
     &     'vertically integrated'
      diagUnits = 'mmol/m^2/d      '
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'NETNAVID'
      diagTitle = 'Nitrogen assimilation minus excretion, diatoms '//
     &     'vertically integrated'
      diagUnits = 'mmol/m^2/d      '
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'GRNAVIS '
      diagTitle = 'Nitrogen assim small Phy, vertically integrated'
      diagUnits = 'mmol/m^2/d      '
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'GRNAVID '
      diagTitle = 'Nitrogen assim diatoms, vertically integrated'
      diagUnits = 'mmol/m^2/d      '
      diagCode  = 'SM      M1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C----------------------------------------
C     surface fluxes
C----------------------------------------
      diagName  = 'PARSURF '
      diagTitle = 'Surface Photosynthetically Available Radiation'
      diagUnits = 'W/m^2           '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'CO2Flx  '
      diagTitle = 'Surface flux of CO2'
      diagUnits = 'mmolCO2/m^2/s^1 '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#ifdef RECOM_WAVEBANDS      
#ifdef RECOM_RADTRANS
      diagName  = 'rmudave '
      diagTitle = 'Inverse cosine zenith solar angle'
      diagUnits = '                '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif /* RECOM_RADTRANS */
      diagName  = 'par2doa '
      diagTitle = 'total PAR below surface      '
      diagUnits = 'uE/m-2/s-1      '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif /* RECOM_WAVEBANDS */
C----------------------------------------
CJH Surface pCO2
C----------------------------------------
      diagName  = 'pCO2surf'
      diagTitle = 'Surface pCO2'
      diagUnits = 'muatm           '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C----------------------------------------
C Surface concentration of H+, for calculation of pH
C----------------------------------------
      diagName  = 'HPlus   '
      diagTitle = 'Concentration of H+'
      diagUnits = 'mmol H+/m^3 '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C----------------------------------------
C     iron fluxes at surface and bottom
C----------------------------------------
      diagName  = 'FeSrfFlx'
      diagTitle = 'Surface flux of iron'
      diagUnits = 'mumolFe/m^2/s^1 '
      diagCode  = 'SM      U1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'FeBtmFlx'
      diagTitle = 'Bottom flux of iron'
      diagUnits = 'mumolFe/m^2/s^1 '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
C----------------------------------------
C     benthic layer variables and fluxes
C----------------------------------------
      diagName  = 'NBENTHOS'
      diagTitle = 'organic N inventory in benthos'
      diagUnits = 'mmolN/m^2       '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'SEDFN   '
      diagTitle = 'organic N flux into sediment'
      diagUnits = 'mmol N/m^2/d   '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'CBENTHOS'
      diagTitle = 'organic C inventory in benthos'
      diagUnits = 'mmolC/m^2       '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'SEDFC   '
      diagTitle = 'organic C flux into sediment'
      diagUnits = 'mmol C/m^2/d   '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'SBENTHOS'
      diagTitle = 'biogenic Si inventory in benthos'
      diagUnits = 'mmolSi/m^2      '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'SEDFSI  '
      diagTitle = 'biogenic Si flux into sediment'
      diagUnits = 'mmol Si/m^2/d   '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'CALCBENT'
      diagTitle = 'CaCO3 inventory in benthos'
      diagUnits = 'mmol CaCO3/m^2  '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'SEDFCALC'
      diagTitle = 'CaCO3 flux into sediment'
      diagUnits = 'mmol CaCO3/m^2/d'
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C-- diagenetic fluxes out of sediment
      diagName  = 'DIAFN   '
      diagTitle = 'diagenetic N flux out of sediment'
      diagUnits = 'mmol N/m^2/d   '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'DIAFC   '
      diagTitle = 'diagenetic C flux out of sediment'
      diagUnits = 'mmol C/m^2/d   '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'DIAFALK'
      diagTitle = 'diagenetic Alk flux out of sediment'
      diagUnits = 'mmol Alk/m^2/d'
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'DIAFSI  '
      diagTitle = 'diagenetic Si flux out of sediment'
      diagUnits = 'mmol Si/m^2/d   '
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'DIAFFE'
      diagTitle = 'diagenetic Fe flux out of sediment'
      diagUnits = 'micromol Fe/m^2/d'
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#ifdef RECOM_EXPORT_DIAGNOSTICS
      diagName  = 'EXPORTN '
      diagTitle = 'sinking flux of particulate organic nitrogen'
      diagUnits = 'micromol N/m^2/d'
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'EXPORTC '
      diagTitle = 'sinking flux of particulate organic carbon'
      diagUnits = 'micromol C/m^2/d'
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'EXPCALC '
      diagTitle = 'sinking flux of calcium carbonate'
      diagUnits = 'micromol C/m^2/d'
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'EXPORTSI'
      diagTitle = 'sinking flux of biogenic silica'
      diagUnits = 'micromol Si/m^2/d'
      diagCode  = 'SM      L1      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif

C----------------------------------------
C     transports of selected variables, example: total carbon
C----------------------------------------
      diagSufx  = 'CARB'
      diagUnits = 'mmolC/m^2/s^2'
C--   Advective flux:
      diagName  = 'ADVr'//diagSufx
      diagTitle = 'Vertical   Advective Flux of Carbon'
      diagCode  = 'WM      LR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid ) 
      diagName  = 'ADVx'//diagSufx
      diagTitle = 'Zonal      Advective Flux of Carbon'
      WRITE(diagCode,'(A,I3.3,A)') 'UU   ',diagNum+2,'MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid ) 
      diagName  = 'ADVy'//diagSufx
      diagTitle = 'Meridional Advective Flux of Carbon'
      WRITE(diagCode,'(A,I3.3,A)') 'VV   ',diagNum,'MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid ) 
C--   Diffusive flux:
      diagName  = 'DFrE'//diagSufx
      diagTitle = 'Vertical Diffusive Flux of Carbon'
     &     //' (Explicit part)'
      diagCode  = 'WM      LR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid ) 
      diagName  = 'DFxE'//diagSufx
      diagTitle = 'Zonal      Diffusive Flux of Carbon'
      WRITE(diagCode,'(A,I3.3,A)') 'UU   ',diagNum+2,'MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid ) 
      diagName  = 'DFyE'//diagSufx
      diagTitle = 'Meridional Diffusive Flux of Carbon'
      WRITE(diagCode,'(A,I3.3,A)') 'VV   ',diagNum,'MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid ) 
      
      diagName  = 'DFrI'//diagSufx
      diagTitle = 'Vertical Diffusive Flux of Carbon'
     &     //' (Implicit part)'
      diagCode  = 'WM      LR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid ) 

#endif /* ALLOW_DIAGNOSTICS */

      RETURN
      END
