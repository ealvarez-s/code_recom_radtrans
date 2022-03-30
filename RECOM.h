C $Header: /csys/software/cvsroot/REcoM/MITgcm/recom/RECOM.h,v 1.13 2008/04/17 14:25:53 mlosch Exp $
C $Name:  $

#ifdef ALLOW_RECOM

CBOP
C    !ROUTINE: RECOM.h
C    !INTERFACE:
C #include RECOM.h
 
C    !DESCRIPTION:
C Contains fields and parameters for biogeochemical model REcoM

C     C2K :: conversion of Celsius into Kelvin
      _RL     C2K
      PARAMETER ( C2K = 273.15 )
C     some time parameter
C     secondsPerHour :: 60*60 seconds
C     secondsPerDay  :: 60*60*24 seconds
C     noonTime :: the 12th hour of the day in seconds
      _RL     recom_secondsPerHour, recom_secondsPerDay
      _RL     noonTime
      PARAMETER ( recom_secondsPerHour =  3600. )
      PARAMETER ( recom_secondsPerDay  = 86400. )
      PARAMETER ( noonTime       = 12.*recom_secondsPerHour )
C     index of refraction for Water
      _RS     nWater
      PARAMETER ( nWater          = 1.33 )
C----------------------------------
C number of 2-dimensional and 3-dimensional REcoM diagnostics 
C (counting only those that need to be transferred between recom_forcing
C and recom_sms, because only those need to get an extra array)
C----------------------------------
      INTEGER     nDiagsREcoM3D, nDiagsREcoM2D
      PARAMETER ( nDiagsREcoM3D   = 9
#ifdef RECOM_MAREMIP
     &  + 28
#endif
#ifdef RECOM_PHOTODAMAGE
     &  + 8
#endif
#ifdef RECOM_MARSHALL
     &  + 4
#endif
#ifdef RECOM_WAVEBANDS
     &  + 3
#endif
     & )
      PARAMETER (nDiagsREcoM2D = 9
#ifdef RECOM_EXPORT_DIAGNOSTICS
     &  + 4
#endif
     & )

C     brockReedInsolation :: switch to compute insolation from atronomical
C                    formulas according to Brock(1981) and Reed(1977),
C                    instead of shortwave heating (the default).
C     diurnalInsolation   :: switch to compute a diurnal cycle of insolation
C                     defaults in .false.
C     computeHalos :: for 2D to 3D runs this must be true (default), because
C                     gchem_forcing is called after the blocking exchanges
C                     (need to convince Steph to change that!)
      LOGICAL brockReedInsolation
      LOGICAL diurnalInsolation
      LOGICAL computeHalos
      COMMON /RECOM_PARM_L/
     &     brockReedInsolation, diurnalInsolation, computeHalos
C     recom_tiny :: lower bound of all passive tracers (passed from
C                   REcoM_para_read)
C     The following parameters are only needed if
C     brockReedInsolation=.true. 
C     solarConstan - solar constant
C     albedo       - planetary albedo
C     parFrac      - photosynthetically reactive fraction of light
C     cloudCover   - could also be read in as a field, but ...
C     daysPerYear  - number of days in a model year, default = 360
C     fracTEP      - TEP (transparent exopolymeric particles) 
C                    fraction of EOM (extra-celluar matter)
C     constantIronSolubility - solubility of iron in water.  
C                    Multiply with iron dust decomposition to get the 
C                    dissolved iron available for bgc-processes.
      _RL recom_tiny 
      _RL solarConstant, parFrac, cloudCover
      _RL daysPerYear
      _RL fracTEP
      _RL constantIronSolubility
C     resuspension parameters 
      _RL recom_sizeFraction
      _RL recom_porosity
C     recomCritShearStress (in m^2/s^2) is really the minimum speed squared
C     that is required to resuspend material
      _RL recomCritShearStress
C     recomShearStressMax (in m^2/s^2) puts a cap on the resuspension flux
      _RL recomShearStressMax
      _RL recom_FeErosionRate 
C     currently not used
      _RL recomDragQuadratic
C     recom_windFile   :: file name of wind speeds
C     recom_pCO2File   :: file name of atmospheric pCO2
C     recom_ironFile   :: file name of aeolian iron flux
C     recom_silicaFile :: file name of surface silica
C     recom_CaCO3File  :: file name of initial benthic layer CaCO3
C     recom_hydroFile  :: file name of hydrothermal iron source

      COMMON /RECOM_FILENAMES/
     &        recom_windFile, recom_pCO2File, recom_iceFile,
     &        recom_CaCO3File,
#ifdef RECOM_WAVEBANDS
     &        darwin_waterabsorbFile, darwin_phytoabsorbFile,
     &        darwin_surfacespecFile, darwin_acdomFile,
     &        darwin_particleabsorbFile,
#ifdef OASIM					 
     &        darwin_oasim_edFile01, darwin_oasim_edFile02,
     &        darwin_oasim_edFile03, darwin_oasim_edFile04,
     &        darwin_oasim_edFile05, darwin_oasim_edFile06,
     &        darwin_oasim_edFile07, darwin_oasim_edFile08,
     &        darwin_oasim_edFile09, darwin_oasim_edFile10,
     &        darwin_oasim_edFile11, darwin_oasim_edFile12,
     &        darwin_oasim_edFile13,
     &        darwin_oasim_esFile01, darwin_oasim_esFile02,
     &        darwin_oasim_esFile03, darwin_oasim_esFile04,
     &        darwin_oasim_esFile05, darwin_oasim_esFile06,
     &        darwin_oasim_esFile07, darwin_oasim_esFile08,
     &        darwin_oasim_esFile09, darwin_oasim_esFile10,
     &        darwin_oasim_esFile11, darwin_oasim_esFile12,
     &        darwin_oasim_esFile13,
#endif
#endif /* RECOM_WAVEBANDS */
     &        recom_ironFile, recom_silicaFile

      CHARACTER*(MAX_LEN_FNAM) recom_windFile
      CHARACTER*(MAX_LEN_FNAM) recom_pCO2File
      CHARACTER*(MAX_LEN_FNAM) recom_iceFile
      CHARACTER*(MAX_LEN_FNAM) recom_ironFile
      CHARACTER*(MAX_LEN_FNAM) recom_silicaFile
      CHARACTER*(MAX_LEN_FNAM) recom_CaCO3File
#ifdef RECOM_WAVEBANDS
      CHARACTER*(MAX_LEN_FNAM) darwin_waterabsorbFile 
      CHARACTER*(MAX_LEN_FNAM) darwin_phytoabsorbFile
      CHARACTER*(MAX_LEN_FNAM) darwin_surfacespecFile
      CHARACTER*(MAX_LEN_FNAM) darwin_acdomFile
      CHARACTER*(MAX_LEN_FNAM) darwin_particleabsorbFile
#ifdef OASIM
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_edFile01
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_edFile02
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_edFile03
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_edFile04
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_edFile05
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_edFile06
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_edFile07
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_edFile08
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_edFile09
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_edFile10
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_edFile11
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_edFile12
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_edFile13
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_esFile01
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_esFile02
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_esFile03
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_esFile04
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_esFile05
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_esFile06
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_esFile07
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_esFile08
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_esFile09
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_esFile10
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_esFile11
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_esFile12
      CHARACTER*(MAX_LEN_FNAM) darwin_oasim_esFile13
#endif
#endif /* RECOM_WAVEBANDS */
  
      COMMON /RECOM_PARM_R/ 
     &     recom_tiny,
     &     solarConstant, parFrac, cloudCover,
     &     daysPerYear, fracTEP,
     &     constantIronSolubility,
     &     recom_sizeFraction, recom_porosity,
     &     recomCritShearStress, recom_FeErosionRate,
     &     recomShearStressMax, recomDragQuadratic
      NAMELIST /RECOM_PARM01/
     &     brockReedInsolation, 
     &     diurnalInsolation,
     &     computeHalos,
     &     recom_tiny,
     &     solarConstant, parFrac, cloudCover,
     &     daysPerYear, fracTEP,
     &     constantIronSolubility,
     &     recom_sizeFraction, recom_porosity,
     &     recomCritShearStress, recom_FeErosionRate,
     &     recomShearStressMax, recomDragQuadratic,
     &     recom_windFile, recom_pCO2File,
     &     recom_CaCO3File,
     &     recom_ironFile, recom_silicaFile
#if defined(RECOM_WAVEBANDS) && defined(OASIM)					 
     &	   ,darwin_oasim_edFile01, darwin_oasim_edFile02
     &     ,darwin_oasim_edFile03, darwin_oasim_edFile04
     &	   ,darwin_oasim_edFile05, darwin_oasim_edFile06
     &     ,darwin_oasim_edFile07, darwin_oasim_edFile08
     &     ,darwin_oasim_edFile09, darwin_oasim_edFile10
     &     ,darwin_oasim_edFile11, darwin_oasim_edFile12
     &     ,darwin_oasim_edFile13
     &	   ,darwin_oasim_esFile01, darwin_oasim_esFile02
     &     ,darwin_oasim_esFile03, darwin_oasim_esFile04
     &	   ,darwin_oasim_esFile05, darwin_oasim_esFile06
     &     ,darwin_oasim_esFile07, darwin_oasim_esFile08
     &     ,darwin_oasim_esFile09, darwin_oasim_esFile10
     &     ,darwin_oasim_esFile11, darwin_oasim_esFile12
     &     ,darwin_oasim_esFile13
#endif /* RECOM_WAVEBANDS */
					 
C     windSpeed     :: scalar wind speed at 10m above the surface
C     atmospCO2     :: atmopheric CO2
C     dicFlux       :: flux of dissolved inorganic surface flux
C     hPlus         :: variable for hPlus
C     feDustDeposit :: iron dust deposition at the surface, in mumol/m^2/s
C     recomC        :: specific sum of tracers (total carbon)
C     recomN        :: specific sum of tracers (total nitrate)
C     recomSi       :: specific sum of tracers (total silica)
C     benthicLayer  :: benthic layer that absorbs sediments;
C                      so far only C, N, and Si
C     feHydrothermal:: hydrothermal iron source in mumol Fe/m^3/s
  
      _RS  windSpeed     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS  atmospCO2     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS  dicFlux       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#ifdef RECOM_O2
      _RS  FluxO2       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
#endif
      _RS  feDustDeposit (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS  recomC        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS  recomN        (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C     unfortunately, recomSi must always be defined, because it is
C     part of a formal parameter list
      _RS  recomSi       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL  hPlus         (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  benthicLayerC (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  benthicLayerN (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CCV#ifdef ALLOW_RECOM_SILICATE
      _RL  benthicLayerSi(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CCV#endif /* ALLOW_RECOM_SILICATE */
      _RL  benthicLayerCaCO3(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
  
CTW increased sinking speed
#ifdef ALLOW_SINK_INCREASE
      _RL Vdetfast(Nr)
#endif /* ALLOW_SINK_INCREASE */

      COMMON /RECOM_FFIELDS/ 
     &     windSpeed, atmospCO2, dicFlux, feDustDeposit,
     &     recomC, recomN, recomSi
     &     , hPlus
     &     , benthicLayerC, benthicLayerN
CCV#ifdef ALLOW_RECOM_SILICATE
     &     , benthicLayerSi
CCV#endif /* ALLOW_RECOM_SILICATE */
#ifdef ALLOW_SINK_INCREASE
CTW vertically increasing sinking speed
     &     , Vdetfast
#endif /* ALLOW_SINK_INCREASE */
     &     , benthicLayerCaCO3
#ifdef RECOM_WAVEBANDS
     &     , alpha_mean
     &     , alpha_mean_dia
#endif /* RECOM_WAVEBANDS */
  
C     Auxilary fields for loading forcing fields
C     windspeed?     :: wind speed at 10m above surface
C     atmosp?        :: atmospheric CO2
C     feDustDeposit? :: iron dust deposition at the surface
      _RS windspeed0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS windspeed1  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS atmosp0     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS atmosp1     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CCV#ifdef ALLOW_RECOM_FEINPUT
      _RS feDustDeposit0 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RS feDustDeposit1 (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
CCV#endif /* ALLOW_RECOM_FEINPUT */
#ifdef RECOM_WAVEBANDS
      _RL  alpha_mean
      _RL  alpha_mean_dia
#endif /* RECOM_WAVEBANDS */
     
      COMMON /RECOM_FFIELDS_AUX/
     &     windspeed0, windspeed1
     &     , atmosp0, atmosp1
CCV#ifdef ALLOW_RECOM_FEINPUT
     &     , feDustDeposit0, feDustDeposit1
CCV#endif /* ALLOW_RECOM_FEINPUT */ 
     
#ifdef ALLOW_TIMEAVE 
      COMMON /TAVE_RECOM/
     &       WCtave,  CDiffRtave
     &     , WNtave,  NDiffRtave
     &     , WSitave, SiDiffRtave
      _RL  WCtave     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL  CDiffRtave (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL  WNtave     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL  NDiffRtave (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL  WSitave    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL  SiDiffRtave(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
#endif /* ALLOW_TIMEAVE */

#ifdef RECOM_ATMOSPCO2_HISTORY
C--------------------------------------------------------------------------
C atmospheric history of prescribed pCO2 (yearly data)
C     co2atmos(i,1): year
C     co2atmos(i,1): pCO2 (uatm)
C--------------------------------------------------------------------------
      COMMON /ATMOS_RECOM/ co2atmos
      _RL co2atmos(5000,2) 
#endif

C--------------------------------------------------------------------------
C   relative tracer indices
C--------------------------------------------------------------------------

      integer nvars_base, nvars_zoo2, nvars_cdom,  
     &        nvars_iron_photo, nvars_iron_lig, nvars_iron_lig2, 
     &        nvars_iron_coll, nvars_iron_lith,
     &        nvars_photo_damage  
      integer idin, idic, ialk, iphyn, iphyc, ipchl, idetn,  
     &        idetc, ihetn, ihetc, idon, idoc, 
     &        idian, idiac, idchl, idiasi, idetsi, isi, ife
#ifdef RECOM_O2
      integer ioxy
#endif
      integer iphycalc, idetcalc
#ifdef RECOM_ZOO2
      integer izoo2n, izoo2c
#endif
#ifdef RECOM_CDOM
      integer icdom
#endif
#ifdef RECOM_MARSHALL
      integer id1, id1d
#endif  

#ifdef RECOM_O2
      parameter (nvars_base=22)
      parameter (idin=1, idic=2, ialk=3, iphyn=4, iphyc=5, ipchl=6,
     &  idetn=7, idetc=8, ihetn=9, ihetc=10, idon=11, idoc=12,
     &  idian=13, idiac=14, idchl=15, idiasi=16, idetsi=17, isi=18,
     &  ife=19, ioxy=20, iphycalc=21, idetcalc=22)
#else
      parameter (nvars_base=21)
      parameter (idin=1, idic=2, ialk=3, iphyn=4, iphyc=5, ipchl=6,
     &  idetn=7, idetc=8, ihetn=9, ihetc=10, idon=11, idoc=12,
     &  idian=13, idiac=14, idchl=15, idiasi=16, idetsi=17, isi=18,
     &  ife=19, iphycalc=20, idetcalc=21)
#endif

#ifdef RECOM_ZOO2
      parameter (nvars_zoo2=2)
      parameter (izoo2c=nvars_base+1,
     &    	 izoo2n=nvars_base+2)
#else
      parameter (nvars_zoo2=0)
#endif
#ifdef RECOM_CDOM
      parameter (nvars_cdom=1)
      parameter (icdom=nvars_base+nvars_zoo2+1)
#else
      parameter (nvars_cdom=0)
#endif
  
      parameter (nvars_iron_photo=0)
      parameter (nvars_iron_lig=0)
      parameter (nvars_iron_lig2=0)
      parameter (nvars_iron_coll=0)
      parameter (nvars_iron_lith=0)
  
#ifdef RECOM_MARSHALL
      parameter (nvars_photo_damage=2)
      parameter (id1 =  nvars_base+nvars_zoo2+nvars_cdom +
     &    nvars_iron_photo+nvars_iron_lig+
     &    nvars_iron_lig2+nvars_iron_coll+
     &    nvars_iron_lith+1)
      parameter (id1d = nvars_base+nvars_zoo2+nvars_cdom +
     &    nvars_iron_photo+nvars_iron_lig+
     &    nvars_iron_lig2+nvars_iron_coll+
     &    nvars_iron_lith+2)
#else
      parameter (nvars_photo_damage=0)
#endif
CEOP
#endif /* ALLOW_RECOM */
