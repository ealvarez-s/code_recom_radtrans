C$Header: /csys/software/cvsroot/REcoM/pkg/recom/RECOM_PARAMS.h,v 1.19 2008/04/17 14:25:51 mlosch Exp $
C$Name:  $
C
C   model run identification: 
C
C     attenuation coefficient [m^{-1}]
C      _RL kappa
C     chlorophyll specific attenuation coefficients [m^{-1} (mg CHL)^{-1}]
      _RL a_chl
 
C-----
C (small) phytoplankton parameters
C-----
C     photosynth_growth
C     CHL-specific intial slope of P-I curve, unit depends on unit of 
C     radiation   
      _RL alpha             
C     Temperature dependent maximum of C-specific rate of photosynthesis [d^{-1}] 
      _RL P_cm
C     CO2 sensitivity of photosynthesis: not used so far!
      _RL k_CO2 
C     photosynth_growth
C
#ifdef RECOM_MARSHALL   
C     photoinhibition
C     relative amount of D1 that keeps QY=QYmax  [rel]
      _RL Drel      
C     chlorophyll absortion cross section [m^{-2} (mg CHL)^{-1}]
      _RL astar
C     Target size for photoinhibition [m^{-2} (J)^{-1}]
      _RL k_deg 
C     Maximum repair rate [d^{-1}] 
      _RL r_max 
C     Half saturation constant for repair [same as DD1] 
      _RL k_rep 
C     photoinhibition
#ifndef RECOM_WAVEBANDS
C     maximum quantum yield of photosynthesis [mmol C (J)^{-1}]   
      _RL QYmax
#endif
#endif  
C  
C     N_utilization
C     Half-saturation constant for nitrate uptake [mmol N m^{-3}] 
      _RL k_din
C     scaling factor for temperature dependent maximum of C-specific 
C     nitrogen uptake 
      _RL V_cm_fact         
C     N_utilization  

C     respi_degrad_exud
C     Maintenance respiration rate constant [d^{-1}]
      _RL res_phy
C     Cost of biosynthesis [mmol C (mmol N)^{-1}]   
      _RL biosynth
C     CHL degradation rate constant [d^{-1}]
      _RL deg_CHL
C     phytoplankton loss of organic N compounds [d^{-1}]
      _RL lossN
C     phytoplankton loss of carbon [d^{-1}]
      _RL lossC
C     excretion of zooplankzon 
      _RL lossN_z
      _RL lossC_z

#ifdef RECOM_ZOO2
C     excretion of zoo2 
      _RL lossN_z2
      _RL lossC_z2
#endif
#ifdef RECOM_CDOM
      _RL fcdom, rho_cdom, phot_cdom, kphot_cdom
#endif
      
C     respi_degrad_exud 

C-----
C Now the same parameters again for Diatoms
C-----
C     photosynth_growth
C     CHL-specific intial slope of P-I curve, unit depends on unit of 
C     radiation   
      _RL alpha_d             
C     Temperature dependent maximum of C-specific rate of photosynthesis [d^{-1}] 
      _RL P_cm_d
C     CO2 sensitivity of photosynthesis: not used so far!
      _RL k_CO2_d 
C     photosynth_growth
C
#ifdef RECOM_MARSHALL  
C     photoinhibition
C     chlorophyll absortion cross section [m^{-2} (mg CHL)^{-1}]
      _RL astar_d
C     target size for photoinhibition, [m^{-2} (J)^{-1}]
      _RL k_deg_d 
C     maximum repair rate, [d^{-1}] 
      _RL r_max_d 
C     half saturation constant for repair, [same as DD1] 
      _RL k_rep_d 
C     photoinhibition
#ifndef RECOM_WAVEBANDS
C     maximum quantum yield of photosynthesis, [mmol C (J)^{-1}]   
      _RL QYmax_d
#endif
#endif  
C
C     N_utilization
C     Half-saturation constant for nitrate uptake [mmol N m^{-3}] 
      _RL k_din_d
C     scaling factor for temperature dependent maximum of C-specific 
C     nitrogen uptake 
      _RL V_cm_fact_d         
C     N_utilization  

C     respi_degrad_exud
C     Maintenance respiration rate constant [d^{-1}]
      _RL res_phy_d
C     Cost of biosynthesis [mmol C (mmol N)^{-1}]   
      _RL biosynth_d
C     CHL degradation rate constant [d^{-1}]
      _RL deg_CHL_d
C     phytoplankton loss of organic N compounds [d^{-1}]
      _RL lossN_d
C     phytoplankton loss of carbon [d^{-1}]
      _RL lossC_d
C     respi_degrad_exud 


C     start remineralisation 
C     temperature dependent remineralisation rate of detritus  [d^{-1}]
      _RL reminN, reminC, reminSi
C     respiration by heterotrophs and mortality (loss to detritus)  [d^{-1}]
      _RL res_het, recip_res_het
      _RL loss_het
#ifdef RECOM_ZOO2
C     respiration by zoo2  and mortality (loss to detritus)  [d^{-1}]
      _RL res_zoo2, recip_res_zoo2
      _RL loss_zoo2
#endif     
C     temperature dependent N degradation of extracellular organic N (EON) [d^{-1}]
      _RL rho_N
C     temperature dependent C degradation of extracellular organic C (EOC) [d^{-1}]
      _RL rho_C1
C     temperature dependent C degradation of TEP-C [d^{-1}]
      _RL rho_C2
C     end remineralisation 

C     start grazing&aggregation
C------------------------------------
C The following 2 Variables are only used when 2 classes of phytoplanton are around. 
C Nevertheless we define them here in any case
C     (maximum) diatom preference [dimensionless, smaller than 1]
      _RL pzdia 
C     half-saturation parameter for grazing preference for diatoms [(mmol N m^{-3})^2]
      _RL sdiasq
C------------------------------------
C     maximum grazing loss parameter [mmol N m^{-3} d^{-1}] 
      _RL graz_max
C     half saturation grazing loss [mmol N^{2} m^{-6}]
      _RL epsilon
C     grazing efficiency (fraction of grazing flux into zooplankton pool)
      _RL grazEff
C     maximum aggregation loss parameters [m^{3} mmol N^{-1} d^{-1}]
      _RL agg_PP, agg_PD
CMLC     half saturation parameter for collision and stickiness probability 
CMLC     (TEP dependend) [mmol C m^{-3}] (TEP-C concentration) 
CML      _RL k_TEP

#ifdef RECOM_ZOO2
C      start grazing  zoo2
C------------------------------------
C The following 2 Variables are only used when 2 classes of phytoplanton are around. 
C Nevertheless we define them here in any case
C     (maximum) diatom preference [dimensionless, smaller than 1]
      _RL pzdia2 
C     half-saturation parameter for grazing preference for diatoms [(mmol N m^{-3})^2]
      _RL sdiasq2
C------------------------------------
C     maximum grazing loss parameter [mmol N m^{-3} d^{-1}] 
      _RL graz_max2
C     half saturation grazing loss [mmol N^{2} m^{-6}]
      _RL epsilon2
C     grazing efficiency (fraction of grazing flux into zooplankton pool)
      _RL grazEff2
#endif

C     phytoplankton sinking velocity [m/d]
      _RL Vphy
      _RL Vdia
C     sinking velocity of detritus [m/d]
      _RL Vdet
#ifdef ALLOW_SINK_INCREASE
C------------------------------------------------------
C The following 2 variables are only used when vertically increasing 
C detritus sinking speed is being used
CTW   increaing  sinking velocity of detritus [d^{-1};m/d]
      _RL Vdet_a
      _RL Vdet_b
C-------------------------------------------------------
#endif

C     decay rate of detritus in the benthic layers [d^{-1}]
      _RL decayRateBenN
      _RL decayRateBenC
      _RL decayRateBenSi
      _RL decayRateBenCa
C     end grazing&aggregation

C     start TEP stuff  
C     TEP fraction (f_TEP-x) found in POM, dimensionless
      _RL f_TEP, x
C     stickiness for PCHO-PCHO 
      _RL agg_PCHO 
C     stickiness for TEP-PCHO 
      _RL agg_TEP
C     end TEP stuff
   
C     start silica parameters
C     Cost of biosynthesis [mmol C (mmol Si)^{-1}]   
      _RL biosynthSi
C     silica saturation constant [mmol Si m^{-3}]
      _RL k_Si
C     minimum Si:N ratio 
      _RL     SiNmin
C     maximum Si:N ratio 
      _RL     SiNmax
C     minimum Si:C ratio 
      _RL     SiCmin
C     maximum Si:C ratio 
      _RL     SiCmax
      _RL     SiCuptakeRatio
C     end silica and iron related parameters
      logical  FeLimit   
#ifdef RECOM_CONSTANT_FE2N
C     iron to carbon ratio [micromol Fe mmol C^{-1}
      _RL     Fe2N
#ifdef RECOM_IRON_BENTHOS
      _RL     Fe2N_benthos
#endif
C     iron to carbon ration [micromol Fe mmol C^{-1}
      _RL     Fe2N_d
#else    
C     iron to carbon ratio [micromol Fe mmol C^{-1}
      _RL     Fe2C
#ifdef RECOM_IRON_BENTHOS
      _RL     Fe2C_benthos
#endif
C     iron to carbon ration [micromol Fe mmol C^{-1}
      _RL     Fe2C_d
#endif
C     iron saturation constant [micromol m^{-3}]
      _RL     k_Fe
C     iron saturation constant [micromol m^{-3}]
      _RL     k_Fe_d
C     start iron chemistry parameters [d^{-1}] [order 0.0005-0.005 per day]
C     rate constant for scavenging on detritus [(mmol C m^{-3})^{-1} d^{-1}]
      _RL     kScavFe
C     rate constant for scaveing on dust particles (kg m^{-3})^{-1} d^{-1})]
C     -RL     kScavFe_dust
C     total free ligand [mumol m^{-3}] [order 1]
      _RL     totalLigand
C     ligand-free iron stability constanty [m^{3}/mumol] [order 100]
      _RL     ligandStabConst
C     end iron chemistry parameters
C     some external parameters (irrelevant to recom_sms and below
C     deviation from Si2N for initial conditions  
      _RL     ini_Si
C     nmol Fe m^{-3} initial in upper 3 layers
      _RL     ini_Fe

      
C     factor for initialisation (only relevant for mesocosm simulations)
      _RL ini  

C     start fixed parameters 
C     [m^{-1}] attenuation due to water 
      _RL k_w
C     slope of the linear part of the Arrhenius function [Kelvin]
      _RL Ae
C     Reference temperature [Kelvin]
      _RL recom_Tref

C-----
C small phytoplanton quota 
C-----
C     Minimum cell quota of nitrogen (N:C) [mmol N (mmol C)^{-1}]
      _RL NCmin
C     Maximum cell quota of nitrogen (N:C) [mmol N (mmol C)^{-1}]
      _RL NCmax
      _RL NCuptakeRatio
C     Maximum CHL a : N ratio [mg CHL (mmol N)^{-1}] = 0.3 gCHL gN^-1
      _RL CHL_N_max
C-----
C diatom quota 
C-----
C     Minimum cell quota of nitrogen (N:C) [mmol N (mmol C)^{-1}]
      _RL NCmin_d
C     Maximum cell quota of nitrogen (N:C) [mmol N (mmol C)^{-1}]
      _RL NCmax_d
      _RL NCuptakeRatio_d
C     Maximum CHL a : N ratio [mg CHL (mmol N)^{-1}] = 0.3 gCHL gN^-1
      _RL CHL_N_max_d

C     constant Redfield C:N ratio 
      _RL redfield
C     slopes for limiting functions
      _RL NMinSlope
      _RL NMaxSlope
      _RL SiMinSlope
      _RL SiMaxSlope
#ifdef RECOM_O2
      _RL redO2C
      _RL sox1, sox2, sox3,sox4
      _RL oA0,oA1,oA2,oA3,oA4,oA5
      _RL oB0,oB1,oB2,oB3
      _RL oC0
#endif
      _RL calc_prod_ratio 
      _RL calc_diss_guts
      _RL calc_diss_rate

C     lower limit for variables (sometimes it is necessary for numerical
C     stability reasons to increase this to 1e-8 or even 1e-6 )
      _RL tiny              

C     some model parameters  
C     ---------------------

CML      _RL    dt             ! time step within REcoM unit 'day'
C     internal biological loop size (the model is stepped bio_step times 
C     during one physical time step
      integer bio_step
CML      integer bgc_num        ! total number of biochemical variables

C     parameters
      _RL c0, c1, p5, c2 
      parameter ( c0 = 0.D0, c1 = 1.D0, p5 = 0.5D0, c2 = 2.D0 )
C     seconds
      _RL one_day           
      parameter ( one_day = 86400. ) 
C     seconds
      _RL one_hour          
      parameter ( one_hour = 3600. )
C     celsius2K:: conversion of Celsius into Kelvin
      _RL recom_celsius2K
      parameter ( recom_celsius2K = 273.15 )

      integer unset_integer
      parameter ( unset_integer = 123456789 )
      _RL    unset_real   
      parameter ( unset_real    = 123456.789)

      logical useReminSiT
      logical TEPaggregation
      logical EvansParslow
      logical SINKadv
      logical benthicLayer
      logical write_flag

#ifndef ALLOW_MODULES
      COMMON /RECOM_INTERNAL_PARM_L/
     &                     useReminSiT,
     &                     TEPaggregation,
     &                     EvansParslow,
     &                     SINKadv,  
     &                     benthicLayer, 
     &                     FeLimit, 
     &                     write_flag
      
      COMMON /RECOM_INTERNAL_PARM_R/
     &                     a_chl,  
     &                     alpha,
     &                     P_cm,
     &                     k_din,
     &                     V_cm_fact,
     &                     res_phy,
     &                     deg_CHL,
     &                     lossN,
     &                     lossC,
     &                     alpha_d,
     &                     P_cm_d,
     &                     k_din_d,
     &                     V_cm_fact_d,
     &                     res_phy_d,
     &                     deg_CHL_d,
     &                     lossN_d,
     &                     lossC_d,
     &                     rho_N,
     &                     rho_C1,
     &                     rho_C2,
     &                     lossN_z,
     &                     lossC_z,
#ifdef RECOM_ZOO2
     &                     lossN_z2,
     &                     lossC_z2,			     
#endif      
     &                     decayRateBenN,
     &                     decayRateBenC,
     &                     decayRateBenSi,
#ifdef RECOM_O2
     &                     oA0,oA1,oA2,oA3,oA4,oA5,
     &                     oB0,oB1,oB2,oB3,
     &                     oC0, sox1, sox2, sox3, sox4,
     &                     redO2C,
#endif
#ifdef RECOM_CDOM
     &                     fcdom, rho_cdom, phot_cdom, kphot_cdom, 
#endif
     &                     decayRateBenCa,
#ifdef RECOM_MARSHALL
     &                     Drel,
     &                     astar,
     &                     k_deg,
     &                     r_max,
     &                     k_rep,
     &                     astar_d,
     &                     k_deg_d,
     &                     r_max_d,
     &                     k_rep_d,
#ifndef RECOM_WAVEBANDS
     &                     QYmax,
     &                     QYmax_d,
#endif
#endif
     &                     f_TEP, 
     &                     x,
     &                     agg_PCHO, 
     &                     agg_TEP,
     &                     agg_PD,
     &                     agg_PP,
     &                     reminN,
     &                     reminC,
     &                     res_het,
     &                     recip_res_het,
     &                     loss_het,
#ifdef RECOM_ZOO2				     
     &                     res_zoo2,
     &                     recip_res_zoo2,
     &                     loss_zoo2,
#endif
     &                     pzdia, 
#ifdef RECOM_GRAZING_VARIABLE_PREFERENCE
     &                     sdiasq,
#endif
#ifdef RECOM_ZOO2				     
     &                     pzdia2, 
#ifdef RECOM_GRAZING_VARIABLE_PREFERENCE
     &                     sdiasq2,
#endif
#endif				     
     &                     graz_max,
     &                     epsilon,
     &                     grazEff,
#ifdef RECOM_ZOO2
     &                     graz_max2,
     &                     epsilon2,
     &                     grazEff2,
#endif				     
     &                     Vphy,
     &                     Vdia,
     &                     Vdet,
#ifdef ALLOW_SINK_INCREASE
     &                     Vdet_a,
     &                     Vdet_b, 
#endif
     &                     k_CO2,
     &                     ini,
     &                     biosynth,
     &                     biosynthSi,
     &                     k_w,
     &                     Ae,
     &                     recom_Tref, 
     &                     redfield,  
     &                     NCuptakeRatio, NCmax, NCmin,
     &                     CHL_N_max, 
     &                     NCuptakeRatio_d, NCmax_d, NCmin_d,
     &                     CHL_N_max_d, 
     &                     k_Si,
     &                     reminSi,
     &                     SiNmin, SiNmax,
     &                     SiCmin, SiCmax, SiCuptakeRatio,
#ifdef RECOM_CONSTANT_FE2N
     &                     Fe2N,
#ifdef RECOM_IRON_BENTHOS
     &                     Fe2N_benthos,
#endif
     &                     Fe2N_d,
#else
     &                     Fe2C,
#ifdef RECOM_IRON_BENTHOS
     &                     Fe2C_benthos,
#endif
     &                     Fe2C_d,
#endif
     &                     k_Fe, 
     &                     k_Fe_d, 
     &                     kScavFe,
     &                     totalLigand,
     &                     ligandStabConst,
     &                     ini_Si,
     &                     ini_Fe,
     &                     NMinSlope, NMaxSlope, SiMinSlope, SiMaxSlope,
     &                     calc_prod_ratio, 
     &                     calc_diss_guts,
     &                     calc_diss_rate,
     &                     tiny
      
      COMMON /RECOM_INTERNAL_PARM_I/
     &     bio_step

#ifdef RECOM_ATMOSPCO2_HISTORY
      integer recom_pco2_int1, recom_pco2_int2, 
     &        recom_pco2_int3, recom_pco2_int4
      COMMON /RECOM_ATMOSPCO2_PARAMS/ 
     &        recom_pco2_int1, recom_pco2_int2, 
     &        recom_pco2_int3, recom_pco2_int4
#endif /* RECOM_ATMOSPCO2_HISTORY */

#endif /* not ALLOW_MODULES */

#ifndef ALLOW_MITGCM
C     parameters that are only needed for support of vectorization
C     do not change in 1D configuration!!!!
      integer  sNx, sNy, Olx, Oly, nSx, nSy
      parameter ( sNx = 1, nSx = 1, Olx = 0 )
      parameter ( sNy = 1, nSy = 1, Oly = 0 )
#endif /* ALLOW_MITGCM */

