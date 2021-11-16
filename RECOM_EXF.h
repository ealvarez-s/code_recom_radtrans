C $Header: /csys/software/cvsroot/REcoM/MITgcm/recom/RECOM_EXF.h,v 1.4 2008/04/17 14:25:52 mlosch Exp $
C $Name:  $

#ifdef ALLOW_EXF
      NAMELIST /RECOM_PARM02/
     &     pco2startdate1, pco2startdate2, pco2period,
     &     pco2repeatperiod,
     &     feduststartdate1, feduststartdate2, fedustperiod,
     &     fedustrepeatperiod
#ifdef OASIM
     &     ,OASIMstartdate1, OASIMstartdate2, OASIMperiod
     &     ,OASIMrepeatperiod
     &     ,exf_inscal_OASIM, exf_outscal_OASIM
#endif
#ifdef USE_EXF_INTERPOLATION
     &     ,pco2_lon0, pco2_lon_inc, pco2_lat0, pco2_lat_inc,
     &     fedust_lon0, fedust_lon_inc, fedust_lat0, fedust_lat_inc,
     &     pco2_nlon, pco2_nlat, fedust_nlon, fedust_nlat
#ifdef OASIM
     &     ,OASIM_lon0, OASIM_lon_inc, OASIM_nlon
     &     ,OASIM_lat0, OASIM_lat_inc, OASIM_nlat
#endif
#endif /* USE_EXF_INTERPOLATION */
      INTEGER pco2startdate1
      INTEGER pco2startdate2
      _RL     pco2startdate
      _RL     pco2period
      _RL     pco2repeatperiod

      INTEGER feduststartdate1
      INTEGER feduststartdate2
      _RL     feduststartdate
      _RL     fedustperiod
      _RL     fedustrepeatperiod
#ifdef OASIM
      INTEGER OASIMstartdate1
      INTEGER OASIMstartdate2  
      _RL     OASIMstartdate
      _RL     OASIMperiod
      _RL     OASIMrepeatperiod
      _RL     exf_inscal_OASIM
      _RL     exf_outscal_OASIM
#endif
      COMMON /RECOM_EXF_PARMS_I/
     &     pco2startdate1, pco2startdate2,
     &     feduststartdate1, feduststartdate2
#ifdef OASIM
     &     ,OASIMstartdate1, OASIMstartdate2
#endif
      COMMON /RECOM_EXF_PARMS_R/
     &     pco2startdate, pco2period, pco2repeatperiod,
     &     feduststartdate, fedustperiod, fedustrepeatperiod
#ifdef OASIM
     &     ,OASIMstartdate, OASIMperiod, OASIMrepeatperiod
     &     ,exf_inscal_OASIM, exf_outscal_OASIM
#endif
#ifdef USE_EXF_INTERPOLATION
      _RL pco2_lon0, pco2_lon_inc
      _RL pco2_lat0, pco2_lat_inc(MAX_LAT_INC)
      INTEGER pco2_nlon, pco2_nlat
      _RL fedust_lon0, fedust_lon_inc
      _RL fedust_lat0, fedust_lat_inc(MAX_LAT_INC)
      INTEGER fedust_nlon, fedust_nlat
#ifdef OASIM
      _RL OASIM_lon0, OASIM_lon_inc
      _RL OASIM_lat0, OASIM_lat_inc(MAX_LAT_INC)
      INTEGER OASIM_nlon
      INTEGER OASIM_nlat
#endif
      COMMON /RECOM_INTERPOLATION/
     & pco2_lon0, pco2_lon_inc,
     & pco2_lat0, pco2_lat_inc,
     & pco2_nlon, pco2_nlat,
     & fedust_lon0, fedust_lon_inc,
     & fedust_lat0, fedust_lat_inc,
     & fedust_nlon, fedust_nlat
#ifdef OASIM
     & ,OASIM_lon0, OASIM_lon_inc,
     & OASIM_lat0, OASIM_lat_inc,
     & OASIM_nlon, OASIM_nlat
#endif
#endif /* USE_EXF_INTERPOLATION */
#endif /* ALLOW_EXF */
