C $Header: /u/gcmpack/MITgcm_contrib/darwin2/pkg/darwin/SPECTRAL.h,v 1.5 2012/10/23 17:55:30 stephd Exp $
C $Name:  $

C SPECTRAL.h 
C description: spectral runtime parameters and fields

C  darwin_waves         :: 'central' wavelengths of wavebands (nm)
C  darwin_wavebands     :: waveband boundaries (nm)
      COMMON/SPECTRAL_PARAMS_R/darwin_waves, darwin_wavebands
      _RL darwin_waves(tlam)
      _RL darwin_wavebands(tlam+1)

C WtouEins :: W to uEin/s conversion factor
      COMMON/darwin_oasim/ WtouEins
      _RL WtouEins(tlam)

C oasim_ed :: spectral direct downwelling irradiance at surface read from file
C oasim_es :: spectral diffuse downwelling irradiance at surface read from file
#ifdef OASIM
      COMMON /SPECTRAL_INPUT/ oasim_ed, oasim_es
      _RL oasim_ed(1-OLx:sNx+OLx,1-OLy:sNy+OLy,tlam,nSx, nSy)
      _RL oasim_es(1-OLx:sNx+OLx,1-OLy:sNy+OLy,tlam,nSx, nSy)

      COMMON /SPECTRAL_INPUT_AUX/ tmp01ed, tmp01ed0, tmp01ed1
     &  , tmp02ed, tmp02ed0, tmp02ed1
     &  , tmp03ed, tmp03ed0, tmp03ed1
     &  , tmp04ed, tmp04ed0, tmp04ed1
     &  , tmp05ed, tmp05ed0, tmp05ed1
     &  , tmp06ed, tmp06ed0, tmp06ed1
     &  , tmp07ed, tmp07ed0, tmp07ed1
     &  , tmp08ed, tmp08ed0, tmp08ed1
     &  , tmp09ed, tmp09ed0, tmp09ed1
     &  , tmp10ed, tmp10ed0, tmp10ed1
     &  , tmp11ed, tmp11ed0, tmp11ed1
     &  , tmp12ed, tmp12ed0, tmp12ed1
     &  , tmp13ed, tmp13ed0, tmp13ed1
     &  , tmp01es, tmp01es0, tmp01es1
     &  , tmp02es, tmp02es0, tmp02es1
     &  , tmp03es, tmp03es0, tmp03es1
     &  , tmp04es, tmp04es0, tmp04es1
     &  , tmp05es, tmp05es0, tmp05es1
     &  , tmp06es, tmp06es0, tmp06es1
     &  , tmp07es, tmp07es0, tmp07es1
     &  , tmp08es, tmp08es0, tmp08es1
     &  , tmp09es, tmp09es0, tmp09es1
     &  , tmp10es, tmp10es0, tmp10es1
     &  , tmp11es, tmp11es0, tmp11es1
     &  , tmp12es, tmp12es0, tmp12es1
     &  , tmp13es, tmp13es0, tmp13es1
#ifndef ALLOW_EXF
     &  , oasim_ed0, oasim_ed1 , oasim_es0,  oasim_es1, 
      _RS oasim_ed0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,tlam,nSx, nSy)
      _RS oasim_ed1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,tlam,nSx, nSy)
      _RS oasim_es0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,tlam,nSx, nSy)
      _RS oasim_es1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,tlam,nSx, nSy)
#endif
C direct
      _RL tmp01ed(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp01ed0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp01ed1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp02ed(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp02ed0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp02ed1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp03ed(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp03ed0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp03ed1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp04ed(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp04ed0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp04ed1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp05ed(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp05ed0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp05ed1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp06ed(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp06ed0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp06ed1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp07ed(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp07ed0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp07ed1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp08ed(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp08ed0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp08ed1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp09ed(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp09ed0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp09ed1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp10ed(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp10ed0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp10ed1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp11ed(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp11ed0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp11ed1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp12ed(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp12ed0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp12ed1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp13ed(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp13ed0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp13ed1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
C diffuse
      _RL tmp01es(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp01es0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp01es1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp02es(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp02es0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp02es1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp03es(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp03es0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp03es1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp04es(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp04es0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp04es1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp05es(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp05es0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp05es1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp06es(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp06es0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp06es1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp07es(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp07es0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp07es1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp08es(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp08es0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp08es1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp09es(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp09es0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp09es1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp10es(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp10es0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp10es1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp11es(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp11es0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp11es1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp12es(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp12es0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp12es1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp13es(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp13es0(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
      _RL tmp13es1(1-OLx:sNx+OLx,1-OLy:sNy+OLy, nSx, nSy)
#endif
