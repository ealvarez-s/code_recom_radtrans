C$Header: /csys/software/cvsroot/REcoM/MITgcm/recom/RECOM_OPTIONS.h,v 1.4 2007/02/09 15:20:14 mlosch Exp $
C$Name:  $
CBOP
C    !ROUTINE: RECOM_OPTIONS.h
C    !INTERFACE:

C    !DESCRIPTION:
C 
C CPP flags controlling which code included in the REcoM files that
C will be compiled.
C

#ifndef RECOM_OPTIONS_H
#define RECOM_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

#ifdef ALLOW_RECOM

#define ALLOW_MITGCM

C---------------------
CCV #define RECOM_2CLASSES
#define RECOM_GRAZING_VARIABLE_PREFERENCE
CCV june 2017: RECOM_CALCIFICATION is retired / always switched on
CCV #define RECOM_CALCIFICATION
#define RECOM_ATMOSPCO2_HISTORY
C---------------------

CTW increased sinking speed 
#define ALLOW_SINK_INCREASE

CCV Include/exclude oxygen as state variable
#undef RECOM_O2

CCV include/exclude a second zooplankton group
#undef RECOM_ZOO2

CCV include/exclude colored dissolved organic matter (CDOM) as
CCV state variable
#define RECOM_CDOM

CCV Short iron recycling loop 
C   Excreted organic iron is directly put into bioavailable iron
#define RECOM_FE_RECYCLING_SHORT
#define RECOM_CONSTANT_FE2N

CCV vertical fluxes of PON, POC, Calc and BSi are defined as diagnostics
#define RECOM_EXPORT_DIAGNOSTICS

CCV many different additional diagnostics for the MAREMIP model intercomparison
#define RECOM_MAREMIP

C- settings for how complex iron chemistry should be
#define RECOM_IRON_BENTHOS

CEA Chlorophyll loss term dependent on PAR and Chl:C		    
#define RECOM_PHOTODAMAGE

CEA Include/exclude non-photosynthetic pigments as state variable
#define RECOM_MARSHALL		    
		    
CEA Radiative Transfer Model
#define RECOM_WAVEBANDS
#define RECOM_RADTRANS
#define OASIM		    
#define RECOM_CALC_ACDOM
#define RECOM_CALC_APART
#define RECOM_2GROUPS		    
#define RECOM_BMASS
#define RECOM_CALC_REFLEC  
CEA undef to keep chl-specific phyto absorption constant
#define RECOM_CALC_APHYT
  
C catch bad CPP-flag combination
#ifdef DIC_BIOTIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
This text is put here deliberately:
RECOM should not be defined together with DIC_BIOTIC     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#endif /* DIC_BIOTIC */
#endif /* ALLOW_RECOM */ 

C NOTE: none of the flags below will be used anywhere, because 
C       RECOM_OPTIONS.h will not be include into the REcoM model. 
C       Use, e.g.
C       set DEFINES='-DEULER -DWRITE -DPARVE -DASTRO_RAD' 
C       in the .genmakerc of REcoM/bin instead.
C#undef EULER
C#undef WRITE
C#undef PARVE 
C#undef ASTRO_RAD
C#undef EP85  
C#undef TEP_AGG
#endif /* RECOM_OPTIONS_H */
