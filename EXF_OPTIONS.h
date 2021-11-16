C $Header: /u/gcmpack/MITgcm/verification/lab_sea/code/ECCO_CPPOPTIONS.h,v 1.6 2004/03/17 23:30:51 dimitri Exp $
C $Name:  $

#ifndef EXF_OPTIONS_H
#define EXF_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#ifdef ALLOW_EXF

#include "CPP_OPTIONS.h"

C CPP flags controlling which code is included in the files that
C will be compiled.
C

C ********************************************************************
C ***             External forcing Package                         ***
C ********************************************************************
C 

C   Do more printout for the protocol file than usual.
#undef EXF_VERBOSE

C   Options that are required to use pkg/exf with pkg/seaice.
#define  ALLOW_ATM_TEMP
#define  ALLOW_ATM_WIND
#define  ALLOW_DOWNWARD_RADIATION
#define  ALLOW_BULKFORMULAE
#define  ALLOW_RUNOFF
#undef  EXF_READ_EVAP
#define ALLOW_EXF_BALANCE_FLUXES

C   Options that control relaxation terms.
#undef   ALLOW_CLIMTEMP_RELAXATION
#undef   ALLOW_CLIMSALT_RELAXATION
#define  ALLOW_CLIMSST_RELAXATION
#define  ALLOW_CLIMSSS_RELAXATION

C   Use spatial interpolation to interpolate
C   forcing files from input grid to model grid.
#define USE_EXF_INTERPOLATION
C   runoff is a special case for which one might want to bypass
C   interpolation from an input grid. This in now done in data.exf, not
C in the Code!
C#ifdef USE_EXF_INTERPOLATION
C# define USE_NO_INTERP_RUNOFF
C#endif

C In order to use Large and Yeager bulk formulae 2004
#define  ALLOW_BULK_LARGEYEAGER04

#endif /* ALLOW_EXF */
#endif /* EXF_OPTIONS_H */
