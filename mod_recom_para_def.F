C$Header: /csys/software/cvsroot/REcoM/pkg/recom/mod_recom_para_def.F,v 1.7 2006/12/09 15:45:44 mlosch Exp $
#include "PTRACERS_OPTIONS.h"
#include "GCHEM_OPTIONS.h"
#include "RECOM_OPTIONS.h"
#ifdef ALLOW_MODULES
      module mod_REcoM_para_def 
!---------------------------------------------------------------------
#include "RECOM_PARAMS.h"
      end module mod_REcoM_para_def
!---------------------------------------------------------------------
#else 
      subroutine mod_recom_para_def_dummy
C     dummy routine
      end
#endif /* ALLOW_MODULES */
