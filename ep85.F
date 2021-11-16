C$Header: /csys/software/cvsroot/REcoM/pkg/recom/ep85.F,v 1.9 2006/12/09 15:45:44 mlosch Exp $
#include "PTRACERS_OPTIONS.h"
#include "GCHEM_OPTIONS.h"
#include "RECOM_OPTIONS.h"
      subroutine ep85(zdepth,daylen,kappa,vp,
     &     alpha,irrad_top,irrad_bot,light_growth)   
C=====================================================================
C
C compute daily averaged light-limited growth rate analytical integration 
C over layer thickness 'zdepth' and when with 'daylen'
C after Evans and Parslow (1985) ==> EP85
C
      implicit none

C     input variables
      _RL zdepth
      _RL daylen
      _RL kappa
      _RL vp, alpha
      _RL irrad_top, irrad_bot
C     output variables
      _RL light_growth
C     local variables
      _RL f(2,2)
      _RL y(2)
C      _RL u1, u2
      _RL tau, top, bot 
      _RL tiny
      PARAMETER( tiny=2.23e-16 )
      integer i

      tau  = 0.5*daylen
    
      top  = max(tiny,irrad_top)
      bot  = max(tiny,irrad_bot)

CML      y(1) = (vp*tau)/(alpha*top)
CML      y(2) = (vp*tau)/(alpha*bot)
      y(1) = max(tiny,(vp*tau)/(max(tiny,alpha)*top))
      y(2) = max(tiny,(vp*tau)/(max(tiny,alpha)*bot))

      do i = 1, 2
  
        f(i,1) = sqrt(y(i)**2.+tau**2.)-tau*log(
     &        (tau+sqrt(y(i)**2.+tau**2.))/y(i))

        f(i,2) = y(i)

      end do
      
       light_growth=abs(2.*vp)/(zdepth*kappa)*
     &     (f(2,1)-f(1,1)-f(2,2)+f(1,2))
      
       light_growth = max(tiny,light_growth)      

! approximation as reported in Garcon and Evans (1997)
!         top=top*daylen*0.5 ! total light during the day 
!         u1=(2.*alpha*top)/(daylen*vp)
!         u2=u1*exp(-kappa*zdepth)
!         F(1,1)=(0.555588*u1+0.004926*u1*u1)/(1.+0.188721*u1)
!         F(2,2)=(0.555588*u2+0.004926*u2*u2)/(1.+0.188721*u2)
!         light_growth=(vp*daylen)/(kappa*zdepth)*(F(1,1)-F(2,2))
      
      end
!=====================================================================
