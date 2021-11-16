C$Header: /csys/software/cvsroot/REcoM/pkg/recom/recom_sinking.F,v 1.11 2007/10/01 14:49:15 mlosch Exp $
C$Name:  $
#include "PTRACERS_OPTIONS.h"
#include "GCHEM_OPTIONS.h"
#include "RECOM_OPTIONS.h"

      subroutine recom_sinking(
     I     dt, recipDz, recipDzC, recip_hFacC, maskC, ws, c,
     O     sink, 
#ifdef RECOM_EXPORT_DIAGNOSTICS
     O     export,
#endif
#ifdef RECOM_MAREMIP
     O     export3d,
#endif
     I     n_max, iMin, iMax, jMin, jMax, bi, bj )
C---------------------------------------------------------------------     
C     this code is adopted from the MITgcm routine: 
C     |==========================================================|
C     | SUBROUTINE GAD_DST3_ADV_R                                |
C     | o Compute Vertical advective Flux of Tracer using        |
C     |   3rd Order DST Scheme with flux limiting               |
C     |==========================================================|
C
C     Author: Martin Losch, July, 2003

#ifdef ALLOW_MODULES
      use mod_REcoM_para_def      
      implicit none
#else /* not ALLOW_MODULES */
      implicit none
#include "RECOM_PARAMS.h"
#ifdef ALLOW_MITGCM
#include "SIZE.h"
#endif 
#endif /* not ALLOW_MODULES */
 
C     time increment [day] 
      _RL     dt       
C     array boundaries
      integer n_max, iMin, iMax, jMin, jMax, bi, bj
C     
C     reciprocal of thickness of grid boxes [1/m] 
      _RS recipDz(n_max)
C     recipDzC :: reciprocal of the distance between cell centers 
C	          (i.e., tracer points)
      _RS recipDzC(n_max)
      _RS recip_hFacC(1-Olx:sNx+Olx,1-Oly:sNy+Oly,n_max,nSx,nSy)
      _RS maskC      (1-Olx:sNx+Olx,1-Oly:sNy+Oly,n_max,nSx,nSy)
C      sinking velocity [m/d]
      _RL ws(1-Olx:sNx+Olx,1-Oly:sNy+Oly,n_max) 
C     concentration [mmol m^{-3}]
      _RL c   (1-Olx:sNx+Olx,1-Oly:sNy+Oly,n_max,nSx,nSy)
C---- sink due to sinking [mmolm^{-3}d^{-1}]
      _RL sink(1-Olx:sNx+Olx,1-Oly:sNy+Oly,n_max)
C----               
C     loop counters and field indices
      integer i, j, k, km2, km1, kp1
C     abbreviations and flux limiters
      _RL     Rjp, Rj, Rjm, wLoc, wP, wM, cfl
      _RL     thetaP, thetaM, psiP, psiM, d0, d1
      _RL     oneSixth
      parameter ( oneSixth = 1./6. )
C     vertical flux (positive upwards) through kth vertical cell interface
      _RL     wFlux(1-Olx:sNx+Olx,1-Oly:sNy+Oly) !(1:n_max)
C     auxillary variable      
      _RL     wFluxkp1(1-Olx:sNx+Olx,1-Oly:sNy+Oly) 
#ifdef RECOM_EXPORT_DIAGNOSTICS
C     auxilliary variable for export diagnostics
      _RL     export(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
#endif
#ifdef RECOM_MAREMIP
      _RL export3d(1-Olx:sNx+Olx,1-Oly:sNy+Oly,n_max) 
#endif
#ifdef RECOM_EXPORT_DIAGNOSTICS
      integer kexp
      kexp = 7
#endif

C---------------------------------------------------------------------
C     some initializations
      do k = 1, n_max
       do j = jMin, jMax
        do i = iMin, iMax
         sink(i,j,k)    = 0.
        enddo
       enddo
      enddo
      do j = jMin, jMax
       do i = iMin, iMax
        wFlux   (i,j) = 0.
        wFluxkp1(i,j) = 0.
       end do
      end do

C     run the k-loop from bottom to top
      do k=n_max,2,-1
C     take care of loop boundaries by defining new field indices
       km2=MAX(1,k-2)
       km1=MAX(1,k-1)
       kp1=MIN(n_max,k+1)
       do j = jMin, jMax
        do i = iMin, iMax
C     use the "right" sign: downward velocity wLoc is negative, if sinking
C     (i.e., downward) velocity ws is positive 
         wLoc = -ws(i,j,k)*maskC(i,j,k,bi,bj)*maskC(i,j,km1,bi,bj)
C     some abbreviations
         wP = wLoc+abs(wLoc)
         wM = wLoc-abs(wLoc)
C     more abbreviations (tracer differences)
         Rjp=(c(i,j,k,  bi,bj)-c(i,j,kp1,bi,bj))*maskC(i,j,kp1,bi,bj)
         Rj =(c(i,j,km1,bi,bj)-c(i,j,k,  bi,bj))
     &        *maskC(i,j,k,bi,bj)*maskC(i,j,km1,bi,bj)
         Rjm=(c(i,j,km2,bi,bj)-c(i,j,km1,bi,bj))*maskC(i,j,km1,bi,bj)
C     compute Courant number cfl 
         cfl=abs(wLoc*dt*recipDzC(k))
C     DST3 parameters
         d0=(2.D0-cfl)*(1.D0-cfl)*oneSixth
         d1=(1.D0-cfl*cfl)*oneSixth
C     compute flux limiters psiP, and psiM
         thetaP=Rjm/(1.D-20+Rj)
c$$$       thetaP=0.
c$$$       IF (Rj.NE.0.) thetaP=Rjm/Rj
         psiP=d0+d1*thetaP
         psiP=MAX(0. _d 0, MIN( MIN(1. _d 0,psiP), 
     &        (1. _d 0-cfl)/(1. _d -20+cfl)*thetaP))
         thetaM=Rjp/(1.D-20+Rj)
c$$$       thetaM=0.
c$$$       IF (Rj.NE.0.) thetaM=Rjp/Rj
         psiM=d0+d1*thetaM
         psiM=MAX(0. _d 0, MIN( MIN(1. _d 0,psiM),
     &      (1. _d 0-cfl)/(1. _d -20+cfl)*thetaM))
C     compute flux at interface between cell k and (k-1)
         wflux(i,j)= 
     &        ( 0.5*wP*( c(i,j,k,  bi,bj) + psiM*Rj )
     &        + 0.5*wM*( c(i,j,km1,bi,bj) - psiP*Rj ) )
C     sink due to sinking for layer/cell k
C     minus sign because this has been moved to the right hand side
         sink(i,j,k) = -(wflux(i,j)-wFluxkp1(i,j))
     &        *recipDz(k)*recip_hFacC(i,j,k,bi,bj)*dt
C     store flux at level kp1 for the next cycle
         wfluxkp1(i,j) = wflux(i,j)       !(k)
        end do
       end do
#ifdef RECOM_EXPORT_DIAGNOSTICS
C vertical flux at the upper edge of layer KEXP is written to diagnostic field
       if (k.eq.kexp) then
          do j = jMin, jMax
             do i = iMin, iMax
                export(i,j) = wflux(i,j)
             end do
          end do
       endif
#endif
#ifdef RECOM_MAREMIP
C 
       do j = jMin, jMax
          do i = iMin, iMax
             export3d(i,j,k) = wflux(i,j)
          end do
       end do
#endif
      end do
C     surface flux is zero (no flux from above the surface)
      k=1
      do j = jMin, jMax
       do i = iMin, iMax
        wflux(i,j) = 0.
        sink(i,j,k) = -(wflux(i,j)-wFluxkp1(i,j))
     &       *recipDz(k)*recip_hFacC(i,j,k,bi,bj)*dt
       end do
      end do

      return
      end
