C$Header: /csys/software/cvsroot/REcoM/pkg/recom/recom_iterate_ph.F,v 1.3 2008/02/11 09:52:59 mlosch Exp $
#include "RECOM_OPTIONS.h"
#ifdef TARGET_NEC_SX
# define RECOM_VECTORIZE_ITERATION
#endif /* TARGET_NEC_SX */
#ifdef ALLOW_MODULES
      module m_iterate_ph
      contains
#endif /* ALLOW_MODULES */
C-----------------------------------------------------------------------
C 'safe' iterative zero-finding routine for minimizing the mismatch
c between measured alkalinity and alkalinity calculated from pH. 
c
c Based on the routine DRTSAFE from Numerical Recipes. 
c Modified by R.M.Key 4/94, Geir Evensen (NERSC), Markus Schartau (AWI),
c Martin Losch (AWI), Christoph Voelker (AWI)
c-----------------------------------------------------------------------         

C  Newton-Raphson method --- rewritten by Tingting Wang
      subroutine recom_iterate_ph(
     I     bt, dic_molal, talk_molal, k1, k2, kw, kb,
     I     X1, X2, XACC,
     U     xsol, 
     O     NITER,
     I     imin, imax, jmin, jmax )
      implicit none
#include "SIZE.h"
C     routine arguments
C     Input:
C     X1, X2: minimum and maximum values expected for result
C     xguess: first guess of result
C     XACC  : required accuracy for convergence
      _RL X1, X2, XACC
C     i/jmin/max: loop boundaries
      integer imin, imax, jmin, jmax
C     Output:
C     xsol  : solution of iteration
C     NITER : actual number of iterations required to obtain xsol
      _RL xsol(1-olx:snx+olx,1-oly:sny+oly)
      integer NITER
C     former common block variables to be passed to recom_talk_difference
      _RL bt(1-olx:snx+olx,1-oly:sny+oly)
      _RL dic_molal(1-olx:snx+olx,1-oly:sny+oly)
      _RL talk_molal(1-olx:snx+olx,1-oly:sny+oly)
      _RL k1(1-olx:snx+olx,1-oly:sny+oly)
      _RL k2(1-olx:snx+olx,1-oly:sny+oly)
      _RL kw(1-olx:snx+olx,1-oly:sny+oly)
      _RL kb(1-olx:snx+olx,1-oly:sny+oly)

C     local variables
      integer MAXIT, i, j
      _RL DF, DFM, SWAP, F, FM, TEMP, xm, xguess      
      _RL FL(1-olx:snx+olx,1-oly:sny+oly)
      _RL FH(1-olx:snx+olx,1-oly:sny+oly)
      _RL XL(1-olx:snx+olx,1-oly:sny+oly)
      _RL XH(1-olx:snx+olx,1-oly:sny+oly)

#ifdef RECOM_VECTORIZE_ITERATION
      integer globalDoIterate
      integer doIterate(1-olx:snx+olx,1-oly:sny+oly)
#endif /* RECOM_VECTORIZE_ITERATION */

      MAXIT=100
C     first guess/preparation
      DO j = jmin, jmax
       DO i = imin, imax

#ifdef RECOM_VECTORIZE_ITERATION
        doIterate(i,j) = 1
#endif /* RECOM_VECTORIZE_ITERATION */

        CALL recom_talk_difference(
     I       bt(i,j), dic_molal(i,j), talk_molal(i,j), 
     I       k1(i,j), k2(i,j), kw(i,j), kb(i,j),
     &       X1,FL(i,j),DF)
        CALL recom_talk_difference(
     I       bt(i,j), dic_molal(i,j), talk_molal(i,j), 
     I       k1(i,j), k2(i,j), kw(i,j), kb(i,j),
     &       X2,FH(i,j),DF)
        IF(FL(i,j) .LT. 0.0) THEN
         XL(i,j)=X1
         XH(i,j)=X2
        ELSE
         XH(i,j)=X1
         XL(i,j)=X2
         SWAP=FL(i,j)
         FL(i,j)=FH(i,j)
         FH(i,j)=SWAP
        END IF
         
C     use [H+] from last timestep for starting iteration 
        xguess = xsol(i,j)
        CALL recom_talk_difference(
     I       bt(i,j), dic_molal(i,j), talk_molal(i,j), 
     I       k1(i,j), k2(i,j), kw(i,j), kb(i,j),
     &       xguess,F,DF)
        IF (F .LT. 0.0) THEN
         XL(i,j)=xguess
        ELSE
         XH(i,j)=xguess
        END IF

C     enddo i/j-loops
       ENDDO
      ENDDO

#ifdef RECOM_VECTORIZE_ITERATION
C     start iteration
      DO NITER=1,MAXIT

       DO j = jmin, jmax
        DO i = imin, imax
         IF ( doIterate(i,j) .EQ. 1 ) THEN
#else
      DO j = jmin, jmax
       DO i = imin, imax
C     start iteration separately for each (i,j)
        DO NITER=1,MAXIT
#endif /* RECOM_VECTORIZE_ITERATION */

          xm=0.5*(XL(i,j)+XH(i,j))
          CALL recom_talk_difference(
     I         bt(i,j), dic_molal(i,j), talk_molal(i,j), 
     I         k1(i,j), k2(i,j), kw(i,j), kb(i,j),
     &         xm,FM,DFM)
          TEMP=xm-FM/DFM
          IF (TEMP .LT. XH(i,j) .AND. TEMP .GT. XL(i,j)) THEN
           CALL recom_talk_difference(
     I          bt(i,j), dic_molal(i,j), talk_molal(i,j), 
     I          k1(i,j), k2(i,j), kw(i,j), kb(i,j),
     &          TEMP,F,DF)
           xsol(i,j)=TEMP
          ELSE
           xsol(i,j)=xm
           F=FM
          END IF
CML         IF(ABS(F) .LT. XACC)RETURN 

          IF(ABS(F) .LT. XACC) THEN

#ifdef RECOM_VECTORIZE_ITERATION
           doIterate(i,j) = 0
#else
C     we can leave the innermost loop and go to the next (i,j)
           EXIT
#endif /* RECOM_VECTORIZE_ITERATION */

          ELSE
          
           IF(F .LT. 0.0) THEN
            XL(i,j)=xsol(i,j)
            FL(i,j)=F
           ELSE
            XH(i,j)=xsol(i,j)
            FH(i,j)=F
           END IF
           
          ENDIF
          
#ifdef RECOM_VECTORIZE_ITERATION
C     endif doIterate
         ENDIF
C     end of i/j-loops
        ENDDO
       ENDDO
C     decide if we need to continue
       globalDoIterate=0
       DO j = jmin, jmax
        DO i = imin, imax
C     collect the convergence status
         globalDoIterate = globalDoIterate+doIterate(i,j)
        ENDDO
       ENDDO
       IF (globalDoIterate.EQ.0) EXIT
C     enddo iter
      ENDDO
#else
C     enddo iter
        ENDDO
C     end of i/j-loops
       ENDDO
      ENDDO
#endif /* RECOM_VECTORIZE_ITERATION */

      RETURN
#ifdef ALLOW_MODULES
      END subroutine recom_iterate_ph
#else
      END 
#endif /* ALLOW_MODULES */

cccc 
cccc      This subroutine is adapted from the OCMIP program 
cccc      (updated and extended by Christoph Voelker)   
cccc
      subroutine recom_talk_difference(
     I     bt, dic_molal, talk_molal, k1, k2, kw, kb,
     U     x,fn,df)
      implicit none

C     routine arguments
      _RL x, fn, df
C     common blocks variables
      _RL bt, dic_molal, talk_molal
C      common /species/ bt,dic_molal,talk_molal
      _RL k1,k2,kw,kb
C      common /equilibrium_constants/ k1,k2,kw,kb,ff

C     local variables
      _RL X2, X3, B, B2, DB
C     reciprocal values for safer division
      _RL rb, rx
      _RL k12
C
C This routine expresses TA as a function of DIC, htotal and constants.
C It also calculates the derivative of this function with respect to 
C htotal. It is used in the iterative solution for htotal. In the call
C "x" is the input value for htotal, "fn" is the calculated value for TA
C and "df" is the value for dTA/dhtotal
C
      x2=x*x
      x3=x2*x
      k12 = k1*k2
      b = x2 + k1*x + k12
      b2=b*b
      db = 2.0*x + k1
      rb = 0.
      if ( b .ne. 0. ) rb = 1./b
      rx = 0.
      if ( x .ne. 0 ) rx = 1./x
C     
C     fn = hco3+2*co3+borate+oh-hfree-ta
C     
      fn = k1*x*dic_molal*rb +
     &     2.0*dic_molal*k12*rb +
     &     bt*kb/(kb + x) +
     &     kw*rx -
     &     x -
     &     talk_molal
c     print components of alkalinity
c     print *, k1*x*dic_molal/b, 2.0*dic_molal*k12/b, bt/(1.0 + x/kb), kw/x, -x 
C     
C     df = dfn/dx
C     
      df = k1*dic_molal*rb - k1*x*dic_molal*db*(rb*rb) -
     &     2.0*dic_molal*k12*db*(rb*rb) -
     &     bt*kb/(kb+x)**2 -
     &     kw*(rx*rx) -
     &     1.0
      
#ifdef ALLOW_MODULES
      END subroutine recom_talk_difference
#else 
      RETURN
      END 
#endif /* ALLOW_MODULES */

#ifdef ALLOW_MODULES
      end module m_iterate_ph
#endif /* ALLOW_MODULES */
