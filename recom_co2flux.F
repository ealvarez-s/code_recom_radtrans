C$Header: /csys/software/cvsroot/REcoM/pkg/recom/recom_co2flux.F,v 1.9 2008/02/15 07:27:33 mlosch Exp $
C$Name:  $
#include "PTRACERS_OPTIONS.h"
#include "GCHEM_OPTIONS.h"
#include "RECOM_OPTIONS.h"

      subroutine recom_co2flux(pco2a, dic, talk, temp, sali, u10,
     &   dflux, pco2surf, co2star, hplus, imin, imax, jmin, jmax)
c==========================================================================  
cccc  Calculates atmosphere-ocean flux of carbon, given atmospheric CO2 
cccc  partial pressure in dry air, dissolved inorganic carbon, total 
cccc  alkalinity, temperature, salinity, and 10m wind speed
cccc 
cccc      This subroutine is adapted from the OCMIP program 
cccc      (updated and extended by Christoph Voelker)   
cccc
c
c pco2a : atmospheric pCO2 in dry air          [microatm]           (input)
c dic   : dissolved inorganic carbon           [mmol/m^3]           (input)
c talk  : total alkalinity                     [mmol/m^3]           (input)
c temp  : temperature                          [deg. centigrade]    (input)
c sali  : salinity                             [1]                  (input)
c u10   : wind speed at 10 meters height       [m/s]                (input)
c dflux : flux of co2                          [mmol C / m^2 / d]  (output)
c pco2surf: seawater pCO2                      [microatm]          (output)
c co2star: seawater [CO2*]                     [mmol/m^3]          (output)
c hplus : concentration of h+                  [mol/kg]            (in/out)
c nit   : number of iterations needed          [1]                 (output)
c==========================================================================  
#ifdef ALLOW_MODULES
      use m_iterate_ph
      implicit none
      _RL piston_velocity, schmidt_number
#else
      implicit none
      _RL piston_velocity, schmidt_number
      external piston_velocity, schmidt_number
#endif /* ALLOW_MODULES */
 
#include "SIZE.h"
C     input variables
      integer imin, imax, jmin, jmax
      _RL dic  (1-olx:snx+olx,1-oly:sny+oly)
      _RL talk (1-olx:snx+olx,1-oly:sny+oly)
      _RL temp (1-olx:snx+olx,1-oly:sny+oly)
      _RL sali (1-olx:snx+olx,1-oly:sny+oly)
      _RL pco2a(1-olx:snx+olx,1-oly:sny+oly)
      _RL u10  (1-olx:snx+olx,1-oly:sny+oly)
C     output variables
      _RL dflux   (1-olx:snx+olx,1-oly:sny+oly)
      _RL pco2surf(1-olx:snx+olx,1-oly:sny+oly)
      _RL co2star (1-olx:snx+olx,1-oly:sny+oly)
C     updated variables
      _RL hplus(1-olx:snx+olx,1-oly:sny+oly)

      integer nit
C     common block variables
      _RL scl
      _RL bt(1-olx:snx+olx,1-oly:sny+oly)
      _RL  dic_molal(1-olx:snx+olx,1-oly:sny+oly)
      _RL talk_molal(1-olx:snx+olx,1-oly:sny+oly)
      _RL k1(1-olx:snx+olx,1-oly:sny+oly)
      _RL k2(1-olx:snx+olx,1-oly:sny+oly)
      _RL kw(1-olx:snx+olx,1-oly:sny+oly)
      _RL kb(1-olx:snx+olx,1-oly:sny+oly)
      _RL ff(1-olx:snx+olx,1-oly:sny+oly)
      
C     local variables
      integer i,j
      _RL piston_vel, schmidt_no
      _RL dco2star, co2starair, hplus2
      _RL phlo, phhi, phacc, phgu 
      _RL permil, permeg, pco2

      
c     initialising some constants
      permil = 1.e-3 / 1024.5
      permeg = 1.e-6
c-- maximum and minimum values for [H+] that are expected
CML      phlo  = 10**(-5.0)
CML      phhi  = 10**(-9.0)
      phlo  = 1.e-5
      phhi  = 1.e-9
c-- precision for [H+] in iteration
      phacc = 1.e-12


      do j = jmin, jmax
       do i = imin, imax

        call recom_calculate_constants(temp(i,j), sali(i,j), 
     O       k1(i,j), k2(i,j), kw(i,j), kb(i,j), ff(i,j))

C------------------------------------------------------------------------
C Calculate concentrations for borate, sulfate, and fluoride
C
C total borate: Uppstrom (1974)
        scl = sali(i,j) / 1.80655    ! Chlorinity
        bt(i,j)  = 0.000232 * scl/10.811 ! total borate (unit: [mol/kg])
#ifdef CHRISTOPH_FULL_CHEMISTRY
C total sulfate: Morris & Riley (1966)
        st = 0.14 * scl/96.062
C total fluoride: Riley (1965)
        ft = 0.000067 * scl/18.9984
#endif

c     input concentrations DIC and Alkalinity are given in 
c     mmol/m^3 = micromol/l, change them to mol/kg
        talk_molal(i,j) = talk(i,j) * permil
        dic_molal (i,j) = dic(i,j)  * permil
c      print *,'DIC,TA',dic_molal,talk_molal

c------------------------------------------------------------------------
c Iteration to find the pH (or rather the [H+] concentration) for given 
c DIC, TALK, 

       enddo
      enddo
c-- this is the iteration routine

      call recom_iterate_ph(
     I     bt, dic_molal, talk_molal, k1, k2, kw, kb,
     I     phlo, phhi, phacc,
     U     hplus,
     O     nit,
     I     imin, imax, jmin, jmax )

      do j = jmin, jmax
       do i = imin, imax
c---------------------------------------------------------------------
c     Change units from uatm to atm. That is, atm is what the body of 
c     this routine uses for calculations.
c---------------------------------------------------------------------
        pco2 = pco2a(i,j) * permeg

C     Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2, 
C     ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
C
        hplus2     = hplus(i,j) * hplus(i,j)
        co2star(i,j) = dic_molal(i,j) * hplus2 / 
     &       (hplus2 + k1(i,j)*hplus(i,j) + k1(i,j)*k2(i,j))
        co2starair = pco2 * ff(i,j)
        dco2star   = co2starair - co2star(i,j)

c     ---------------------------------------------------------------
c     Add output argument for storing pco2surf
c     Should we be using K0 or ff for the solubility here?
c     ---------------------------------------------------------------

        pco2surf(i,j) = co2star(i,j) / ff(i,j)

C     Convert units of output arguments
c     Note: co2star and dco2star are calculated in mol/kg within this routine 
c     Thus convert now from mol/kg -> mmol/m^3

        co2star(i,j)  = co2star(i,j) / permil
        dco2star      = dco2star / permil

c     Note: pCO2surf and dpCO2 are calculated in atm above. 
c     Thus convert now to uatm

        pco2surf(i,j) = pco2surf(i,j) / permeg

c-- piston velocity at Schmidt number = 660 in [cm/h]
        piston_vel = piston_velocity(temp(i,j),u10(i,j))

c-- actual Schmidt number
        schmidt_no = schmidt_number(temp(i,j))

c-- piston velocity in [m/d]
        piston_vel = piston_vel * sqrt(660.0/schmidt_no)*0.24

c-- finally, calculate CO2 flux
        dflux(i,j) = piston_vel * dco2star ! [mmol C m^{-2} d^{-1}]
        
       enddo
      enddo

      return
      end

      subroutine recom_calculate_constants(temp, sali,
     O     k1, k2, kw, kb, ff)
c-- definition equilibrium constants
      implicit none
C     input variables
      _RL temp, sali
C     common block variables
      _RL k1, k2, kw, kb, ff
C      common /equilibrium_constants/ k1, k2, kw, kb, ff

c-- definition auxiliary variables
      _RL tk, tk100, tk1002, invtk, dlogtk, s2, sqrts, s15

c-- auxiliary variables
      tk = 273.15 + temp
      tk100 = tk/100.0
      tk1002=tk100*tk100
      invtk=1.0/tk
      dlogtk=log(tk)
      s2=sali*sali
      sqrts=sqrt(sali)
CML      s15=sali**1.5
      s15=sali*sqrts

C------------------------------------------------------------------------
C ff = k0*(1-pH2O)*correction term for non-ideality
C
C Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)

      ff = exp(-162.8301 + 218.2968/tk100  +
     &     90.9241*log(tk100) - 1.47696*tk1002 +
     &     sali * (.025695 - .025225*tk100 + 
     &     0.0049867*tk1002))
      
C------------------------------------------------------------------------
C k1 = [H][HCO3]/[H2CO3]
C k2 = [H][CO3]/[HCO3]
C
C Millero p.664 (1995) using Mehrbach et al. data on seawater scale 

      k1=10**(-1*(3670.7*invtk - 62.008 + 9.7944*dlogtk -
     &     0.0118 * sali + 0.000116*s2))
C     
      k2=10**(-1*(1394.7*invtk + 4.777 - 
     &     0.0184*sali + 0.000118*s2))

C------------------------------------------------------------------------
C kb = [H][BO2]/[HBO2]
C
C Millero p.669 (1995) using data from Dickson (1990)

      kb=exp((-8966.90 - 2890.53*sqrts - 77.942*sali +
     &     1.728*s15 - 0.0996*s2)*invtk +
     &     (148.0248 + 137.1942*sqrts + 1.62142*sali) +
     &     (-24.4344 - 25.085*sqrts - 0.2474*sali) *
     &     dlogtk + 0.053105*sqrts*tk)

C------------------------------------------------------------------------
C kw = [H][OH]
C
C Millero p.670 (1995) using composite data

      kw = exp(-13847.26*invtk + 148.9652 - 23.6521 * dlogtk +
     &     (118.67*invtk - 5.977 + 1.0495 * dlogtk) *
     &     sqrts - 0.01615 * sali)

#ifdef CHRISTOPH_FULL_CHEMISTRY

C k1p = [H][H2PO4]/[H3PO4]
C
C DOE(1994) eq 7.2.20 with footnote using data from Millero (1974)
C
      k1p = exp(-4576.752*invtk + 115.525 - 18.453 * dlogtk +
     &     (-106.736*invtk + 0.69171) * sqrts +
     &     (-0.65643*invtk - 0.01844) * s)
      
C k2p = [H][HPO4]/[H2PO4]
C
C DOE(1994) eq 7.2.23 with footnote using data from Millero (1974))
C
      k2p = exp(-8814.715*invtk + 172.0883 - 27.927 * dlogtk +
     &     (-160.340*invtk + 1.3566) * sqrts +
     &     (0.37335*invtk - 0.05778) * s)
      
C------------------------------------------------------------------------
C k3p = [H][PO4]/[HPO4]
C
C DOE(1994) eq 7.2.26 with footnote using data from Millero (1974)
C
      k3p = exp(-3070.75*invtk - 18.141 +
     &     (17.27039*invtk + 2.81197) *
     &     sqrts + (-44.99486*invtk - 0.09984) * s)
      
C------------------------------------------------------------------------
C ksi = [H][SiO(OH)3]/[Si(OH)4]
C
C Millero p.671 (1995) using data from Yao and Millero (1995)
C
      ksi = exp(-8904.2*invtk + 117.385 - 19.334 * dlogtk +
     &     (-458.79*invtk + 3.5913) * sqrtis +
     &     (188.74*invtk - 1.5998) * is +
     &     (-12.1652*invtk + 0.07871) * is2 +
     &     log(1.0-0.001005*s))

C------------------------------------------------------------------------
C ks = [H][SO4]/[HSO4]
C
C Dickson (1990, J. chem. Thermodynamics 22, 113)
C
      ks=exp(-4276.1*invtk + 141.328 - 23.093*dlogtk +
     &     (-13856*invtk + 324.57 - 47.986*dlogtk) * sqrtis +
     &     (35474*invtk - 771.54 + 114.723*dlogtk) * is -
     &     2698*invtk*is**1.5 + 1776*invtk*is2 +
     &     log(1.0 - 0.001005*s))

C------------------------------------------------------------------------
C kf = [H][F]/[HF]
C
C Dickson and Riley (1979) -- change pH scale to total
C
      kf=exp(1590.2*invtk - 12.641 + 1.525*sqrtis +
     &     log(1.0 - 0.001005*s) +
     &     log(1.0 + (0.1400/96.062)*(scl)/ks))
#endif
      
      return
      end

C$Header: /csys/software/cvsroot/REcoM/pkg/recom/recom_co2flux.F,v 1.9 2008/02/15 07:27:33 mlosch Exp $
c#include "CPP_OPTIONS.h"
c#include "PTRACERS_OPTIONS.h"
c#include "GCHEM_OPTIONS.h"
c#include "RECOM_OPTIONS.h"

cccc 
cccc      This subroutine is adapted from the OCMIP program 
cccc      (updated and extended by Christoph Voelker)   
cccc
      _RL function piston_velocity(tc,u)
c--------------------------------------------------------------------------
      _RL tc,u
c chemical_enhancement
        piston_velocity = (2.5*(0.5246 + 1.6256e-2*tc + 
     &                     4.9946e-4*tc*tc) + 0.3*u*u)
c else
c  ifdef cubic_wann
c        u2 = u*u
c        u3 = u2*u
c        wann = (1.09*u - 0.333*u2 + 0.078*u3)
c  else
c	wann = 0.39*u*u
c  endif
c endif
c      print*, 'wann = ',wann
      end 

cccc 
cccc      This subroutine is adapted from the OCMIP program 
cccc      (updated and extended by Christoph Voelker)   
cccc
      _RL function schmidt_number(tc)
c--------------------------------------------------------------------------
c calculates polynomial approximation to Schmidt number
c for seawater (at salinity = 35 psu) after Wanninkhof (1992),
c table A1
c
      _RL tc
      schmidt_number = 2073.1 - 125.62*tc + 3.6276*tc*tc - 
     &                 4.3219e-2*tc*tc*tc
      end 
