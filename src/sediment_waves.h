/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Architect: Hans Bihs
--------------------------------------------------------------------*/

#ifndef SEDIMENT_WAVES_H_
#define SEDIMENT_WAVES_H_

#include"increment.h"
#include"slice4.h"

class lexer;
class fdm2D;
class fdm_nhf;
class ghostcell;
class sediment_fdm;
class seastate_coupling;

using namespace std;

/*--------------------------------------------------------------------
Waves in the sediment transport of SFLOW and NHFLOW (S 110 > 0)

Wave parameters of the cells
  S 110 1  from REEF3D::SEASTATE (A 750 1, seastate_coupling::bed_waves),
           after every wave step: Hs, Tm-1,0, Tp, the rms orbital velocity
           at the bed Urms = sqrt(sum sig^2/sinh^2(kd) E), its mean period
           Tbot and direction (ex,ey), the roller dissipation (A 748 1)
  S 110 2  prescribed uniform waves Hs = S 111, Tp = S 112, direction
           S 113 [deg, ccw from +x], depth-limited Hs <= S 114 h, linear
           theory: Urms = sig Hs/(4 sinh kh), Tm-1,0 = 0.9 Tp
Derived
  Hrms = Hs/sqrt(2), Tw = Tbot (S 110 1) or Tp (S 110 2)
  breaking dissipation per area and density Dbr [m^3/s^3]: the roller
           dissipation (SEASTATE surfbeat with roller), otherwise
           Baldock et al. (1998) with Hb = 0.73 h:
             Dbr = 1/4 g/Tp exp(-(Hb/Hrms)^2) (Hb^2 + Hrms^2)
  S 118 1  stirring by breaking (XBeach, Roelvink et al. 2009):
             kb = (Dbr)^(2/3)/(exp(h/Hrms) - 1),
             Urms_eff^2 = Urms^2 + 1.45 kb
  Uw = sqrt(2) Urms_eff   amplitude of the equivalent monochromatic wave
  Skewness Sk and asymmetry As of the near-bed orbital velocity after
  Ruessink, Ramaekers and van Rijn (2012) from the Ursell number
             Ur = 3/8 Hs k/(kh)^3, k of Tm-1,0
             B   = 0.857/(1 + exp((-0.471 - log10 Ur)/0.297))
             psi = -90 deg + 90 deg tanh(0.815/Ur^0.672)
             Sk = B cos(psi), As = B sin(psi)
  asymmetry transport velocity (XBeach, van Thiel de Vries 2009)
             ua = S 117 (Sk - As) Urms_eff, in the direction of the waves
Wave cycle (orbital): the near-bed orbital velocity u_w(t) along (ex,ey),
  zero mean and rms Urms_eff, S 116 0 sinusoidal, S 116 1 Abreu et al.
  (2010) with r and phi from B and psi (Ruessink et al. 2012):
             b = sqrt(2) B/sqrt(9 + 2 B^2), r = 2b/(1 + b^2), phi = -psi - pi/2
             u ~ [sin(wt) + r sin(phi)/(1 + sqrt(1-r^2))]/[1 - r cos(wt + phi)]
Wave and current bed shear stress, Soulsby (1997, DATA13), after the
current bed shear stress tau_c of S 16 (bedshear::taubed):
  fw = 1.39 (A/z0)^-0.52 (<= 0.3), A = Uw Tw/(2 pi), z0 = ks_eff/30
  tau_w   = 1/2 rho fw Uw^2
  tau_m   = tau_c (1 + 1.2 (tau_w/(tau_c + tau_w))^3.2)
  tau_max = sqrt((tau_m + tau_w |cos phi|)^2 + (tau_w sin phi)^2),
            phi the angle between current and waves
  the current is the Eulerian one (S 119 0): with the radiation stress
  (A 751 1) the host velocity minus the Stokes drift M/h of SEASTATE,
  with the vortex force (A 751 2) the host velocity;
  the bed load formulas S 11 1-5 see tau_m (S 115 1) or tau_max (S 115 2),
  the wave formulas S 11 6-9 (bedload_waves) use the wave fields directly.
NHFLOW (A 10 5) as the current model of phase-averaged waves (SEASTATE
A 750 1 or prescribed S 110 2; no wave generation in NHFLOW, B 90 0, B 98 < 2):
  the current is the velocity of the bottom sigma layer (s->P, s->Q);
  with the radiation stress (A 751 1) the SEASTATE force is uniform over
  the layers and the NHFLOW velocity Lagrangian, the Eulerian current is
  the bottom velocity minus M/h (S 119 0), and tau_c of S 16 (from the
  NHFLOW velocity) is corrected with the rough log law of the bottom
  layer, tau_c + rho (kappa/ln(30 z_b/ks_eff))^2 (|u_E|^2 - |u_L|^2) >= 0
  (eulerian_nhflow);
  the wave formulas S 11 6-9 take the depth-averaged current of the log
  profile through the bottom velocity (depth_mean_nhflow, as Delft3D for
  the depth-averaged formulas in 3D):
    U = u_b (ln(h/z0) - 1)/ln(z_b/z0),  z_b = half the bottom layer,
    z0 = ks_eff/30
v1: Cartesian grid, one grid (no mesh refinement, curvilinear grid or
blocks).
--------------------------------------------------------------------*/

class sediment_waves : public increment
{
public:
    sediment_waves(lexer*);
    virtual ~sediment_waves() = default;
    
    // wave parameters of the cells (every call; SEASTATE spectra only after a new wave step)
    void start(lexer*, seastate_coupling*, ghostcell*, sediment_fdm*);
    
    // S 119 0 with SEASTATE and the radiation stress (A 751 1, host velocity = Lagrangian): the Eulerian current
    // of the sediment transport, s->P and s->Q minus the Stokes drift M/h of the faces (ramped as the host)
    void eulerian(lexer*, fdm2D*, ghostcell*, sediment_fdm*);
    
    // NHFLOW, S 119 0 with SEASTATE and the radiation stress: the Eulerian bottom velocity s->P, s->Q (minus M/h)
    // and the correction dtc of tau_c for the Eulerian current (bedshear)
    void eulerian_nhflow(lexer*, fdm_nhf*, ghostcell*, sediment_fdm*);
    
    // NHFLOW, S 11 6-9: s->P, s->Q from the bottom layer to the depth-averaged current of the log profile
    static void depth_mean_nhflow(lexer*, ghostcell*, sediment_fdm*);
    
    // wave and current bed shear stress, after bedshear::taubed (tau_c in s->tau_eff)
    void bedshear(lexer*, ghostcell*, sediment_fdm*);
    
    // near-bed orbital velocity of the wave cycle of cell (ii,jj) at nph phases, along (ex,ey); zero mean, rms Urms_eff
    void cycle(lexer*, int ii, int jj, int nph, double *uw);
    
    // unit-rms shape of the Abreu et al. (2010) velocity for B, psi (radians) at nph phases
    static void abreu(double B, double psi, int nph, double *u);
    
    // Ruessink et al. (2012): B and psi (radians) of the Ursell number
    static void ruessink(double Ur, double &B, double &psi);
    
    // wave cycle phases (S 11 7, 8)
    static const int nphase = 36;
    
    // SEASTATE / prescribed
    slice4 Hs,Tm10,Tp,Urms,Tbot,ex,ey,Dr;
    // derived
    slice4 Tw,Dbr,kb,Ure,Uw,Sk,As,ua;
    // bed shear stresses [N/m^2] and wave friction factor
    slice4 tauc,tauw,taum,taumax,fw;
    // NHFLOW: correction of the current bed shear stress for the Eulerian current [N/m^2] (0 without)
    slice4 dtc;
    
private:
    void seastate(lexer*, seastate_coupling*);
    void prescribed(lexer*, sediment_fdm*);
    void derived(lexer*, ghostcell*, sediment_fdm*);
    
    int stamp;
};

#endif
