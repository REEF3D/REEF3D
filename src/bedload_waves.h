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

#ifndef BEDLOAD_WAVES_H_
#define BEDLOAD_WAVES_H_

#include"bedload.h"
#include"increment.h"
#include"slice4.h"

/*--------------------------------------------------------------------
Wave-related transport formulas for SFLOW (A 10 2), S 11 6-9: depth-averaged
current U, waves of sediment_waves (S 110 > 0; without waves the formulas
reduce to their current-only forms). Volume transport [m^2/s] as a vector
(magnitude s->qbe, unit direction tx, ty); the Exner step takes the
direction of the transport (sediment_f::transport_direction).
d = S 20, Delta = S 22/W 1 - 1, D* = d (Delta g/nu^2)^(1/3).

S 11 6  Soulsby-van Rijn (Soulsby 1997), total load, transport velocity
        Ut = U + ua (ex,ey) (asymmetry, S 117):
          q = (Asb + Ass) |Ut| [ (|Ut|^2 + 0.018/CD Urms^2)^0.5 - Ucr ]^2.4
          Asb = 0.005 h (d/h)^1.2/(Delta g d)^1.2
          Ass = 0.012 d D*^-0.6/(Delta g d)^1.2
          CD = (0.40/(ln(h/z0) - 1))^2, z0 = 0.006 m
          Ucr = 0.19 d^0.1 log10(4h/d90) (d <= 0.5 mm), 8.5 d^0.6 log10(4h/d90)
          (0.5-2 mm), d90 = 2 d
        bed part Asb, suspended part Ass; equilibrium suspended load per
        area Ls = Ass [...]^2.4 (S 120 2); with the suspended load model
        (S 120 > 0) only the bed part is bed load.
        The slope term (1 - 1.6 tan(beta)) is left to S 80/S 85.
S 11 7  van Rijn (2007a) bed load over the wave cycle (quasi-steady, 36 phases),
        friction velocity of current and orbital velocity u_w(t) (S 116):
          u*(t) = sqrt(fc/8) U + sqrt(fw/2) u_w(t) (ex,ey)
          fc = 0.24 log10(12 h/ks)^-2, fw = exp(-6 + 5.2 (A/ks)^-0.19) <= 0.3,
          ks = S 21 d (grain roughness), A = Uw Tw/(2 pi)
          q(t) = 0.5 d D*^-0.3 |u*| ((rho |u*|^2 - tau_cr)/tau_cr)^1.5  along u*(t)
          tau_cr of taucritbed (Shields S 30, slope reduction S 80)
S 11 8  Bailard (1981) energetics, total load over the wave cycle,
        u(t) = U + u_w(t) (ex,ey):
          qb = cf eb/(Delta g tan(phi)) <|u|^2 u>, eb = 0.1, tan(phi) = 0.63
          qs = cf es/(Delta g ws) <|u|^3 u>,       es = 0.02
          cf = a fc/8 + (1-a) fw/2, a = |U|/(|U| + Uw)
        Ls = cf es/(Delta g ws) <|u|^3>; the slope terms are left to S 80/S 85.
S 11 9  Bijker (1971), bed load in the direction of the current:
          q = b d sqrt(g) |U|/C exp(-0.27 Delta d C^2/(mu (U^2 + 0.5 xi^2 Uw^2)))
          b = 2, C = 18 log10(12 h/ks_eff), mu = min((C/C90)^1.5, 1),
          C90 = 18 log10(12 h/d90), xi = C sqrt(fw/(2g)),
          fw = exp(-5.977 + 5.213 (A/ks_eff)^-0.194) <= 0.3 (Swart 1974)
v1: single fraction (S 51 0), Cartesian grid, one grid.
--------------------------------------------------------------------*/

class bedload_waves : public bedload, public increment
{
public:
    bedload_waves(lexer*);
    virtual ~bedload_waves() = default;
    
    void start(lexer*, ghostcell*, sediment_fdm*) override;
    
    // unit direction of the transport, equilibrium suspended load per area [m] (S 11 6, 8)
    slice4 tx,ty,Ls;
    
private:
    void soulsby_vanrijn(lexer*, sediment_fdm*, double&, double&);
    void vanrijn2007(lexer*, sediment_fdm*, double&, double&);
    void bailard(lexer*, sediment_fdm*, double&, double&);
    void bijker(lexer*, sediment_fdm*, double&, double&);
    
    double h,u,v,d,Delta,g,Dstar;
};

#endif
