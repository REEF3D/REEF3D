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

#ifndef SEDIMENT_SUSP2D_H_
#define SEDIMENT_SUSP2D_H_

#include"increment.h"
#include"slice1.h"
#include"slice2.h"
#include"slice4.h"

class lexer;
class fdm2D;
class ghostcell;
class sediment_fdm;
class sediment_inflow;
class sediment_waves;
class bedload_waves;

using namespace std;

/*--------------------------------------------------------------------
Depth-averaged suspended load for SFLOW (S 120 > 0)

  d(hC)/dt + div(q C) = div(h K grad C) + E - D,   E - D = ws F (C_eq - C)

C depth-averaged volume concentration, M = hC the transported variable.
q    the mass flux of the SFLOW continuity: the RK-weighted HLL face
     fluxes of the time step (sflow_eta::update, fdm2D::FEacc; other
     momentum schemes: face velocity x face depth), plus the Stokes
     transport M with the vortex force (A 751 2, ramped as the host).
     The residual of the continuity, (h^n+1 - h^n) + dt div(q) (eta
     boundary conditions, relaxation zones, wetting and drying), is put
     back with the cell concentration, so that a uniform C stays uniform
     with S 119 1.
     The sediment moves with the Eulerian flux q - M (S 119 0, default;
     the Stokes drift is near the surface, the sediment near the bed, as
     XBeach and Delft3D) or with the Lagrangian flux q (S 119 1); with
     waves (S 110) also with h ua (ex,ey), the asymmetry velocity of
     sediment_waves (S 117, XBeach). The sediment volume is conserved in
     both cases.
     Advection: first-order upwind (S 125 1) or MUSCL with the minmod
     limiter (S 125 2), explicit with the flow time step.
K    = S 122 + nu_t/S 123, explicit, sub-steps when dt K (1/dx^2+1/dy^2) > 0.45
E-D  implicit: M^n+1 = (M* + dt ws F C_eq)/(1 + dt ws F/h): unconditionally
     stable, C between C* and C_eq; ws (1-C)^4.7 with S 127 1 (Richardson-Zaki)
C_eq S 120 1: van Rijn (2007b) reference concentration at a = max(ks_eff/2, 0.01 h)
              c_a = 0.015 d D*^-0.3 T^1.5/a,  T = (tau'_cw - tau_cr)/tau_cr
              tau'_cw = 1/8 rho fc U^2 + mu_w 1/4 rho fw Uw^2 (waves: sediment_waves),
              mu_w = 0.7/D* (0.14 .. 0.35), the wave efficiency factor
              fc = 0.24 log10(12h/ks)^-2, fw = exp(-6 + 5.2 (A/ks)^-0.19), ks = S 21 d
              C_eq = c_a/F
     S 120 2: C_eq = Ls/h, the suspended load of the transport formula (S 11 6, 8)
     limited to S 126
F    = c_a/C of the Rouse profile, Z = ws/(kappa u*), u* = sqrt(tau_b/rho)
     with the bed shear stress of the bed (with waves tau_max), 1 <= F <= 200;
     S 121 > 0: F = S 121
Inflow: C = C_eq of the first cell (S 124 < 0) or S 124 at the faces of the
domain edges where q enters; outflow upwind; walls and solids closed.
Cells that fall dry (h < S 128) deposit their sediment on an erodible bed,
cells without one keep it until they are wet again.
Bed: the exchange of a step, M^n+1 - M*, is summed over the flow steps
(ED [m]) and handed to the Exner equation as the mean rate ED/T over the
flow time T since the last bed update (sediment_exner::susp_ED), so the
bed takes S 35 dtsed/T of the exchange, as the morphological factor acts
on the bed load. Bed and water column conserve the sediment exactly when
the sediment time step is the flow time since the last bed update
(S 15 4, or S 15 3 with S 42 1 and S 44 1) with S 17 1 and S 35 1, outside
the bed relaxation zones (S 73, S 75); a note is printed otherwise.
Output: ST_cb = C, ST_cbe = C_eq, ST_qbs = |q C| (cells).
v1: Cartesian grid, one grid, one suspended class (with S 51 > 0 the
exchange is distributed with the active-layer composition).
--------------------------------------------------------------------*/

class sediment_susp2D : public increment
{
public:
    sediment_susp2D(lexer*);
    virtual ~sediment_susp2D();
    
    // allocate the mass flux accumulator of the continuity
    void ini(lexer*, fdm2D*, ghostcell*, sediment_fdm*);
    
    // one flow step: C_eq (bed shear, waves and transport formula already evaluated), then transport and exchange
    void start(lexer*, fdm2D*, ghostcell*, sediment_fdm*, bedload_waves*);
    
    // the mean exchange rate E - D since the last bed update into s->edrate; reset after the bed update
    void bed_rate(lexer*, ghostcell*, sediment_fdm*);
    void bed_done(lexer*);
    
    double total(lexer*, ghostcell*);   // suspended volume in the domain [m^3]
    
    slice4 M,C,Ceq,ED,hn,rate,beta;
    
private:
    void equilibrium(lexer*, fdm2D*, ghostcell*, sediment_fdm*, bedload_waves*);
    void advection(lexer*, fdm2D*, ghostcell*, sediment_fdm*);
    void diffusion(lexer*, fdm2D*, ghostcell*);
    void exchange(lexer*, fdm2D*, ghostcell*, sediment_fdm*);
    double rouse(double h, double a, double Z) const;
    double cin(lexer*, int, int);
    
    // the inflow concentration (S 124 or the time series of S 129 2)
    const sediment_inflow *pin = nullptr;
    static double minmod(double, double);
    
    slice1 Fx,Gx;
    slice2 Fy,Gy;
    
    double tacc[2];             // flow time since the last bed update, initialised flag (checkpoints)
    double clip;                // volume added by clipping negative M (log)
    unsigned long long ckid;
};

#endif
