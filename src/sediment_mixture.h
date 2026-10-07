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
Author: Hans Bihs
--------------------------------------------------------------------*/

#ifndef SEDIMENT_MIXTURE_H_
#define SEDIMENT_MIXTURE_H_

#include"increment.h"
#include"slice4.h"
#include<fstream>

class lexer;
class ghostcell;
class sediment_fdm;
class bedload;
class slice;

using namespace std;

/*--------------------------------------------------------------------
Multi-fraction (graded) sediment bed, Hirano active-layer model.

The bed is split into an active (mixing) layer of constant thickness La
on top of one well-mixed substrate layer of thickness Hs. Each layer
stores the volume fractions F_k (active) and Fs_k (substrate) of the
nf grain size classes d_k, with sum_k F_k = sum_k Fs_k = 1.

Bedload per fraction
    qb_k = F_k * qb( d_k, theta_k, theta_cr,k )
    theta_k    = tau_b / ((rho_s - rho) g d_k)
    theta_cr,k = S30 * reduce * xi_k
with the hiding/exposure factor of Wu, Wang & Jia (2000)
    p_hk = sum_j F_j d_j/(d_k + d_j),   p_ek = 1 - p_hk
    xi_k = (p_ek/p_hk)^(-m),            m = S55 (0.6)
The existing bedload formulas (S11) are evaluated once per fraction
with d_k. Bedload direction, relaxation zones and the non-equilibrium
closure act on the fractions exactly as on the single-fraction flux.

Exner per fraction (same discretisation as the single-fraction solver)
    dz_k/dt = -1/(1-n) div(qb_k)           dz/dt = sum_k dz_k/dt
Active-layer sorting (Hirano 1971), discrete form per sediment step
    La F_k^new = La F_k + dz_k - dz FI_k
    FI_k = F_k   (aggradation, dz>0: material passed to the substrate)
    FI_k = Fs_k  (degradation, dz<0: material taken from the substrate)
Substrate: Hs += dz, and on aggradation Fs_k is mixed with F_k.
Per-fraction volume is conserved exactly by this update. If a fraction
is eroded beyond its content, its erosion is limited locally and the
bed level is corrected; the limited volume is reported in the mixture
log (REEF3D_Log/REEF3D_sediment_mixture.dat).

Sandslide: material moved by the slide algorithms carries the
composition of the eroding cell (active layer, and substrate for the
part of the loss exceeding La), so the slide is also conservative
per fraction.

Control
  S 51 d_k Fa_k Fs_k   one line per fraction: diameter [m], initial
                       active-layer fraction, initial substrate fraction
  S 52 La              active layer thickness [m], <=0: 2*d90 (initial)
  S 53 Hs              initial substrate thickness [m]
  S 54 0/1             hiding/exposure off / Wu, Wang & Jia (2000)
  S 55 m               hiding/exposure exponent
  S 56 0/1/2/3         bed roughness ks = S21*d with d = S20/d50/dm/d90
--------------------------------------------------------------------*/

class sediment_mixture : public increment
{
public:
    sediment_mixture(lexer*);
    virtual ~sediment_mixture();

    void ini(lexer*, ghostcell*, sediment_fdm*);

    // bedload per fraction, s->qbe returns the sum over all fractions
    void bedload_fractions(lexer*, ghostcell*, sediment_fdm*, bedload*);

    // Exner helpers
    void exner_begin(lexer*, ghostcell*, sediment_fdm*);
    void load_fraction(lexer*, ghostcell*, sediment_fdm*, int);
    void exner_end(lexer*, ghostcell*, sediment_fdm*);
    void bedchange(lexer*, ghostcell*, sediment_fdm*, slice4**);

    // sandslide hooks
    void slide_zero(lexer*, ghostcell*);
    void slide_transfer(int,int,int,int,double);
    void slide_pde(lexer*, sediment_fdm*, slice&, int, int, double*);
    void slide_finish(lexer*, ghostcell*, sediment_fdm*);

    // grain size statistics and roughness
    void grain_stats(lexer*, ghostcell*);
    double ks_diameter(lexer*, int, int);

    // output
    void print_log(lexer*, ghostcell*);
    int nfields();

    int nf;
    double *d;
    double La;

    slice4 **F, **Fs;
    slice4 **qbe_k, **qb_k, **qbn_k, **vz_k, **dh_k, **fh_k;
    int *noneq_ini_k;

    slice4 Hs, d50, d90, dm;
    slice4 qbe_raw, qbe_tot, fac;

private:
    double hiding(int,int,int);
    void save_base(lexer*, sediment_fdm*);
    void restore_base(lexer*, sediment_fdm*);
    void composition_eroded(int,int,double,double*);

    slice4 shields_eff0, shields_crit0, tau_crit0, shearvel_crit0;
    slice4 **out;
    int *order;
    double *Fe, *V, *FI, *dhc;
    int *inS;
    double *clip_k, *clip_sum;
    double rho_s, S20;
    int hide_type;
    double hide_m;
    ofstream mixlog;
};

#endif
