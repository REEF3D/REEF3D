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
Author: Ahmet Soydan
--------------------------------------------------------------------*/

#ifndef ACOUSTICS_FWH_KERNEL_H_
#define ACOUSTICS_FWH_KERNEL_H_

#include<vector>

//  Permeable-surface Ffowcs Williams-Hawkings (FW-H) acoustic analogy, time domain.
//
//  Solver independent: no lexer, fdm or ghostcell. The coupling delivers p' and u on the panels
//  of a closed control surface at each source time; the class returns the acoustic pressure at
//  the observers. With MPI each rank holds the panels of its subdomain; the observer signals of
//  the ranks are summed by the coupling (data() starts at the same index on all ranks).
//
//  Formulation: Farassat 1A with the permeable surface of di Francescantonio (1997), surface and
//  observers at rest in a medium at rest, viscous stresses and the quadrupoles outside the
//  surface neglected:
//
//   4 pi p'(x,t) = int [ Qdot/r + Ldot_r/(c0 r) + L_r/r^2 ]_ret dS
//
//   Q   = rho u_n                    ret: source time tau = t - r/c0
//   L_i = p' n_i + rho u_i u_n       n:   unit normal, out of the enclosed region (to the observers)
//   rho = rho0 + p'/c0^2             (pressure-based density: the input is incompressible)
//
//  Time discretisation (advanced time, Casalino 2003): the integrand of a panel is evaluated at
//  the source times, with central differences on the non-uniform source time steps, and spread
//  with the piecewise linear (hat) interpolant to the uniform observer time grid t_k = k dto,
//  shifted by the constant delay r/c0. The sample of time level n-1 is spread at level n, when
//  both neighbours are known, so only three time levels of the panel data are stored and
//  adaptive time steps need no special treatment.
//
//  Mean flow: with set_medium_velocity(U) the medium moves uniformly with U through the surface and
//  the observers at rest (CFD frame of a ship or a propeller in a current; incompressible input).
//  The formulation is applied in the frame of the medium at rest, where surface and observers move
//  with v = -U (Galilean transformation, Farassat 1A for a moving permeable surface):
//
//   4 pi p'(x,t) = int [ Qdot/r + Ldot_r/(c0 r) + L_r/r^2 + Q v_r/r^2 ]_ret dS
//
//   Q   = rho u_n - rho0 U_n,   L_i = p' n_i + rho (u_i - U_i) u_n,   v_r = -U.rhat
//
//  u is the fluid velocity of the CFD frame, so the uniform stream itself gives Q = L = 0. The
//  Doppler factors 1/(1-M_r) and the other terms of order M = |U|/c0 are neglected (M ~ 1e-3 in
//  water); the term Q v_r/r^2 is the part of order 1 (c0 M_r). U = 0 gives the formulation above.
//
//  Images: an observer can carry image points with a weight, e.g. the mirror image across a
//  pressure-release free surface with weight -1 (Lloyd's mirror): evaluating the real sources at
//  the image point equals evaluating the image sources at the observer.

struct fwh_panel
{
    double x[3];    // panel centre
    double n[3];    // unit normal, out of the enclosed region
    double dS;      // panel area
};

class fwh_permeable
{
public:
    fwh_permeable(double c0, double rho0, double dto);

    void add_panel(const fwh_panel&);
    int add_observer(const double *x);                              // returns the observer id
    void add_image(int obs, const double *x, double weight);        // extra point of observer obs
    void set_medium_velocity(const double *U);                      // uniform mean flow, default 0

    // one source time level: p'[np] and u[3*np] on the panels, in the order of add_panel
    void step(double tau, const double *p, const double *u);

    int panels() const {return int(panel.size());}
    int observers() const {return int(obs.size());}
    double dt_obs() const {return dto;}

    // accumulated signal of observer o of this rank: data(o)[m] at t = (first_index()+m)*dto;
    // first_index() is the same on all ranks (from the first source time), the length may differ
    long first_index() const {return k0;}
    const std::vector<double>& data(int o) const {return obs[size_t(o)].sig;}

    // delay bounds r/c0 of observer o over the panels of this rank and all points of the observer
    void delays(int o, double &dmin, double &dmax) const;

    // observer times whose signal is complete, given the delay bounds over all ranks;
    // false if there is none yet
    bool complete_range(double dmin, double dmax, double &t0, double &t1) const;

private:
    struct point
    {
        double x[3];
        double w;
    };

    struct observer
    {
        std::vector<point> pt;
        std::vector<double> sig;
    };

    void spread(observer&, double f, double d, double ta, double tm, double tb);
    void integrand();

    const double c0, rho0, dto;
    double Um[3];

    std::vector<fwh_panel> panel;
    std::vector<observer> obs;

    // p', ux, uy, uz per panel at the last three source times (lev[2] newest)
    std::vector<double> lev[3];
    double tau[3];
    int nlev;

    long k0;
    double tau_first_mid, tau_last_mid;
    bool spread_any;
};

#endif
