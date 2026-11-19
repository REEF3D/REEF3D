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

#ifndef RHEOLOGY_CFD_3PHASE_H_
#define RHEOLOGY_CFD_3PHASE_H_

#include"increment.h"
#include"strain.h"
#include"field4.h"
#include<ostream>
#include<vector>

class lexer;
class fdm;
class ghostcell;
class solver;

/*--------------------------------------------------------------------
Three-phase slide model of REEF3D::CFD (W 140 1): a slide (granular or debris material) with the
rheology of the core (W 130), water and air.

Phases
- water (W 1, W 2) and air (W 3, W 4): the level set phi as in the two-phase model, advected with the
  mixture velocity. Inside the slide phi only marks which fluid fills the space the slide leaves.
- slide: the volume fraction alpha_s (0..1) of the bulk slide material (grains and pore fluid, bulk
  density W 141), transported conservatively with the mixture velocity:
    split in the directions (x y z, the next step z y x), flux form with the face velocities of the
    staggered grid, THINC face values (Xiao et al. 2005, beta 3.5) and the divergence correction of
    Weymouth & Yue (2010): sum(alpha dV) is conserved to round-off and the solver tolerance of the
    projection, alpha stays in [0, 1] for a Courant number up to 1/2 per direction (sub-steps above
    0.4). Closed boundaries (no slide flux through walls, inflow and outflow).

Mixture (fluid update, density of the faces), with alpha smoothed by two passes of the filter (1/4 1/2 1/4)
in each direction (the band of about +-2 cells of the level set, F 45): with the sharp alpha of the transport
the density jump of 1000 between slide and air in one cell drove interface currents of 1 m/s at the front
of a viscous current (0.1 m/s with the level set)
    rho = alpha rho_s + (1 - alpha) (H rho_w + (1 - H) rho_a)
    mu  = alpha mu_s + (1 - alpha) mu_f      (W 149 1; W 149 2: 1/mu = alpha/mu_s + (1 - alpha)/mu_f)
    mu_f = H rho_w nu_w + (1 - H) rho_a nu_a,  mu_s = mu_eff(gamma, p') of the core
with the dynamic form of the momentum diffusion (D 24, set with W 140).

Effective pressure p' = p - p_gage - p_hyd - p_e in the slide
- W 143 0: p_hyd = p_e = 0
- W 143 1: p_hyd = rho_w g (eta_w - z) below the water surface eta_w of the column (the highest water
  cell outside the slide, z + phi), 0 above: the pore pressure of a drained slide
- W 143 2: and the excess pore pressure p_e, advected with the slide and diffusing,
    beta (dp_e/dt + u.grad p_e) = div(k/eta_w grad p_e) - gamma tan(psi)
  (Iverson & Denlinger 2001, Iverson & George 2014) with p_e = 0 outside the slide (drained surface),
  no flux through walls, implicit in time. Initial p_e = lambda_0 (p - p_hyd) (W 144) after the first
  step (the solved pressure of the slide at rest). Dilatancy (W 145 phi_0 > 0): the solid fraction phi
  advected with the slide, tan(psi) = K3 (phi - phi_c/(1 + J^1/2)) as the source of p_e and as extra
  friction of mu(I) / mu(J), and dphi/dt = -phi D_v with the dilation D_v = div(k/eta_w grad p_e), the
  inflow of pore fluid (the volume change the drainage allows).

Interface band (W 131 3, as in the two-phase model): the smoothed band (alpha_s > 1e-3 after smoothing) gets the full
viscosity of the slide, max(mu_s, mu_f), instead of alpha mu_s + (1 - alpha) mu_f, and p' at least the buoyant weight
of the band above the cell, g sum (rho_s - rho_f) alpha dz over the band cells (alpha_s < 1) of the column. With the share alpha mu_s
the light band cells on the inclined surface of a deposit at rest were weaker than their weight: they flowed down the
surface, drove a jet of air along the bed at the toe (0.5 - 0.9 m/s) and carried the thin tail with them.
--------------------------------------------------------------------*/

class rheology_cfd_3phase final : public strain
{
public:
    rheology_cfd_3phase(lexer*, fdm*, ghostcell*);
    ~rheology_cfd_3phase();

    // the instance of the run (W 140 > 0), else nullptr
    static rheology_cfd_3phase *self;

    // initial slide fraction (boxes W 147, triangles W 148)
    void ini(lexer*, fdm*, ghostcell*);

    // one time step, before the level set: transport of alpha_s (p_e, phi), hydrostatic and excess pore pressure
    void start(lexer*, fdm*, ghostcell*, solver*);

    // mixture density of the face between (i,j,k) and (i+aa,j+bb,k+cc) for the Heaviside H of the face
    double rho_face(lexer*, int aa, int bb, int cc, double H);

    // effective pressure p' and the dilatancy tan(psi) in cell (i,j,k)
    double peff(lexer*, fdm*);
    double tanpsi() { return dil ? tpsi(i,j,k) : 0.0; }

    // VTU output (printer_CFD): alpha_s, p_eff, p_hyd (W 143 >= 1), p_e (W 143 2), phi_s (dilatancy)
    void offset_vtu(lexer*, int*, int&);
    void name_vtu(lexer*, std::ostream&, int*, int&);
    void name_pvtu(lexer*, std::ostream&);
    void print_vtu(lexer*, fdm*, ghostcell*, std::vector<char>&, size_t&);

    field4 alpha;
    field4 alphas;   // alpha smoothed over about two cells: density and viscosity of the mixture
    field4 pband;    // W 131 3: the weight of the smoothed band above the cell, the floor of p'
    double rho_s, rho_w, rho_a, eta_w;
    double volume;

private:
    void sweep(lexer*, fdm*, ghostcell*, int dir, double dt);
    double thinc(double aback, double aup, double aahead, double c) const;
    bool face_open(lexer*, int dir, int ii, int jj, int kk) const;
    double face_vel(lexer*, fdm*, int dir, int ii, int jj, int kk) const;
    void hydrostatic(lexer*, fdm*, ghostcell*);
    void excess(lexer*, fdm*, ghostcell*, solver*);
    double permeability(double phi) const;
    void diagnostics(lexer*, fdm*, ghostcell*, double clipped, int nsub);
    void bc(lexer*, ghostcell*, field&);
    void smooth(lexer*, ghostcell*);
    void band(lexer*, fdm*, ghostcell*);
    int nvtu(lexer*) const;

    field4 alpha0;   // alpha at the start of the step (the divergence correction)
    field4 flux;     // flux of alpha through the face i+1/2 (j+1/2, k+1/2) of the sweep [m]
    field4 phyd;     // hydrostatic pore pressure
    field4 pe;       // excess pore pressure (W 143 2)
    field4 phis;     // solid volume fraction (dilatancy)
    field4 tpsi;     // tan(psi) of the step
    field4 qpe, qphi;// alpha p_e, alpha phi during the transport
    field4 fpe, fphi;// their fluxes

    const int pore;  // W 143
    const bool dil;  // W 145 phi_0 > 0
    const double beta_thinc;
    double lambda0, kperm, beta, K3, phic, phi0, dgrain;
    bool pe_ini;
    int sweeps;
    double volume0;
    double clipsum;
    std::vector<double> eta;   // water surface of the columns (global i, j)
};

#endif
