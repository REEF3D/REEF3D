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

#ifndef CAVITATION_H_
#define CAVITATION_H_

#include"increment.h"
#include"field4.h"
#include<fstream>
#include<map>
#include<vector>

class lexer;
class fdm;
class ghostcell;
class convection;

//  Cavitation for REEF3D::CFD (W 60 1): homogeneous mixture of liquid and vapour with the vapour
//  volume fraction alpha, Schnerr & Sauer (2001) mass transfer.
//
//  Mixture: rho = (1-alpha) rho_l + alpha rho_v, mu = (1-alpha) rho_l nu_l + alpha rho_v nu_v,
//  rho_l = W 1, nu_l = W 2, rho_v = W 62, nu_v = W 63.
//
//  Mass transfer from the Rayleigh growth rate dR/dt = sign(p_v - p) sqrt(2/3 |p_v - p|/rho_l) of
//  bubbles of radius R_B = (3 alpha/(4 pi n0 (1-alpha)))^(1/3), n0 = W 64 bubbles per volume of liquid:
//
//      S = alpha (1-alpha) 3/R_B dR/dt          (vapour production, volume fraction per time)
//      mdot = rho_l rho_v/rho S                 (mass per volume and time)
//
//  For evaporation alpha is increased by the nuclei, alpha_nuc = n0 pi d^3/6/(1 + n0 pi d^3/6),
//  d = W 65, so that pure liquid can cavitate; W 66 C_v C_c scale evaporation and condensation.
//  With the mixture continuity and the vapour equation
//
//      div u = mdot (1/rho_v - 1/rho_l) = S (rho_l - rho_v)/rho
//      d alpha/dt + div(alpha u) - alpha div u = S
//
//  (the vapour equation d alpha/dt + div(alpha u) = mdot/rho_v = S rho_l/rho minus alpha times the
//  continuity). The term alpha div u uses the divergence of the convection scheme itself (the
//  scheme applied to alpha = 1), so a uniform alpha stays uniform whatever the discrete div u; the
//  conservative form would multiply any mismatch between the divergence of the scheme and that of
//  the projection by rho_l/rho, up to 6e4 in vapour.
//
//  Linearised in the pressure, S = C (p_v - p) with C = alpha(1-alpha) 3/R_B sqrt(2/(3 rho_l))/
//  sqrt(|p - p_v| + eps), from alpha and p at the start of the time step (coefficients); for
//  condensation C <= alpha/(dt (p - p_v + eps)), at most the vapour there is in one step (else
//  traces of vapour behind a cavity pull the pressure to p_v every other step). In each
//  projection of the time step the divergence source is implicit in the new pressure:
//  (C (rho_l - rho_v)/rho)/(a dt) is added to the diagonal of the Poisson matrix and
//  (C (rho_l - rho_v)/rho)(p_v - p)/(a dt) to its right-hand side (projection), so that the
//  pressure in a cavity is held near p_v instead of the vapour growing without limit. After the
//  momentum step alpha is transported (TVD, van Leer, RK3) with the source from the final pressure,
//  the source whose divergence the last projection imposed, and bounded to [0, 1] (start).
//
//  The pressure level must be fixed by a Dirichlet boundary: an outlet patch with a pressure
//  (B 441 on a C 1x 21 wall, B 412 id p). p_v (W 61) is in the same scale: for the cavitation
//  number sigma = (p_inf - p_v)/(rho U^2/2) with p_inf = 0 at the outlet, W 61 = -sigma rho U^2/2. Single phase only (F 30 0, F 80 0).
//  Cavitation number control W 69 sigma U x y z tau: from t = tau on p_v = p_ref - sigma rho_l U^2/2,
//  p_ref the pressure at the point (x, y, z) (upstream, as the reference of an experiment) averaged
//  with the time constant tau; before p_v = W 61, a good estimate. A slow correction of the level
//  (outlet patch, outflow B 77 10, the drag of the body): a fast one couples with the cavity (its
//  growth pushes the water out and raises the upstream pressure). p_v is in the output file.
//  Initial vapour (alpha = 1): the cells with the centre in the box W 67 xs xe ys ye zs ze or the
//  sphere W 68 x y z r.
//  Output: alpha_v in the VTU files (the mixture density with P 24 1), REEF3D_CFD_Cavitation/REEF3D-CFD-Cavitation.dat
//  (time, vapour volume, equivalent radius (3V/4pi)^(1/3), vapour production rate, minimum pressure).
//  Restart: the checkpoints (P 48, I 45) hold alpha and the other fields of this class (every field
//  of the run, ckpt_registry), the restart continues bitwise (Regression/validation/cavitation 01);
//  the state files (P 40, I 40) do not hold alpha.
//  Mesh refinement (A 10 6, G 1 > 0, G 7 0): one object per grid (of(lexer)), alpha in the registry of
//  the refinement (filled around the patches, restricted, new patches from their parent), the density
//  of the patch edges from it (cfd_amr::rmix); the phase-change rows of every grid enter the composite
//  projection through its pjm_corr. After the step of level 0 the patches level by level
//  (cfd_amr::cav_step), then the restriction; alpha is not flux-matched at the patch edges.

class cavitation : virtual public increment
{
public:
    // primary: level 0 (checks, messages, the output file); a patch of the mesh refinement: false
    cavitation(lexer*, fdm*, ghostcell*, bool primary=true);
    virtual ~cavitation();

    // the cavitation model of the grid of a lexer (level 0, or a patch of the mesh refinement),
    // used by its density, its projection and the printer; nullptr without one
    static cavitation* of(const lexer *p)
    {
        auto it = grids.find(p);
        return it==grids.end() ? nullptr : it->second;
    }

    void coefficients(lexer*, fdm*, ghostcell*);    // before the momentum step
    void projection(lexer*, fdm*, double);          // after the Poisson matrix is assembled
    void start(lexer*, fdm*, ghostcell*);           // after the momentum step
    void properties(lexer*, fdm*, ghostcell*);      // density and viscosity of the cells from alpha

    double rho_mix(lexer*, double) const;            // mixture density of a vapour fraction
    double vapour_pressure() const { return pv; }
    void vapour_pressure(double v) { pv = v; }        // a patch of the mesh refinement: p_v of level 0 (W 69)
    double val(int ii, int jj, int kk) { return alpha(ii,jj,kk); }

    void name_ParaView_parallel(lexer*, std::ofstream&);
    void name_ParaView(lexer*, std::ostream&, int*, int&);
    void offset_ParaView(lexer*, int*, int&);
    void print_3D(lexer*, fdm*, ghostcell*, std::vector<char>&, size_t&);

    field4 alpha;

private:
    inline static std::map<const lexer*, cavitation*> grids;
    const lexer *owner;
    bool primary;

    void transport_step(lexer*, fdm*, ghostcell*, field&, field&, double, double, double);
    void log(lexer*, fdm*, ghostcell*);
    void sigma_control(lexer*, fdm*, ghostcell*);   // W 69: p_v from the reference pressure

    field4 Cs;          // dS/d(p_v - p), the linearised vapour production
    field4 Sv;          // vapour source S of the alpha equation, from the final pressure
    field4 Dv;          // div u of the convection scheme (minus the scheme applied to alpha = 1)
    convection *pconvec;
    std::ofstream out;

    double rhol, rhov, nul, nuv, pv, n0, alnuc, Cv, Cc, eps;
    double Smax;
    double pref_avg;    // W 69: averaged reference pressure
    bool pref_set;
    int gcval_alpha;
};

#endif
