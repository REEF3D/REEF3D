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

#ifndef RHEOLOGY_CFD_CORE_H_
#define RHEOLOGY_CFD_CORE_H_

#include"rheology.h"
#include"strain.h"
#include"rheology_core.h"

// CFD adapter of the rheology core (W 130 > 0): the kinematic viscosity of the rheological phase
// (phase 1 of the level set) from the local laws of the core,
//   nu_1 = mu_eff(gamma, p')/rho_1,  gamma = sqrt(2 D:D) of the velocities of the fluid update,
// with the effective pressure p' = (1 - lambda) p (W 126 1, W 127 lambda) and p either
//   W 131 2: the solved pressure minus the gauge pressure (default)
//   W 131 3: as 2, at least the hydrostatic pressure of the smoothed interface band above the cell
//   W 131 1: lithostatic, rho_1 (g . n) phi with n = grad(phi)/|grad(phi)|: the distance to the nearest
//            surface times the gravity normal to it; right under a flat surface (also on an inclined frame
//            W 20, W 21), too low next to steep flanks (a column, the front of a slide), where phi is the
//            horizontal distance to the flank and not the overburden.
// No momentum source terms: the stress acts through the (implicit) diffusion of the momentum
// equations with the variable viscosity (D 20 2), in the dynamic form (D 24 1, set with W 130).
// W 130 5: mu(J) / mu(K) of a submerged granular material (Boyer et al. 2011, Trulsson et al. 2012), the
// interstitial fluid phase 2 (W 3 W 4) or, in the three-phase model, the water (W 1 W 2).
// Three-phase model (W 140, rheology_cfd_3phase): the slide is not phase 1 of the level set; mu_eff()
// gives its dynamic viscosity with the effective pressure and the dilatancy of the three-phase model.

class rheology_cfd_core final : public rheology, public strain
{
public:
    rheology_cfd_core(lexer*);
    virtual ~rheology_cfd_core() = default;

    double viscosity(lexer*,fdm*,ghostcell*, field&, field&, field&) override final;

    // dynamic viscosity [Pa s] of the material in cell (i,j,k)
    double mu_eff(lexer*,fdm*, field&, field&, field&);
    
    void u_source(lexer*,fdm*) override final {}
    void v_source(lexer*,fdm*) override final {}
    void w_source(lexer*,fdm*) override final {}
    
    void filltau(lexer*,fdm*,ghostcell*) override final {}

private:
    double pressure(lexer*,fdm*);
    double gravity_normal(lexer*,fdm*);
    double local(lexer*,fdm*,field&,field&,field&);
    bool in_solid(lexer*,fdm*,int,int,int);
    
    rheo::law3D_params q;
    const double rho1;
    double gx,gy,gz;
};

#endif
