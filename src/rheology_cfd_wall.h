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


#ifndef RHEOLOGY_CFD_WALL_H_
#define RHEOLOGY_CFD_WALL_H_

#include"increment.h"

class lexer;
class fdm;
class field;

// Side-wall friction of a channel that the grid does not resolve (2D x-z or a quasi-2D slice; W 135 delta_w b):
// two smooth side walls at distance b, Coulomb friction tau_w = tan(delta_w) p' on each (Jop, Forterre & Pouliquen
// 2005, JFM 541: the wall term -mu_w p 2/b of their depth-averaged balance, here per cell). Per unit mass, against
// the in-plane velocity u = (u, w):
//   f = - 2 tan(delta_w) p'/(rho_1 b) u/|u|,
// p' = (1 - lambda) (p - p_gauge) >= 0, only in phase 1 of the level set (smoothed Heaviside over psi).
// Applied after the explicit update of the stage like the VRANS drag (B 268 1), as a Coulomb step with the stage
// weight w: u* <- u* max(0, 1 - w dt F/|u*|), F = 2 tan(delta_w) p'/(rho_1 b); the friction stops the velocity instead
// of reversing it and holds the material where |u*| <= w dt F. |u*| from the updated component and the stage value of
// the other in-plane component. Only the components u and w (the walls are normal to y).

class rheology_cfd_wall : public increment
{
public:
    rheology_cfd_wall(lexer*);
    
    void implicit(lexer*, fdm*, int, double, field&, field&, field&, field&);
    
private:
    double rate(lexer*, fdm*, int, field&, field&, field&, double&);
    
    const double coef;      // 2 tan(delta_w)/(rho_1 b)
    const double lambda;
};

#endif
