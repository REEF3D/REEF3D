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

#ifndef SOLID_APERTURE_H_
#define SOLID_APERTURE_H_

#include<utility>

class lexer;
class fdm;
class ghostcell;
class field;

// D 1 2, cut cells: the static solids (S, solid level set) and the bed (topoforcing) of
// REEF3D::CFD as impermeable walls of the pressure projection, instead of the direct forcing
// (D 1 1, the default; the moving bodies keep the direct forcing with both).
//
// With the direct forcing the solid is a part of the fluid whose velocity is forced to 0 after the
// predictor; the projection solves through it and moves fluid across it again (a medium of the
// permeability dt/rho, 830 times more open for air than for water), and the solid holds a free
// surface and a hydrostatic column of its own.  D 32 1/2 close the faces inside the solid; the wall
// then lies on the cell faces, half a cell from the surface of the solid at most.
//
// Here the wall is where the surface of the solid cuts the cell faces.  a->apx, apy, apz hold the
// open fraction (aperture) A_f of the faces at the positions of u, v, w:
//
//   continuity          sum over the faces of  A_f u_f S_f = 0            (pjm_corr::rhs)
//   pressure correction coefficient of a face  A_f /(rho_f dx_P dx_N)     (poisson_pcorr)
//   velocity correction on the faces with A_f > 0, unchanged (the aperture is in the divergence,
//                       not in the gradient)
//   A_f = 0             a closed face: Neumann, the velocity 0, a wall for the diffusion and the
//                       wall law (the flags DF, ghostcell::solid_forcing_flag_update)
//   all A_f = 0         a cell of the solid: no row; its level set from the fluid (X 48), its
//                       pressure from the fluid (solid_press_extend)
//
// The momentum equation stays on the full cells (no small-cell limit of the time step), and there
// is no forcing at these solids (gc_solid_forcing.cpp): a face is closed or a face of the fluid.
//
// A_f: the signed distance to the solid at the corners of the face, interpolated from the cell
// centres, and the part of the face where it is positive, cut linearly (a line in 2D; in 3D the
// four triangles to the centre of the face).  The solids: the exact distance to the triangles of
// the grid file (the solid level set of the grid is the distance along the grid lines, and not
// below half a cell inside); the bed: its level set (a moving bed is followed: calculated again
// when the level sets changed).  Below solid_aperture_min the face is closed.  A solid of one cell
// across a direction (its two neighbours in the fluid) cannot be seen at the corners of its faces:
// these faces are closed, as D 32 2.

const double solid_aperture_min = 0.01;

void solid_aperture_update(lexer*, fdm*, ghostcell*);

// the pressure of the cells with all faces closed, from the fluid (forces, probes, print)
void solid_press_extend(lexer*, fdm*, ghostcell*, field&, int);

// the part of a line with the level set a, b at its ends where the level set is positive
inline double solid_cut_line(double a, double b)
{
    if(a>0.0 && b>0.0)
    return 1.0;
    if(a<=0.0 && b<=0.0)
    return 0.0;
    return (a>0.0) ? a/(a-b) : b/(b-a);
}

// the part of a triangle with the level set a, b, c at its corners where the level set is positive
inline double solid_cut_triangle(double a, double b, double c)
{
    const int np = (a>0.0 ? 1 : 0) + (b>0.0 ? 1 : 0) + (c>0.0 ? 1 : 0);

    if(np==3)
    return 1.0;
    if(np==0)
    return 0.0;

    if(np==1)
    {
        // the positive corner first
        if(b>0.0) std::swap(a,b);
        if(c>0.0) std::swap(a,c);
        return (a/(a-b))*(a/(a-c));
    }

    // the corner not positive first
    if(b<=0.0) std::swap(a,b);
    if(c<=0.0) std::swap(a,c);
    return 1.0 - (a/(a-b))*(a/(a-c));
}

// the part of a rectangle with the level set at its corners (in order around it) where the level
// set is positive: the four triangles to the centre
inline double solid_cut_quad(double a, double b, double c, double d)
{
    if(a>0.0 && b>0.0 && c>0.0 && d>0.0)
    return 1.0;
    if(a<=0.0 && b<=0.0 && c<=0.0 && d<=0.0)
    return 0.0;

    const double m = 0.25*(a+b+c+d);

    return 0.25*(solid_cut_triangle(a,b,m) + solid_cut_triangle(b,c,m) + solid_cut_triangle(c,d,m) + solid_cut_triangle(d,a,m));
}


#endif
