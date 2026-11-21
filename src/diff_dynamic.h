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

#ifndef DIFF_DYNAMIC_H_
#define DIFF_DYNAMIC_H_

class lexer;
class fdm;
class field;
class field4;

/*--------------------------------------------------------------------
Dynamic form of the implicit momentum diffusion (D 24 > 0, idiff2_FS, idiff2_FS_2D):

    (1/rho) div( mu (grad u + grad u^T) ),   mu = rho (nu + nu_t)

instead of div( nu (grad u + grad u^T) ). Across a density jump (a granular or viscous material
under air or water) the kinematic form gets the stress wrong by the density ratio: air driven by the
gravity of a tilted frame dragged a film with the stress per unit mass of the air (Validation/
rheology 06). The cell values of mu are rho (nu + nu_t); the viscosities at the faces between the
velocity nodes are the cell values (normal stresses) or the mean of four cells (shear stresses):
arithmetic (D 24 1) or harmonic (D 24 2). rho at the velocity node: the mean of its two cells.

With D 24 0 the assembly reads a->visc and a->eddyv and the factor 1/rho is 1 (bitwise the
kinematic form).
--------------------------------------------------------------------*/

class diff_dynamic
{
public:
    diff_dynamic(lexer*);
    ~diff_dynamic();

    // mu = rho (nu + nu_t) on all cells (with ghost cells), at the start of each assembly
    void update(lexer*, fdm*);

    // the viscosities the assembly reads: mu and zero (D 24 > 0), or a->visc and a->eddyv
    field& visc(fdm*);
    field& eddyv(fdm*);

    // 1/rho at the velocity node between cell (i,j,k) and the next cell in direction c (0 x, 1 y, 2 z); 1 with D 24 0
    double irho(fdm*, int c, int i, int j, int k) const;

    // the shear-stress viscosity from the four cells around an edge: arithmetic mean (D 24 0, 1), harmonic (D 24 2)
    double corner(double m1, double m2, double m3, double m4) const;

    const int mode;

private:
    field4 *mu;
    field4 *zero;
    int size;
};

#endif
