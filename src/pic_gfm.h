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

#ifndef PIC_GFM_H_
#define PIC_GFM_H_

class lexer;
class fdm;
class field;

// PIC single phase (F 81 1, PICIN): the pressure only in the water (phi >= 0), p = 0 at the free
// surface by the ghost fluid method (Gibou et al. 2002): across a face to a cell in the air the
// pressure of the ghost is p_i (theta - 1)/theta, theta = phi_i/(phi_i - phi_air) the distance
// to the surface in the spacing, at least PIC_GFM_THETA

const double PIC_GFM_THETA = 0.05;

// the Poisson matrix (poisson_pcorr): rows of the cells in the air p = 0, the faces of the water
// to them with the diagonal 1/theta
void pic_gfm_matrix(lexer*, fdm*);

// the gradient of the pressure across the face from (i,j,k) in the direction (di,dj,dk), spacing dx
double pic_gfm_grad(lexer*, fdm*, field&, int, int, int, int, int, int, double);

#endif
