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

#include<functional>

class lexer;
class fdm;
class field;

// PIC single phase (F 81 1, PICIN): the pressure only in the water (phi >= 0), p = 0 at the free
// surface by the ghost fluid method (Gibou et al. 2002): across a face to a cell in the air the
// pressure of the ghost is p_i (theta - 1)/theta, theta = phi_i/(phi_i - phi_air) the distance
// to the surface in the spacing, at least PIC_GFM_THETA

const double PIC_GFM_THETA = 0.05;

// the single phase with p = 0 at the free surface by the ghost fluid method: PIC single phase (F 81 1)
// or the single-phase level set (F 41 1, cfd_singlephase.h)
bool gfm_single(lexer*);

// the Poisson matrix (poisson_pcorr): rows of the cells in the air p = 0, the faces of the water
// to them with the diagonal 1/theta
void pic_gfm_matrix(lexer*, fdm*);

// theta of the ghost fluid method: the distance from the cell centre in the water (level set fin)
// to the surface, in the spacing to the centre in the air (fout), at least PIC_GFM_THETA
double pic_gfm_theta(double fin, double fout);

// the composite projection of the mesh refinement (cfd_amr, F 41 1): the faces at the patch edges
// are assembled by the composite operator, with the ghost fluid method across them there.
// pic_gfm_matrix leaves the faces for which keep(p,i,j,k,s) is true (s: 0 +i, 1 -i, 2 +j, 3 -j,
// 4 +k, 5 -k); nullptr: none
void pic_gfm_keep(const std::function<bool(lexer*,int,int,int,int)> *keep);

// F 41 2, the closed air pockets (cfd_pocket.h): the pressure of the air at the surface is that of the
// pocket of the cell in the air (0 in the atmosphere). pic_gfm_air_field registers the field of these
// pressures of a grid (cfd_singlephase::pa; lexer nullptr: removes the field); pic_gfm_air the value of
// a cell (0 without a field)
void pic_gfm_air_field(lexer*, field*);
double pic_gfm_air(lexer*, int, int, int);

// the gradient of the pressure across the face from (i,j,k) in the direction (di,dj,dk), spacing dx.
// For the pressure (a->press) the surface has the pressure of the air pocket (F 41 2, pic_gfm_air);
// for the correction (pcorr) 0: the projection keeps it
double pic_gfm_grad(lexer*, fdm*, field&, int, int, int, int, int, int, double);

#endif
