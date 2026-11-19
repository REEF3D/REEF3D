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

#ifndef SFLOW_MEMORY_H_
#define SFLOW_MEMORY_H_

#include"lexer.h"

// Which face states of SFLOW (fdm2D) are stored. Each of these is reconstructed and read only with the options
// named; without them the slice has no array (V = nullptr, slice_base) - 8 bytes per cell each. The checkpoints
// skip a slice without an array (sflow_f_ckpt, sflow_amr_ckpt).

// QUs, QUn, QUe, QUw, QVs, QVn, QVe, QVw: the conserved momentum V of the Boussinesq model (A 220 4) at the faces,
// for the HLL dissipation (sflow_momentum_func::reconstruct, sflow_HLL::aij_U/aij_V)
inline bool sflow_bous_faces(lexer *p)
{
    return p->A220==4;
}

// WHs, WHn, WHe, WHw: the vertical momentum at the faces, convected by the non-hydrostatic models
// (A 220 1, 2, 3 with A 214 1; sflow_flux_build_f::start_W, sflow_curv)
inline bool sflow_wh_faces(lexer *p)
{
    return p->A220>=1 && p->A220<=3;
}

// Vs, Vn (x-faces) and Ue, Uw (y-faces), the tangential velocities of the faces: the HLLC flux (A 225 2), the
// face-normal velocities of the curvilinear grid (sflow_curv) and the advection of the semi-implicit step
// (A 210 11/12/13, sflow_momentum_SI::advection). The HLL flux (A 225 1) reads the normal velocities only.
inline bool sflow_tangential_faces(lexer *p)
{
    return p->A225==2 || p->geo_curv==1 || p->A210>=11 || p->A226==1;
}

// Ws, Wn, We, Ww, the vertical velocity at the faces: the HLLC flux of the vertical momentum (A 225 2 with a
// non-hydrostatic model)
inline bool sflow_w_faces(lexer *p)
{
    return p->A225==2 && sflow_wh_faces(p);
}

// Fs, Fn, Fe, Fw, the physical fluxes of the two states of a face: stored for the HLLC flux (A 225 2), which
// reads them after the HLL flux, and on the curvilinear grid (sflow_curv::flux). The HLL flux of the Cartesian
// grid takes them on the fly (sflow_HLL::flux_HLL).
inline bool sflow_face_fluxes(lexer *p)
{
    return p->A225==2 || p->geo_curv==1;
}

#endif
