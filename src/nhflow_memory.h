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

#ifndef NHFLOW_MEMORY_H_
#define NHFLOW_MEMORY_H_

#include"lexer.h"

// Which 3D fields of NHFLOW (fdm_nhf) are stored. Each of these is read only by the options named; without
// them the array is not allocated (nullptr) - 8 bytes per cell each, about 0.8 GB per array and 1e8 cells.
// The checkpoints skip a nullptr field (driver_ckpt, nhflow_amr_ckpt).

// test: the output P 23
inline bool nhflow_test_field(lexer *p)
{
    return p->P23==1;
}

// KIN: the turbulence models (A 560); the bed shear stress (bedshear_nhflow) takes k = 0 without one
inline bool nhflow_kin_field(lexer *p)
{
    return p->A560!=0;
}

// CONC: suspended sediment (A 10 5 with S 12 > 0) and the output P 26
inline bool nhflow_conc_field(lexer *p)
{
    return (p->A10==5 && p->S12>0) || p->P26==1;
}

// PORPART and PORSTRUC: the porous structures of VRANS (B 200 1, vrans_nhflow_f)
inline bool nhflow_vrans_fields(lexer *p)
{
    return p->B200==1;
}

// Fext, Gext, Hext: the external momentum sources of VRANS (B 200 > 0, vrans_nhflow_f::u/v/w_source) and of the nets
// (X 320 > 0, net_interface::dlm_nhflow), added to F, G, H by irhs/jrhs/krhs; nullptr without them
inline bool nhflow_ext_fields(lexer *p)
{
    return p->B200>0 || p->X320>0;
}

// PORDEM: the porosity of the REEF3D::DEM particles (E 10 > 0); 1 elsewhere, so VRANS skips the factor
inline bool nhflow_dem_field(lexer *p)
{
    return p->E10>0;
}

// The tangential face states Vs, Vn, Ws, Wn (x faces) and Ue, Uw, We, Ww (y faces): only the HLLC flux
// (A 511 2), the curvilinear grid (face-normal velocities, nhflow_curv) and the semi-implicit momentum
// (A 510 12/13) reconstruct and read them. The HLL flux (A 511 1) reconstructs the normal velocity and the
// conserved variables only (nhflow_momentum_func::reconstruct).
inline bool nhflow_tangential_faces(lexer *p)
{
    return p->A511==2 || p->geo_curv==1 || p->A534==1;
}

// The vertical face states Ub, Vb, Wb, Wt and the face states of the wet-dry faces (wetdry_fluxes: ETAs..ETAw,
// Ds..Dw, dfx, dfy): read only by the HLL/HLLC fluxes of the explicit NHFLOW, set by its kinematic conditions
// (kinematic_fsf, kinematic_bed) and by wetdry_fluxes. The semi-implicit momentum (A 534 1) reads none of them on
// one Cartesian grid; mesh refinement (G 1), multi-block and curvilinear grids keep them.
inline bool nhflow_hll_faces(lexer *p)
{
    return !(p->A534==1 && p->G1==0 && p->geo_mb==0 && p->geo_curv==0);
}

class fdm_nhf;
class ghostcell;

// The semi-implicit NHFLOW after its initialisation (driver::nhflow_driver): the pages of the 3D fields, slices and
// solver arrays that hold zeros only are given back to the system; they read as zeros again (nhflow_memory.cpp).
void nhflow_memory_trim(lexer*, fdm_nhf*, ghostcell*);

#endif
