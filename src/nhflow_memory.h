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

// FSW, DWDT and the stored face fluxes Fs, Fn, Fe, Fw: the HLLYL flux and the Yamazaki pressure (A 520 3).
// The HLL and HLLC fluxes of nhflow_HLL evaluate the face fluxes on the fly (nhflow_flux_face.h).
inline bool nhflow_yl_fields(lexer *p)
{
    return p->A520==3;
}

// test: the output P 23 and the work array of the Yamazaki pressure (A 520 3)
inline bool nhflow_test_field(lexer *p)
{
    return p->P23==1 || p->A520==3;
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

// PORDEM: the porosity of the REEF3D::DEM particles (E 10 > 0); 1 elsewhere, so VRANS skips the factor
inline bool nhflow_dem_field(lexer *p)
{
    return p->E10>0;
}

#endif
