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

#ifndef NHFLOW_PRINT_SCOUR_H_
#define NHFLOW_PRINT_SCOUR_H_

#include"increment.h"
#include<vector>

class lexer;
class fdm_nhf;
class ghostcell;
class scour_gauge;

// NHFLOW scour gauges P 83 (Sumer-Fredsoe incremental) and P 84 (sediment balance), see scour_gauge.h:
// the bed velocity is the velocity at the centre of the lowest cell, the current the velocity at D/2
// above the bed (at most the depth, at least the lowest cell centre), both interpolated at the gauge

class nhflow_print_scour : public increment
{
public:
    nhflow_print_scour(lexer*,fdm_nhf*);
    virtual ~nhflow_print_scour();

    void start(lexer*,fdm_nhf*,ghostcell*);

    void ckpt_state(std::vector<double>&) const;
    void ckpt_restore(const std::vector<double>&, size_t&);

private:
    void sample(lexer*,fdm_nhf*,ghostcell*,scour_gauge*);

    scour_gauge *psf = nullptr;
    scour_gauge *psb = nullptr;
    std::vector<double> smp;
};

#endif
