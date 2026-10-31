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

#ifndef FNPF_PRINT_SCOUR_H_
#define FNPF_PRINT_SCOUR_H_

#include"increment.h"
#include<vector>

class lexer;
class fdm_fnpf;
class ghostcell;
class scour_gauge;

// FNPF scour gauges P 83 (Sumer-Fredsoe incremental) and P 84 (sediment balance), see scour_gauge.h:
// the bed velocity from Fi at the bed node of the column that holds the gauge (as velcalc_sig, no 3D
// velocity arrays needed); potential flow: the current is the mean of the bed velocity

class fnpf_print_scour : public increment
{
public:
    fnpf_print_scour(lexer*,fdm_fnpf*);
    virtual ~fnpf_print_scour();

    void start(lexer*,fdm_fnpf*,ghostcell*);

    void ckpt_state(std::vector<double>&) const;
    void ckpt_restore(const std::vector<double>&, size_t&);

private:
    void sample(lexer*,fdm_fnpf*,ghostcell*,scour_gauge*);

    scour_gauge *psf = nullptr;
    scour_gauge *psb = nullptr;
    std::vector<double> smp;
};

#endif
