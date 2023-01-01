/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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
Author: Weizhi Wang
--------------------------------------------------------------------*/

#include "wave_lib_spectrum.h"
#include "lexer.h"
#include "ghostcell.h"

double wave_lib_spectrum::Goda_JONSWAP(lexer* p, double w) 
{
    beta_J = 0.06238 / (0.230 + 0.0336 * p->B88 - 0.185 * pow((1.9 + p->B88), -1.0)) * (1.094 - 0.01915 * log(p->B88));

    if(w <= p->wwp)
        sigma = 0.07;

    if(w > p->wwp)
        sigma = 0.09;

    // Goda-1999-version JONSWAP
    Sval = beta_J * pow(p->wHs, 2.0) * pow(p->wwp,4.0) * pow(w,-5.0) * exp(-1.25 * pow(w / p->wwp, -4.0)) * pow(p->B88, exp(-pow(w / p->wwp - 1.0, 2.0) / (2.0 * pow(sigma, 2.0))));

    return Sval;
}
