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

#ifndef PROBE_CORE_H_
#define PROBE_CORE_H_

#include<fstream>

class lexer;
class ghostcell;

using namespace std;

// shared part of the point probes and gauges of all solvers

// probes (one file per probe): val[n*ncomp+c] holds component c of probe n on the rank that owns
// the probe and -1e20 elsewhere; one reduction for all probes, then rank 0 writes the row
// time \t val_0 \t val_1 ... to pout[n]
void probe_rows(lexer*, ghostcell*, ofstream *pout, double *val, int num, int ncomp);

// gauges (one file for all gauges): the row time \t val_0 \t val_1 ... (rank 0, values already reduced),
// flushed every fileFlushMaxCount gauges to limit data loss for many gauges
void gauge_row(lexer*, ofstream &out, const double *val, int num, int fileFlushMaxCount);

#endif
