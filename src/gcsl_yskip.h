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

#ifndef GCSL_YSKIP_H_
#define GCSL_YSKIP_H_

#include"lexer.h"
#include"slice.h"

//  SFLOW on a grid with one cell in y (j_dir 0): every cell has two y-side boundary entries, and their
//  ghost cells took half the time of a run.  They are not written: the y sides of the slice ghost
//  cells (gcsl_start1/4) and of the SFLOW velocity and conserved variables (sflow_momentum_func::
//  vel_bc, ghostcells) are skipped.  The readers of these ghost rows read the cell row instead (the
//  y derivatives, 0 on such a grid, as 0; the node values of the output from row 0), or have the
//  rows copied from row 0 first (gcsl_yrows), which is what their Neumann ghost cells held.  Found
//  with a run that filled the skipped ghost cells with signalling NaN and trapped the first read;
//  all SFLOW and SEASTATE SFLOW cases bit for bit.
inline bool gcsl_yskip(const lexer *p) { return p->j_dir==0 && p->A10==2; }

// gcsl_yskip: the y ghost rows of f copied from the cell row, as the Neumann ghost cells, for a reader
// of them
inline void gcsl_yrows(const lexer *p, slice &f)
{
    if(!gcsl_yskip(p))
    return;
    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    if(jj<0 || jj>=p->knoy)
    f(ii,jj) = f(ii,(jj<0) ? 0 : p->knoy-1);
}

#endif
