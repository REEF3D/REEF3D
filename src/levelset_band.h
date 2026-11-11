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

#ifndef LEVELSET_BAND_H_
#define LEVELSET_BAND_H_

#include<vector>
#include<cmath>
#include"lexer.h"
#include"field.h"

// F 48: the narrow band of the level set.  The mask of the cells of the array of f (the cells around
// the grid included) whose level set is below w cells (of the smallest cell size of the grid), and
// two layers of neighbours around them: a cell next to the band whose value is not a distance (it
// was left alone while it was outside) is taken in and corrected from inside the band before the
// interface gets there.  work: a scratch array
inline void levelset_band_mask(lexer *p, field &f, double w, std::vector<unsigned char> &band, std::vector<unsigned char> &work)
{
    const size_t nn = (size_t)p->imax*p->jmax*p->kmax;
    band.resize(nn);
    const double bw = w*p->DXM;
    for(size_t q=0; q<nn; ++q)
    band[q] = std::fabs(f.data()[q])<bw ? 1 : 0;

    const long sI = (long)p->jmax*p->kmax, sJ = p->kmax;
    const bool y = (p->j_dir==1);
    for(int it=0; it<2; ++it)
    {
        work = band;
        for(int ii=1; ii<p->imax-1; ++ii)
        for(int jj=(y ? 1 : 0); jj<(y ? p->jmax-1 : p->jmax); ++jj)
        for(int kk=1; kk<p->kmax-1; ++kk)
        {
            const long q = ii*sI + jj*sJ + kk;
            if(work[q])
            continue;
            if(work[q-sI] || work[q+sI] || work[q-1] || work[q+1] || (y && (work[q-sJ] || work[q+sJ])))
            band[q] = 1;
        }
    }
}

// the mask the convection of the level set (ipol 4) keeps to while it is set (momentum_rk sets it
// around its transport of the level set, per grid; nullptr: all cells)
extern const unsigned char *levelset_band_active;

#endif
