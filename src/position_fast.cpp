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

#include"position.h"
#include"lexer.h"

/*--------------------------------------------------------------------
cell index lookup without the bisection (posc_*, posf_*)

Inside the grid the index is the unique interval of a strictly increasing coordinate array,
A[ii] <= x < A[ii+1]. A table over [A[lo], A[hi]] in buckets no wider than the smallest
interval (at most 16 buckets per interval on strongly stretched grids) gives the interval at
the bucket start; a short walk corrects it. The table is built at the first call and rebuilt
when the array or its end points change. The callers use it only where their bisection cannot
leave through an out-of-bounds branch, so the result is the same as before.
--------------------------------------------------------------------*/

int position::fast_index(fastfind &t, const double *A, int lo, int hi, double x)
{
    const int m = marge;

    if(t.A!=A || t.lo!=lo || t.hi!=hi || t.X0!=A[lo+m] || t.X1!=A[hi+m])
    {
        t.A = A;
        t.lo = lo;
        t.hi = hi;
        t.X0 = A[lo+m];
        t.X1 = A[hi+m];

        double dmin = t.X1 - t.X0;

        for(int q=lo; q<hi; ++q)
        dmin = MIN(dmin, A[q+1+m] - A[q+m]);

        int n = MAX(hi-lo,1);
        int nb = 1;

        if(dmin>0.0)
        nb = int(MIN((t.X1-t.X0)/dmin + 1.0, 16.0*double(n) + 1.0));

        nb = MAX(nb,1);
        t.w = (t.X1-t.X0)/double(nb);
        t.tbl.assign(nb+1, lo);

        int q = lo;

        for(int b=0; b<=nb; ++b)
        {
            double xb = t.X0 + double(b)*t.w;

            while(q+1<hi && A[q+1+m]<=xb)
            ++q;

            t.tbl[b] = q;
        }
    }

    int b = t.w>0.0 ? int((x-t.X0)/t.w) : 0;
    b = MAX(0, MIN(int(t.tbl.size())-1, b));

    int q = t.tbl[b];

    while(q+1<hi && A[q+1+m]<=x)
    ++q;

    while(q>lo && A[q+m]>x)
    --q;

    return q;
}
