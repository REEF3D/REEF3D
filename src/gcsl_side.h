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

#ifndef GCSL_SIDE_H_
#define GCSL_SIDE_H_

// The side of a 2D boundary cell of the lists gcbsl4, gcslin, gcslout (entry [3]): 1 x-, 2 y+, 3 y-, 4 x+.
// The boundary cell (i,j) is inside the domain; its ghost cells lie q = 1, 2, 3 cells across the side,
// (i + q di, j + q dj). The boundary face is an x-face (slice1, P) for the sides 1 and 4 and a y-face
// (slice2, Q) for 2 and 3, with the index (fi(i), fj(j)). inward() is +1 on the - sides and -1 on the
// + sides: the axis velocity times inward() is the velocity into the domain.
struct gcsl_side
{
    int cs, di, dj;

    explicit gcsl_side(int c) : cs(c), di(c==1 ? -1 : (c==4 ? 1 : 0)), dj(c==3 ? -1 : (c==2 ? 1 : 0)) {}

    bool xdir() const { return cs!=2 && cs!=3; }
    double inward() const { return (cs==4 || cs==2) ? -1.0 : 1.0; }

    int gi(int i, int q=1) const { return i + q*di; }
    int gj(int j, int q=1) const { return j + q*dj; }

    int fi(int i) const { return cs==1 ? i-1 : i; }
    int fj(int j) const { return cs==3 ? j-1 : j; }
};

#endif
