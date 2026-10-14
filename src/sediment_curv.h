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

#ifndef SEDIMENT_CURV_H_
#define SEDIMENT_CURV_H_

#include"increment.h"
#include"lexer.h"
#include<cmath>

class sediment_fdm;
class ghostcell;
class slice;

using namespace std;

// Sediment transport (S 10) on a curvilinear horizontal grid (G 50, p->geo_curv==1), NHFLOW and SFLOW.
//
// The sediment state stays cell-centred; the near-bed velocity of a cell is the Cartesian vector
// (ub,vb) of the flow (NHFLOW U,V of the bottom layer, SFLOW U,V). The face velocities of the
// Cartesian code become vectors:
//   - i-face i+1/2: Pvec = mean of the two cells, its normal component s->P = Pvec.n1 (the slice
//     boundary conditions of the Cartesian P act on it: walls 0); j-faces likewise (s->Q = Qvec.n2)
//   - the transport direction of a face (Exner S 31 1): the normal component of the face and the
//     tangential component of the mean of the four crossing faces (as the Cartesian code takes
//     P(i-1) and the mean of four Q): sg = u.n/|u|
//   - the velocity of a cell: the i-component from the mean of its two i-faces, the j-component
//     from the mean of its two j-faces (contravariant, with the cell's grid vectors g1, g2)
//   - bed gradients by the chain rule from index differences of the cell centres (one-sided next
//     to cells without a bed, as the Cartesian code)
//   - Exner: sum of face fluxes times face lengths over the cell area; filter, sandslide, mixture
//     with cell areas and centre distances
// On Cartesian nodes all of this is the Cartesian scheme (to round-off); the Cartesian path does not
// use it.

class sediment_curv : public increment
{
public:
    // unsupported options are refused
    static void check(lexer*);

    static inline int q(lexer *p, int ii, int jj)
    {
        return (ii-p->imin)*p->jmax + (jj-p->jmin);
    }

    // cell area (Cartesian: DXN*DYN)
    static inline double area(lexer *p, int ii, int jj)
    {
        if(p->geo_curv==1)
        return p->AREA2D[q(p,ii,jj)];

        return p->DXN[ii+marge]*p->DYN[jj+marge];
    }

    // distance between the centres of cell (ii,jj) and (ii+di,jj+dj)
    static inline double dist(lexer *p, int ii, int jj, int di, int dj)
    {
        const int a = q(p,ii,jj), b = q(p,ii+di,jj+dj);
        const double dx = p->XC2D[b]-p->XC2D[a], dy = p->YC2D[b]-p->YC2D[a];
        return sqrt(dx*dx + dy*dy);
    }

    // area ratio cell / neighbour (sandslide: a slid height on the receiver)
    static inline double aratio(lexer *p, int ii, int jj, int di, int dj)
    {
        return p->AREA2D[q(p,ii,jj)]/p->AREA2D[q(p,ii+di,jj+dj)];
    }

    // unit normal (towards +i / +j) and length of face i+1/2 (dir 1) or j+1/2 (dir 2) of cell (ii,jj)
    static inline void normal(lexer *p, int dir, int ii, int jj, double &nx, double &ny, double &len)
    {
        const int c = q(p,ii,jj);
        len = (dir==1) ? p->LEN1[c] : p->LEN2[c];

        if(len>0.0)
        {
        nx = ((dir==1) ? p->SX1[c] : p->SX2[c])/len;
        ny = ((dir==1) ? p->SY1[c] : p->SY2[c])/len;
        }

        else
        {
        nx = (dir==1) ? 1.0 : 0.0;
        ny = (dir==1) ? 0.0 : 1.0;
        }
    }

    // s->P, s->Q: normal velocities of the faces from the cell velocities (ub,vb), with the slice
    // boundary conditions (gcv 10, 11); wet: SFLOW (both cells wet), 0 none
    static void face_normals(lexer*, ghostcell*, sediment_fdm*, int wetcheck);

    // face vector with the normal component of s->P / s->Q (dir 1: i-face i+1/2, dir 2: j-face j+1/2)
    static void facevec(lexer*, sediment_fdm*, int dir, int ii, int jj, double &u, double &v);

    // velocity of cell (ii,jj): i-part from the mean of its i-faces, j-part from its j-faces
    static void cellvel(lexer*, sediment_fdm*, int ii, int jj, double &u, double &v);

    // transport direction u.n/|u| of face i+1/2 (dir 1) or j+1/2 (dir 2) of cell (ii,jj)
    static double face_sign(lexer*, sediment_fdm*, int dir, int ii, int jj);

    // bed gradient of cell (ii,jj) (bedslope, bedload direction): central index differences,
    // one-sided next to cells without a bed (DFBED<0); zero=1: no gradient along a line closed on
    // both sides (bedload_direction_f), 0: the central difference stays
    static void bedgrad(lexer*, slice &f, int ii, int jj, int zero, double &fx, double &fy);

    // contravariant component s.g^dir of a vector at cell (ii,jj) (per index), cell grid vectors
    // from the centre differences over two cells
    static double contra(lexer*, int dir, int ii, int jj, double u, double v);

    // width of cell (ii,jj) across its i-faces (dir 1) or j-faces (dir 2): area / mean face length
    static double width(lexer*, int dir, int ii, int jj);
};

#endif
