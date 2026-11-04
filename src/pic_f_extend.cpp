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

#include"pic_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include<vector>

// The velocity of the water for the particles (F 81 1, 2; Kelly et al. 2015, section 3.3):
// a velocity face is in the water when phi at the face (the mean of its two cells) is >= 0;
// the faces in the air get, layer by layer (3 layers), the mean over the directions in which the
// neighbour face of the same component is in the water or got a value in an earlier layer: the
// linear extrapolation from that face and the next one, if that has a value too, else the
// neighbour's value (a constant extension moved the surface of a rotating disk by a cell per turn). The particles
// close to the surface then move with the water, not with the light air of the two-phase
// velocity, whose spurious velocities at the interface would otherwise move the surface.
// The layers run over the faces of the subdomain; after each the values (gcval 10, 11, 12) and the
// marks of the faces of the ghost cells come from the neighbours, so the layers do not depend on the
// decomposition (marked from phi in the ghost cells alone, the third layer at a subdomain boundary
// lacked the cells beyond the margin).

namespace
{
inline int fid(lexer *p, int i, int j, int k)
{
    return (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + k-p->kmin;
}
}

void pic_f::water_velocity(lexer *p, fdm *a, ghostcell *pgc)
{
    const int layers = 3;
    const int fm = fmarg;
    const int jlo = p->j_dir==1 ? -fm : 0;
    const int jhi = p->j_dir==1 ? p->knoy+fm-1 : 0;
    const int size = p->imax*p->jmax*p->kmax;

    for(int c=0; c<3; ++c)
    {
        if(c==1 && p->j_dir==0)
        {
        for(int n=0; n<size; ++n)
        vw->V[n] = 0.0;
        continue;
        }

        field &src = c==0 ? static_cast<field&>(a->u) : (c==1 ? static_cast<field&>(a->v) : static_cast<field&>(a->w));
        field &dst = c==0 ? static_cast<field&>(*uw) : (c==1 ? static_cast<field&>(*vw) : static_cast<field&>(*ww));
        const int di = c==0, dj = c==1, dk = c==2;

        for(int n=0; n<size; ++n)
        dst.V[n] = src.V[n];

        // the faces in the water, with the margin (phi is known there)
        std::vector<char> known(size,0);

        for(int i=-fm; i<p->knox+fm-di; ++i)
        for(int j=jlo; j<=jhi-(p->j_dir==1 ? dj : 0); ++j)
        for(int k=-fm; k<p->knoz+fm-dk; ++k)
        if(0.5*(a->phi(i,j,k) + a->phi(i+di,j+dj,k+dk)) >= 0.0)
        known[fid(p,i,j,k)] = 1;

        field &km = c==0 ? static_cast<field&>(*km1) : (c==1 ? static_cast<field&>(*km2) : static_cast<field&>(*km3));

        for(int L=0; L<layers; ++L)
        {
            std::vector<char> next = known;

            // the faces of this subdomain; those of the ghost cells, value and mark, from the neighbours
            for(int i=0; i<p->knox; ++i)
            for(int j=0; j<p->knoy; ++j)
            for(int k=0; k<p->knoz; ++k)
            {
                const int f = fid(p,i,j,k);

                if(known[f])
                continue;

                double sum=0.0;
                int cnt=0;

                // from the neighbour in a direction: linear from the two faces there, else its value
                auto add = [&](int si, int sj, int sk)
                {
                    const int g1 = fid(p,i+si,j+sj,k+sk);

                    if(!known[g1])
                    return;

                    const int i2 = i+2*si, j2 = j+2*sj, k2 = k+2*sk;
                    const bool inside = i2>=-fm && i2<p->knox+fm && k2>=-fm && k2<p->knoz+fm
                                     && (p->j_dir==0 || (j2>=jlo && j2<=jhi));

                    if(inside && known[fid(p,i2,j2,k2)])
                    sum += 2.0*dst.V[g1] - dst.V[fid(p,i2,j2,k2)];
                    else
                    sum += dst.V[g1];

                    ++cnt;
                };

                add(-1,0,0); add(1,0,0);
                add(0,0,-1); add(0,0,1);

                if(p->j_dir==1)
                {
                add(0,-1,0); add(0,1,0);
                }

                if(cnt>0)
                {
                dst.V[f] = sum/double(cnt);
                next[f] = 1;
                }
            }

            known.swap(next);

            for(int n=0; n<size; ++n)
            km.V[n] = known[n] ? 1.0 : 0.0;

            pgc->gcparax(p,km,c+1);

            for(int n=0; n<size; ++n)
            known[n] = km.V[n]>0.5 ? 1 : 0;

            if(c==0)
            pgc->start1(p,dst,10);
            if(c==1)
            pgc->start2(p,dst,11);
            if(c==2)
            pgc->start3(p,dst,12);
        }
    }
}
