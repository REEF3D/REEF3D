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
#include<algorithm>
#include<cmath>
#include<mpi.h>

// The particles of the inflow (B 60; F 81 1, 2): the ghost cells of the inflow boundary (x-min,
// the cells of lexer::gcin, 3 layers) hold particles too, in rows at the heights (and y) of the
// sub-cells of the seeding. Every step each row gets new particles one spacing of the seeding
// (h / F 83) behind its last one, up to the far end of the third layer, where the inflow is
// water (phi of the ghost cells, ioflow::fsfinflow); an empty row starts half a spacing outside
// the boundary. They move with the velocity of the inflow there (ioflow::inflow) and enter the
// domain at the spacing of the seeding, so they bring the discharge in. (Refilling every sub-cell
// whose particle had moved on, the first version, gave them half the spacing and twice the
// discharge.) They stand in for the mirror images at that boundary in the water fraction
// (pic_f::phi_build), and leave (pic_f::domain) only beyond the third layer. Their ids are
// negative, so they do not meet those of the seeding and of pic_f::reseed.

void pic_f::inflow_seed(lexer *p, fdm *a, ghostcell *pgc)
{
    long long mine = 0;
    const int nl = 3;
    const int ny = p->j_dir==1 ? npd : 1;

    if(xloc_s==xglo_s && p->gcin_count>0)
    {
        const double xend = p->XN[marge-nl];

        // the last (smallest x) particle of each row outside the boundary: (j, k, sub j, sub k)
        const int nj = p->knoy, nk = p->knoz;
        std::vector<double> last(size_t(nj)*nk*ny*npd, 1.0e30);

        auto row = [&](int cj, int ck, int qj, int qk) { return ((size_t(cj)*nk + ck)*ny + qj)*npd + qk; };

        for(const auto &q : P)
        {
            if(q.x>=xglo_s)
            continue;

            const int cj = p->j_dir==1 ? cell(p->YN,p->knoy,q.y) : 0;
            const int ck = cell(p->ZN,p->knoz,q.z);

            if(cj<0 || cj>=nj || ck<0 || ck>=nk)
            continue;

            const int qj = p->j_dir==1 ? std::min(npd-1, std::max(0, int(double(npd)*(q.y - p->YN[cj+marge])/p->DYN[cj+marge]))) : 0;
            const int qk = std::min(npd-1, std::max(0, int(double(npd)*(q.z - p->ZN[ck+marge])/p->DZN[ck+marge])));

            double &l = last[row(cj,ck,qj,qk)];
            l = std::min(l,q.x);
        }

        const double gin = double(p->gknoy)*double(p->gknoz)*double(ny*npd);

        for(int n=0; n<p->gcin_count; ++n)
        {
            if(p->gcin[n][3]!=1)
            continue;

            const int i0 = p->gcin[n][0];
            const int cj = p->gcin[n][1];
            const int ck = p->gcin[n][2];
            const double sx = p->DXN[i0-1+marge]/double(npd);

            for(int qj=0; qj<ny; ++qj)
            for(int qk=0; qk<npd; ++qk)
            {
                double &l = last[row(cj,ck,qj,qk)];

                if(l>1.0e29)
                l = xglo_s + 0.5*sx;

                const double y = p->j_dir==1 ? p->YN[cj+marge] + (double(qj)+0.5)/double(npd)*p->DYN[cj+marge] : p->YP[cj+marge];
                const double z = p->ZN[ck+marge] + (double(qk)+0.5)/double(npd)*p->DZN[ck+marge];
                const double gidx = ((double(cj+p->origin_j)*double(p->gknoz) + double(ck+p->origin_k))*double(ny) + double(qj))*double(npd) + double(qk);
                int m = 0;

                for(double x=l-sx; x>=xend; x-=sx, ++m)
                {
                    if(ipol(p,a->phi,4,x,y,z)<=0.0 || ipol(p,a->topo,4,x,y,z)<0.0)
                    continue;

                    pic_data q;
                    q.x = x;
                    q.y = y;
                    q.z = z;
                    q.id = -1.0 - ((double(p->count)*gin + gidx)*double(4*nl*npd) + double(m));
                    velocity(p,a,q.x,q.y,q.z,q.u,q.v,q.w);
                    if(apic)
                    affine(p,q);
                    P.push_back(q);
                    ++mine;
                }
            }
        }
    }

    const long long all = (long long)pgc->globalsum(double(mine));
    added += all;
    added_total += all;
}
