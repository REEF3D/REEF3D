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
#include"ioflow.h"
#include<algorithm>
#include<cmath>
#include<mpi.h>

// Particles kept in the water (every F 86 steps; Edwards & Bridson 2012): the particles drift
// apart in places and gather in others, and with the free surface from the particles an empty
// cell in the water is a bubble. Cells inside the water (phi > h, so the surface is not touched)
// with fewer than a quarter of F 83^dim particles (empty ones in 2D with F 83 2; a rotating
// lattice leaves cells with one) get new ones in their empty sub-cells, up to F 83^dim; cells
// with more than twice as many lose those with the largest ids. The choice
// depends only on the particles in the cell, so it is the same for any decomposition; the new
// particles sit at the centres of their sub-cells and take the velocity of the water and ids
// after those of the seeding.
// In the relaxation zones of the waves the columns of the particles follow the relaxed surface
// (pic_f::zone_stretch); no particles are removed there: removing those above the relaxed surface
// (patch 0005) and only stretching the columns where it is higher drained the zone of the wave
// generation once the waves reflected from the beach of the spilling breaker came back into it.
// interior: the particles in the water (F 86).

// in a relaxation zone: the cell column of the position (in the subdomain)
bool pic_f::inzone(lexer *p, double x, double y)
{
    if(!relaxzones)
    return false;

    if(zcol.empty())
    {
        zcol.assign(size_t(p->knox)*p->knoy,0);

        for(int i=0; i<p->knox; ++i)
        for(int j=0; j<p->knoy; ++j)
        zcol[size_t(i)*p->knoy + j] = pflow->relax_zone(p->XP[i+marge],p->YP[j+marge]) ? 1 : 0;
    }

    const int ci = std::min(std::max(cell(p->XN,p->knox,x),0),p->knox-1);
    const int cj = p->j_dir==1 ? std::min(std::max(cell(p->YN,p->knoy,y),0),p->knoy-1) : 0;

    return zcol[size_t(ci)*p->knoy + cj]!=0;
}

void pic_f::reseed(lexer *p, fdm *a, ghostcell *pgc, bool interior)
{
    const int n = p->knox*p->knoy*p->knoz;
    const int ny = p->j_dir==1 ? npd : 1;

    std::vector<int> count(n,0);
    std::vector<int> where(P.size(),-1);
    std::vector<std::vector<char> > used;
    std::vector<int> slot(n,-1);

    auto deep = [&](int ci, int cj, int ck)
    {
        if(p->flag4[(ci-p->imin)*p->jmax*p->kmax + (cj-p->jmin)*p->kmax + ck-p->kmin]<0)
        return false;

        const double h = p->j_dir==1 ? std::cbrt(p->DXN[ci+marge]*p->DYN[cj+marge]*p->DZN[ck+marge])
                                     : std::sqrt(p->DXN[ci+marge]*p->DZN[ck+marge]);
        return a->phi(ci,cj,ck) > h;
    };

    // the particles of each cell and the sub-cells they are in
    for(size_t q=0; q<P.size(); ++q)
    {
        const int ci = cell(p->XN,p->knox,P[q].x);
        const int cj = p->j_dir==1 ? cell(p->YN,p->knoy,P[q].y) : 0;
        const int ck = cell(p->ZN,p->knoz,P[q].z);

        if(ci<0 || ci>=p->knox || cj<0 || cj>=p->knoy || ck<0 || ck>=p->knoz)
        continue;

        const int l = lid(p,ci,cj,ck);
        where[q] = l;
        ++count[l];

        if(slot[l]<0)
        {
        slot[l] = int(used.size());
        used.emplace_back(npc,0);
        }

        const int qi = std::min(npd-1, int(double(npd)*(P[q].x - p->XN[ci+marge])/p->DXN[ci+marge]));
        const int qj = p->j_dir==1 ? std::min(npd-1, int(double(npd)*(P[q].y - p->YN[cj+marge])/p->DYN[cj+marge])) : 0;
        const int qk = std::min(npd-1, int(double(npd)*(P[q].z - p->ZN[ck+marge])/p->DZN[ck+marge]));
        used[slot[l]][(std::max(qi,0)*ny + std::max(qj,0))*npd + std::max(qk,0)] = 1;
    }

    // crowded cells in the water: the particles with the largest ids go
    std::vector<char> drop(P.size(),0);
    long long mine_removed = 0;

    if(interior)
    {
        std::vector<std::vector<size_t> > members(n);

        for(size_t q=0; q<P.size(); ++q)
        if(where[q]>=0 && count[where[q]]>2*npc)
        members[where[q]].push_back(q);

        for(int ci=0; ci<p->knox; ++ci)
        for(int cj=0; cj<p->knoy; ++cj)
        for(int ck=0; ck<p->knoz; ++ck)
        {
            auto &m = members[lid(p,ci,cj,ck)];

            if(m.empty() || !deep(ci,cj,ck))
            continue;

            // not twice
            m.erase(std::remove_if(m.begin(),m.end(),[&](size_t s){ return drop[s]!=0; }),m.end());

            if(m.size()<=size_t(2*npc))
            continue;

            std::sort(m.begin(), m.end(), [&](size_t s, size_t t){ return P[s].id > P[t].id; });

            for(size_t r=0; r<m.size()-size_t(2*npc); ++r)
            {
            drop[m[r]] = 1;
            ++mine_removed;
            }
        }
    }

    // sparse cells in the water: new particles in the empty sub-cells
    const double idbase = double(p->gknox)*double(p->gknoy)*double(p->gknoz)*double(npc);
    const double gcells = double(p->gknox)*double(p->gknoy)*double(p->gknoz);
    std::vector<pic_data> fresh;

    for(int ci=0; ci<p->knox; ++ci)
    for(int cj=0; cj<p->knoy; ++cj)
    for(int ck=0; ck<p->knoz; ++ck)
    {
        const int l = lid(p,ci,cj,ck);
        if(!interior || 4*count[l]>=npc || !deep(ci,cj,ck))
        continue;

        const long long gcell = ((long long)(ck+p->origin_k)*(long long)p->gknoy + (long long)(cj+p->origin_j))*(long long)p->gknox
                              + (long long)(ci+p->origin_i);
        int have = count[l];

        for(int qi=0; qi<npd && have<npc; ++qi)
        for(int qj=0; qj<ny && have<npc; ++qj)
        for(int qk=0; qk<npd && have<npc; ++qk)
        {
            const int s = (qi*ny + qj)*npd + qk;

            if(slot[l]>=0 && used[slot[l]][s])
            continue;

            pic_data q;
            q.id = idbase + (double(p->count)*gcells + double(gcell))*double(npc) + double(s);

            q.x = p->XN[ci+marge] + (double(qi)+0.5)/double(npd)*p->DXN[ci+marge];
            q.y = p->j_dir==1 ? p->YN[cj+marge] + (double(qj)+0.5)/double(npd)*p->DYN[cj+marge] : p->YP[cj+marge];
            q.z = p->ZN[ck+marge] + (double(qk)+0.5)/double(npd)*p->DZN[ck+marge];

            if(blocked(p,a,q.x,q.y,q.z))
            continue;

            velocity(p,a,q.x,q.y,q.z,q.u,q.v,q.w);
            if(apic)
            affine(p,q);
            fresh.push_back(q);
            ++have;
        }
    }

    if(mine_removed>0)
    {
        size_t w=0;
        for(size_t q=0; q<P.size(); ++q)
        if(!drop[q])
        P[w++] = P[q];
        P.resize(w);
    }

    const long long mine_added = (long long)fresh.size();
    P.insert(P.end(), fresh.begin(), fresh.end());

    long long loc[2] = {mine_added, mine_removed}, glob[2] = {0,0};
    MPI_Allreduce(loc,glob,2,MPI_LONG_LONG,MPI_SUM,pgc->mpi_comm);
    added = glob[0];
    removed = glob[1];
    added_total += added;
    removed_total += removed;
}
