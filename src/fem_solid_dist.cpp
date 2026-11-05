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

#include"fem_solid.h"
#include<cstring>
#include<algorithm>

// Distribution of the solid over the ranks (MPI ranks of the flow solver). Every rank keeps
// the whole solid; the deformable elements are split into contiguous blocks, one per rank
// (owner), and a node belongs to the owner of its first deformable element. A rank computes
// the elements with a node it owns (its own and the neighbours at the block border) and moves
// the nodes it owns; after every substep the owners send the positions and velocities of
// their nodes to all ranks, so contact, rigid bodies and supports are computed from the same
// state everywhere. The node sums run over the same elements in the same order as on one
// rank: the results do not depend on the number of ranks (or threads). Failed elements are
// exchanged every substep, the integration points at the end of the fluid step.

double fem_solid::dissipated_energy() const
{
    double w = wdiss;
    for(double q : wde)
    w += q;
    return w;
}

void fem_solid::dist_partition()
{
    nranks = comm->size();
    myrank = comm->rank();
    const int ne = nelem(), nn = nnode(), P = nranks;

    std::vector<int> D;
    for(int e=0; e<ne; ++e)
    if(elems[e].alive && !elems[e].rigid)
    D.push_back(e);
    const long nd = (long)D.size();

    eown.assign(ne,-1);
    relems.assign(P,std::vector<int>());
    for(int r=0; r<P; ++r)
    for(long k=r*nd/P; k<(r+1)*nd/P; ++k)
    {
        eown[D[k]] = r;
        relems[r].push_back(D[k]);
    }

    nown.assign(nn,-1);
    rnodes.assign(P,std::vector<int>());
    for(int i=0; i<nn; ++i)
    {
        if(!rnode.empty() && rnode[i]>=0)
        continue;
        int o = -1;
        for(int q=node_elem_start[i]; q<node_elem_start[i+1] && o<0; ++q)
        o = eown[node_elem[q]];
        if(o<0)
        o = (int)(((long)i*P)/nn);       // debris, nodes of failed elements
        nown[i] = o;
        rnodes[o].push_back(i);
    }

    ecomp.assign(ne,0);
    for(int e : D)
    for(int a=0; a<8; ++a)
    if(nown[elems[e].n[a]]==myrank)
    {
        ecomp[e] = 1;
        break;
    }
    act_dirty = true;
}

void fem_solid::dist_exchange_nodes(bool end)
{
    // 6 doubles per owned node: x and v (substep), vbar and fint (end of the fluid step)
    const int P = nranks;
    const int rec = 6*sizeof(double);
    cnt.resize(P);
    dsp.resize(P);
    int tot = 0;
    for(int r=0; r<P; ++r)
    {
        cnt[r] = (int)rnodes[r].size()*rec;
        dsp[r] = tot;
        tot += cnt[r];
    }
    sbuf.resize(std::max(1,cnt[myrank]));
    rbuf.resize(std::max(1,tot));
    double *s = reinterpret_cast<double*>(sbuf.data());
    for(int i : rnodes[myrank])
    {
        const Vec3& a = end ? vbar[i] : x[i];
        const Vec3& b = end ? fint[i] : v[i];
        *s++ = a(0); *s++ = a(1); *s++ = a(2);
        *s++ = b(0); *s++ = b(1); *s++ = b(2);
    }
    comm->allgatherv(sbuf.data(),cnt[myrank],rbuf.data(),cnt.data(),dsp.data());
    for(int r=0; r<P; ++r)
    {
        const double *q = reinterpret_cast<const double*>(rbuf.data()+dsp[r]);
        for(int i : rnodes[r])
        {
            Vec3& a = end ? vbar[i] : x[i];
            Vec3& b = end ? fint[i] : v[i];
            a = Vec3(q[0],q[1],q[2]);
            b = Vec3(q[3],q[4],q[5]);
            q += 6;
        }
    }
}

void fem_solid::dist_exchange_elements()
{
    // per owned element: the stored integration points, von Mises, J, the dissipated
    // energy and the integration flag
    const int P = nranks;
    const size_t gsz = (size_t)ngp*sizeof(gpstate);
    const int rec = (int)(gsz + 4*sizeof(double));
    cnt.resize(P);
    dsp.resize(P);
    int tot = 0;
    for(int r=0; r<P; ++r)
    {
        cnt[r] = (int)relems[r].size()*rec;
        dsp[r] = tot;
        tot += cnt[r];
    }
    sbuf.resize(std::max(1,cnt[myrank]));
    rbuf.resize(std::max(1,tot));
    char *s = sbuf.data();
    for(int e : relems[myrank])
    {
        std::memcpy(s,&gps[(size_t)e*ngp],gsz);
        const double d[4] = {elems[e].svm, elems[e].J, wde[e], elems[e].full ? 1.0 : 0.0};
        std::memcpy(s+gsz,d,sizeof(d));
        s += rec;
    }
    comm->allgatherv(sbuf.data(),cnt[myrank],rbuf.data(),cnt.data(),dsp.data());
    for(int r=0; r<P; ++r)
    {
        const char *q = rbuf.data()+dsp[r];
        for(int e : relems[r])
        {
            std::memcpy(&gps[(size_t)e*ngp],q,gsz);
            double d[4];
            std::memcpy(d,q+gsz,sizeof(d));
            elems[e].svm = d[0];
            elems[e].J = d[1];
            wde[e] = d[2];
            elems[e].full = d[3]>0.5;
            q += rec;
        }
    }
}

void fem_solid::dist_failures(std::vector<int>& failed)
{
    // failed: the elements that failed here in this substep; on return: on any rank,
    // ascending, without repeats
    const int P = nranks;
    std::vector<int> n(P);
    comm->allgather_int((int)failed.size(),n.data());
    int tot = 0;
    for(int r=0; r<P; ++r) tot += n[r];
    if(tot==0)
    return;
    cnt.resize(P);
    dsp.resize(P);
    int off = 0;
    for(int r=0; r<P; ++r)
    {
        cnt[r] = n[r]*(int)sizeof(int);
        dsp[r] = off;
        off += cnt[r];
    }
    std::vector<int> all(tot);
    comm->allgatherv(failed.data(),cnt[myrank],all.data(),cnt.data(),dsp.data());
    std::sort(all.begin(),all.end());
    all.erase(std::unique(all.begin(),all.end()),all.end());
    failed.swap(all);
}
