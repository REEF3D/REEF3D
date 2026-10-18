/*--------------------------------------------------------------------
DIVEMesh
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D (port of the DIVEMesh decomposition, dmpart).

DIVEMesh is free software; you can redistribute it and/or modify it
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

#include "dmpart_decomp.h"
#include <algorithm>
#include <vector>
#include <cmath>

namespace dmpart
{

void decomp::partition_correction(lexer* p, dive* a)
{
    cout<<"maindir: "<<maindir<<endl;

    if(maindir==1)
    partition_correct_x(p,a);
    else if(maindir==2)
    partition_correct_y(p,a);
    else if(maindir==3)
    partition_correct_z(p,a);
}

namespace
{
    // Split planes 0..N-1 with weights w into m contiguous, non-empty parts
    // (node[0]=0 < node[1] < ... < node[m]=N) such that the largest part is
    // as small as possible:
    // 1) binary search for the smallest cap B that m parts can meet
    // 2) greedy packing up to B, leaving at least one plane per remaining part
    // 3) move each boundary to where it best balances its two neighbours;
    //    a move only happens if it lowers the larger of the two, so the
    //    largest part never grows
    void balance_1d(const std::vector<long long> &w, int m, int *node)
    {
        const int N = int(w.size());

        std::vector<long long> P(N+1,0);
        for(int i=0;i<N;++i)
        P[i+1] = P[i] + w[i];

        // last split b in [lo,hi] with P[b] <= limit (at least lo)
        auto last_below = [&](int lo, int hi, long long limit)
        {
            int b = int(std::upper_bound(P.begin()+lo, P.begin()+hi+1, limit) - P.begin()) - 1;
            return std::max(b,lo);
        };

        auto pack = [&](long long cap, int *nd)
        {
            nd[0]=0;
            for(int a=1;a<m;++a)
            nd[a] = last_below(nd[a-1]+1, N-(m-a), P[nd[a-1]]+cap);
            nd[m]=N;

            return P[N]-P[nd[m-1]] <= cap;
        };

        long long lo = *std::max_element(w.begin(),w.end());
        long long hi = P[N];
        std::vector<int> tmp(m+1);

        while(lo<hi)
        {
            const long long mid = lo + (hi-lo)/2;

            if(pack(mid,tmp.data()))
            hi = mid;
            else
            lo = mid+1;
        }

        pack(lo,node);

        for(int sweep=0; sweep<100*m; ++sweep)
        {
            bool moved = false;

            for(int a=1;a<m;++a)
            {
                const int l = node[a-1], h = node[a+1];
                const long long half = (P[l]+P[h])/2;

                int b = last_below(l+1, h-1, half);
                if(b+1<=h-1 && std::max(P[b+1]-P[l],P[h]-P[b+1]) < std::max(P[b]-P[l],P[h]-P[b]))
                ++b;

                const long long now  = std::max(P[node[a]]-P[l], P[h]-P[node[a]]);
                const long long then = std::max(P[b]-P[l], P[h]-P[b]);

                if(then<now)
                {
                    node[a] = b;
                    moved = true;
                }
            }

            if(!moved)
            break;
        }
    }
}

void decomp::partition_balance(lexer* p, dive* a, int dir)
{
    const char name = dir==1 ? 'x' : (dir==2 ? 'y' : 'z');
    const int N = dir==1 ? a->knox : (dir==2 ? a->knoy : a->knoz);
    const int m = dir==1 ? a->mx : (dir==2 ? a->my : a->mz);
    int *node = dir==1 ? a->xnode : (dir==2 ? a->ynode : a->znode);
    double *orig = dir==1 ? a->xorig : (dir==2 ? a->yorig : a->zorig);
    double *coord = dir==1 ? p->XN : (dir==2 ? p->YN : p->ZN);
    int *cnt = dir==1 ? xcount : (dir==2 ? ycount : zcount);
    int *cross = dir==1 ? xcross : (dir==2 ? ycross : zcross);

    cout<<"starting "<<name<<"-dir partition correction"<<endl;

    // active cells per plane normal to the correction direction
    std::vector<long long> w(N,0);

    LOOP
    if(a->flag(i,j,k)>0 && a->solid(i,j,k)>0)
    ++w[dir==1 ? i : (dir==2 ? j : k)];

    long long total=0;
    for(int n=0;n<N;++n)
    {
        cross[n] = int(w[n]);
        total += w[n];
    }

    auto slab_counts = [&]()
    {
        cnt[0]=0;
        for(int s=1;s<=m;++s)
        {
            long long c=0;
            for(int n=node[s-1];n<node[s];++n)
            c+=w[n];
            cnt[s]=int(c);
        }
    };

    slab_counts();
    for(int s=0;s<=m;++s)
    ddout<<"old "<<name<<"count"<<s<<" : "<<cnt[s]<<"  "<<name<<"node: "<<node[s]<<"  "<<name<<"orig: "<<orig[s]<<endl;
    ddout<<name<<"average: "<<double(total)/double(m)<<endl;

    std::vector<int> node_in(node,node+m+1);

    balance_1d(w,m,node);

    auto largest = [&](const int *nd)
    {
        long long mx=0, c;
        for(int s=1;s<=m;++s)
        {
            c=0;
            for(int n=nd[s-1];n<nd[s];++n)
            c+=w[n];
            mx=std::max(mx,c);
        }
        return mx;
    };

    if(largest(node)>largest(node_in.data()))
    std::copy(node_in.begin(),node_in.end(),node);

    slab_counts();
    for(int s=0;s<=m;++s)
    orig[s] = coord[node[s]+marge];

    for(int s=0;s<=m;++s)
    ddout<<"new "<<name<<"count"<<s<<" : "<<cnt[s]<<"  "<<name<<"node: "<<node[s]<<"  "<<name<<"orig: "<<orig[s]<<endl;

    // subdomain ids
    MALOOP
    a->subgrid(i,j,k)=-1;

    NLOOP
    SUBLOOP
    {
        a->subgrid(i,j,k)=PARANUM;
        a->subslice(i,j)=PARANUM;
    }

    // report
    for(int q=0;q<p->M10;++q)
    subcell[q]=0;

    int q=0;
    NLOOP
    {
        SUBLOOP
        if(a->flag(i,j,k)>0 && a->solid(i,j,k)>0)
        ++subcell[q];

        ++q;
    }

    int maxcell=0;
    for(q=0;q<p->M10;++q)
    {
        maxcell = std::max(maxcell,subcell[q]);
        ddout<<q<<" new subcell_count: "<<subcell[q]<<endl;
    }

    cout<<name<<"-dir correction: largest subdomain "<<maxcell<<" active cells, average "
        <<double(total)/double(p->M10)<<" (ratio "<<double(maxcell)*double(p->M10)/std::max(1.0,double(total))<<")"<<endl;
}
}
