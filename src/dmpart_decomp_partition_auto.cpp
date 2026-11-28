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
Architect: Hans Bihs
--------------------------------------------------------------------*/

//  M 20 1: the partition, chosen automatically.
//
//  The subdomains stay a tensor product (one neighbour per face).  Of all layouts mx x my x mz of the M 10
//  ranks (a direction split only where M 11-13 allow it, at least 3 cells per subdomain) the one with the
//  lightest heaviest subdomain is taken, with the cut planes placed so that the heaviest subdomain is as
//  light as possible, and among the layouts within 1 % of it the one with the least cut area.  It replaces
//  the layout of the least cut area with equal widths (the partition before) only when its heaviest
//  subdomain is at least 5 % lighter, so grids without a load imbalance stay as they were.
//
//  The work of a cell: 1 for an active cell, 0 for the others (solid).  For the column models of REEF3D
//  (ctrl.txt in the case folder: A 10 2 SFLOW, 3 FNPF, 5 NHFLOW) the work is that of the water columns:
//  1 per column, 0.1 for a dry column (geodat bed at or above the still water level F 60 of ctrl.txt, or
//  no active cell), spread over the cells of the column.  The solvers skip most of the work of the dry
//  columns; the 0.1 stands for the loops that still visit them.
//
//  The same code as DIVEMesh decomp_partition_auto.cpp, so that REEF3D M 40 decomposes the global grid
//  file into the same subdomains as DIVEMesh (the geodat bed from the section GEOB of the grid file).

#include "dmpart_decomp.h"
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

namespace dmpart
{

namespace
{
// least bottleneck cuts of [0,N) into m intervals for the profiles P[s][i] of ns slabs; wmin: least width
double cuts_1d(const std::vector<std::vector<double>> &P, int N, int m, int wmin, std::vector<int> &cut)
{
    const int ns = (int)P.size();

    auto greedy = [&](double B, std::vector<int> *out) -> bool
    {
        std::vector<double> acc(ns);
        int pos=0, nint=0;
        if(out)
        out->assign(1,0);
        while(pos<N)
        {
            std::fill(acc.begin(),acc.end(),0.0);
            int len=0;
            // cells left for the intervals after this one
            const int rest = (m-nint-1)*wmin;
            while(pos<N)
            {
                if(N-pos<=rest && len>=wmin)
                break;
                bool fits=true;
                for(int s=0; s<ns; ++s)
                if(acc[s]+P[s][pos]>B)
                {
                    fits=false;
                    break;
                }
                if(!fits && len>=wmin)
                break;
                if(!fits && len<wmin)
                return false;
                for(int s=0; s<ns; ++s)
                acc[s] += P[s][pos];
                ++pos;
                ++len;
            }
            ++nint;
            if(nint>m)
            return false;
            if(out)
            out->push_back(pos);
        }
        return true;
    };

    double lo=0.0, hi=0.0;
    for(int s=0; s<ns; ++s)
    {
        double sum=0.0;
        for(int i=0; i<N; ++i)
        {
            sum += P[s][i];
            lo = std::max(lo,P[s][i]);
        }
        hi = std::max(hi,sum);
    }
    hi *= 1.0000001;
    hi += 1.0e-12;

    if(!greedy(hi,nullptr))
    {
        cut.clear();
        return -1.0;
    }

    for(int it=0; it<80 && hi-lo>1.0e-12*hi; ++it)
    {
        const double mid = 0.5*(lo+hi);
        if(greedy(mid,nullptr))
        hi = mid;
        else
        lo = mid;
    }

    greedy(hi,&cut);

    // fewer intervals than m: split the widest ones (never heavier)
    while((int)cut.size()-1<m)
    {
        int wbest=-1, kbest=-1;
        for(int k=0; k+1<(int)cut.size(); ++k)
        if(cut[k+1]-cut[k]>wbest)
        {
            wbest = cut[k+1]-cut[k];
            kbest = k;
        }
        if(wbest<2)
        break;
        cut.insert(cut.begin()+kbest+1,cut[kbest]+wbest/2);
    }

    // the bottleneck of the cuts
    double mx=0.0;
    for(int k=0; k+1<(int)cut.size(); ++k)
    for(int s=0; s<ns; ++s)
    {
        double sum=0.0;
        for(int i=cut[k]; i<cut[k+1]; ++i)
        sum += P[s][i];
        mx = std::max(mx,sum);
    }
    return mx;
}

// the cuts C of the subdomains M[d] per direction for the cell weights w (i, j, k; Nd cells): the equal
// widths (h0 their heaviest subdomain), then the directions swept in turn, each with the least bottleneck
// cuts for the cuts of the others, until the heaviest does not get lighter; returns the heaviest
// subdomain, mean the mean
double optimise(const std::vector<double> &w, const int *Nd, const int *M, std::vector<int> *C, double &h0, double &mean)
{
    const int NX=Nd[0], NY=Nd[1], NZ=Nd[2];
    auto W = [&](int ii, int jj, int kk) -> double { return w[((size_t)ii*NY+jj)*NZ+kk]; };

    for(int d=0; d<3; ++d)
    {
        C[d].resize(M[d]+1);
        for(int n=0; n<=M[d]; ++n)
        C[d][n] = (int)(((long)n*Nd[d])/M[d]);
    }

    auto heaviest = [&](double &mn)
    {
        std::vector<double> sub((size_t)M[0]*M[1]*M[2],0.0);
        std::vector<int> ix(NX), jy(NY), kz(NZ);
        for(int n=0; n<M[0]; ++n) for(int q=C[0][n]; q<C[0][n+1]; ++q) ix[q]=n;
        for(int n=0; n<M[1]; ++n) for(int q=C[1][n]; q<C[1][n+1]; ++q) jy[q]=n;
        for(int n=0; n<M[2]; ++n) for(int q=C[2][n]; q<C[2][n+1]; ++q) kz[q]=n;
        double tot=0.0;
        for(int i=0; i<NX; ++i)
        for(int j=0; j<NY; ++j)
        for(int k=0; k<NZ; ++k)
        {
            sub[((size_t)ix[i]*M[1]+jy[j])*M[2]+kz[k]] += W(i,j,k);
            tot += W(i,j,k);
        }
        mn = tot/double(sub.size());
        return *std::max_element(sub.begin(),sub.end());
    };

    h0 = heaviest(mean);
    double hcur = h0;

    for(int sweep=0; sweep<20; ++sweep)
    {
        const double hstart = hcur;

        for(int d=0; d<3; ++d)
        {
            if(M[d]<=1)
            continue;

            const int e = (d+1)%3, f = (d+2)%3;
            const int ns = M[e]*M[f];
            std::vector<std::vector<double>> P(ns,std::vector<double>(Nd[d],0.0));
            std::vector<int> ie(Nd[e]), jf(Nd[f]);
            for(int n=0; n<M[e]; ++n) for(int q=C[e][n]; q<C[e][n+1]; ++q) ie[q]=n;
            for(int n=0; n<M[f]; ++n) for(int q=C[f][n]; q<C[f][n+1]; ++q) jf[q]=n;

            for(int i=0; i<NX; ++i)
            for(int j=0; j<NY; ++j)
            for(int k=0; k<NZ; ++k)
            {
                const int idx[3] = {i,j,k};
                P[ie[idx[e]]*M[f]+jf[idx[f]]][idx[d]] += W(i,j,k);
            }

            const int wmin = (Nd[d]>=3*M[d]) ? 3 : std::max(1,Nd[d]/M[d]);
            std::vector<int> cut;
            const double hd = cuts_1d(P,Nd[d],M[d],wmin,cut);

            if(hd>=0.0 && (int)cut.size()==M[d]+1 && cut[0]==0 && cut[M[d]]==Nd[d])
            {
                std::vector<int> old = C[d];
                C[d] = cut;
                double m2;
                const double hnew = heaviest(m2);
                if(hnew<=hcur*(1.0+1.0e-12))
                hcur = hnew;
                else
                C[d] = old;
            }
        }

        if(hcur>=hstart*(1.0-1.0e-6))
        break;
    }
    return hcur;
}

struct layout
{
    int T[3];
    double h, area;
};

// the lightest layout of m ranks: the least heaviest subdomain, among those within 1 % of it the least
// cut area; T[0] = 0: none
layout lightest(const std::vector<double> &w, const int *Nd, int m, const int *allow)
{
    std::vector<layout> ls;
    for(int a=1; a<=m; ++a)
    for(int b=1; a*b<=m; ++b)
    {
        if(m%(a*b)!=0)
        continue;
        layout L;
        L.T[0]=a; L.T[1]=b; L.T[2]=m/(a*b);
        bool ok = true;
        for(int d=0; d<3; ++d)
        if(L.T[d]>1 && (allow[d]==0 || Nd[d]<3*L.T[d]))
        ok = false;
        if(!ok)
        continue;
        std::vector<int> CC[3];
        double h0, mn;
        L.h = optimise(w,Nd,L.T,CC,h0,mn);
        L.area = 0.0;
        for(int d=0; d<3; ++d)
        L.area += double(L.T[d]-1)*double(Nd[(d+1)%3])*double(Nd[(d+2)%3]);
        ls.push_back(L);
    }

    layout b;
    b.T[0]=b.T[1]=b.T[2]=0;
    b.h=b.area=0.0;
    bool found=false;
    for(const layout &L : ls)
    if(!found || L.h<b.h)
    {
        b = L;
        found = true;
    }
    if(!found)
    return b;
    const double hb = b.h;
    for(const layout &L : ls)
    if(L.h<=1.01*hb && (L.area<b.area || (L.area==b.area && L.h<b.h)))
    b = L;
    return b;
}

// A 10 and F 60 of the REEF3D control file in the case folder (absent: 0 and false)
void reef3d_ctrl(int &a10, double &f60, bool &hasf60)
{
    a10 = 0;
    f60 = 0.0;
    hasf60 = false;
    std::ifstream in("ctrl.txt");
    std::string line;
    while(std::getline(in,line))
    {
        if(line.size()<2 || (line[0]!='A' && line[0]!='F'))
        continue;
        std::istringstream ls(line.substr(1));
        int num;
        if(!(ls>>num))
        continue;
        if(line[0]=='A' && num==10)
        ls>>a10;
        if(line[0]=='F' && num==60 && (ls>>f60))
        hasf60 = true;
    }
}
}

void decomp::partition_auto(lexer* p, dive* a)
{
    // the layout of the least cut area, equal widths: the partition before
    nodecalc(p,a);
    costfunc(p,a);

    const int NX=a->knox, NY=a->knoy, NZ=a->knoz;
    int M[3] = {a->mx,a->my,a->mz};
    // ---- the work of the cells
    int a10;
    double f60;
    bool hasf60;
    reef3d_ctrl(a10,f60,hasf60);
    const bool columns = (a10==2 || a10==3 || a10==5);

    std::vector<double> w((size_t)NX*NY*NZ,0.0);
    auto W = [&](int ii, int jj, int kk) -> double& { return w[((size_t)ii*NY+jj)*NZ+kk]; };

    int ndry=0;
    if(columns)
    {
        ILOOP
        JLOOP
        {
            int nact=0;
            KLOOP
            if(a->flag(i,j,k)>0)
            ++nact;
            const bool dry = (nact==0) || (hasf60 && p->Geobed!=nullptr && p->Geobed->size()==(size_t)NX*NY && (*p->Geobed)[(size_t)i*NY+j]>=f60);
            if(dry)
            ++ndry;
            KLOOP
            W(i,j,k) = (dry ? 0.1 : 1.0)/double(NZ);
        }
    }
    else
    LOOP
    W(i,j,k) = (a->flag(i,j,k)>0) ? 1.0 : 0.0;

    // without splits in z (M 13 0) the columns are summed: the same subdomains, much less work
    const int allow[3] = {p->M11, p->M12, p->M13};
    const bool col = (allow[2]==0 && M[2]==1);
    const int Ns[3] = {NX, NY, col ? 1 : NZ};
    if(col)
    {
        std::vector<double> wc((size_t)NX*NY,0.0);
        LOOP
        wc[(size_t)i*NY+j] += W(i,j,k);
        w.swap(wc);
    }

    // ---- the lightest layout, taken when it is at least 5 % lighter than that of the least cut area
    std::vector<int> C[3];
    double mean=0.0, h0=0.0;
    optimise(w,Ns,M,C,h0,mean);

    const layout b = lightest(w,Ns,M[0]*M[1]*M[2],allow);

    if(b.T[0]==0 || b.h>=0.95*h0)
    {
        partition(p,a);
        cout<<"M 20 1: partition "<<a->mx<<" x "<<a->my<<" x "<<a->mz<<", heaviest subdomain / mean "<<h0/mean;
        if(columns)
        cout<<" (columns of REEF3D A 10 "<<a10<<", "<<ndry<<" of "<<NX*NY<<" dry)";
        cout<<endl;
        ddout<<"M 20 1: partition "<<a->mx<<" "<<a->my<<" "<<a->mz<<", heaviest subdomain / mean "<<h0/mean<<endl;
        return;
    }

    for(int d=0; d<3; ++d)
    M[d] = b.T[d];
    double h1;
    const double hcur = optimise(w,Ns,M,C,h1,mean);
    if(col)
    C[2] = {0, NZ};
    a->mx=M[0];
    a->my=M[1];
    a->mz=M[2];

    cout<<"M 20 1: partition by the load "<<M[0]<<" x "<<M[1]<<" x "<<M[2]<<", heaviest subdomain / mean "<<h0/mean<<" -> "<<hcur/mean;
    if(columns)
    cout<<" (columns of REEF3D A 10 "<<a10<<", "<<ndry<<" of "<<NX*NY<<" dry)";
    cout<<endl;
    ddout<<"M 20 1: partition by the load "<<M[0]<<" "<<M[1]<<" "<<M[2]<<", heaviest subdomain / mean "<<hcur/mean<<endl;

    // ---- nodes and origins
    for(n=0; n<=a->mx; ++n)
    {
        a->xnode[n] = C[0][n];
        a->xorig[n] = (n==0) ? p->xmin : (n==a->mx ? p->xmax : p->XN[a->xnode[n]+marge]);
    }
    for(n=0; n<=a->my; ++n)
    {
        a->ynode[n] = C[1][n];
        a->yorig[n] = (n==0) ? p->ymin : (n==a->my ? p->ymax : p->YN[a->ynode[n]+marge]);
    }
    for(n=0; n<=a->mz; ++n)
    {
        a->znode[n] = C[2][n];
        a->zorig[n] = (n==0) ? p->zmin : (n==a->mz ? p->zmax : p->ZN[a->znode[n]+marge]);
    }

    for(aa=0;aa<=a->mx;++aa)
    ddout<<"xnode: "<<a->xnode[aa]<<"  xorig: "<<a->xorig[aa]<<endl;
    for(bb=0;bb<=a->my;++bb)
    ddout<<"ynode: "<<a->ynode[bb]<<"  yorig: "<<a->yorig[bb]<<endl;

    MALOOP
    {
        a->subgrid(i,j,k)=-1;
        a->subslice(i,j)=-1;
    }

    NLOOP
    SUBLOOP
    {
        a->subgrid(i,j,k)=PARANUM;
        a->subslice(i,j)=PARANUM;
    }
}
}
