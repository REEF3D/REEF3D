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

//  M 20 5: weighted rectilinear partition.
//
//  The subdomains stay a tensor product (one neighbour per face, as for M 20 1-4); the cut planes
//  are placed so that the heaviest subdomain is as light as possible.  The weight of a cell:
//    1 for an active cell, M 21 for the others (solid, default 0),
//    with M 35 1: the load map of a REEF3D run (REEF3D G 45: the work per step of every level-0
//      column, the patches of the mesh refinement included), spread over the active cells of the
//      column, read from ./loadmap.dat or ./REEF3D_<NHFLOW|FNPF|SFLOW>_AMR/REEF3D_<..>_AMR_loadmap.dat;
//      a map of the cells (header nx ny nz: CFD, REEF3D_CFD_AMR/REEF3D_CFD_AMR_loadmap.dat) gives
//      the work of every active cell,
//    times the factors of the weight boxes M 34 xs xe ys ye zs ze f (cell centre in the box).
//  The subdomains per direction: M 30 if given, else as for M 20 1 (M 11-13).
//  The cuts of one direction for fixed cuts of the others are optimal: the least bottleneck B for
//  which the greedy sweep (each slab of the other directions at most B per subdomain) needs at most
//  m subdomains, by bisection.  The directions are swept in turn until the heaviest subdomain does
//  not get lighter.  Subdomains are at least 3 cells wide where the grid allows it.

#include "dmpart_decomp.h"
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>

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
}

void decomp::partition_weighted(lexer* p, dive* a)
{
    // ---- subdomains per direction
    if(p->M30_x>0 && p->M30_y>0 && p->M30_z>0)
    {
        a->mx=p->M30_x;
        a->my=p->M30_y;
        a->mz=p->M30_z;
    }
    else
    {
        nodecalc(p,a);
        costfunc(p,a);
    }

    const int NX=a->knox, NY=a->knoy, NZ=a->knoz;
    const int M[3] = {a->mx,a->my,a->mz};
    const int Nd[3] = {NX,NY,NZ};

    // ---- weights
    std::vector<double> w((size_t)NX*NY*NZ);
    auto W = [&](int ii, int jj, int kk) -> double& { return w[((size_t)ii*NY+jj)*NZ+kk]; };

    LOOP
    W(i,j,k) = (a->flag(i,j,k)>0) ? 1.0 : p->M21;

    if(p->M35==1)
    {
        const char *names[5] = {"./loadmap.dat",
                                "./REEF3D_NHFLOW_AMR/REEF3D_NHFLOW_AMR_loadmap.dat",
                                "./REEF3D_FNPF_AMR/REEF3D_FNPF_AMR_loadmap.dat",
                                "./REEF3D_SFLOW_AMR/REEF3D_SFLOW_AMR_loadmap.dat",
                                "./REEF3D_CFD_AMR/REEF3D_CFD_AMR_loadmap.dat"};
        ifstream in;
        std::string used;
        for(int f=0; f<5 && !in.is_open(); ++f)
        {
            in.open(names[f]);
            if(in.is_open())
            used = names[f];
        }

        if(!in.is_open())
        {
            cout<<"M 35 1: no load map found (loadmap.dat or REEF3D_<M>_AMR/REEF3D_<M>_AMR_loadmap.dat) -- the active cells only"<<endl;
            ddout<<"M 35 1: no load map found"<<endl;
        }
        else
        {
            std::string line;
            int gnx=-1, gny=-1, gnz=-1;     // gnz > 0: a map of the cells (CFD), else of the columns
            std::vector<double> map;
            while(std::getline(in,line))
            {
                if(line.empty() || line[0]=='#')
                continue;
                std::istringstream ls(line);
                if(gnx<0)
                {
                    ls>>gnx>>gny;
                    if(!(ls>>gnz))
                    gnz=-1;
                    continue;
                }
                double v;
                while(ls>>v)
                map.push_back(v);
            }

            // a map of the cells: the work of every cell (the inactive ones keep M 21)
            if(gnz>0)
            {
                if(gnx!=NX || gny!=NY || gnz!=NZ || (long)map.size()!=(long)NX*NY*NZ)
                {
                    cout<<"M 35 1: the load map "<<used<<" is "<<gnx<<" x "<<gny<<" x "<<gnz<<" cells, the grid "<<NX<<" x "<<NY<<" x "<<NZ
                        <<" -- the active cells only"<<endl;
                    ddout<<"M 35 1: load map does not fit the grid"<<endl;
                }
                else
                {
                    LOOP
                    if(a->flag(i,j,k)>0)
                    W(i,j,k) = map[((size_t)i*NY+j)*NZ+k];
                    cout<<"M 35 1: load map "<<used<<" ("<<NX<<" x "<<NY<<" x "<<NZ<<" cells)"<<endl;
                    ddout<<"M 35 1: load map "<<used<<endl;
                }
            }
            else if(gnx!=NX || gny!=NY || (long)map.size()!=(long)NX*NY)
            {
                cout<<"M 35 1: the load map "<<used<<" is "<<gnx<<" x "<<gny<<" columns ("<<map.size()<<" values), the grid "<<NX<<" x "<<NY
                    <<" -- the active cells only"<<endl;
                ddout<<"M 35 1: load map does not fit the grid"<<endl;
            }
            else
            {
                double mean=0.0;
                for(double v : map)
                mean += v;
                mean /= double(map.size());

                ILOOP
                JLOOP
                {
                    int nact=0;
                    KLOOP
                    if(a->flag(i,j,k)>0)
                    ++nact;
                    const double v = map[(size_t)i*NY+j];
                    if(nact>0)
                    KLOOP
                    if(a->flag(i,j,k)>0)
                    W(i,j,k) = v/double(nact);
                }
                cout<<"M 35 1: load map "<<used<<" ("<<NX<<" x "<<NY<<" columns, mean "<<mean<<" per column)"<<endl;
                ddout<<"M 35 1: load map "<<used<<endl;
            }
        }
    }

    for(int b=0; b<p->M34; ++b)
    LOOP
    {
        const double xc = 0.5*(p->XN[i+marge]+p->XN[i+1+marge]);
        const double yc = 0.5*(p->YN[j+marge]+p->YN[j+1+marge]);
        const double zc = 0.5*(p->ZN[k+marge]+p->ZN[k+1+marge]);
        if(xc>=p->M34_xs[b] && xc<p->M34_xe[b] && yc>=p->M34_ys[b] && yc<p->M34_ye[b] && zc>=p->M34_zs[b] && zc<p->M34_ze[b])
        W(i,j,k) *= p->M34_f[b];
    }

    // ---- cuts: start from the equal widths of M 20 1
    std::vector<int> C[3];
    for(int d=0; d<3; ++d)
    {
        C[d].resize(M[d]+1);
        for(int n=0; n<=M[d]; ++n)
        C[d][n] = (int)(((long)n*Nd[d])/M[d]);
    }

    // the heaviest subdomain and the mean
    auto heaviest = [&](double &mean)
    {
        std::vector<double> sub((size_t)M[0]*M[1]*M[2],0.0);
        std::vector<int> ix(NX), jy(NY), kz(NZ);
        for(int n=0; n<M[0]; ++n) for(int q=C[0][n]; q<C[0][n+1]; ++q) ix[q]=n;
        for(int n=0; n<M[1]; ++n) for(int q=C[1][n]; q<C[1][n+1]; ++q) jy[q]=n;
        for(int n=0; n<M[2]; ++n) for(int q=C[2][n]; q<C[2][n+1]; ++q) kz[q]=n;
        double tot=0.0;
        LOOP
        {
            sub[((size_t)ix[i]*M[1]+jy[j])*M[2]+kz[k]] += W(i,j,k);
            tot += W(i,j,k);
        }
        mean = tot/double(sub.size());
        return *std::max_element(sub.begin(),sub.end());
    };

    double mean=0.0;
    const double h0 = heaviest(mean);
    double hcur = h0;

    for(int sweep=0; sweep<20; ++sweep)
    {
        const double hstart = hcur;

        for(int d=0; d<3; ++d)
        {
            if(M[d]<=1)
            continue;

            // the slabs: the subdomains of the other two directions
            const int e = (d+1)%3, f = (d+2)%3;
            const int ns = M[e]*M[f];
            std::vector<std::vector<double>> P(ns,std::vector<double>(Nd[d],0.0));
            std::vector<int> ie(Nd[e]), jf(Nd[f]);
            for(int n=0; n<M[e]; ++n) for(int q=C[e][n]; q<C[e][n+1]; ++q) ie[q]=n;
            for(int n=0; n<M[f]; ++n) for(int q=C[f][n]; q<C[f][n+1]; ++q) jf[q]=n;

            LOOP
            {
                const int idx[3] = {i,j,k};
                P[ie[idx[e]]*M[f]+jf[idx[f]]][idx[d]] += W(i,j,k);
            }

            const int wmin = (Nd[d]>=3*M[d]) ? 3 : std::max(1,Nd[d]/M[d]);
            std::vector<int> cut;
            const double hd = cuts_1d(P,Nd[d],M[d],wmin,cut);

            // take the new cuts when they are complete and not heavier
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

    cout<<"M 20 5: weighted partition "<<M[0]<<" x "<<M[1]<<" x "<<M[2]<<", heaviest subdomain / mean "<<h0/mean<<" (equal widths) -> "<<hcur/mean<<endl;
    ddout<<"M 20 5: heaviest subdomain / mean "<<h0/mean<<" (equal widths) -> "<<hcur/mean<<endl;

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
