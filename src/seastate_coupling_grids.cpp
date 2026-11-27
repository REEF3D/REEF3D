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

#include"seastate_coupling.h"
#include"seastate_remote.h"
#include"lexer.h"
#include"ghostcell.h"
#include"met_forcing.h"
#include<algorithm>
#include<iomanip>
#include<mpi.h>
#include<cmath>

// Remote waves (A 750 2) with a curvilinear or multi-block host (G 50 2, G 50 1 with blocks): the host maps both ways
// itself, at the physical cell centres XC2D, YC2D of every grid it has on the rank (XP, YP are the logical box of
// DIVEMesh there, and the grids of a multi-block host are no lattice). The waves send their fields as for a
// rectilinear host; the host sends per wave cell the sums of its cells instead of its own lattice
// (seastate_remote_wave::map_host_grids divides them):
//   the cells whose centre lies in the wave cell:   area, area WL, wet area, wet area eta, U, V
//   the cell that holds the wave cell centre:       count, WL, wet count, eta, U, V
// The waves take the area mean where host centres fall into the wave cell (a host finer than the waves), else the
// value of the host cell at the wave centre (a coarser host), else the own still water (outside the host). Wet: the
// wet flag of the host cell. The host cells of a grid: its own cells (no ghost cells) with flagslice4 > 0; a wave
// centre on a face between two grids is counted by both and averaged (the counts).
// Waves -> host: fields_at at the physical centres, ghost cells included, for every grid (the wave fields of the
// grids of the blocks are held by their own seastate objects, w of seastate_host_grid).

namespace
{
// faces of the cells of a lattice of centres (one centre: the whole line)
void faces(const std::vector<double> &c, std::vector<double> &f)
{
    const int n = int(c.size());
    f.assign(n+1,0.0);

    if(n==1)
    {
    f[0] = -1.0e300;
    f[1] = 1.0e300;
    return;
    }

    f[0] = c[0] - 0.5*(c[1]-c[0]);
    for(int k=1; k<n; ++k)
    f[k] = 0.5*(c[k-1]+c[k]);
    f[n] = c[n-1] + 0.5*(c[n-1]-c[n-2]);
}

// the cell of the lattice that holds x (the last face closed), -1 outside
int cell_of(const std::vector<double> &f, double x)
{
    const int n = int(f.size())-1;

    if(n<1 || x<f[0] || x>f[n])
    return -1;

    const int a = int(std::upper_bound(f.begin(),f.end(),x) - f.begin()) - 1;
    return std::min(std::max(a,0),n-1);
}

void centre(lexer *p, int i, int j, double &x, double &y)
{
    if(p->geo_curv==1)
    {
    x = p->XC2D[IJ];
    y = p->YC2D[IJ];
    }
    else
    {
    x = p->XP[i+increment::marge];
    y = p->YP[j+increment::marge];
    }
}
}

void seastate_coupling::remote_ini_grids(lexer *p, ghostcell *pgc, const std::vector<seastate_host_grid> &G)
{
    hgrids = true;

    remote_setup(p,pgc);

    // the wave cells of the own host cells and the own host cells at the wave centres (the grids do not move)
    faces(rx,wfx);
    faces(ry,wfy);

    const int nyw = int(ry.size());
    long nown=0, nin=0;

    hbin.assign(G.size(),std::vector<int>());
    hsamp.clear();

    for(size_t g=0; g<G.size(); ++g)
    {
    lexer *q = G[g].p;
    std::vector<int> &B = hbin[g];

        for(int ii=0; ii<q->knox; ++ii)
        for(int jj=0; jj<q->knoy; ++jj)
        {
        const int s = (ii-q->imin)*q->jmax + (jj-q->jmin);

            if(q->flagslice4[s]<=0)
            continue;

        ++nown;

        double x, y;
        centre(q,ii,jj,x,y);

        const int a = cell_of(wfx,x), b = cell_of(wfy,y);

            if(a<0 || b<0)
            continue;

        B.push_back(ii);
        B.push_back(jj);
        B.push_back(a*nyw + b);
        ++nin;
        }

        // the wave centres in the box of the grid's nodes: the cell of the grid that holds them
        double x0, x1, y0, y1;
        q->geo_box(x0,x1,y0,y1);

        const int a0 = int(std::lower_bound(rx.begin(),rx.end(),x0) - rx.begin());
        const int a1 = int(std::upper_bound(rx.begin(),rx.end(),x1) - rx.begin());
        int b0 = 0, b1 = nyw;

        if(nyw>1)
        {
        b0 = int(std::lower_bound(ry.begin(),ry.end(),y0) - ry.begin());
        b1 = int(std::upper_bound(ry.begin(),ry.end(),y1) - ry.begin());
        }

        for(int a=a0; a<a1; ++a)
        for(int b=b0; b<b1; ++b)
        {
        int ii, jj;

            // a single row of waves: the row of the host at its own y (the lattice is one line)
            const double yw = (nyw>1) ? ry[b] : 0.5*(y0+y1);

            if(!q->geo_cell(rx[a],yw,ii,jj))
            continue;

            if(ii<0 || ii>=q->knox || jj<0 || jj>=q->knoy)
            continue;

            if(q->flagslice4[(ii-q->imin)*q->jmax + (jj-q->jmin)]<=0)
            continue;

        hsamp.push_back(int(g));
        hsamp.push_back(ii);
        hsamp.push_back(jj);
        hsamp.push_back(a*nyw + b);
        }
    }

    long cnt[3] = {nown, nin, long(hsamp.size()/4)}, tot[3] = {0,0,0};
    MPI_Allreduce(cnt,tot,3,MPI_LONG,MPI_SUM,pgc->mpi_comm);

    long ng = long(G.size()), ngt = 0;
    MPI_Allreduce(&ng,&ngt,1,MPI_LONG,MPI_SUM,pgc->mpi_comm);

    if(p->mpirank==0)
    cout<<"SEASTATE remote: host "<<ngt<<" grid piece(s), "<<tot[0]<<" cells, "<<tot[1]<<" of them in the wave grid; "
        <<tot[2]<<" wave cell centres in host cells"<<endl;

    remote_step_grids(p,pgc,G,true);
}

void seastate_coupling::remote_step_grids(lexer *p, ghostcell *pgc, const std::vector<seastate_host_grid> &G, bool force)
{
    if(!force && p->simtime < tnext - 1.0e-9*p->A706)
    return;

    const size_t Nw = rx.size()*ry.size();
    std::vector<double> loc(2+size_t(hg_nsc)*Nw,0.0);
    double *S[hg_nsc];
    for(int n=0; n<hg_nsc; ++n)
    S[n] = &loc[2+size_t(n)*Nw];

    // the host cells inside the wave cells: area-weighted sums
    for(size_t g=0; g<G.size(); ++g)
    {
    lexer *q = G[g].p;
    const std::vector<int> &B = hbin[g];
    const seastate_host_grid &H = G[g];

        for(size_t m=0; m<B.size(); m+=3)
        {
        const int ii = B[m], jj = B[m+1];
        const size_t w = size_t(B[m+2]);
        const double A = q->cellarea(ii,jj), wl = (*H.WL)(ii,jj);

        S[0][w] += A;
        S[1][w] += A*wl;

            if(q->wet[(ii-q->imin)*q->jmax + (jj-q->jmin)]==1)
            {
            S[2][w] += A;
            S[3][w] += A*(*H.eta)(ii,jj);
            S[4][w] += A*(*H.U)(ii,jj);
            S[5][w] += A*(*H.V)(ii,jj);
            }
        }
    }

    // the host cells at the wave centres
    for(size_t m=0; m<hsamp.size(); m+=4)
    {
    const seastate_host_grid &H = G[size_t(hsamp[m])];
    lexer *q = H.p;
    const int ii = hsamp[m+1], jj = hsamp[m+2];
    const size_t w = size_t(hsamp[m+3]);

    S[6][w] += 1.0;
    S[7][w] += (*H.WL)(ii,jj);

        if(q->wet[(ii-q->imin)*q->jmax + (jj-q->jmin)]==1)
        {
        S[8][w] += 1.0;
        S[9][w] += (*H.eta)(ii,jj);
        S[10][w] += (*H.U)(ii,jj);
        S[11][w] += (*H.V)(ii,jj);
        }
    }

    sbuf.assign(loc.size(),0.0);
    MPI_Reduce(loc.data(),sbuf.data(),int(loc.size()),MPI_DOUBLE,MPI_SUM,0,pgc->mpi_comm);
    sbuf[0] = 1.0;
    sbuf[1] = p->simtime;

    exchange(p,pgc,sbuf,rbuf,4961,4962);

    // the wave fields at the cells of every grid
    for(const seastate_host_grid &H : G)
    {
    fields_at(H.p,H.w,true);

        if(H.w!=this)
        {
        ++H.w->nupd;

            if(H.w->pmetw!=nullptr)
            H.w->pmetw->waves_update(H.p);
        }
    }

    ++nupd;

    if(pmetw!=nullptr)
    pmetw->waves_update(p);

    tlast = p->simtime;
    tnext = p->simtime + p->A706;
    ++nstep;

    if(p->mpirank==0)
    {
    std::ios fmt(nullptr);
    fmt.copyfmt(cout);
    cout<<"SEASTATE remote exchange "<<nstep<<" at simtime "<<std::fixed<<std::setprecision(1)<<p->simtime<<" s (wave time "<<rbuf[1]
        <<" s), forcing ramp "<<std::setprecision(3)<<ramp(p)<<endl;
    cout.copyfmt(fmt);
    }
}
