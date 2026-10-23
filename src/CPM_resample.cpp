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

#include"CPM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<algorithm>
#include<cmath>

/*--------------------------------------------------------------------
Parcel merge and split (Q 78 N)

All seeded parcels have the same volume V_1 = ParcelFactor V_grain, set by the smallest seeded cell
(Q 24 parcels per smallest cell). A cell of volume V then holds theta_0 V/V_1 parcels when packed: on a
graded grid the large cells hold many times Q 24 parcels, which cost but carry no information, and on a
refined grid (AMR) the parcels of a coarse cell would be too few for the fine cells after a refinement.

A parcel carries a weight w (P.Wt): its volume is w V_1, w = 1 for a seeded parcel. The target weight of a
cell is the largest power of two with

    w_c theta_0 V/V_1 >= ... :   w_c = 2^m <= theta_0 V/(V_1 N),   m >= 0

so that a packed cell holds about N parcels (N = Q 78), and at least one: w_c = 1 in the cells of the size
the parcels were seeded for.

Once per step:
  - split: a parcel with w > 1 that is awake (not held at rest by the jammed bed: moving, in a hop, in the
    suspension) or heavier than the target weight of its cell is split into two parcels of weight w/2, until
    it is neither. The two halves keep the velocity and all other values; they sit at -+ a quarter of the
    parcel spacing (w V_1/theta_0)^(1/dim) along one axis (z, x, y in turn with the generation), within the
    cell. The parcels the grid-limited step and the layer see near the bed (erosion, pickup, the active
    columns) are thus seeded parcels.
  - merge: resting grains of the jammed bed (the criterion of the sleeping grains, Q 58, S 10 1, not with
    Q 67 or Q 69) of equal weight w below the target weight of their cell, equal diameter and density, are
    merged in pairs into one parcel of weight 2 w at the middle of the two (in the cell, ordered by height),
    with the mean of their velocities. Repeated until the cell's target weight is reached.

The volume of the sediment is conserved exactly (weights are whole numbers, the volumes of the cells and of
the bed follow the weights). The momentum is conserved by both. The positions of the merged grains inside
their cell are averaged: only the distribution within a resting cell changes, by less than a cell.
The weights are whole numbers, so the arrival counts of the grid-limited step remain exact.
Q 78 0: no weights other than 1; all results are those without merge and split.
--------------------------------------------------------------------*/

double CPM::wcell(lexer *p, int ii, int jj, int kk)
{
    if(p->Q78<=0)
    return 1.0;

    const double vpar = P.ParcelFactor*Vp;
    const int ny = p->knoy+2, nz = p->knoz+2;

    auto calc = [&](int qi, int qj, int qk) -> double
    {
        const double V = p->DXN[qi+marge]*p->DYN[qj+marge]*p->DZN[qk+marge];
        const double r = theta_0*V/(vpar*double(p->Q78));
        double w = 1.0;

        while(2.0*w <= r*(1.0+1.0e-9) && w<1048576.0)
        w *= 2.0;

        return w;
    };

    if(wcc.size()!=size_t(p->knox+2)*ny*nz || wcc_vpar!=vpar)
    {
        wcc_vpar = vpar;
        wcc.resize(size_t(p->knox+2)*ny*nz);

        for(int qi=-1;qi<=p->knox;++qi)
        for(int qj=-1;qj<=p->knoy;++qj)
        for(int qk=-1;qk<=p->knoz;++qk)
        wcc[(size_t(qi+1)*ny + qj+1)*nz + qk+1] = calc(qi,qj,qk);
    }

    if(ii>=-1 && ii<=p->knox && jj>=-1 && jj<=p->knoy && kk>=-1 && kk<=p->knoz)
    return wcc[(size_t(ii+1)*ny + jj+1)*nz + kk+1];

    return calc(ii,jj,kk);
}

// a grain held at rest by the jammed bed (the criterion of bed_asleep, without its counters)
bool CPM::jammed(lexer *p, fdm *a, int q)
{
    if(p->Q58<=0 || p->S10!=1 || p->Q67>0 || p->Q69>0 || P.Hop[q]!=0.0 || P.Test[q]<0.0)
    return false;

    const int kk = pck(p,q,P.Z[q]);
    const double psi = p->Q64*p->DZN[kk+marge];

    return ptopo_p(p,a,q,P.X[q],P.Y[q],P.Z[q]) + 1.5*psi <= -0.5*psi;
}

void CPM::parcel_resample(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->Q78<=0)
    return;

    const double vpar = P.ParcelFactor*Vp;
    const int dim = p->j_dir==1 ? 3 : 2;

    index_cache(p);

    auto cell = [&](int q, int &ci0, int &cj0, int &ck0)
    {
        ci0 = MAX(0, MIN(p->knox-1, pci(p,q,P.X[q])));
        cj0 = p->j_dir==1 ? MAX(0, MIN(p->knoy-1, pcj(p,q,P.Y[q]))) : 0;
        ck0 = MAX(0, MIN(p->knoz-1, pck(p,q,P.Z[q])));
    };

    // ---- split
    std::vector<int> &hv = rs_list;
    hv.clear();

    for(n=0;n<P.index;++n)
    if(P.Flag[n]==ACTIVE && P.Wt[n]>1.0)
    hv.push_back(n);

    while(!hv.empty())
    {
        const int q = hv.back();
        hv.pop_back();

        int ci0,cj0,ck0;
        cell(q,ci0,cj0,ck0);

        if(P.Wt[q]<=1.0 || (jammed(p,a,q) && P.Wt[q]<=wcell(p,ci0,cj0,ck0)))
        continue;

        if(P.index_empty<=0)
        P.resize(p,P.capacity+1000);

        --P.index_empty;
        const int c = P.Empty[P.index_empty];

        const double w = 0.5*P.Wt[q];

        // the axis by the generation: z, x, y in turn
        const int gen = int(std::lround(std::log2(P.Wt[q]))) - 1;
        const int axes3[3] = {2,0,1}, axes2[2] = {2,0};
        const int ax = dim==3 ? axes3[gen%3] : axes2[gen%2];

        const double s = dim==3 ? std::cbrt(P.Wt[q]*vpar/theta_0) : std::sqrt(P.Wt[q]*vpar/(theta_0*p->DYN[cj0+marge]));
        double lo, hi;

        if(ax==0) { lo = p->XN[ci0+marge]; hi = p->XN[ci0+1+marge]; }
        else if(ax==1) { lo = p->YN[cj0+marge]; hi = p->YN[cj0+1+marge]; }
        else { lo = p->ZN[ck0+marge]; hi = p->ZN[ck0+1+marge]; }

        const double eps = 1.0e-3*(hi-lo);
        double *xa = ax==0 ? P.X : (ax==1 ? P.Y : P.Z);
        double *xr = ax==0 ? P.XRK1 : (ax==1 ? P.YRK1 : P.ZRK1);
        const double x0 = xa[q];
        const double d = 0.25*s;

        // copy of the parcel into the new slot
        P.X[c]=P.X[q]; P.Y[c]=P.Y[q]; P.Z[c]=P.Z[q];
        P.XRK1[c]=P.XRK1[q]; P.YRK1[c]=P.YRK1[q]; P.ZRK1[c]=P.ZRK1[q];
        P.U[c]=P.U[q]; P.V[c]=P.V[q]; P.W[c]=P.W[q];
        P.URK1[c]=P.URK1[q]; P.VRK1[c]=P.VRK1[q]; P.WRK1[c]=P.WRK1[q];
        P.Uf[c]=P.Uf[q]; P.Vf[c]=P.Vf[q]; P.Wf[c]=P.Wf[q];
        P.D[c]=P.D[q]; P.RO[c]=P.RO[q]; P.Test[c]=P.Test[q]; P.Hop[c]=P.Hop[q];
        P.Flag[c]=P.Flag[q];

        xa[q] = xr[q] = MAX(lo+eps, MIN(hi-eps, x0 - d));
        xa[c] = xr[c] = MAX(lo+eps, MIN(hi-eps, x0 + d));
        P.Wt[q] = P.Wt[c] = w;

        ++rs_nsplit;

        if(w>1.0)
        {
            hv.push_back(q);
            hv.push_back(c);
        }
    }

    // ---- merge: the resting grains below the target weight of their cell, per cell and weight
    struct cand { size_t cell; double w, d, ro, z, x, y; int q; };
    std::vector<cand> cs;
    const int ny = p->knoy, nz = p->knoz;

    for(int pass=0; pass<40; ++pass)
    {
        cs.clear();
        index_cache(p);

        for(n=0;n<P.index;++n)
        if(P.Flag[n]==ACTIVE)
        {
            int ci0,cj0,ck0;
            cell(n,ci0,cj0,ck0);

            if(P.Wt[n] >= wcell(p,ci0,cj0,ck0) || !jammed(p,a,n))
            continue;

            cs.push_back({(size_t(ci0)*ny + cj0)*nz + ck0, P.Wt[n], P.D[n], P.RO[n], P.Z[n], P.X[n], P.Y[n], int(n)});
        }

        std::sort(cs.begin(), cs.end(), [](const cand &a1, const cand &a2)
        {
            if(a1.cell!=a2.cell) return a1.cell<a2.cell;
            if(a1.w!=a2.w) return a1.w<a2.w;
            if(a1.d!=a2.d) return a1.d<a2.d;
            if(a1.ro!=a2.ro) return a1.ro<a2.ro;
            if(a1.z!=a2.z) return a1.z<a2.z;
            if(a1.x!=a2.x) return a1.x<a2.x;
            if(a1.y!=a2.y) return a1.y<a2.y;
            return a1.q<a2.q;
        });

        int nm=0;

        for(size_t m=0; m+1<cs.size(); )
        {
            const cand &A = cs[m], &B = cs[m+1];

            if(A.cell!=B.cell || A.w!=B.w || A.d!=B.d || A.ro!=B.ro)
            {
                ++m;
                continue;
            }

            const int qa = A.q, qb = B.q;

            P.X[qa] = P.XRK1[qa] = 0.5*(P.X[qa] + P.X[qb]);
            P.Y[qa] = P.YRK1[qa] = 0.5*(P.Y[qa] + P.Y[qb]);
            P.Z[qa] = P.ZRK1[qa] = 0.5*(P.Z[qa] + P.Z[qb]);
            P.U[qa] = P.URK1[qa] = 0.5*(P.U[qa] + P.U[qb]);
            P.V[qa] = P.VRK1[qa] = 0.5*(P.V[qa] + P.V[qb]);
            P.W[qa] = P.WRK1[qa] = 0.5*(P.W[qa] + P.W[qb]);
            P.Uf[qa] = 0.5*(P.Uf[qa] + P.Uf[qb]);
            P.Vf[qa] = 0.5*(P.Vf[qa] + P.Vf[qb]);
            P.Wf[qa] = 0.5*(P.Wf[qa] + P.Wf[qb]);
            P.Wt[qa] = 2.0*A.w;

            P.remove(qb);

            ++nm;
            m += 2;
        }

        rs_nmerge += nm;

        if(pgc->globalisum(nm)==0)
        break;
    }
}
