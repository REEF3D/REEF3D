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

#include"net_membrane.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"
#include<mpi.h>
#include<cmath>
#include<iostream>
#include<algorithm>

// Link mode (membrane.dat 'mobility link').
//
// In the default layer mode the implicit porous factor beta = 1/(1 + a K_n H) is the mobility of every cell of the
// smeared layer in the pressure Poisson equation. It is isotropic, so the pressure can not move the layer fluid along
// the membrane either, and that fluid has to be held with the membrane (R_t = R_n): the flow sees a body that is
// thicker by the layer, with a no-slip surface smeared over ~1.5 delta (rigid cage of Strand et al. 2013 in current:
// C_D 1.26 at delta = 1.5 dx, 0.96 at 1.0 dx, measured 0.66).
//
// Link mode keeps the normal resistance of the layer in the momentum forcing, but in the projection only the links
// whose two end points lie on opposite sides of the membrane are restricted:
//   horizontal links  between the cell centres (i,j,k) and (i+1,j,k) / (i,j+1,k)      d->MBX, d->MBY
//   vertical links    between the nodes k and k+1 of a column (through cell k)        d->MBZ
// with the porous-jump mobility beta = 1/(1 + a R_n / l) (l: link length), i.e. a flux (dp/rho)/R_n across the link.
// All other links keep mobility 1, so the pressure acts on the tangential velocity in the layer as anywhere else and
// the layer fluid slides along the membrane (R_t = 0 by default). Every path from inside to outside crosses a blocked
// link, so the sealing does not depend on the layer width. Sides: closest point on the membrane and the angle-weighted
// pseudo normal (Baerentzen & Aanaes 2005), robust at the floor edge and other corners.
//
// Loads: the momentum the forcing takes out of the layer cells (normal resistance) plus the pressure difference across
// the blocked links times the link face area, F = (P_a - P_b) A e on the membrane.

void net_membrane::link_ini(lexer *p)
{
    // edge -> triangles (the topology of the mesh does not change)
    etri_.assign(edge_.size(),{-1,-1});

    for(size_t t=0; t<tri_.size(); ++t)
    for(int q=0; q<3; ++q)
    {
        const int e = tedge_[t][q];

        if(etri_[e][0]<0)
        etri_[e][0]=t;
        else
        etri_[e][1]=t;
    }

    p->Darray(sideC_,p->imax*p->jmax*(p->kmax+2));
    sideN_.assign(p->imax*p->jmax*p->kmaxF,0);

    pseudonormals();
}

void net_membrane::pseudonormals()
{
    vpn_.assign(x_.size(),Eigen::Vector3d::Zero());
    epn_.assign(edge_.size(),Eigen::Vector3d::Zero());

    for(size_t t=0; t<tri_.size(); ++t)
    for(int q=0; q<3; ++q)
    {
        const int a = tri_[t][q], b = tri_[t][(q+1)%3], c = tri_[t][(q+2)%3];
        const Eigen::Vector3d e1 = x_[b]-x_[a], e2 = x_[c]-x_[a];
        const double n1 = e1.norm(), n2 = e2.norm();

        if(n1<1.0e-20 || n2<1.0e-20)
        continue;

        const double ang = acos(MAX(-1.0,MIN(1.0,e1.dot(e2)/(n1*n2))));
        vpn_[a] += ang*tn_[t];

        epn_[tedge_[t][q]] += tn_[t];
    }

    for(auto &v : vpn_)
    if(v.norm()>1.0e-20)
    v.normalize();

    for(auto &v : epn_)
    if(v.norm()>1.0e-20)
    v.normalize();
}

int net_membrane::side_of(const Eigen::Vector3d &P, int t, double u, double v, double w) const
{
    // closest point u a + v b + w c on triangle t; normal of the feature it lies on
    const Eigen::Vector3d C = u*x_[tri_[t][0]] + v*x_[tri_[t][1]] + w*x_[tri_[t][2]];
    const double eps = 1.0e-9;

    Eigen::Vector3d N = tn_[t];

    if(v<eps && w<eps)      N = vpn_[tri_[t][0]];
    else if(u<eps && w<eps) N = vpn_[tri_[t][1]];
    else if(u<eps && v<eps) N = vpn_[tri_[t][2]];
    else if(w<eps)          N = epn_[tedge_[t][0]];
    else if(u<eps)          N = epn_[tedge_[t][1]];
    else if(v<eps)          N = epn_[tedge_[t][2]];

    return (P-C).dot(N)>=0.0 ? 1 : -1;
}

int net_membrane::side_near(const Eigen::Vector3d &P, const cellentry &e) const
{
    // a point next to a layer cell (its nodes): closest among the triangles closest to the cell
    double dbest=1.0e20, ub=1.0, vb=0.0, wb=0.0, u, v, w;
    int tb=e.tc;

    auto test = [&](int t)
    {
        const Eigen::Vector3d C = closest_point(P,x_[tri_[t][0]],x_[tri_[t][1]],x_[tri_[t][2]],u,v,w);
        const double dd = (P-C).squaredNorm();

        if(dd<dbest)
        {
            dbest=dd; tb=t; ub=u; vb=v; wb=w;
        }
    };

    test(e.tc);

    for(int r=0; r<e.ns; ++r)
    test(e.t[r]);

    return side_of(P,tb,ub,vb,wb);
}

void net_membrane::link_mobility(lexer *p, fdm_nhf *d, ghostcell *pgc, double a)
{
    if(moving())
    pseudonormals();

    // sides of the cell centres and nodes of the layer cells
    for(int q : sideCq_)
    sideC_[q]=0.0;

    for(int q : sideNq_)
    sideN_[q]=0;

    sideCq_.clear();
    sideNq_.clear();

    for(const auto &e : cells_)
    {
        i=e.i; j=e.j; k=e.k;

        const Eigen::Vector3d P(p->XP[IP], p->YP[JP], p->ZSP[IJK]);
        sideC_[IJK] = side_of(P,e.tc,e.w0,e.w1,e.w2);
        sideCq_.push_back(IJK);

        for(int kk=k; kk<=k+1; ++kk)
        {
            const int qf = (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + kk-p->kmin;

            if(sideN_[qf]!=0)
            continue;

            const Eigen::Vector3d N(p->XP[IP], p->YP[JP], p->ZSN[qf]);
            sideN_[qf] = side_near(N,e);
            sideNq_.push_back(qf);
        }
    }

    // neighbours across subdomain borders
    pgc->start4V(p,sideC_,1);

    // blocked links: porous-jump mobility 1/(1 + a R_n/l)
    blocked_.clear();

    auto beta = [&](double l) {return 1.0/(1.0 + a*prm.Rn/MAX(l,1.0e-20));};

    for(const auto &e : cells_)
    {
        i=e.i; j=e.j; k=e.k;

        if(p->wet[IJ]==0)
        continue;

        const double s = sideC_[IJK];
        double bmin=1.0;

        // +x, -x
        if(sideC_[Ip1JK]*s<0.0 && p->flag4[Ip1JK]>0)
        {
            const double b = beta(p->DXP[IP]);
            d->MBX[IJK] = MIN(d->MBX[IJK],b);
            bmin = MIN(bmin,b);
            blocked_.push_back({i,j,k,0,e.tc,{e.w0,e.w1,e.w2}});
        }

        if(sideC_[Im1JK]*s<0.0 && p->flag4[Im1JK]>0)
        {
            const double b = beta(p->DXP[IM1]);
            bmin = MIN(bmin,b);

            if(i-1>=0)
            d->MBX[Im1JK] = MIN(d->MBX[Im1JK],b);
        }

        // +y, -y
        if(p->j_dir==1)
        {
            if(sideC_[IJp1K]*s<0.0 && p->flag4[IJp1K]>0)
            {
                const double b = beta(p->DYP[JP]);
                d->MBY[IJK] = MIN(d->MBY[IJK],b);
                bmin = MIN(bmin,b);
                blocked_.push_back({i,j,k,1,e.tc,{e.w0,e.w1,e.w2}});
            }

            if(sideC_[IJm1K]*s<0.0 && p->flag4[IJm1K]>0)
            {
                const double b = beta(p->DYP[JM1]);
                bmin = MIN(bmin,b);

                if(j-1>=0)
                d->MBY[IJm1K] = MIN(d->MBY[IJm1K],b);
            }
        }

        // vertical link through the cell: nodes k and k+1
        {
            const int qf  = (i-p->imin)*p->jmax*p->kmaxF + (j-p->jmin)*p->kmaxF + k-p->kmin;

            if(sideN_[qf]*sideN_[qf+1]<0)
            {
                const double b = beta(p->DZN[KP]*d->WL(i,j));
                d->MBZ[IJK] = MIN(d->MBZ[IJK],b);
                bmin = MIN(bmin,b);
                blocked_.push_back({i,j,k,2,e.tc,{e.w0,e.w1,e.w2}});
            }
        }

        // cells next to a blocked link: face corrections, Rhie-Chow continuity flux (nhflow_membrane_beta.h)
        d->MBETA[IJK] = MIN(d->MBETA[IJK],bmin);
    }
}

void net_membrane::link_loads(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL,
                              const function<void(int,const double*,const Eigen::Vector3d&)> &add)
{
    // 1. momentum taken out of the layer cells by the forcing in this stage
    for(size_t s=0; s<cells_.size() && s<fimp_.size(); ++s)
    {
        const cellentry &e = cells_[s];
        const double w[3] = {e.w0,e.w1,e.w2};
        add(e.tc,w,fimp_[s]);
    }

    // 2. pressure difference across the blocked links (cell level of the nodes, as in the velocity correction)
    const double *P = d->P;

    for(const auto &b : blocked_)
    {
        i=b.i; j=b.j; k=b.k;

        Eigen::Vector3d f = Eigen::Vector3d::Zero();

        if(b.dir==0)
        {
            const double Pa = 0.5*(P[FIJK]+P[FIJKp1]);
            const double Pb = 0.5*(P[FIp1JK]+P[FIp1JKp1]);
            f(0) = (Pa-Pb)*p->DYN[JP]*p->DZN[KP]*WL(i,j);
        }

        if(b.dir==1)
        {
            const double Pa = 0.5*(P[FIJK]+P[FIJKp1]);
            const double Pb = 0.5*(P[FIJp1K]+P[FIJp1Kp1]);
            f(1) = (Pa-Pb)*p->DXN[IP]*p->DZN[KP]*WL(i,j);
        }

        if(b.dir==2)
        f(2) = (P[FIJK]-P[FIJKp1])*p->DXN[IP]*p->DYN[JP];

        add(b.t,b.w,f);
    }
}
