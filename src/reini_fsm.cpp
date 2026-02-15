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

#include"reini_fsm.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"field4.h"
#include"increment.h"
#include<algorithm>
#include<cmath>
#include<iostream>
#include<vector>

//  The signed distance of the level set at the start: fast sweeping with closest-point propagation
//  (vector distance transform, Danielsson 1980, the sweep orders of Zhao 2005), as the coastline of
//  FNPF (fnpf_coastline_fsm.cpp), in 3D.
//
//  1. The cells next to the interface (one of the 26 (2D: 8) neighbours of the other sign) keep
//     their distance and get the point of the interface closest to them (the foot point).  With the
//     diagonal neighbours the foot points cover the interface densely (with the face neighbours
//     alone the cells along a diagonal of the interface are left out, and the distance next to it
//     comes from foot points beside it: 0.35 dx too large next to a sphere).  The distance is that of
//     the subcell fix of Russo & Smereka (2000), |f|/D with D the largest of the gradient and the
//     one-sided differences of the face neighbours, and not larger than the distance to the crossing
//     of the zero level set on the line to a neighbour of the other sign (linear between the two
//     cells): a signed distance (the regions of iniphi, F 60, the wave theory of the relaxation
//     zones) stays as it is (D = 1), a step +1/-1 gives the crossing half-way between the cells.  The
//     foot point lies at that distance against the gradient; the crossing itself is the distance and
//     the foot point if the subcell fix reaches beyond it, or if the gradient does not point to its
//     side (within 60 degrees: inside a union of overlapping regions of iniphi a face of one region
//     within the other is no interface, and the gradient there points to that face).
//  2. Gauss-Seidel sweeps in the 8 (2D: 4) orders of i, j, k: a cell takes the foot point of one of
//     its 26 (2D: 8) neighbours if it is closer than its own.  The foot points of the ghost cells
//     across the partition edges are exchanged after every sweep, the sweeps are repeated until no
//     cell on any rank changes.
//  3. f = sign(f) |x - foot point|.
//
//  The cells next to the interface keep their value: the zero level set does not move.  The
//  sweeping runs over the cells of LOOP (flag4 > 0), as the cells that took part in the PDE of the
//  start before: the solids of the direct forcing (G 505 0) and the 6DOF bodies are such cells, the
//  free surface of the start continues through them and is passed on (levelset_solid_extend works
//  on them in the time steps).  Masked out, the air above a dry beach of a direct-forcing wedge
//  took the distance to the shore line instead of that to the still water level through the wedge,
//  and still water drove 0.025 m/s along the beach (cfd_2d_beach_still, 3e-4 m/s with the level
//  set of the start).  The solid cells (flag4 <= 0: topography, solids of the ghost-cell method)
//  keep their value and pass nothing on, the ghost cells beyond the edges of the domain are not
//  used (a periodic edge is not continued across).  Cells that no foot point reaches (no interface
//  in the domain, or in a part of the fluid closed off by solids) keep their value.
namespace
{
const double nofoot = 1.0e30;
}

int reini_fsm(lexer *p, fdm *a, ghostcell *pgc, field &f)
{
    int i,j,k;
    const int marge = increment::marge;
    const bool jd = p->j_dir==1;
    const int si = p->jmax*p->kmax, sj = p->kmax;
    const int n = p->imax*p->jmax*p->kmax;

    auto id = [&](int ii, int jj, int kk) {return (ii-p->imin)*si + (jj-p->jmin)*sj + (kk-p->kmin);};

    // ------------------------------------------------------------
    // 0. the cells of the sweeping (m = 1): the cells of LOOP, on the partition and on the ghost
    //    cells across a partition edge (exchange); 0 elsewhere
    // ------------------------------------------------------------
    field4 m(p), cx(p), cy(p), cz(p);

    for(int q=0; q<n; ++q)
    {
    m.data()[q] = 0.0;
    cx.data()[q] = cy.data()[q] = cz.data()[q] = nofoot;
    }

    LOOP
    m(i,j,k) = 1.0;

    pgc->gcparax(p,m,4);

    std::vector<double> d2(n,nofoot);    // squared distance to the foot point
    std::vector<char> frozen(n,0);

    // ------------------------------------------------------------
    // 1. the cells next to the interface
    // ------------------------------------------------------------
    int nseed = 0;

    LOOP
    if(m(i,j,k)>0.5)
    {
        const double fc = f(i,j,k);
        const bool water = fc>=0.0;
        const double xc[3] = {p->XP[IP], p->YP[JP], p->ZP[KP]};

        double grad[3] = {0.0,0.0,0.0};
        double dmax = 0.0;                  // largest one-sided difference
        double dcross = nofoot;             // distance to the nearest crossing
        double xcross[3] = {xc[0],xc[1],xc[2]};
        bool iface = false;

        // the crossings of the zero level set on the lines to the 26 (2D: 8) neighbours
        for(int oi=-1; oi<=1; ++oi)
        for(int oj=(jd?-1:0); oj<=(jd?1:0); ++oj)
        for(int ok=-1; ok<=1; ++ok)
        {
            if((oi==0 && oj==0 && ok==0) || m(i+oi,j+oj,k+ok)<0.5)
            continue;

            const double fn = f(i+oi,j+oj,k+ok);

            if((fn>=0.0)==water)
            continue;

            iface = true;

            const double xn[3] = {p->XP[IP+oi], p->YP[JP+oj], p->ZP[KP+ok]};
            const double t = fc/(fc-fn);
            const double dist = t*sqrt((xn[0]-xc[0])*(xn[0]-xc[0]) + (xn[1]-xc[1])*(xn[1]-xc[1]) + (xn[2]-xc[2])*(xn[2]-xc[2]));

            if(dist<dcross)
            {
                dcross = dist;

                for(int ax=0; ax<3; ++ax)
                xcross[ax] = xc[ax] + t*(xn[ax]-xc[ax]);
            }
        }

        if(!iface)
        continue;

        // gradient and one-sided differences from the face neighbours
        for(int ax=0; ax<3; ++ax)
        {
            if(ax==1 && !jd)
            continue;

            double slope[2] = {0.0,0.0};
            double xn[2] = {0.0,0.0};
            double fn[2] = {0.0,0.0};
            bool valid[2] = {false,false};

            for(int s=0; s<2; ++s)
            {
                const int o = s==0 ? -1 : 1;
                const int in = i + (ax==0 ? o : 0);
                const int jn = j + (ax==1 ? o : 0);
                const int kn = k + (ax==2 ? o : 0);

                if(m(in,jn,kn)<0.5)
                continue;

                valid[s] = true;
                fn[s] = f(in,jn,kn);
                xn[s] = ax==0 ? p->XP[in+marge] : (ax==1 ? p->YP[jn+marge] : p->ZP[kn+marge]);
                slope[s] = (fn[s]-fc)/(xn[s]-xc[ax]);
                dmax = std::max(dmax,fabs(slope[s]));
            }

            if(valid[0] && valid[1])
            grad[ax] = (fn[1]-fn[0])/(xn[1]-xn[0]);
            else if(valid[0])
            grad[ax] = slope[0];
            else if(valid[1])
            grad[ax] = slope[1];
        }

        const double g = sqrt(grad[0]*grad[0] + grad[1]*grad[1] + grad[2]*grad[2]);
        const double D = std::max(g,dmax);

        // the distance: the subcell fix if it lies before the nearest crossing, the foot point
        // against the gradient if that points to the side of the crossing (within 60 degrees);
        // otherwise the nearest crossing itself.  Inside a union of overlapping regions (iniphi) a
        // face of one region within the other is no interface, the level set there not a distance,
        // and its gradient points to that face: the crossing is the interface
        double d0 = dcross;
        double fp[3] = {xcross[0],xcross[1],xcross[2]};

        if(D>0.0 && fabs(fc)/D<dcross && g>1.0e-12*D)
        {
            const double sg = water ? 1.0 : -1.0;
            const double dn = fabs(fc)/D;
            double tp[3], dot = 0.0;

            for(int ax=0; ax<3; ++ax)
            {
            tp[ax] = xc[ax] - sg*dn*grad[ax]/g;
            dot += (tp[ax]-xc[ax])*(xcross[ax]-xc[ax]);
            }

            if(dot>0.5*dn*dcross)
            {
                d0 = dn;

                for(int ax=0; ax<3; ++ax)
                fp[ax] = tp[ax];
            }
        }

        const int q = id(i,j,k);
        d2[q] = d0*d0;
        frozen[q] = 1;
        cx.data()[q] = fp[0];
        cy.data()[q] = fp[1];
        cz.data()[q] = fp[2];
        ++nseed;
    }

    nseed = pgc->globalisum(nseed);

    if(nseed==0)
    {
        if(p->mpirank==0)
        cout<<"level set fast sweeping: no interface"<<endl;

        return 0;
    }

    pgc->gcparax(p,cx,4);
    pgc->gcparax(p,cy,4);
    pgc->gcparax(p,cz,4);

    // ------------------------------------------------------------
    // 2. sweeping
    // ------------------------------------------------------------
    const int nsweep = jd ? 8 : 4;
    const int maxiter = 1000;
    const int oj0 = jd ? -1 : 0, oj1 = jd ? 1 : 0;
    int iter;

    for(iter=0; iter<maxiter; ++iter)
    {
        int changed = 0;

        for(int sweep=0; sweep<nsweep; ++sweep)
        {
            const bool ri = (sweep&1)!=0;
            const bool rk = (sweep&2)!=0;
            const bool rj = (sweep&4)!=0;

            for(int ii=0; ii<p->knox; ++ii)
            {
            i = ri ? p->knox-1-ii : ii;

                for(int jj=0; jj<p->knoy; ++jj)
                {
                j = rj ? p->knoy-1-jj : jj;

                    for(int kk=0; kk<p->knoz; ++kk)
                    {
                    k = rk ? p->knoz-1-kk : kk;

                    const int q = id(i,j,k);

                    if(m.data()[q]<0.5 || frozen[q])
                    continue;

                    const double x = p->XP[IP];
                    const double y = p->YP[JP];
                    const double z = p->ZP[KP];

                    double best = d2[q];
                    int qbest = -1;

                    for(int oi=-1; oi<=1; ++oi)
                    for(int oj=oj0; oj<=oj1; ++oj)
                    for(int ok=-1; ok<=1; ++ok)
                    {
                        const int qn = q + oi*si + oj*sj + ok;

                        if(qn==q || m.data()[qn]<0.5 || cx.data()[qn]>=0.5*nofoot)
                        continue;

                        const double ddx = x-cx.data()[qn];
                        const double ddy = y-cy.data()[qn];
                        const double ddz = z-cz.data()[qn];
                        const double dd = ddx*ddx + ddy*ddy + ddz*ddz;

                        if(dd<best)
                        {
                            best = dd;
                            qbest = qn;
                        }
                    }

                    if(qbest>=0)
                    {
                        d2[q] = best;
                        cx.data()[q] = cx.data()[qbest];
                        cy.data()[q] = cy.data()[qbest];
                        cz.data()[q] = cz.data()[qbest];
                        ++changed;
                    }
                    }
                }
            }

            pgc->gcparax(p,cx,4);
            pgc->gcparax(p,cy,4);
            pgc->gcparax(p,cz,4);
        }

        changed = pgc->globalisum(changed);

        if(changed==0)
        break;
    }

    if(p->mpirank==0)
    cout<<"level set fast sweeping: "<<nseed<<" interface cells, "<<MIN(iter+1,maxiter)<<" iterations"<<endl;

    if(p->mpirank==0 && iter==maxiter)
    cout<<"WARNING: level set fast sweeping not converged"<<endl;

    // ------------------------------------------------------------
    // 3. signed distance
    // ------------------------------------------------------------
    LOOP
    {
        const int q = id(i,j,k);

        if(m.data()[q]>0.5 && d2[q]<0.5*nofoot)
        f(i,j,k) = (f(i,j,k)>=0.0 ? 1.0 : -1.0)*sqrt(d2[q]);
    }

    return MIN(iter+1,maxiter);
}
