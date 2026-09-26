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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"ice_contact_nscd.h"
#include<cmath>
#include<algorithm>
#include<numeric>

ice_contact_nscd::ice_contact_nscd(double mu_, double e_, int iter_, int walls_, double xmin_, double xmax_, double ymin_, double ymax_)
              : mu(mu_),e(e_),iter(iter_),wallflag(walls_),xmin(xmin_),xmax(xmax_),ymin(ymin_),ymax(ymax_),
                beta(0.2),slop(0.0),vrest(1.0e-3)
{
    if(iter<1)
    iter=1;
}

void ice_contact_nscd::solve(vector<ice_body2D> &b, double dt, int is2D)
{
    const int nb = int(b.size());

    for(auto &bd : b)
    {
    bd.fx=bd.fy=0.0;
    bd.ncontact=0;
    }

    rec.clear();
    
    if(nb==0 || dt<=0.0)
    return;

    // penetration slop: 1% of the smallest body
    double rmin=1.0e20;
    for(auto &bd : b)
    if(bd.rb>0.0)
    rmin = min(rmin,bd.rb);
    slop = 0.01*rmin;

    // detection
    con.clear();
    broadphase(b,dt);

    if(wallflag==1)
    for(int n=0; n<nb; ++n)
    if(b[n].im>0.0)
    walls(b,n,dt);

    const int nc = int(con.size());

    // prepare, warm start
    map<tuple<int,int,long>,pair<double,double>> newcache;

    for(auto &c : con)
    {
        ice_body2D &A = b[c.a];
        const int wall = (c.b<0);

        double imb=0.0, iIb=0.0, vbx=0.0, vby=0.0, wb=0.0;

        c.rax = c.px - A.x;
        c.ray = c.py - A.y;
        c.rbx = c.rby = 0.0;

        if(!wall)
        {
        ice_body2D &B = b[c.b];
        c.rbx = c.px - B.x;
        c.rby = c.py - B.y;
        imb = B.im; iIb = B.iI;
        vbx = B.vx; vby = B.vy; wb = B.w;
        }

        const double tx = -c.ny, ty = c.nx;
        const double rna = c.rax*c.ny - c.ray*c.nx;
        const double rnb = c.rbx*c.ny - c.rby*c.nx;
        const double rta = c.rax*ty - c.ray*tx;
        const double rtb = c.rbx*ty - c.rby*tx;

        const double kn = A.im + imb + A.iI*rna*rna + iIb*rnb*rnb;
        const double kt = A.im + imb + A.iI*rta*rta + iIb*rtb*rtb;
        c.mn = kn>0.0 ? 1.0/kn : 0.0;
        c.mt = kt>0.0 ? 1.0/kt : 0.0;

        // relative normal velocity before the solve
        const double dvx = (vbx - wb*c.rby) - (A.vx - A.w*c.ray);
        const double dvy = (vby + wb*c.rbx) - (A.vy + A.w*c.rax);
        const double vn0 = dvx*c.nx + dvy*c.ny;

        if(c.sep>0.0)
        c.bias = -c.sep/dt;        // speculative: may close the gap, not more

        else
        {
        double brest = (vn0 < -vrest) ? -e*vn0 : 0.0;
        double bpos  = (-c.sep > slop) ? beta*(-c.sep - slop)/dt : 0.0;
        c.bias = max(brest,bpos);
        }

        // warm start
        auto it = cache.find(make_tuple(b[c.a].id, c.b<0 ? c.b : b[c.b].id, c.key));
        if(it!=cache.end() && c.sep<=0.0)
        {
        c.Pn = it->second.first;
        c.Pt = it->second.second;

        const double Px = c.Pn*c.nx + c.Pt*tx;
        const double Py = c.Pn*c.ny + c.Pt*ty;

        A.vx -= A.im*Px;
        A.vy -= A.im*Py;
        A.w  -= A.iI*(c.rax*Py - c.ray*Px);

        if(!wall)
        {
        ice_body2D &B = b[c.b];
        B.vx += B.im*Px;
        B.vy += B.im*Py;
        B.w  += B.iI*(c.rbx*Py - c.rby*Px);
        }
        }
    }

    // projected Gauss-Seidel on the impulses
    for(int it=0; it<iter; ++it)
    for(auto &c : con)
    {
        ice_body2D &A = b[c.a];
        const int wall = (c.b<0);
        ice_body2D *B = wall ? nullptr : &b[c.b];

        const double tx = -c.ny, ty = c.nx;

        auto relvel = [&](double &dvx, double &dvy)
        {
            double vbx=0.0, vby=0.0;
            if(!wall)
            {
            vbx = B->vx - B->w*c.rby;
            vby = B->vy + B->w*c.rbx;
            }
            dvx = vbx - (A.vx - A.w*c.ray);
            dvy = vby - (A.vy + A.w*c.rax);
        };

        auto apply = [&](double Px, double Py)
        {
            A.vx -= A.im*Px;
            A.vy -= A.im*Py;
            A.w  -= A.iI*(c.rax*Py - c.ray*Px);
            if(!wall)
            {
            B->vx += B->im*Px;
            B->vy += B->im*Py;
            B->w  += B->iI*(c.rbx*Py - c.rby*Px);
            }
        };

        double dvx,dvy;

        // normal: Signorini with restitution / speculative gap
        relvel(dvx,dvy);
        const double vn = dvx*c.nx + dvy*c.ny;
        const double lam = c.mn*(c.bias - vn);
        const double Pn_new = max(c.Pn + lam, 0.0);
        const double dPn = Pn_new - c.Pn;
        c.Pn = Pn_new;
        apply(dPn*c.nx, dPn*c.ny);

        // tangential: Coulomb
        relvel(dvx,dvy);
        const double vt = dvx*tx + dvy*ty;
        const double lamt = -c.mt*vt;
        const double Pt_max = mu*c.Pn;
        const double Pt_new = max(-Pt_max, min(Pt_max, c.Pt + lamt));
        const double dPt = Pt_new - c.Pt;
        c.Pt = Pt_new;
        apply(dPt*tx, dPt*ty);
    }

    // loads and cache
    for(auto &c : con)
    {
        const double tx = -c.ny, ty = c.nx;
        const double Px = c.Pn*c.nx + c.Pt*tx;
        const double Py = c.Pn*c.ny + c.Pt*ty;

        if(c.Pn>0.0)
        {
        b[c.a].fx -= Px/dt;
        b[c.a].fy -= Py/dt;
        ++b[c.a].ncontact;

        if(c.b>=0)
        {
        b[c.b].fx += Px/dt;
        b[c.b].fy += Py/dt;
        ++b[c.b].ncontact;
        }
        }

        if(c.Pn>0.0)
        newcache[make_tuple(b[c.a].id, c.b<0 ? c.b : b[c.b].id, c.key)] = make_pair(c.Pn,c.Pt);
    }

    cache.swap(newcache);
    
    // per pair records: force-weighted contact point, summed normal force
    rec.clear();
    map<pair<int,int>,int> pairidx;
    for(auto &c : con)
    if(c.Pn>0.0)
    {
        auto key = make_pair(c.a,c.b);
        auto it = pairidx.find(key);
        const double Fn = c.Pn/dt;
        
        if(it==pairidx.end())
        {
        ice_contact_record r;
        r.a=c.a; r.b=c.b;
        r.px=c.px*Fn; r.py=c.py*Fn;
        r.nx=c.nx; r.ny=c.ny;
        r.Fn=Fn;
        pairidx[key]=int(rec.size());
        rec.push_back(r);
        }
        else
        {
        ice_contact_record &r = rec[it->second];
        r.px += c.px*Fn;
        r.py += c.py*Fn;
        r.Fn += Fn;
        }
    }
    for(auto &r : rec)
    {
    r.px /= r.Fn;
    r.py /= r.Fn;
    }

    if(is2D==1)
    for(auto &bd : b)
    {
    bd.vy = 0.0;
    bd.w = 0.0;
    }

    (void)nc;
}

void ice_contact_nscd::broadphase(vector<ice_body2D> &b, double dt)
{
    const int nb = int(b.size());
    vector<int> order(nb);
    iota(order.begin(),order.end(),0);
    sort(order.begin(),order.end(),[&](int p, int q){return b[p].x-b[p].rb < b[q].x-b[q].rb;});

    double vmax=0.0;
    for(auto &bd : b)
    vmax = max(vmax, sqrt(bd.vx*bd.vx + bd.vy*bd.vy) + fabs(bd.w)*bd.rb);

    for(int n=0; n<nb; ++n)
    {
        const int ia = order[n];
        const ice_body2D &A = b[ia];
        const double xa_max = A.x + A.rb + 2.0*vmax*dt + slop;

        for(int m=n+1; m<nb; ++m)
        {
            const int ib = order[m];
            const ice_body2D &B = b[ib];

            if(B.x - B.rb > xa_max)
            break;

            if(A.im==0.0 && B.im==0.0)
            continue;

            // speculative margin from the closing speed bound
            const double va = sqrt(A.vx*A.vx + A.vy*A.vy) + fabs(A.w)*A.rb;
            const double vb = sqrt(B.vx*B.vx + B.vy*B.vy) + fabs(B.w)*B.rb;
            const double margin = (va+vb)*dt + slop;

            const double dx = B.x - A.x, dy = B.y - A.y;
            const double rr = A.rb + B.rb + margin;

            if(dx*dx + dy*dy > rr*rr)
            continue;

            // keep a < b so the warm start keys are stable
            if(ia<ib)
            collide(b,ia,ib,margin);
            else
            collide(b,ib,ia,margin);
        }
    }
}

double ice_contact_nscd::max_separation(const ice_body2D &A, const ice_body2D &B, int &edge) const
{
    const int na = int(A.px.size());
    const int nbv = int(B.px.size());
    double best = -1.0e20;
    edge = 0;

    for(int i=0; i<na; ++i)
    {
        const int i2 = (i+1)%na;
        double ex = A.px[i2]-A.px[i];
        double ey = A.py[i2]-A.py[i];
        const double len = sqrt(ex*ex+ey*ey);
        if(len<=0.0)
        continue;
        const double nx = ey/len, ny = -ex/len;   // outward normal, CCW polygon

        double smin = 1.0e20;
        for(int j=0; j<nbv; ++j)
        smin = min(smin, nx*(B.px[j]-A.px[i]) + ny*(B.py[j]-A.py[i]));

        if(smin>best)
        {
        best = smin;
        edge = i;
        }
    }
    return best;
}

void ice_contact_nscd::collide(vector<ice_body2D> &b, int ia, int ib, double margin)
{
    const ice_body2D &A = b[ia];
    const ice_body2D &B = b[ib];

    int ea,eb;
    const double sepA = max_separation(A,B,ea);
    if(sepA>margin)
    return;

    const double sepB = max_separation(B,A,eb);
    if(sepB>margin)
    return;

    // reference face: the one with the larger separation, with a bias for A to avoid flip-flopping
    const double tol = 0.1*slop;
    int flip;
    const ice_body2D *R, *I;
    int er;

    if(sepB > sepA + tol)
    {
    R=&B; I=&A; er=eb; flip=1;
    }
    else
    {
    R=&A; I=&B; er=ea; flip=0;
    }

    const int nr = int(R->px.size());
    const int ni = int(I->px.size());

    const double v1x = R->px[er], v1y = R->py[er];
    const double v2x = R->px[(er+1)%nr], v2y = R->py[(er+1)%nr];
    double tx = v2x-v1x, ty = v2y-v1y;
    const double len = sqrt(tx*tx+ty*ty);
    tx/=len; ty/=len;
    const double nx = ty, ny = -tx;

    // incident edge: most anti-parallel normal
    int ei=0;
    double dmin=1.0e20;
    for(int i=0; i<ni; ++i)
    {
        const int i2=(i+1)%ni;
        double ex = I->px[i2]-I->px[i];
        double ey = I->py[i2]-I->py[i];
        const double l = sqrt(ex*ex+ey*ey);
        if(l<=0.0)
        continue;
        const double d = (ey/l)*nx + (-ex/l)*ny;
        if(d<dmin)
        {
        dmin=d;
        ei=i;
        }
    }

    double cx[2] = {I->px[ei], I->px[(ei+1)%ni]};
    double cy[2] = {I->py[ei], I->py[(ei+1)%ni]};
    int cid[2] = {0,1};

    // clip against the side planes of the reference edge
    auto clip = [&](double mx, double my, double off) -> int
    {
        const double d0 = mx*cx[0] + my*cy[0] - off;
        const double d1 = mx*cx[1] + my*cy[1] - off;

        if(d0<=0.0 && d1<=0.0)
        return 2;

        if(d0>0.0 && d1>0.0)
        return 0;

        const double s = d0/(d0-d1);
        const double xi = cx[0] + s*(cx[1]-cx[0]);
        const double yi = cy[0] + s*(cy[1]-cy[0]);

        if(d0>0.0)
        {
        cx[0]=xi; cy[0]=yi; cid[0]=2;
        }
        else
        {
        cx[1]=xi; cy[1]=yi; cid[1]=3;
        }
        return 2;
    };

    if(clip(-tx,-ty,-(tx*v1x+ty*v1y))<2)
    return;
    if(clip(tx,ty,tx*v2x+ty*v2y)<2)
    return;

    for(int q=0; q<2; ++q)
    {
        const double sep = nx*(cx[q]-v1x) + ny*(cy[q]-v1y);

        if(sep<=margin)
        {
        contact c;
        c.a = ia;
        c.b = ib;
        c.nx = flip ? -nx : nx;
        c.ny = flip ? -ny : ny;
        c.px = cx[q] - 0.5*sep*nx;
        c.py = cy[q] - 0.5*sep*ny;
        c.sep = sep;
        c.key = long(er) + 1000L*long(ei) + 1000000L*long(cid[q]) + 10000000L*long(flip);
        con.push_back(c);
        }
    }
}

void ice_contact_nscd::walls(vector<ice_body2D> &b, int n, double dt)
{
    const ice_body2D &A = b[n];
    const double va = sqrt(A.vx*A.vx + A.vy*A.vy) + fabs(A.w)*A.rb;
    const double margin = va*dt + slop;

    if(A.x - A.rb > xmin + margin && A.x + A.rb < xmax - margin
    && A.y - A.rb > ymin + margin && A.y + A.rb < ymax - margin)
    return;

    const int nv = int(A.px.size());
    const int is2D = (ymax-ymin<=0.0);

    for(int v=0; v<nv; ++v)
    {
        const double wsep[4] = {A.px[v]-xmin, xmax-A.px[v], A.py[v]-ymin, ymax-A.py[v]};
        const double wnx[4] = {-1.0, 1.0, 0.0, 0.0};
        const double wny[4] = { 0.0, 0.0,-1.0, 1.0};

        for(int q=0; q<(is2D?2:4); ++q)
        if(wsep[q]<=margin)
        {
        contact c;
        c.a = n;
        c.b = -1-q;
        c.nx = wnx[q];
        c.ny = wny[q];
        c.px = A.px[v];
        c.py = A.py[v];
        c.sep = wsep[q];
        c.key = long(v);
        con.push_back(c);
        }
    }
}
