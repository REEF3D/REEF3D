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

#include"dem_core.h"
#include<cmath>
#include<algorithm>
#include<iostream>

namespace
{
    inline dem_mat skew(const dem_vec &r)
    {
        dem_mat S;
        S <<    0.0, -r(2),  r(1),
               r(2),   0.0, -r(0),
              -r(1),  r(0),   0.0;
        return S;
    }

    inline uint64_t mix64(uint64_t z)
    {
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31);
    }
}

uint64_t dem_core::make_key(int a, int b, int feature) const
{
    uint64_t h = mix64(uint64_t(a) + 0x9E3779B97F4A7C15ULL);
    h = mix64(h ^ (uint64_t(int64_t(b)+1024) * 0xC2B2AE3D27D4EB4FULL));
    h = mix64(h ^ (uint64_t(feature) * 0x165667B19E3779F9ULL));
    return h;
}

void dem_core::initialize()
{
    for(auto &B : bodies)
    {
        const dem_shape &S = shapes[B.shape];
        B.m  = mats[B.mat].rho*S.volume;
        B.Ib = mats[B.mat].rho*S.inertia_unit;
        B.q.normalize();
        B.R = B.q.toRotationMatrix();
        B.vold = B.v;
        B.wold = B.w;
    }
}

double dem_core::kinetic_energy() const
{
    double E=0.0;
    for(auto &B : bodies)
    if(B.active && !B.fixed && !B.ghost)
    {
        dem_vec wb = B.R.transpose()*B.w;
        E += 0.5*B.m*B.v.squaredNorm() + 0.5*wb.dot(B.Ib.cwiseProduct(wb));
    }
    return E;
}

double dem_core::max_velocity() const
{
    double vmax=0.0;
    for(auto &B : bodies)
    if(B.active && !B.fixed)
    vmax = std::max(vmax, B.v.norm() + shapes[B.shape].rbound*B.w.norm());
    return vmax;
}

double dem_core::min_rbound() const
{
    double r=1.0e20;
    for(auto &B : bodies)
    if(B.active)
    r = std::min(r,shapes[B.shape].rbound);
    return r;
}

// ---------------------------------------------------------------------------------------------
// time step
// ---------------------------------------------------------------------------------------------

void dem_core::step(double dt, const wallfunc &walls, const dem_hooks *hooks)
{
    hk = hooks;
    update_mass(dt);

    // contact margin from the pre-step velocities
    double vm = max_velocity();
    if(hk)
    vm = std::max(vm,vmax_ext);     // distributed: particles arriving from other ranks
    double margin = 1.5*vm*dt + gravity.norm()*dt*dt + 2.0*slop;

    // free velocities: gravity, external loads, implicit drag and added mass, implicit gyroscopic term
    for(auto &B : bodies)
    {
        if(B.ghost)
        continue;

        B.vold = B.v;
        B.wold = B.w;

        if(!B.active || B.fixed)
        continue;

        double meff = B.m + B.madd + dt*B.K;
        B.v = ((B.m + B.madd)*B.v + dt*(B.F + B.m*gravity + B.K*B.uf + B.madd*B.af))/meff;

        dem_vec wg = gyro_free(B,dt);
        dem_mat Iw = B.R*B.Ib.asDiagonal()*B.R.transpose();
        B.w = B.invI*(Iw*wg + dt*B.T);

        if(plane2D)
        {
            B.v(1) = 0.0;
            B.w(0) = 0.0;
            B.w(2) = 0.0;
        }
    }

    // ghosts receive the free velocities, masses and positions of their owners
    if(hk && hk->ghosts)
    hk->ghosts(*this);

    contacts.clear();
    detect(margin);
    planes_contacts(margin);
    if(walls)
    walls(*this,margin,contacts);

    reduce();

    // mass splitting for particles whose contacts are solved on several ranks
    if(hk && hk->split)
    hk->split(*this);

    prepare(dt);
    solve(dt);
    integrate(dt);

    hk = nullptr;
}

void dem_core::update_mass(double dt)
{
    for(auto &B : bodies)
    {
        if(B.ghost)
        continue;

        B.R = B.q.toRotationMatrix();

        if(!B.active || B.fixed)
        {
            B.invm.setZero();
            B.invI.setZero();
            continue;
        }

        double meff = B.m + B.madd + dt*B.K;
        B.invm = dem_vec::Constant(1.0/meff);

        dem_mat Iw = B.R*B.Ib.asDiagonal()*B.R.transpose() + dt*B.Kr*dem_mat::Identity();
        B.invI = Iw.inverse();

        if(plane2D)
        {
            B.invm(1) = 0.0;
            B.invI.setZero();
            B.invI(1,1) = 1.0/Iw(1,1);
        }
    }
}

dem_vec dem_core::gyro_free(const dem_body &B, double dt) const
{
    // implicit gyroscopic term, one Newton step (Catto 2015)
    dem_vec wb = B.R.transpose()*B.w;
    dem_mat Ib = B.Ib.asDiagonal();
    dem_vec f = dt*wb.cross(Ib*wb);
    dem_mat J = Ib + dt*(skew(wb)*Ib - skew(Ib*wb));
    wb = wb - J.inverse()*f;
    return B.R*wb;
}

// ---------------------------------------------------------------------------------------------
// contact detection
// ---------------------------------------------------------------------------------------------

void dem_core::detect(double margin)
{
    int N = bodies.size();
    if(N<2)
    return;

    double rmax=0.0;
    for(auto &B : bodies)
    if(B.active)
    rmax = std::max(rmax,shapes[B.shape].rbound);

    auto test = [&](int a, int b)
    {
        const dem_body &A = bodies[a], &Bb = bodies[b];
        if(!A.active || !Bb.active)
        return;
        if(A.fixed && Bb.fixed)
        return;
        if(!mine(a,b))
        return;
        double rs = shapes[A.shape].rbound + shapes[Bb.shape].rbound + margin;
        if((A.x-Bb.x).squaredNorm() < rs*rs)
        narrow(a,b,margin);
    };

    if(N<=64)
    {
        for(int a=0; a<N; ++a)
        for(int b=a+1; b<N; ++b)
        test(a,b);
        return;
    }

    // uniform hash grid
    double h = 2.0*rmax + margin;
    unordered_map<uint64_t,vector<int>> cells;
    auto ckey = [](long long ix, long long iy, long long iz)
    {
        return uint64_t((ix & 0x1FFFFF) | ((iy & 0x1FFFFF)<<21) | ((iz & 0x1FFFFF)<<42));
    };

    vector<long long> cx(N),cy(N),cz(N);
    for(int a=0; a<N; ++a)
    {
        cx[a] = (long long)floor(bodies[a].x(0)/h);
        cy[a] = (long long)floor(bodies[a].x(1)/h);
        cz[a] = (long long)floor(bodies[a].x(2)/h);
        cells[ckey(cx[a],cy[a],cz[a])].push_back(a);
    }

    for(int a=0; a<N; ++a)
    for(int di=-1; di<=1; ++di)
    for(int dj=-1; dj<=1; ++dj)
    for(int dk=-1; dk<=1; ++dk)
    {
        auto it = cells.find(ckey(cx[a]+di,cy[a]+dj,cz[a]+dk));
        if(it==cells.end())
        continue;
        for(int b : it->second)
        if(b>a)
        test(a,b);
    }
}

void dem_core::narrow(int a, int b, double margin)
{
    const dem_body &A = bodies[a], &B = bodies[b];
    const dem_shape &SA = shapes[A.shape], &SB = shapes[B.shape];

    double mu = 0.5*(mats[A.mat].friction + mats[B.mat].friction);
    double e  = 0.5*(mats[A.mat].restitution + mats[B.mat].restitution);

    auto add = [&](int ca, int cb, int feature, const dem_vec &x, const dem_vec &n, double gap)
    {
        dem_contact C;
        C.a = ca;
        C.b = cb;
        C.key = make_key(bodies[ca].id,bodies[cb].id,feature);
        C.x = x;
        C.n = n;
        C.gap = gap;
        C.mu = mu;
        C.e = e;
        contacts.push_back(C);
    };

    // sphere - sphere
    if(SA.type==DEM_SPHERE && SB.type==DEM_SPHERE)
    {
        dem_vec d = A.x - B.x;
        double dist = d.norm();
        double gap = dist - SA.dim(0) - SB.dim(0);
        if(gap<margin && dist>1.0e-14)
        {
            dem_vec n = d/dist;
            add(a,b,0,B.x + n*(SB.dim(0)+0.5*gap),n,gap);
        }
        return;
    }

    // sphere - general: sphere centre against the level set of the other body
    if(SA.type==DEM_SPHERE || SB.type==DEM_SPHERE)
    {
        int s = SA.type==DEM_SPHERE ? a : b;
        int g = SA.type==DEM_SPHERE ? b : a;
        const dem_body &Bs = bodies[s], &Bg = bodies[g];
        const dem_shape &Ss = shapes[Bs.shape], &Sg = shapes[Bg.shape];
        double r = Ss.dim(0);

        dem_vec y = Bg.R.transpose()*(Bs.x - Bg.x);
        if(y.norm() > Sg.rbound + r + margin)
        return;

        double gap = Sg.sdf(y) - r;
        if(gap<margin)
        {
            dem_vec n = Bg.R*Sg.sdf_grad(y);
            add(s,g,0,Bs.x - n*(r+0.5*gap),n,gap);
        }
        return;
    }

    // general - general: surface nodes of one body against the level set of the other, both ways
    for(int side=0; side<2; ++side)
    {
        int na = side==0 ? a : b;     // node body, receives +P
        int nb = side==0 ? b : a;     // level-set body
        const dem_body &Bn = bodies[na], &Bl = bodies[nb];
        const dem_shape &Sn = shapes[Bn.shape], &Sl = shapes[Bl.shape];

        dem_mat Rrel = Bl.R.transpose()*Bn.R;
        dem_vec trel = Bl.R.transpose()*(Bn.x - Bl.x);
        double rlim = Sl.rbound + margin;

        for(size_t q=0; q<Sn.nodes.size(); ++q)
        {
            dem_vec y = trel + Rrel*Sn.nodes[q];
            if(y.squaredNorm() > rlim*rlim)
            continue;

            double phi = Sl.sdf(y);
            if(phi<margin)
            {
                dem_vec n = Bl.R*Sl.sdf_grad(y);
                dem_vec xw = Bn.x + Bn.R*Sn.nodes[q];
                add(na,nb,int(q) | (side<<23),xw - 0.5*phi*n,n,phi);
            }
        }
    }
}

void dem_core::planes_contacts(double margin)
{
    for(size_t pl=0; pl<planes.size(); ++pl)
    {
        const dem_plane &P = planes[pl];

        for(size_t nb=0; nb<bodies.size(); ++nb)
        {
            const dem_body &B = bodies[nb];
            if(!B.active || B.fixed || B.ghost)
            continue;
            if(B.tier==1 && myrank!=0)
            continue;

            const dem_shape &S = shapes[B.shape];
            double dc = P.n.dot(B.x) - P.d;
            if(dc > S.rbound + margin)
            continue;

            double mu = 0.5*(mats[B.mat].friction + wallmat.friction);
            double e  = 0.5*(mats[B.mat].restitution + wallmat.restitution);

            auto add = [&](int feature, const dem_vec &x, double gap)
            {
                dem_contact C;
                C.a = nb;
                C.b = -2-int(pl);
                C.key = make_key(B.id,C.b,feature);
                C.x = x;
                C.n = P.n;
                C.gap = gap;
                C.mu = mu;
                C.e = e;
                contacts.push_back(C);
            };

            if(S.type==DEM_SPHERE)
            {
                double gap = dc - S.dim(0);
                if(gap<margin)
                add(0,B.x - P.n*(S.dim(0)+0.5*gap),gap);
                continue;
            }

            for(size_t q=0; q<S.nodes.size(); ++q)
            {
                dem_vec xw = B.x + B.R*S.nodes[q];
                double gap = P.n.dot(xw) - P.d;
                if(gap<margin)
                add(q,xw - 0.5*gap*P.n,gap);
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------
// non-smooth contact dynamics: Moreau-Jean with projected Gauss-Seidel on the Coulomb cone
// ---------------------------------------------------------------------------------------------

dem_vec dem_core::relvel(const dem_contact &C) const
{
    const dem_body &A = bodies[C.a];
    dem_vec u = A.v + A.w.cross(C.ra);

    if(C.b>=0)
    {
        const dem_body &B = bodies[C.b];
        u -= B.v + B.w.cross(C.rb);
    }
    else
    u -= C.vwall;

    return u;
}

void dem_core::apply(int c, const dem_vec &Pw)
{
    dem_contact &C = contacts[c];
    dem_body &A = bodies[C.a];
    A.v += A.invm.cwiseProduct(Pw);
    A.w += A.invI*(C.ra.cross(Pw));

    if(C.b>=0)
    {
        dem_body &B = bodies[C.b];
        B.v -= B.invm.cwiseProduct(Pw);
        B.w -= B.invI*(C.rb.cross(Pw));
    }
}

void dem_core::reduce()
{
    // contact manifold reduction: node-based detection gives many redundant points per patch.
    // Per body pair (or body-wall), contacts are clustered by normal and each cluster is reduced to
    // the deepest point plus farthest-point samples. This keeps the patch extent (and so the
    // moment arm) but removes most of the static indeterminacy that slows projected Gauss-Seidel.
    if(manifold<=0 || contacts.size()<2)
    return;

    std::stable_sort(contacts.begin(),contacts.end(),[](const dem_contact &c1, const dem_contact &c2)
    {
        if(c1.a!=c2.a) return c1.a<c2.a;
        return c1.b<c2.b;
    });

    vector<dem_contact> out;
    out.reserve(contacts.size());

    size_t s=0;
    while(s<contacts.size())
    {
        size_t e=s;
        while(e<contacts.size() && contacts[e].a==contacts[s].a && contacts[e].b==contacts[s].b)
        ++e;

        if(int(e-s)<=manifold)
        {
            for(size_t c=s; c<e; ++c)
            out.push_back(contacts[c]);
            s=e;
            continue;
        }

        // cluster by normal
        vector<int> cl(e-s,-1);
        vector<dem_vec> cn;
        for(size_t c=s; c<e; ++c)
        {
            int found=-1;
            for(size_t q=0; q<cn.size(); ++q)
            if(cn[q].dot(contacts[c].n)>0.95)
            {
                found=q;
                break;
            }
            if(found<0)
            {
                found=cn.size();
                cn.push_back(contacts[c].n);
            }
            cl[c-s]=found;
        }

        for(size_t q=0; q<cn.size(); ++q)
        {
            vector<size_t> mem;
            for(size_t c=s; c<e; ++c)
            if(cl[c-s]==int(q))
            mem.push_back(c);

            if(int(mem.size())<=manifold)
            {
                for(size_t c : mem)
                out.push_back(contacts[c]);
                continue;
            }

            // farthest point sampling started from the point farthest from the patch centroid
            // (geometric, hence stable between steps for warm starting), plus the deepest point
            dem_vec cen = dem_vec::Zero();
            for(size_t c : mem)
            cen += contacts[c].x;
            cen /= double(mem.size());

            vector<size_t> sel;
            size_t first = mem[0];
            double dfar=-1.0;
            for(size_t c : mem)
            {
                double d = (contacts[c].x-cen).squaredNorm();
                if(d>dfar+1.0e-14*(1.0+d))
                {
                    dfar=d;
                    first=c;
                }
            }
            sel.push_back(first);

            vector<double> dmin(mem.size(),1.0e30);
            while(int(sel.size())<manifold)
            {
                size_t last = sel.back();
                double best=-1.0;
                size_t bi=mem[0];
                for(size_t m=0; m<mem.size(); ++m)
                {
                    dmin[m] = std::min(dmin[m],(contacts[mem[m]].x - contacts[last].x).squaredNorm());
                    if(dmin[m]>best+1.0e-14*(1.0+dmin[m]))
                    {
                        best=dmin[m];
                        bi=mem[m];
                    }
                }
                if(best<=0.0)
                break;
                sel.push_back(bi);
            }

            size_t deep = mem[0];
            for(size_t c : mem)
            if(contacts[c].gap<contacts[deep].gap)
            deep=c;
            if(std::find(sel.begin(),sel.end(),deep)==sel.end())
            sel.push_back(deep);

            for(size_t c : sel)
            out.push_back(contacts[c]);
        }
        s=e;
    }

    contacts.swap(out);
}

void dem_core::prepare(double dt)
{
    maxpen = 0.0;

    // drop contacts that cannot carry impulse (both sides immovable)
    vector<dem_contact> keep;
    keep.reserve(contacts.size());

    for(auto &C : contacts)
    {
        const dem_body &A = bodies[C.a];
        bool amov = A.active && !A.fixed;
        bool bmov = C.b>=0 && bodies[C.b].active && !bodies[C.b].fixed;
        if(amov || bmov)
        keep.push_back(C);
    }
    contacts.swap(keep);


    for(size_t c=0; c<contacts.size(); ++c)
    {
        dem_contact &C = contacts[c];
        const dem_body &A = bodies[C.a];

        C.ra = C.x - A.x;
        C.rb = C.b>=0 ? dem_vec(C.x - bodies[C.b].x) : dem_vec(dem_vec::Zero());

        // tangent basis
        dem_vec aux = fabs(C.n(0))<0.57 ? dem_vec(dem_vec::UnitX()) : dem_vec(dem_vec::UnitY());
        C.t1 = C.n.cross(aux).normalized();
        C.t2 = C.n.cross(C.t1);

        dem_mat W = dem_mat(A.invm.asDiagonal()) - skew(C.ra)*A.invI*skew(C.ra);
        if(C.b>=0)
        {
            const dem_body &B = bodies[C.b];
            W += dem_mat(B.invm.asDiagonal()) - skew(C.rb)*B.invI*skew(C.rb);
        }

        dem_mat Cm;
        Cm.row(0) = C.n.transpose();
        Cm.row(1) = C.t1.transpose();
        Cm.row(2) = C.t2.transpose();
        C.Wl = Cm*W*Cm.transpose();

        // pre-step normal velocity
        dem_vec u0 = A.vold + A.wold.cross(C.ra);
        if(C.b>=0)
        u0 -= bodies[C.b].vold + bodies[C.b].wold.cross(C.rb);
        else
        u0 -= C.vwall;

        C.un0 = C.n.dot(u0);

        double uimp = C.un0;
        auto itu = cacheU.find(C.key);
        if(itu!=cacheU.end())
        uimp = std::min(uimp,itu->second);

        double pred = C.gap + theta*dt*C.un0;

        if(C.gap<=0.0 || pred<=0.0)
        {
            // active: Newton impact law plus penetration stabilisation
            C.speculative = false;
            C.target = (C.e>0.0 && uimp<-vrest) ? -C.e*uimp : -std::max(C.gap,0.0)/dt;

            // penetration is removed by split impulses, which change positions but not velocities
            C.stab = std::min(beta*std::max(0.0,-(C.gap+slop))/dt, vstab_max);
        }
        else
        {
            // speculative: no penetration within this step
            C.speculative = true;
            C.target = -C.gap/dt;
            C.stab = 0.0;
        }

        maxpen = std::max(maxpen,-C.gap);

        // warm start
        C.P.setZero();
        auto itp = cacheP.find(C.key);
        if(itp!=cacheP.end() && !C.speculative)
        {
            dem_vec Pl = Cm*(warmstart*itp->second);
            Pl(0) = std::max(Pl(0),0.0);
            double pt = sqrt(Pl(1)*Pl(1)+Pl(2)*Pl(2));
            if(pt > C.mu*Pl(0) && pt>0.0)
            {
                Pl(1) *= C.mu*Pl(0)/pt;
                Pl(2) *= C.mu*Pl(0)/pt;
            }
            C.P = Pl;
            apply(c,C.n*Pl(0) + C.t1*Pl(1) + C.t2*Pl(2));
        }
    }

    ncontacts = contacts.size();
}

double dem_core::sweep(const vector<int> &idx, bool reverse)
{
    double res = 0.0;
    int nc = idx.size();

    for(int cc=0; cc<nc; ++cc)
    {
        int c = reverse ? idx[nc-1-cc] : idx[cc];
        dem_contact &C = contacts[c];
        dem_vec u = relvel(C);
        double un = C.n.dot(u);
        Eigen::Vector2d ut(C.t1.dot(u),C.t2.dot(u));

        // normal
        double Wnn = C.Wl(0,0);
        if(Wnn<=1.0e-20)
        continue;

        double Pn = std::max(0.0, C.P(0) + sor*(C.target - un)/Wnn);
        double dPn = Pn - C.P(0);

        // tangential, projected on the friction disc
        ut += C.Wl.block<2,1>(1,0)*dPn;
        Eigen::Matrix2d Wt = C.Wl.block<2,2>(1,1);
        Eigen::Vector2d Pt = C.P.tail<2>();
        Eigen::Vector2d Ptn = Pt;

        // isotropic tangential effective mass: the fixed point then satisfies maximum dissipation
        // (friction impulse anti-parallel to the slip velocity); a full 2x2 solve would not
        double wt = std::max(Wt(0,0),Wt(1,1));
        if(wt>1.0e-20)
        Ptn = Pt - sor*ut/wt;

        double pt = Ptn.norm();
        double lim = C.mu*Pn;
        if(pt>lim && pt>0.0)
        Ptn *= lim/pt;

        Eigen::Vector2d dPt = Ptn - Pt;

        dem_vec dP = C.n*dPn + C.t1*dPt(0) + C.t2*dPt(1);
        apply(c,dP);

        C.P(0) = Pn;
        C.P(1) = Ptn(0);
        C.P(2) = Ptn(1);

        res = std::max(res, fabs(dPn)*Wnn + dPt.norm()*std::max(Wt(0,0),Wt(1,1)));
    }
    return res;
}

double dem_core::sweep_pseudo(const vector<int> &idx, bool reverse)
{
    double res=0.0;
    int nc = idx.size();

    for(int cc=0; cc<nc; ++cc)
    {
        int c = reverse ? idx[nc-1-cc] : idx[cc];
        dem_contact &C = contacts[c];
        double Wnn = C.Wl(0,0);
        if(Wnn<=1.0e-20)
        continue;

        const dem_body &A = bodies[C.a];
        dem_vec u = A.vp + A.wp.cross(C.ra);
        if(C.b>=0)
        u -= bodies[C.b].vp + bodies[C.b].wp.cross(C.rb);

        double Pn = std::max(0.0, C.Pp + (C.stab - C.n.dot(u))/Wnn);
        double dPn = Pn - C.Pp;
        C.Pp = Pn;
        apply_pseudo(c,dPn);
        res = std::max(res,fabs(dPn)*Wnn);
    }
    return res;
}

void dem_core::solve(double dt)
{
    iterations = 0;
    residual = 0.0;

    bool dist = hk && hk->sync;

    for(auto &B : bodies)
    {
        B.vp.setZero();
        B.wp.setZero();
        B.vps.setZero();
        B.wps.setZero();
    }

    if(contacts.empty() && !dist)
    return;

    // contacts of particles solved on one rank only are swept freely, contacts of shared particles
    // in the colour phases (distributed runs only)
    vector<int> freec, sharedc;
    for(size_t c=0; c<contacts.size(); ++c)
    {
        const dem_contact &C = contacts[c];
        bool shared = dist && (bodies[C.a].split>1.5 || (C.b>=0 && bodies[C.b].split>1.5));
        if(shared)
        sharedc.push_back(c);
        else
        freec.push_back(c);
    }

    int ncol = dist ? std::max(1,hk->ncolors) : 1;

    // residual scale: current velocities with a floor of 1 mm/s (resting stacks need a strict tolerance)
    double vref = vref_ext>0.0 ? vref_ext : std::max(max_velocity(),1.0e-3);

    // distributed: make the warm-start impulses of all ranks visible before the first sweep
    if(dist)
    hk->sync(*this,0,0.0,false);

    for(int it=0; it<maxiter; ++it)
    {
        bool rev = (it%2==1);   // symmetric Gauss-Seidel: alternate the order

        if(!dist)
        {
            double res = sweep(freec,rev);
            iterations = it+1;
            residual = res/vref;
            if(residual<tol)
            break;
            continue;
        }

        double res = 0.0, gres;
        if(!rev)
        res = std::max(res,sweep(freec,false));

        for(int cc=0; cc<ncol; ++cc)
        {
            int col = rev ? ncol-1-cc : cc;
            if(col==hk->mycolor)
            res = std::max(res,sweep(sharedc,rev));
            if(cc<ncol-1)
            hk->sync(*this,0,res,false);
        }

        if(rev)
        res = std::max(res,sweep(freec,true));

        gres = hk->sync(*this,0,res,true);

        iterations = it+1;
        residual = gres/vref;
        if(residual<tol)
        break;
    }

    // split impulse: normal-only projected Gauss-Seidel on the pseudo velocities
    bool anypen=false;
    for(auto &C : contacts)
    {
        C.Pp = 0.0;
        if(C.stab>0.0)
        anypen=true;
    }

    if(!anypen && !dist)
    return;

    for(int it=0; it<pseudoiter; ++it)
    {
        bool rev = (it%2==1);

        if(!dist)
        {
            double res = sweep_pseudo(freec,rev);
            if(res<1.0e-3*vstab_max)
            break;
            continue;
        }

        double res = sweep_pseudo(freec,rev);
        for(int cc=0; cc<ncol; ++cc)
        {
            int col = rev ? ncol-1-cc : cc;
            if(col==hk->mycolor)
            res = std::max(res,sweep_pseudo(sharedc,rev));
            if(cc<ncol-1)
            hk->sync(*this,1,res,false);
        }

        double g = hk->sync(*this,1,res,true);
        if(g<1.0e-3*vstab_max)
        break;
    }
}

void dem_core::apply_pseudo(int c, double Pn)
{
    dem_contact &C = contacts[c];
    dem_vec Pw = C.n*Pn;
    dem_body &A = bodies[C.a];
    A.vp += A.invm.cwiseProduct(Pw);
    A.wp += A.invI*(C.ra.cross(Pw));

    if(C.b>=0)
    {
        dem_body &B = bodies[C.b];
        B.vp -= B.invm.cwiseProduct(Pw);
        B.wp -= B.invI*(C.rb.cross(Pw));
    }
}

void dem_core::integrate(double dt)
{
    for(auto &B : bodies)
    {
        if(!B.active || B.fixed || B.ghost)
        continue;

        B.x += dt*(B.v + B.vp);

        dem_vec wt = B.w + B.wp;
        double wn = wt.norm();
        if(wn*dt>1.0e-14)
        {
            dem_quat dq(Eigen::AngleAxisd(wn*dt,wt/wn));
            B.q = dq*B.q;
            B.q.normalize();
        }
        B.R = B.q.toRotationMatrix();
    }

    // impulse cache for warm starting, approach velocity of binding speculative contacts
    cacheP.clear();
    cacheU.clear();
    for(auto &C : contacts)
    {
        if(C.P(0)>0.0)
        cacheP[C.key] = C.n*C.P(0) + C.t1*C.P(1) + C.t2*C.P(2);

        if(C.speculative && C.P(0)>0.0)
        cacheU[C.key] = C.un0;
    }
}

bool dem_core::mine(int a, int b) const
{
    const dem_body &A = bodies[a], &B = bodies[b];

    if(A.ghost && B.ghost)
    return false;

    // replicated - replicated: rank 0
    if(A.tier==1 && B.tier==1)
    return myrank==0;

    // replicated - distributed: owner of the distributed particle
    if(A.tier==1)
    return B.owner==myrank;
    if(B.tier==1)
    return A.owner==myrank;

    // distributed - distributed: owner of the lower global id
    return (A.id<B.id ? A.owner : B.owner)==myrank;
}
