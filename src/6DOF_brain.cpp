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

#include"6DOF_brain.h"
#include<cmath>
#include<algorithm>

namespace
{
    const double brain_tau = 1.5;          // time to turn onto the bearing of the target [s]
    const double brain_kturn = 1.41;       // turn rate of the gait: omega = kturn b U/L
    const double brain_slow = 1.5;         // slower within this distance of the target [L]
    const double brain_arrive = 0.25;      // at the target within this distance [L]
    const double brain_front = 0.35;       // the head stations of the lateral line: s/L below
}

sixdof_brain::sixdof_brain(double U_des, double b_max, double D_vis, double th_esc, double f0, double L, double rho,
                           double t_start, int lateral, const std::vector<target_point> &tp)
             : Udes_(U_des), bmax_(b_max>0.0 ? b_max : 1.0), Dvis_(D_vis), the_(th_esc), f0_(f0), L_(L), rho_(rho),
               tstart_(t_start), lat_(lateral), tp_(tp)
{
    std::stable_sort(tp_.begin(), tp_.end(), [](const target_point &a, const target_point &b) {return a.t<b.t;});
    t_start_ = tstart_;
}

const std::vector<double> &sixdof_brain::rays()
{
    static const double d = 3.14159265358979323846/180.0;
    static const std::vector<double> r = {0.0, 25.0*d, -25.0*d, 50.0*d, -50.0*d, 75.0*d, -75.0*d};
    return r;
}

Eigen::Vector3d sixdof_brain::target(double t) const
{
    // linear in t between the points, constant before the first and after the last
    if(tp_.empty())
    return c_;

    if(t<=tp_.front().t)
    return Eigen::Vector3d(tp_.front().x, tp_.front().y, tp_.front().z);

    for(size_t k=1; k<tp_.size(); ++k)
    if(t<=tp_[k].t)
    {
        const target_point &a = tp_[k-1], &b = tp_[k];
        const double w = (b.t>a.t) ? (t - a.t)/(b.t - a.t) : 1.0;
        return Eigen::Vector3d(a.x + w*(b.x - a.x), a.y + w*(b.y - a.y), a.z + w*(b.z - a.z));
    }

    return Eigen::Vector3d(tp_.back().x, tp_.back().y, tp_.back().z);
}

void sixdof_brain::sense(double dt, const Eigen::Vector3d &ax, const Eigen::Vector3d &lat, const Eigen::Vector3d &c,
                         const Eigen::Vector3d &U, const std::vector<double> *ll_p, const std::vector<double> &ll_s,
                         const std::vector<int> &ll_side)
{
    T_ += dt;
    ax_ += dt*ax;
    lt_ += dt*lat;
    Uax_ += dt*U.dot(ax);
    c_ = c;

    if(ll_p!=nullptr && !ll_p->empty())
    {
        if(llp_.size()!=ll_p->size())
        llp_.assign(ll_p->size(),0.0);

        for(size_t q=0; q<ll_p->size(); ++q)
        llp_[q] += dt*(*ll_p)[q];

        llT_ += dt;
        lls_ = ll_s;
        llside_ = ll_side;
    }
}

sixdof_brain::decision sixdof_brain::decide(double t, double f, double a, double b, const std::vector<double> *vis)
{
    decision d;
    d.t = t;

    if(T_<=0.0)
    {
        next_ = t + std::min(0.5/std::max(f,1.0e-3), 0.5);
        return d;
    }

    // the plane of the undulation: its normal (3D: z, 2D: y); the mean axes in it
    const int nrm = (lat_==1) ? 2 : 1;
    Eigen::Vector3d ex = ax_/T_, el = lt_/T_;
    ex(nrm) = 0.0;
    el(nrm) = 0.0;
    ex.normalize();
    el = el - el.dot(ex)*ex;
    el.normalize();

    const double U = Uax_/T_;
    const double Uf = std::max(U, 0.1*L_);
    const double q = 0.5*rho_*Uf*Uf;

    // the stride U/f while swimming straight at the full amplitude, after the start
    if(std::fabs(b)<0.15 && std::fabs(a - 1.0)<0.05 && t>tstart_ + 2.0/f0_ && U>0.0 && f>0.0)
    {
        const double s = U/(f*L_);
        stride_ = (nstride_==0) ? s : 0.8*stride_ + 0.2*s;
        ++nstride_;
    }

    // the target: distance and bearing to the mean axis (> 0: towards the lateral axis)
    const Eigen::Vector3d T = target(t);
    Eigen::Vector3d r = T - c_;
    r(nrm) = 0.0;
    const double dist = tp_.empty() ? 0.0 : r.norm();
    const double psi = tp_.empty() ? 0.0 : std::atan2(r.dot(el), r.dot(ex));

    // the lateral line at the head: the pressure difference across (p+ - p-) and the excess of the
    // head pressure over the mean of all probes, in units of q
    double dpl = 0.0, phd = 0.0;

    if(llT_>0.0 && !llp_.empty())
    {
        double dsum = 0.0, hsum = 0.0, asum = 0.0;
        int nd = 0, nh = 0;

        for(size_t k=0; k<llp_.size(); ++k)
        {
            const double pk = llp_[k]/llT_;
            asum += pk;

            if(lls_[k]<brain_front)
            {
                dsum += llside_[k]*pk;
                hsum += pk;
                ++nh;

                if(llside_[k]>0)
                ++nd;
            }
        }

        if(nd>0 && nh>0)
        {
            dpl = (dsum/nd)/q;
            phd = (hsum/nh - asum/double(llp_.size()))/q;
        }
    }

    // vision: the closest ray, the free space on both sides
    double dvis = -1.0, threat = 0.0, side = 0.0;
    
    if(vis!=nullptr && Dvis_>0.0 && vis->size()==rays().size())
    {
        const double D = Dvis_*L_;
        double fp = 0.0, fm = 0.0;
        dvis = D;
        
        for(size_t k=0; k<vis->size(); ++k)
        {
            const double dk = std::min((*vis)[k], D);
            dvis = std::min(dvis, dk);
            threat = std::max(threat, (1.0 - dk/D)*(1.0 - dk/D));
            
            if(rays()[k]>0.0)
            fp += dk;
            if(rays()[k]<0.0)
            fm += dk;
        }
        
        // the side with the longer rays (+: the lateral axis); equal: away from the target's side; kept
        // while the threat lasts (no dithering in front of an obstacle)
        side = (fp>fm) ? 1.0 : ((fm>fp) ? -1.0 : (psi>=0.0 ? -1.0 : 1.0));

        if(threat>0.0 && side_!=0.0)
        side = side_;

        side_ = (threat>0.0) ? side : 0.0;
        dvis /= L_;
    }
    
    rec_ = {t, T(0), T(1), T(2), dist, psi, U, f, a, b, dpl, phd, dvis, 0};

    // reset the senses
    T_ = Uax_ = llT_ = 0.0;
    ax_.setZero();
    lt_.setZero();
    std::fill(llp_.begin(), llp_.end(), 0.0);

    d.on = true;

    // 1. escape: a C-start away from the side of the higher pressure
    if(the_>0.0 && t>=refr_ && (std::fabs(dpl)>the_ || phd>the_))
    {
        const double sg = (dpl>=0.0) ? 1.0 : -1.0;
        d.f = f; d.a = 0.5; d.b = -2.5*sg; d.ts = 0.08;
        d.second = true;
        d.t2 = t + 0.15; d.f2 = std::min(1.5*f, fmax()); d.a2 = 1.0; d.b2 = 0.0; d.ts2 = 0.2;
        refr_ = t + 2.0;
        next_ = t + 0.4;
        rec_.Ff = d.f2; rec_.Fa = d.a; rec_.Fb = d.b; rec_.mode = 2;
        return d;
    }

    // 3. the target: the turn rate psi/tau, the speed from the stride
    double ft = f, bt = 0.0;

    if(!tp_.empty())
    {
        bt = (psi/brain_tau)*L_/(brain_kturn*Uf);

        if(Udes_>0.0)
        {
            const double Ud = Udes_*std::min(1.0, dist/(brain_slow*L_));
            const double st = std::max(stride_,0.1)*L_;

            // the closing speed since the last decision (a moving target: the speed in it)
            if(dprev_>=0.0 && t>tprev_ && fctl_>0.0)
            ft = fctl_ + 0.3*(Ud - (dprev_ - dist)/(t - tprev_))/st;
            else
            ft = Ud/st;

            if(std::fabs(psi)>1.0472)
            ft = std::max(ft, 0.8*f0_);

            ft = std::max(0.3*f0_, std::min(fmax(), ft));
            fctl_ = ft;
        }

        dprev_ = dist;
        tprev_ = t;

        rec_.mode = (dist<brain_arrive*L_) ? 3 : 0;
    }

    bt = std::max(-bmax_, std::min(bmax_, bt));

    // 2. avoidance: a turn to the free side, blended in with the threat, slower
    if(threat>0.0)
    {
        const double w = std::min(1.0, threat/0.25);
        bt = (1.0 - w)*bt + w*side*bmax_;
        ft = std::max(0.3*f0_, ft*(1.0 - 0.3*w));

        if(w>0.2)
        rec_.mode = 1;
    }

    d.f = ft; d.a = 1.0; d.b = bt; d.ts = 0.0;
    next_ = t + std::min(0.5/ft, 0.5);
    rec_.Ff = ft; rec_.Fa = 1.0; rec_.Fb = bt;

    return d;
}

void sixdof_brain::ckpt_state(std::vector<double> &v) const
{
    v.push_back(1.0);      // version
    for(double x : {T_, Uax_, llT_, t_start_, next_, stride_, refr_, double(nstride_), dprev_, tprev_, side_, fctl_})
    v.push_back(x);
    for(int k=0; k<3; ++k)
    {
        v.push_back(ax_(k));
        v.push_back(lt_(k));
        v.push_back(c_(k));
    }
    v.push_back(double(llp_.size()));
    v.insert(v.end(), llp_.begin(), llp_.end());
}

bool sixdof_brain::ckpt_restore(const std::vector<double> &v)
{
    if(v.size()<22 || int(v[0])!=1)
    return false;

    size_t k = 1;
    T_ = v[k++]; Uax_ = v[k++]; llT_ = v[k++]; t_start_ = v[k++]; next_ = v[k++]; stride_ = v[k++]; refr_ = v[k++];
    nstride_ = int(v[k++]);
    dprev_ = v[k++]; tprev_ = v[k++]; side_ = v[k++]; fctl_ = v[k++];
    for(int j=0; j<3; ++j)
    {
        ax_(j) = v[k++];
        lt_(j) = v[k++];
        c_(j) = v[k++];
    }
    const size_t n = size_t(v[k++]);
    if(v.size()<k+n)
    return false;
    llp_.assign(v.begin()+k, v.begin()+k+n);
    return true;
}
