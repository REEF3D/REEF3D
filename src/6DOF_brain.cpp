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
    const double brain_tau = 1.5;          // time to turn onto the wanted heading [s]
    const double brain_kturn = 1.3;        // default turn rate of the gait: omega = k b U_w/L (U_w through the water)
    const double brain_slow = 1.5;         // slower within this distance of the target [L]
    const double brain_arrive = 0.25;      // at the target within this distance [L]
    const double brain_front = 0.35;       // the head stations of the lateral line: s/L below
    const double brain_tauu = 0.8;         // default response time of the swimming speed to the gait (X 201) [s]
    const double brain_s0 = 0.5;           // default prior of the stride per tail beat at full amplitude [L]
    const double brain_tfloor = 2.5;       // the speed for a turn: the turn in this time at b_max [s]
    const double brain_psifloor = 1.0;     // ... beyond this angle to the wanted heading [rad]
    const double brain_reaff = 0.5;        // head pressure difference of a turn: -reaff b q_w
    const double brain_burst = 0.5;        // a gliding fish beats again beyond this distance of the target [L]
    const double brain_tauc0 = 2.0;        // prior of the coasting time: the speed of a glide decays with it [s]
    const double brain_cknown = 0.2;       // a glide: the current known within this (sqrt of its variance) [L/s]
    const int brain_nknown = 3;            // a glide at the target: the filter updated this many times at least
    const double brain_cdrift = 0.1;       // ... and its drift over the coasting time below this [L]
    const double brain_brake = 2.5;        // braking while gliding: the body bent with this bias (0: no braking)
    const double brain_uesc = 0.3;         // an escape: the pressures in units of q at this speed at least [L/s]
}

sixdof_brain::sixdof_brain(double U_des, double b_max, double D_vis, double th_esc, double f0, double L, double rho,
                           double t_start, int lateral, const std::vector<target_point> &tp, double s0, double k_turn,
                           double tau_u, double tau_d)
             : Udes_(U_des), bmax_(b_max>0.0 ? b_max : 1.0), Dvis_(D_vis), the_(th_esc), f0_(f0), L_(L), rho_(rho),
               tstart_(t_start), s0_(s0>0.0 ? s0 : brain_s0), kturn_(k_turn!=0.0 ? k_turn : brain_kturn),
               tauu_(tau_u>0.0 ? tau_u : brain_tauu), taud_(tau_d>0.0 ? tau_d : brain_tauu), lat_(lateral), tp_(tp)
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

    const double U = Uax_/T_;                                   // the speed along the axis over the ground
    const int p0 = 0, p1 = (nrm==1) ? 2 : 1;                    // the axes of the plane of the undulation
    const double fmx = fmax();

    // 1. the speed of the gait over the time since the last decision; from it and the track over the
    // ground the estimate of the stride and of the current (Kalman filter on s, c_p0, c_p1):
    //    ground velocity = s U_g ex + c,   dU_g/dt = (f a L - U_g)/tau_u
    if(!kf_on_)
    {
        kf_x_[0] = s0_; kf_x_[1] = kf_x_[2] = 0.0;
        for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
        kf_P_[i][j] = 0.0;
        kf_P_[0][0] = 0.09*s0_*s0_; kf_P_[1][1] = kf_P_[2][2] = 0.3*0.3;
        // the speed of the gait: from the track over the first interval (a body that starts from rest may
        // still be accelerating), not yet known
        ug_ = -1.0;
        tdec_ = t;
        pdec_ = c_;
        bdec_ = b;
        adec_ = a;
        Tset_ = false;
        vT_.setZero();
        tauc_ = brain_tauc0;
        nkf_ = 0;
        kf_on_ = true;
    }

    const double dt = t - tdec_;

    if(dt>0.0)
    {
        const double y[2] = {(c_(p0) - pdec_(p0))/dt, (c_(p1) - pdec_(p1))/dt};

        if(ug_<0.0)
        {
            // the first interval: the speed of the gait from the track (no current assumed), no update
            const double uw = y[0]*ex(p0) + y[1]*ex(p1);
            ug_ = std::min(f*a*L_, std::max(0.0, uw)/kf_x_[0]);
        }
        else if(std::min(a, adec_)<0.95)
        {
            // gliding (or just after it): the body coasts, the model of the speed does not hold; no update of
            // the estimate, the speed of the gait from the track (the velocity through the water along the axis)
            const double uw = (y[0] - kf_x_[1])*ex(p0) + (y[1] - kf_x_[2])*ex(p1);
            const double u0 = ug_*kf_x_[0];
            ug_ = std::max(0.0, uw)/kf_x_[0];

            // the coasting time from the slowing down over a whole glide (gliding at both decisions), well
            // above the uncertainty of the current
            const double cm = std::sqrt(kf_x_[1]*kf_x_[1] + kf_x_[2]*kf_x_[2]);
            if(std::max(a, adec_)<0.05 && std::max(std::fabs(b), std::fabs(bdec_))<0.05 && uw>3.0*cm + 0.05*L_ && uw<u0)
            tauc_ = 0.7*tauc_ + 0.3*std::min(30.0, std::max(0.5, dt/std::log(u0/uw)));
        }
        else
        {
            const double fal = f*a*L_;
            // the speed follows a faster gait with tau_u, a slower one with tau_d (X 201)
            const double tu = (fal>=ug_) ? tauu_ : taud_;
            const double e = std::exp(-dt/tu);
            const double ug1 = fal + (ug_ - fal)*e;
            const double ugm = fal + (ug_ - fal)*(tu/dt)*(1.0 - e);
            ug_ = ug1;
            const double bm = 0.5*(b + bdec_);
            const double R = 0.07*0.07*(1.0 + 4.0*bm*bm);
            const double H[2][3] = {{ugm*ex(p0), 1.0, 0.0}, {ugm*ex(p1), 0.0, 1.0}};
            kf_update(H, y, R, dt);
            ++nkf_;
        }
    }

    const double s = kf_x_[0];
    Eigen::Vector3d cur = Eigen::Vector3d::Zero();
    cur(p0) = kf_x_[1];
    cur(p1) = kf_x_[2];
    const double Uw = std::max((ug_<0.0 ? U : ug_*s), 0.1*L_);
    const double q = 0.5*rho_*Uw*Uw;

    // 2. the target: the wanted velocity over the ground (towards the target, U_des, slower within
    // 1.5 L; plus the velocity of the target as seen), and the velocity through the water that gives it
    const bool has_t = !tp_.empty();
    const Eigen::Vector3d T = target(t);
    Eigen::Vector3d r = T - c_;
    r(nrm) = 0.0;
    const double dist = has_t ? std::sqrt(r(0)*r(0) + r(1)*r(1) + r(2)*r(2)) : 0.0;

    if(has_t && Tset_ && dt>0.0)
    {
        const double al = std::min(1.0, dt/1.0);
        for(int k=0; k<3; ++k)
        vT_(k) = (1.0 - al)*vT_(k) + al*(T(k) - Tprev_(k))/dt;
    }
    Tprev_ = T;
    Tset_ = has_t;

    Eigen::Vector3d vw = Eigen::Vector3d::Zero();

    if(has_t)
    {
        double sp = (Udes_>0.0 ? Udes_ : f0_*s*L_)*std::min(1.0, dist/(brain_slow*L_));

        // cruising: the closing speed also on a detour (the target up to 60 deg off the axis)
        if(dist>brain_slow*L_)
        {
            const double ca = (r(0)*ex(0) + r(1)*ex(1) + r(2)*ex(2))/dist;
            if(ca>0.0)
            sp = sp/std::max(ca, 0.5);
        }

        for(int k=0; k<3; ++k)
        vw(k) = vT_(k) + (dist>1.0e-9 ? sp*r(k)/dist : 0.0);
    }

    Eigen::Vector3d w = vw - cur;
    w(nrm) = 0.0;
    const double Ureq = std::sqrt(w(0)*w(0) + w(1)*w(1) + w(2)*w(2));
    const double psi = (Ureq>0.02*L_) ? std::atan2(w(0)*el(0) + w(1)*el(1) + w(2)*el(2), w(0)*ex(0) + w(1)*ex(1) + w(2)*ex(2)) : 0.0;

    // 3. the frequency from the speed through the water (no integrator), at least the speed a turn-around
    // needs; the bias from the turn rate psi/tau at the speed through the water
    double ft = std::max(0.3*f0_, std::min(fmx, f)), bt = 0.0;
    int mode = 0;

    if(has_t)
    {
        const double Uturn = (std::fabs(psi)>brain_psifloor) ? std::min(std::fabs(psi)/brain_tfloor*L_/(std::fabs(kturn_)*bmax_), fmx*s*L_) : 0.0;
        const double Ucmd = std::max(Ureq, Uturn);
        ft = std::max(0.3*f0_, std::min(fmx, Ucmd/(s*L_)));
        bt = (psi/brain_tau)*L_/(kturn_*Uw);
        mode = (dist<brain_arrive*L_) ? 3 : 0;
    }

    bt = std::max(-bmax_, std::min(bmax_, bt));

    // the lateral line at the head: the pressure difference across (p+ - p-) and the excess of the
    // head pressure over the mean of all probes, in units of q through the water; the signal of the
    // own turn (-reaff b q) removed
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
            dpl = (dsum/nd)/q + brain_reaff*0.5*(b + bdec_);
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
        
        // the side with the longer rays (+: the lateral axis); equal: away from the wanted heading; kept
        // while the threat lasts (no dithering in front of an obstacle)
        side = (fp>fm) ? 1.0 : ((fm>fp) ? -1.0 : (psi>=0.0 ? -1.0 : 1.0));

        if(threat>0.0 && side_!=0.0)
        side = side_;

        side_ = (threat>0.0) ? side : 0.0;
        dvis /= L_;
    }
    
    // glide: in a current and with a target slower than its slowest swimming speed through the water and
    // seeing nothing close, the fish stops beating its tail and coasts (amplitude 0, no bias)
    //  - early: the target ahead (within 30 deg of the axis) and the coasting distance (its speed through
    //    the water times the coasting time, learnt from its glides) reaching it, the current known
    //  - at the target (brain_arrive), needing less than that speed, the current learnt (within brain_cknown,
    //    after brain_nknown updates of the filter: a glide stops the estimate, and a fish that arrives before
    //    it has learnt the current drifts off)
    // It beats again beyond brain_burst of the target when the coast no longer reaches the target ahead, and
    // with the target behind only beyond twice brain_burst, while the bent body turns it back (burst and
    // coast), for at least a tail beat
    double at = 1.0;
    {
        const double Umin = 0.3*f0_*s*L_;
        const double cm = std::sqrt(cur(p0)*cur(p0) + cur(p1)*cur(p1));
        const double dstop = std::max(0.0, ug_)*s*tauc_;
        const double rx = r(0)*ex(0) + r(1)*ex(1) + r(2)*ex(2), ry = r(0)*el(0) + r(1)*el(1) + r(2)*el(2);
        const bool onaxis = has_t && rx>0.0 && std::fabs(std::atan2(ry, rx))<0.5236;
        const bool calm = cm<Umin && vT_.norm()<Umin;
        const bool cvar = kf_P_[1][1] + kf_P_[2][2]<std::pow(brain_cknown*L_, 2.0);
        const bool known = cvar && cm*tauc_<brain_cdrift*L_;
        const bool learnt = cvar && nkf_>=brain_nknown;

        if(!has_t || threat>0.0 || !calm)
        glide_ = false;
        else if(glide_)
        {
            // beyond brain_burst of the target it beats again: with the target ahead when the coast no longer
            // reaches it; with the target behind only beyond twice brain_burst, until then it keeps braking (the
            // bent body turns it towards the target; beating with the target behind costs a turn-round of 1.5 L)
            glide_ = !(dist>brain_burst*L_ && (rx<=0.0 ? dist>2.0*brain_burst*L_ : dstop<dist - brain_burst*L_));
            if(!glide_)
            tburst_ = t + 1.0/f0_;
        }
        else if(t>=tburst_)
        glide_ = (known && onaxis && dstop>=dist - brain_arrive*L_) || (learnt && dist<brain_arrive*L_ && Ureq<Umin);

        if(glide_)
        {
            // the tail still: the slowest frequency, so that the bend of the body is released slowly (released
            // at the frequency of a turn, the bend of the brake throws the fish forward); braking keeps the
            // frequency of the turn, so that the bend builds quickly
            const double fturn = ft;
            ft = 0.3*f0_;
            at = 0.0;
            bt = 0.0;
            mode = 5;

            // braking: the coast would carry the fish past the target, carries it away from it, or past it
            // to the side (the target ahead but off the axis): the body bent, to the side that turns it towards
            // the wanted heading at the start of the braking, kept while it brakes (bending to and fro would
            // swim); so the fish turns towards the target while it coasts, and beats again with it ahead
            if(brain_brake>0.0 && (rx>0.0 ? (dstop>dist + brain_arrive*L_ || (!onaxis && dist>brain_arrive*L_)) : dist>brain_arrive*L_))
            {
                if(bside_==0.0)
                bside_ = (psi*kturn_>=0.0) ? 1.0 : -1.0;
                bt = bside_*brain_brake;
                ft = fturn;
                mode = 6;
            }
            else
            bside_ = 0.0;
        }
        else
        bside_ = 0.0;
    }

    rec_ = {t, T(0), T(1), T(2), dist, psi, U, f, a, b, dpl, phd, dvis, s, cur(0), cur(1), cur(2), Uw, 0};
    
    perc_.ex = ex;
    perc_.el = el;
    perc_.c = c_;
    perc_.T = T;
    perc_.U = U;
    perc_.stride = s;
    perc_.target = has_t;
    perc_.llp.assign(llp_.size(), 0.0);
    for(size_t k=0; k<llp_.size() && llT_>0.0; ++k)
    perc_.llp[k] = llp_[k]/llT_;
    perc_.rays = (vis!=nullptr) ? *vis : std::vector<double>();

    // reset the senses
    T_ = Uax_ = llT_ = 0.0;
    ax_.setZero();
    lt_.setZero();
    std::fill(llp_.begin(), llp_.end(), 0.0);

    tdec_ = t;
    pdec_ = c_;
    bdec_ = bt;
    adec_ = at;
    d.on = true;

    // 1. escape: a C-start away from the side of the higher pressure; not while the fish coasts or brakes (the
    // bent body of the braking gives a pressure across the head far above the -reaff b q of a swimming turn);
    // the pressures against q at brain_uesc at least: at a low speed (the start, a station) the pressures of
    // the own turns and of the start-up are large against q of that speed, and they are no stimulus
    const double qesc = q/std::max(q, 0.5*rho_*std::pow(brain_uesc*L_, 2.0));
    if(the_>0.0 && t>=refr_ && !(glide_ || a<0.95) && (std::fabs(dpl)*qesc>the_ || phd*qesc>the_))
    {
        const double sg = (dpl>=0.0) ? 1.0 : -1.0;
        d.f = f; d.a = 0.5; d.b = -2.5*sg; d.ts = 0.08;
        d.second = true;
        d.t2 = t + 0.15; d.f2 = std::min(1.5*f, fmx); d.a2 = 1.0; d.b2 = 0.0; d.ts2 = 0.2;
        refr_ = t + 2.0;
        next_ = t + 0.4;
        bdec_ = d.b;
        adec_ = d.a;
        glide_ = false;
        bside_ = 0.0;
        rec_.Ff = d.f2; rec_.Fa = d.a; rec_.Fb = d.b; rec_.mode = 2;
        return d;
    }

    // 2. avoidance: a turn to the free side, blended in with the threat, slower
    if(threat>0.0)
    {
        const double wv = std::min(1.0, threat/0.25);
        bt = (1.0 - wv)*bt + wv*side*bmax_;
        ft = std::max(0.3*f0_, ft*(1.0 - 0.3*wv));
        bdec_ = bt;

        if(wv>0.2)
        mode = 1;
    }

    d.f = ft; d.a = at; d.b = bt; d.ts = 0.0;
    next_ = t + std::min(0.5/ft, 0.5);
    rec_.Ff = ft; rec_.Fa = at; rec_.Fb = bt; rec_.mode = mode;

    return d;
}

void sixdof_brain::kf_update(const double H[2][3], const double y[2], double R, double dt)
{
    // random walk of the stride and of the current, then the two components of the ground velocity
    kf_P_[0][0] += 1.0e-4*dt;
    kf_P_[1][1] += 2.5e-3*dt;
    kf_P_[2][2] += 2.5e-3*dt;

    for(int m=0; m<2; ++m)
    {
        double Ph[3], K[3];
        for(int i=0; i<3; ++i)
        Ph[i] = kf_P_[i][0]*H[m][0] + kf_P_[i][1]*H[m][1] + kf_P_[i][2]*H[m][2];
        const double S = H[m][0]*Ph[0] + H[m][1]*Ph[1] + H[m][2]*Ph[2] + R;
        for(int i=0; i<3; ++i)
        K[i] = Ph[i]/S;
        const double e = y[m] - (H[m][0]*kf_x_[0] + H[m][1]*kf_x_[1] + H[m][2]*kf_x_[2]);
        for(int i=0; i<3; ++i)
        kf_x_[i] += K[i]*e;
        for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
        kf_P_[i][j] -= K[i]*Ph[j];
    }

    kf_x_[0] = std::min(1.6*s0_, std::max(0.5*s0_, kf_x_[0]));
}

void sixdof_brain::ckpt_state(std::vector<double> &v) const
{
    v.push_back(5.0);      // version
    for(double x : {T_, Uax_, llT_, t_start_, next_, refr_, side_, kf_on_ ? 1.0 : 0.0, ug_, tdec_, bdec_, Tset_ ? 1.0 : 0.0, glide_ ? 1.0 : 0.0, adec_, tauc_, tburst_, bside_, double(nkf_)})
    v.push_back(x);
    for(int k=0; k<3; ++k)
    {
        v.push_back(ax_(k));
        v.push_back(lt_(k));
        v.push_back(c_(k));
        v.push_back(pdec_(k));
        v.push_back(Tprev_(k));
        v.push_back(vT_(k));
        v.push_back(kf_x_[k]);
        for(int j=0; j<3; ++j)
        v.push_back(kf_P_[k][j]);
    }
    v.push_back(double(llp_.size()));
    v.insert(v.end(), llp_.begin(), llp_.end());
}

bool sixdof_brain::ckpt_restore(const std::vector<double> &v)
{
    // version 4 (patches 0016-0022): without the count of the updates (taken as learnt)
    const int ver = v.empty() ? 0 : int(v[0]);
    if(v.size()<55 || (ver!=4 && ver!=5))
    return false;

    size_t k = 1;
    T_ = v[k++]; Uax_ = v[k++]; llT_ = v[k++]; t_start_ = v[k++]; next_ = v[k++]; refr_ = v[k++]; side_ = v[k++];
    kf_on_ = (v[k++]!=0.0);
    ug_ = v[k++]; tdec_ = v[k++]; bdec_ = v[k++];
    Tset_ = (v[k++]!=0.0);
    glide_ = (v[k++]!=0.0);
    adec_ = v[k++];
    tauc_ = v[k++];
    tburst_ = v[k++];
    bside_ = v[k++];
    nkf_ = (ver>=5) ? int(v[k++]) : brain_nknown;
    for(int j=0; j<3; ++j)
    {
        ax_(j) = v[k++];
        lt_(j) = v[k++];
        c_(j) = v[k++];
        pdec_(j) = v[k++];
        Tprev_(j) = v[k++];
        vT_(j) = v[k++];
        kf_x_[j] = v[k++];
        for(int i=0; i<3; ++i)
        kf_P_[j][i] = v[k++];
    }
    const size_t n = size_t(v[k++]);
    if(v.size()<k+n)
    return false;
    llp_.assign(v.begin()+k, v.begin()+k+n);
    return true;
}

void sixdof_brain::external(const decision &d, double next)
{
    rec_.Ff = d.second ? d.f2 : d.f;
    rec_.Fa = d.a;
    rec_.Fb = d.b;
    rec_.mode = 4;
    next_ = next;
}
