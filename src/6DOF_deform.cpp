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

#include"6DOF_deform.h"
#include"6DOF_geometry.h"
#include<cmath>
#include<algorithm>
#include<utility>
#include<tuple>

namespace
{
    const double fish_pi = 3.14159265358979323846;
    const int fish_bins = 256;
    const int fish_stations = 2049;
}

sixdof_deform::sixdof_deform(int type, double a, double lambda, double f, double tramp, int lateral, bool recoil)
              : type_(type>10 ? type-10 : type), inext_(type>10), lat_(lateral), a_(a), lambda_(lambda>0.0 ? lambda : 1.0), f_(f), tramp_(tramp), recoil_(recoil)
{
}

void sixdof_deform::setup(const sixdof_geometry &g)
{
    // body length and head (largest x of the body frame)
    xhead_ = -1.0e20;
    xtail_ =  1.0e20;

    for(int n=0; n<g.tricount; ++n)
    for(int q=0; q<3; ++q)
    {
        xhead_ = std::max(xhead_, g.tri_x0[n][q]);
        xtail_ = std::min(xtail_, g.tri_x0[n][q]);
    }

    L_ = (xhead_>xtail_) ? xhead_ - xtail_ : 1.0;

    // volume left of the plane x = xi:  V(xi) = int_V H(xi - x) dV = oint min(x,xi) n_x dS
    // (divergence theorem with the field (min(x,xi),0,0)), with a midpoint rule on a uniform
    // subdivision of every triangle into M^2 triangles of an x-extent of at most half a bin
    const double dxb = L_/double(fish_bins);
    std::vector<std::tuple<double,double,double> > qp;     // x of the quadrature point, n_x dS, lateral coordinate
    double **tl = (lat_==1) ? g.tri_y0 : g.tri_z0;

    for(int n=0; n<g.tricount; ++n)
    {
        const double ax = g.tri_x0[n][0], ay = g.tri_y0[n][0], az = g.tri_z0[n][0];
        const double bx = g.tri_x0[n][1], by = g.tri_y0[n][1], bz = g.tri_z0[n][1];
        const double cx = g.tri_x0[n][2], cy = g.tri_y0[n][2], cz = g.tri_z0[n][2];
        const double al = tl[n][0], bl = tl[n][1], cl = tl[n][2];

        // x component of the area vector (outward orientation of the hull)
        const double Sx = 0.5*((by-ay)*(cz-az) - (bz-az)*(cy-ay));

        if(Sx==0.0)
        continue;

        const double ext = std::max(ax,std::max(bx,cx)) - std::min(ax,std::min(bx,cx));
        const int M = std::min(64, std::max(1, int(std::ceil(ext/(0.5*dxb)))));
        const double w = Sx/double(M*M);

        for(int i=0; i<M; ++i)
        for(int j=0; j<M-i; ++j)
        {
            // upward triangle, centroid at barycentric (i+1/3, j+1/3)/M
            double u = (i + 1.0/3.0)/double(M);
            double v = (j + 1.0/3.0)/double(M);
            qp.push_back(std::make_tuple(ax + u*(bx-ax) + v*(cx-ax), w, al + u*(bl-al) + v*(cl-al)));

            // downward triangle
            if(i+j<M-1)
            {
            u = (i + 2.0/3.0)/double(M);
            v = (j + 2.0/3.0)/double(M);
            qp.push_back(std::make_tuple(ax + u*(bx-ax) + v*(cx-ax), w, al + u*(bl-al) + v*(cl-al)));
            }
        }
    }

    std::sort(qp.begin(), qp.end());

    const int nq = int(qp.size());
    std::vector<double> sw(nq+1,0.0), swx(nq+1,0.0);    // prefix sums of w and w x
    // ... and of w l, w l x, w l^2, w l^2 x (the moments of the lateral coordinate l, field (min(x,xi) l^k,0,0))
    std::vector<double> s1(nq+1,0.0), s1x(nq+1,0.0), s2(nq+1,0.0), s2x(nq+1,0.0);

    for(int q=0; q<nq; ++q)
    {
        const double x = std::get<0>(qp[q]), w = std::get<1>(qp[q]), l = std::get<2>(qp[q]);
        sw[q+1]  = sw[q]  + w;
        swx[q+1] = swx[q] + w*x;
        s1[q+1]  = s1[q]  + w*l;
        s1x[q+1] = s1x[q] + w*l*x;
        s2[q+1]  = s2[q]  + w*l*l;
        s2x[q+1] = s2x[q] + w*l*l*x;
    }

    auto first = [&](double xi)
    {
        return int(std::lower_bound(qp.begin(), qp.end(), std::make_tuple(xi,-1.0e300,-1.0e300)) - qp.begin());
    };

    // V(xi) = sum_{x<xi} w x + xi sum_{x>=xi} w
    auto Vleft = [&](double xi)
    {
        const int m = first(xi);
        return swx[m] + xi*(sw[nq] - sw[m]);
    };
    
    // int_V H(xi - x) l dV and int_V H(xi - x) l^2 dV
    auto L1left = [&](double xi)
    {
        const int m = first(xi);
        return s1x[m] + xi*(s1[nq] - s1[m]);
    };
    auto L2left = [&](double xi)
    {
        const int m = first(xi);
        return s2x[m] + xi*(s2[nq] - s2[m]);
    };

    // closed hull: V(x_head) = volume; an inward oriented hull gives a negative volume
    vol_ = Vleft(xhead_) - Vleft(xtail_);
    const double sgn = (vol_<0.0) ? -1.0 : 1.0;
    vol_ = std::fabs(vol_);

    xb_.assign(fish_bins,0.0);
    mb_.assign(fish_bins,0.0);
    S0_ = S1_ = S2_ = 0.0;

    n1b_.assign(fish_bins,0.0);
    n2b_.assign(fish_bins,0.0);
    In_ = 0.0;

    double Vl = Vleft(xtail_);
    double L1l = L1left(xtail_), L2l = L2left(xtail_);

    for(int b=0; b<fish_bins; ++b)
    {
        const double xr = xtail_ + (b+1)*dxb;
        const double Vr = Vleft(xr);

        xb_[b] = xtail_ + (b+0.5)*dxb;
        mb_[b] = std::max(0.0, sgn*(Vr - Vl));
        Vl = Vr;

        S0_ += mb_[b];
        S1_ += mb_[b]*xb_[b];
        S2_ += mb_[b]*xb_[b]*xb_[b];
        
        const double L1r = L1left(xr), L2r = L2left(xr);
        n1b_[b] = sgn*(L1r - L1l);
        n2b_[b] = std::max(0.0, sgn*(L2r - L2l));
        L1l = L1r;
        L2l = L2r;
        
        In_ += mb_[b]*xb_[b]*xb_[b] + n2b_[b];
    }
}

double sixdof_deform::envelope(double s) const
{
    if(type_==2)
    return (a_/0.1)*(0.02 - 0.0825*s + 0.1625*s*s);

    if(type_==3)
    return a_*std::exp(s - 1.0);

    return a_*(s + 0.03125)/1.03125;
}

double sixdof_deform::envelope_s(double s) const
{
    // dA/ds
    if(type_==2)
    return (a_/0.1)*(-0.0825 + 0.325*s);

    if(type_==3)
    return a_*std::exp(s - 1.0);

    return a_/1.03125;
}

void sixdof_deform::raw(double x, double &h, double &ht) const
{
    double s = (xhead_ - x)/L_;
    s = std::min(1.0, std::max(0.0, s));

    const double A = L_*envelope(s);
    const double om = 2.0*fish_pi*f_;
    const double ph = 2.0*fish_pi*s/lambda_ - om*t_;

    h  = r_*A*std::sin(ph);
    ht = -r_*A*om*std::cos(ph) + rdot_*A*std::sin(ph);
}

void sixdof_deform::fit(const std::vector<double> &g, double &k0, double &k1) const
{
    // k0 + k1 x: least squares fit of g with the weights of the mass distribution
    double G0 = 0.0, G1 = 0.0;

    for(int b=0; b<int(xb_.size()); ++b)
    {
        G0 += mb_[b]*g[b];
        G1 += mb_[b]*xb_[b]*g[b];
    }

    const double det = S0_*S2_ - S1_*S1_;

    if(S0_<=0.0 || std::fabs(det)<=1.0e-30*S0_*S2_)
    {
    k0 = k1 = 0.0;
    return;
    }

    k1 = (S0_*G1 - S1_*G0)/det;
    k0 = (G0 - k1*S1_)/S0_;
}

void sixdof_deform::set_time(double t)
{
    t_ = t;

    // ramp r(t) = 0.5 (1 - cos(pi t/t_ramp)) and its rate
    r_ = 1.0;
    rdot_ = 0.0;

    if(tramp_>0.0 && t_<tramp_)
    {
        r_    = (t_>0.0) ? 0.5*(1.0 - std::cos(fish_pi*t_/tramp_)) : 0.0;
        rdot_ = (t_>0.0) ? 0.5*(fish_pi/tramp_)*std::sin(fish_pi*t_/tramp_) : 0.0;
    }

    c0_ = c1_ = d0_ = d1_ = 0.0;
    
    if(inext_)
    {
    midline(t_,r_,rdot_,tab_);
    cx_ = tab_.cx;
    cl_ = tab_.cl;
    vx_ = tab_.vx;
    vl_ = tab_.vl;
    phit_ = tab_.phit;
    phi_ = (recoil_ && !xb_.empty()) ? phi_at(t_) : 0.0;
    return;
    }

    if(!recoil_ || xb_.empty())
    return;

    // recoil: no net lateral momentum and no angular momentum of the deformation about the CoG
    std::vector<double> hb(xb_.size()), htb(xb_.size());

    for(int b=0; b<int(xb_.size()); ++b)
    raw(xb_[b],hb[b],htb[b]);

    fit(hb,c0_,c1_);
    fit(htb,d0_,d1_);
}

double sixdof_deform::disp(double x) const
{
    double h, ht;
    raw(x,h,ht);

    return h - c0_ - c1_*x;
}

double sixdof_deform::vel(double x) const
{
    double h, ht;
    raw(x,h,ht);

    return ht - d0_ - d1_*x;
}

// ---------------------------------------------------------------------------------------------
// inextensible midline (types 11-13): the lateral displacement y(s,t) along the arc length s of the
// midline, the axial position from the inextensibility, the cross-sections normal to the midline
// ---------------------------------------------------------------------------------------------

void sixdof_deform::midline(double t, double r, double rdot, mtab &m) const
{
    // midline tables over s in [0,L] (s from the head) at time t: angle theta (sin theta = dy/ds),
    // its rate, the lateral displacement y and its rate, the axial displacement ax = int_0^s
    // (1 - cos theta) (towards the head) and its rate; the recoil translation
    const int N = fish_stations;
    const double ds = L_/double(N-1);
    const double om = 2.0*fish_pi*f_;
    const double k = 2.0*fish_pi/lambda_;
    
    m.th.resize(N); m.tht.resize(N); m.ym.resize(N); m.ymt.resize(N); m.ax.resize(N); m.axt.resize(N);
    
    for(int q=0; q<N; ++q)
    {
        const double sg = double(q)/double(N-1);     // s/L
        const double A = envelope(sg), As = envelope_s(sg);
        const double ph = k*sg - om*t;
        const double sp = std::sin(ph), cp = std::cos(ph);
        
        // y = r L A sin ph, dy/ds = r (A_s sin ph + A k cos ph)
        const double ys  = r*(As*sp + A*k*cp);
        const double yst = -r*om*(As*cp - A*k*sp) + rdot*(As*sp + A*k*cp);
        
        const double sn = std::max(-0.999, std::min(0.999, ys));
        const double cs = std::sqrt(1.0 - sn*sn);
        
        m.th[q]  = std::asin(sn);
        m.tht[q] = (std::fabs(ys)<0.999) ? yst/cs : 0.0;
        m.ym[q]  = r*L_*A*sp;
        m.ymt[q] = -r*L_*A*om*cp + rdot*L_*A*sp;
    }
    
    // axial displacement (trapezoidal rule)
    m.ax[0] = m.axt[0] = 0.0;
    
    for(int q=1; q<N; ++q)
    {
        m.ax[q]  = m.ax[q-1]  + 0.5*ds*((1.0 - std::cos(m.th[q-1])) + (1.0 - std::cos(m.th[q])));
        m.axt[q] = m.axt[q-1] + 0.5*ds*(std::sin(m.th[q-1])*m.tht[q-1] + std::sin(m.th[q])*m.tht[q]);
    }
    
    // recoil translation: the mean displacement and velocity (axial, lateral) of the volume
    m.cx = m.cl = m.vx = m.vl = m.phit = 0.0;
    
    if(!recoil_ || xb_.empty() || S0_<=0.0)
    return;
    
    double Dx = 0.0, Dl = 0.0, Vx = 0.0, Vl = 0.0;
    
    for(int b=0; b<int(xb_.size()); ++b)
    {
        station st;
        at(m, xhead_ - xb_[b], st);
        const double sn = std::sin(st.th), cs = std::cos(st.th);
        
        // d = (ax + l sin th, y + l (cos th - 1)), v = (axt + l cos th tht, yt - l sin th tht)
        Dx += st.ax*mb_[b] + sn*n1b_[b];
        Dl += st.y*mb_[b] + (cs - 1.0)*n1b_[b];
        Vx += st.axt*mb_[b] + cs*st.tht*n1b_[b];
        Vl += st.yt*mb_[b] - sn*st.tht*n1b_[b];
    }
    
    m.cx = Dx/S0_;
    m.cl = Dl/S0_;
    m.vx = Vx/S0_;
    m.vl = Vl/S0_;
    
    // recoil rotation rate: the angular momentum of the translated deformation, X1 x v1, over the
    // polar moment of the deformed volume, |X1|^2 (both invariant under the rotation itself); per
    // bin X1 = A0 + l A1, v1 = B0 + l B1 with the moments of l
    double Lt = 0.0, Ip = 0.0;
    
    for(int b=0; b<int(xb_.size()); ++b)
    {
        station st;
        at(m, xhead_ - xb_[b], st);
        const double sn = std::sin(st.th), cs = std::cos(st.th);
        const double x = xb_[b];
        
        const double A0a = x + st.ax - m.cx, A0l = st.y - m.cl, A1a = sn, A1l = cs;
        const double B0a = st.axt - m.vx, B0l = st.yt - m.vl, B1a = cs*st.tht, B1l = -sn*st.tht;
        
        Lt += mb_[b]*(A0a*B0l - A0l*B0a) + n1b_[b]*(A0a*B1l - A0l*B1a + A1a*B0l - A1l*B0a) + n2b_[b]*(A1a*B1l - A1l*B1a);
        Ip += mb_[b]*(A0a*A0a + A0l*A0l) + 2.0*n1b_[b]*(A0a*A1a + A0l*A1l) + n2b_[b];
    }
    
    m.phit = (Ip>0.0) ? Lt/Ip : 0.0;
}

void sixdof_deform::ramp(double t, double &r, double &rdot) const
{
    // r(t) = 0.5 (1 - cos(pi t/t_ramp)) and its rate
    r = 1.0;
    rdot = 0.0;
    
    if(tramp_>0.0 && t<tramp_)
    {
        r    = (t>0.0) ? 0.5*(1.0 - std::cos(fish_pi*t/tramp_)) : 0.0;
        rdot = (t>0.0) ? 0.5*(fish_pi/tramp_)*std::sin(fish_pi*t/tramp_) : 0.0;
    }
}

double sixdof_deform::phit_at(double t) const
{
    mtab m;
    double r, rdot;
    ramp(t,r,rdot);
    midline(t,r,rdot,m);
    return m.phit;
}

double sixdof_deform::phi_at(double t)
{
    // recoil rotation phi(t) = int_0^t dphi/dt: the rate depends on t only (not on phi), so phi is a
    // quadrature, kept at the multiples of dtau = T/64 (extended as later times are asked for; a
    // pure function of t: stages, mesh refinement and restarts give the same angle), 4-point Gauss
    // on each interval
    if(t<=0.0)
    return 0.0;
    
    const double dtau = (f_>0.0) ? 1.0/(64.0*f_) : 0.01;
    
    static const double gx[4] = {-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053};
    static const double gw[4] = { 0.347854845137454,  0.652145154862546, 0.652145154862546, 0.347854845137454};
    
    auto gauss = [&](double a, double b)
    {
        double sum = 0.0;
        for(int q=0; q<4; ++q)
        sum += gw[q]*phit_at(0.5*(a+b) + 0.5*(b-a)*gx[q]);
        return 0.5*(b-a)*sum;
    };
    
    const int K = int(t/dtau);
    
    if(phik_.empty())
    phik_.push_back(0.0);
    
    while(int(phik_.size())<=K)
    {
        const int k = int(phik_.size());
        phik_.push_back(phik_[k-1] + gauss((k-1)*dtau, k*dtau));
    }
    
    return phik_[K] + gauss(K*dtau, t);
}

void sixdof_deform::at(double s, station &st) const
{
    at(tab_,s,st);
}

void sixdof_deform::at(const mtab &m, double s, station &st) const
{
    // linear interpolation of the midline tables at the arc length s (clamped to [0,L])
    const int N = int(m.th.size());
    double u = std::min(1.0, std::max(0.0, s/L_))*double(N-1);
    int q = std::min(N-2, int(u));
    const double w = u - q;
    
    st.th  = (1.0-w)*m.th[q]  + w*m.th[q+1];
    st.tht = (1.0-w)*m.tht[q] + w*m.tht[q+1];
    st.y   = (1.0-w)*m.ym[q]  + w*m.ym[q+1];
    st.yt  = (1.0-w)*m.ymt[q] + w*m.ymt[q+1];
    st.ax  = (1.0-w)*m.ax[q]  + w*m.ax[q+1];
    st.axt = (1.0-w)*m.axt[q] + w*m.axt[q+1];
    st.ths = (m.th[q+1] - m.th[q])/(L_/double(N-1));
}

void sixdof_deform::displacement_raw(double x, double l, double &dx, double &dl) const
{
    // displacement of the material point at (x, l) of the undeformed body frame (axial, lateral):
    // the midline map, the recoil translation and the recoil rotation (a rigid rotation by -phi
    // about the CoG, which keeps the length of the midline)
    if(!inext_)
    {
    dx = 0.0;
    dl = disp(x);
    return;
    }
    
    station st;
    at(xhead_ - x, st);
    const double sn = std::sin(st.th), cs = std::cos(st.th);
    
    const double a1 = x + st.ax + l*sn - cx_;
    const double l1 = st.y + l*cs - cl_;
    const double cp = std::cos(phi_), sp = std::sin(phi_);
    
    dx =  cp*a1 + sp*l1 - x;
    dl = -sp*a1 + cp*l1 - l;
}

void sixdof_deform::velocity_raw(double x, double l, double &vx, double &vl) const
{
    // velocity of the material point at (x, l) of the undeformed body frame: the time derivative
    // of displacement
    if(!inext_)
    {
    vx = 0.0;
    vl = vel(x);
    return;
    }
    
    station st;
    at(xhead_ - x, st);
    const double sn = std::sin(st.th), cs = std::cos(st.th);
    
    const double a1 = x + st.ax + l*sn - cx_;
    const double l1 = st.y + l*cs - cl_;
    const double u1 = st.axt + l*cs*st.tht - vx_;
    const double w1 = st.yt - l*sn*st.tht - vl_;
    const double cp = std::cos(phi_), sp = std::sin(phi_);
    
    // d/dt (R(-phi) X1) = R(-phi) dX1/dt - phi_t (-l_f, a_f)... with (a_f, l_f) = R(-phi) X1
    const double af =  cp*a1 + sp*l1;
    const double lf = -sp*a1 + cp*l1;
    
    vx =  cp*u1 + sp*w1 + phit_*lf;
    vl = -sp*u1 + cp*w1 - phit_*af;
}

void sixdof_deform::velocity_at_raw(double px, double pl, double &vx, double &vl) const
{
    // velocity of the deformation at the point (px, pl) of the body frame (a point of the deformed
    // body): its material coordinates (x, l) from Newton iterations on x + d(x,l) = p, with the
    // Jacobian of the midline map (tangent t = (-cos th, sin th) per arc length s = x_head - x,
    // normal (sin th, cos th) per lateral coordinate) rotated by the recoil rotation
    if(!inext_)
    {
    vx = 0.0;
    vl = vel(px);
    return;
    }
    
    const double cp = std::cos(phi_), sp = std::sin(phi_);
    
    // start: the point without the recoil
    double x = cp*px - sp*pl + cx_;
    double l = sp*px + cp*pl + cl_;
    
    for(int it=0; it<10; ++it)
    {
        double dx, dl;
        displacement_raw(x,l,dx,dl);
        
        // residual, rotated back into the frame of the midline map
        const double ra = x + dx - px, rb = l + dl - pl;
        
        if(std::fabs(ra) + std::fabs(rb) < 1.0e-13*L_)
        break;
        
        const double rx = cp*ra - sp*rb;
        const double rl = sp*ra + cp*rb;
        
        station st;
        at(xhead_ - x, st);
        const double sn = std::sin(st.th), cs = std::cos(st.th);
        
        const double J = 1.0 - l*st.ths;
        const double dsv = -(rx*(-cs) + rl*sn)/(std::fabs(J)>0.1 ? J : (J<0.0 ? -0.1 : 0.1));
        const double dlv = -(rx*sn + rl*cs);
        
        x -= dsv;
        l += dlv;
        
        // beyond the head or the tail: the end section
        x = std::min(xhead_ + 0.5*L_, std::max(xtail_ - 0.5*L_, x));
    }
    
    velocity_raw(x,l,vx,vl);
}

// ---------------------------------------------------------------------------------------------
// discrete recoil correction (X 198): a rigid motion in the plane of the undulation after the
// recoil, translation (tx, tl) and rotation psi about the CoG, and their rates; the 6DOF body
// integrates them so that the deformation carries no momentum and angular momentum in the hull
// region of the grid either
// ---------------------------------------------------------------------------------------------

void sixdof_deform::displacement(double x, double l, double &dx, double &dl) const
{
    displacement_raw(x,l,dx,dl);
    
    if(!corr_on_)
    return;
    
    const double a = x + dx, b = l + dl;
    const double cp = std::cos(cr_.psi), sp = std::sin(cr_.psi);
    
    dx = cp*a - sp*b + cr_.tx - x;
    dl = sp*a + cp*b + cr_.tl - l;
}

void sixdof_deform::velocity(double x, double l, double &vx, double &vl) const
{
    velocity_raw(x,l,vx,vl);
    
    if(!corr_on_)
    return;
    
    double dx, dl;
    displacement_raw(x,l,dx,dl);
    correct_velocity(x+dx,l+dl,vx,vl);
}

void sixdof_deform::velocity_at(double px, double pl, double &vx, double &vl) const
{
    if(!corr_on_)
    {
    velocity_at_raw(px,pl,vx,vl);
    return;
    }
    
    // the point without the correction: X = Rot(-psi)(p - t)
    const double cp = std::cos(cr_.psi), sp = std::sin(cr_.psi);
    const double qa = px - cr_.tx, qb = pl - cr_.tl;
    const double a =  cp*qa + sp*qb;
    const double b = -sp*qa + cp*qb;
    
    velocity_at_raw(a,b,vx,vl);
    correct_velocity(a,b,vx,vl);
}

void sixdof_deform::correct_velocity(double a, double b, double &vx, double &vl) const
{
    // d/dt (Rot(psi) X + t) = Rot(psi) dX/dt + psi_t e_n x Rot(psi) X + dt/dt, X = (a, b) the
    // position before the correction
    const double cp = std::cos(cr_.psi), sp = std::sin(cr_.psi);
    const double ra = cp*a - sp*b, rb = sp*a + cp*b;
    const double ua = cp*vx - sp*vl, ub = sp*vx + cp*vl;
    
    vx = ua - cr_.psit*rb + cr_.vx;
    vl = ub + cr_.psit*ra + cr_.vl;
}
