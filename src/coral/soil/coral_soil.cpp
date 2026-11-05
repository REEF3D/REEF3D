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

#include"coral/soil/coral_soil.h"
#include"coral/core/coral_model.h"
#include"coral/core/coral_rotation.h"
#include<algorithm>
#include<cmath>
#include<memory>
#include<stdexcept>

namespace coral
{

double soil_curve::force(double yy) const
{
    if(empty() || (compression_only && yy <= 0.0))
        return 0.0;
    const double s = yy < 0.0 ? -1.0 : 1.0, ya = std::fabs(yy);
    if(ya >= y.back())
        return s*p.back();
    const size_t i = std::upper_bound(y.begin(), y.end(), ya) - y.begin();     // y[i-1] <= ya < y[i]
    const double t = (ya - y[i-1])/(y[i] - y[i-1]);
    return s*(p[i-1] + t*(p[i] - p[i-1]));
}

double soil_curve::slope(double yy) const
{
    if(empty() || (compression_only && yy < 0.0))
        return 0.0;
    const double ya = std::fabs(yy);
    if(ya >= y.back())
        return 0.0;
    size_t i = std::upper_bound(y.begin(), y.end(), ya) - y.begin();
    if(i < 1) i = 1;
    return (p[i] - p[i-1])/(y[i] - y[i-1]);
}

soil_curve soil_curve::scaled(double s) const
{
    soil_curve c = *this;
    for(double &v : c.p) v *= s;
    return c;
}

soil_curve soil_curve::linear(double k, double ymax)
{
    soil_curve c;
    c.y = {0.0, ymax};
    c.p = {0.0, k*ymax};
    return c;
}

double soil_layer::su(double depth) const
{
    if(su_bottom < 0.0 || bottom <= top)
        return su_top;
    const double t = std::min(1.0, std::max(0.0, (depth - top)/(bottom - top)));
    return su_top + t*(su_bottom - su_top);
}

const soil_layer *soil_profile::layer_at(double depth) const
{
    for(const auto &l : layers)
        if(depth >= l.top && depth <= l.bottom)
            return &l;
    // below the last layer: the deepest one
    const soil_layer *d = nullptr;
    for(const auto &l : layers)
        if(!d || l.bottom > d->bottom) d = &l;
    return depth > 0.0 ? d : nullptr;
}

double soil_profile::effective_stress(double depth) const
{
    double s = 0.0;
    for(const auto &l : layers)
    {
        const double a = std::max(0.0, l.top), b = std::min(depth, l.bottom);
        if(b > a) s += l.gamma*(b - a);
    }
    return s;
}

soil_curve soil_profile::py(double X, double D) const
{
    soil_curve c;
    const soil_layer *l = layer_at(X);
    if(!l)
        return c;
    if(l->type==soil_layer::kind::linear)
        return l->kh > 0.0 ? soil_curve::linear(l->kh) : c;
    if(X <= 0.0)
        return c;
    const double p0 = effective_stress(X);
    if(l->type==soil_layer::kind::sand)
    {
        if(l->k <= 0.0)
            throw std::runtime_error("coral::soil " + name + ": sand layer needs the modulus k for p-y curves");
        const double phi = l->phi*M_PI/180.0, al = 0.5*phi, be = 0.25*M_PI + 0.5*phi;
        const double K0 = 0.4, Ka = std::pow(std::tan(0.25*M_PI - 0.5*phi), 2);
        const double tb = std::tan(be), tp = std::tan(phi), tbp = std::tan(be - phi), ta = std::tan(al);
        const double pst = p0*(K0*X*tp*std::sin(be)/(tbp*std::cos(al)) + tb/tbp*(D + X*tb*ta) + K0*X*tb*(tp*std::sin(be) - ta) - Ka*D);
        const double psd = Ka*D*p0*(std::pow(tb, 8) - 1.0) + K0*D*p0*tp*std::pow(tb, 4);
        const double pu = std::max(0.0, std::min(pst, psd));
        const double A = l->cyclic ? 0.9 : std::max(0.9, 3.0 - 0.8*X/D);
        const double k0 = l->k*X;
        if(pu <= 0.0 || k0 <= 0.0)
            return c;
        const double yr = A*pu/k0;
        const double ts[] = {0.0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8, 2.2, 2.6, 3.0, 3.5, 4.0, 5.0};
        for(double t : ts)
        {
            c.y.push_back(t*yr);
            c.p.push_back(A*pu*std::tanh(t));
        }
        return c;
    }
    // soft clay (Matlock)
    const double su = l->su(X), yc = 2.5*l->eps50*D;
    if(su <= 0.0)
        return c;
    const double pu = std::min(3.0*su + p0 + l->J*su*X/D, 9.0*su)*D;
    const double XR = 6.0*D/(l->gamma*D/su + l->J);
    const double r[] = {0.01, 0.02, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    c.y.push_back(0.0);
    c.p.push_back(0.0);
    for(double t : r)
    {
        if(l->cyclic && t > 3.0) break;
        c.y.push_back(t*yc);
        c.p.push_back(std::min(1.0, 0.5*std::cbrt(t))*pu);
    }
    if(l->cyclic)
    {
        if(X >= XR)
        {
            c.y.push_back(15.0*yc);
            c.p.push_back(0.72*pu);
        }
        else
        {
            c.y.push_back(15.0*yc);
            c.p.push_back(0.72*pu*X/XR);
        }
    }
    else
    {
        c.y.push_back(15.0*yc);
        c.p.push_back(pu);
    }
    return c;
}

soil_curve soil_profile::tz(double X, double D) const
{
    soil_curve c;
    const soil_layer *l = layer_at(X);
    if(!l)
        return c;
    if(l->type==soil_layer::kind::linear)
        return l->kv > 0.0 ? soil_curve::linear(l->kv) : c;
    if(X <= 0.0)
        return c;
    const double p0 = effective_stress(X), per = M_PI*D;
    if(l->type==soil_layer::kind::sand)
    {
        const double delta = (l->delta < 0.0 ? l->phi - 5.0 : l->delta)*M_PI/180.0;
        const double tmax = l->K*p0*std::tan(delta)*per;
        if(tmax <= 0.0) return c;
        c.y = {0.0, 0.00254};
        c.p = {0.0, tmax};
        return c;
    }
    const double su = l->su(X);
    const double psi = p0 > 0.0 ? su/p0 : 1.0e30;
    const double alpha = std::min(1.0, psi <= 1.0 ? 0.5/std::sqrt(psi) : 0.5*std::pow(psi, -0.25));
    const double tmax = alpha*su*per;
    if(tmax <= 0.0) return c;
    const double z[] = {0.0, 0.0016, 0.0031, 0.0057, 0.0080, 0.0100, 0.0200};
    const double t[] = {0.0, 0.30, 0.50, 0.75, 0.90, 1.00, l->residual};
    for(int i=0; i<7; ++i)
    {
        c.y.push_back(z[i]*D);
        c.p.push_back(t[i]*tmax);
    }
    return c;
}

soil_curve soil_profile::qz(double X, double D) const
{
    soil_curve c;
    c.compression_only = true;
    const soil_layer *l = layer_at(X);
    if(!l || X <= 0.0)
        return c;
    if(l->type==soil_layer::kind::linear)
    {
        if(l->kq > 0.0) {c = soil_curve::linear(l->kq); c.compression_only = true;}
        return c;
    }
    const double q = l->type==soil_layer::kind::sand ? l->Nq*effective_stress(X) : 9.0*l->su(X);
    const double Qp = q*0.25*M_PI*D*D;
    if(Qp <= 0.0) return c;
    const double z[] = {0.0, 0.002, 0.013, 0.042, 0.073, 0.100};
    const double r[] = {0.0, 0.25, 0.50, 0.75, 0.90, 1.00};
    for(int i=0; i<6; ++i)
    {
        c.y.push_back(z[i]*D);
        c.p.push_back(r[i]*Qp);
    }
    return c;
}

soil_spring::soil_spring(const model &m, int b, const Vec3 &axis, const soil_curve &lateral, const soil_curve &axial,
                         const soil_curve &tip, const std::string &name, double radius)
    : name_(name), b_(b), a_(axis.normalized()), x0_(m.body(b).x), R0_(m.body(b).R), r_(radius), lat_(lateral), ax_(axial), tip_(tip)
{
}

Vec3 soil_spring::force(const model &m) const
{
    const Vec3 u = m.body(b_).x - x0_;
    const double ua = a_.dot(u);
    const Vec3 ul = u - ua*a_;
    const double y = ul.norm();
    Vec3 F = -(ax_.force(ua) + tip_.force(ua))*a_;
    if(y > 0.0)
        F -= lat_.force(y)/y*ul;
    return F;
}

void soil_spring::eval(const model &m, double, VecX &f) const
{
    f.segment<3>(6*b_) += force(m);
    if(r_ > 0.0 && !ax_.empty())
    {
        // twist about the axis (small rotations): torque -t(r phi) r about the axis
        const rigid_body &B = m.body(b_);
        const double phi = a_.dot(log_so3(B.R*R0_.transpose()));
        m.add_moment(f, b_, -ax_.force(r_*phi)*r_*a_);
    }
}

void soil_spring::tangent(const model &m, double, double cK, double, double, std::vector<triplet> &T) const
{
    if(cK==0.0)
        return;
    const Vec3 u = m.body(b_).x - x0_;
    const double ua = a_.dot(u);
    const Vec3 ul = u - ua*a_;
    const double y = ul.norm();
    const Mat3 P = Mat3::Identity() - a_*a_.transpose();
    Mat3 K;
    if(y > 1.0e-12)
    {
        const Vec3 n = ul/y;
        const Mat3 N = n*n.transpose();
        K = -(lat_.slope(y)*N + lat_.force(y)/y*(P - N));
    }
    else
        K = -lat_.slope(0.0)*P;
    K -= (ax_.slope(ua) + tip_.slope(ua))*a_*a_.transpose();
    for(int j=0; j<3; ++j)
        for(int i=0; i<3; ++i)
            if(K(i, j) != 0.0)
                T.emplace_back(6*b_ + i, 6*b_ + j, cK*K(i, j));
    if(r_ > 0.0 && !ax_.empty())
    {
        // moment (body frame) R^T (-t(r phi) r a), phi = a.log(R R0^T): central differences in the rotation
        const rigid_body &B = m.body(b_);
        auto torque = [&](const Mat3 &R)
        {
            const double phi = a_.dot(log_so3(R*R0_.transpose()));
            return Vec3(R.transpose()*(-ax_.force(r_*phi)*r_*a_));
        };
        Mat3 Kr;
        const double e = 1.0e-7;
        for(int j=0; j<3; ++j)
            Kr.col(j) = (torque(B.R*exp_so3(e*Vec3::Unit(j))) - torque(B.R*exp_so3(-e*Vec3::Unit(j))))/(2.0*e);
        for(int j=0; j<3; ++j)
            for(int i=0; i<3; ++i)
                if(Kr(i, j) != 0.0)
                    T.emplace_back(6*b_ + 3 + i, 6*b_ + 3 + j, cK*Kr(i, j));
    }
}

void soil_spring::channel_names(std::vector<std::string> &n) const
{
    n.push_back(name_ + ".y");
    n.push_back(name_ + ".p");
    n.push_back(name_ + ".z");
    n.push_back(name_ + ".t");
}

void soil_spring::channel_values(const model &m, double, std::vector<double> &v) const
{
    const Vec3 u = m.body(b_).x - x0_;
    const double ua = a_.dot(u);
    const double y = (u - ua*a_).norm();
    v.push_back(y);
    v.push_back(lat_.force(y));
    v.push_back(ua);
    v.push_back(ax_.force(ua) + tip_.force(ua));
}

int add_pile_springs(model &m, const soil_profile &s, const std::vector<int> &nodes, double D, bool py, bool tz,
                     bool qz, bool linear, const std::string &name)
{
    const int n = (int)nodes.size();
    if(n < 2)
        return 0;
    std::vector<Vec3> x(n);
    for(int i=0; i<n; ++i) x[i] = m.body(nodes[i]).x;
    // the tip is the deeper end; the axis points from the head to the tip
    const bool rev = x[0](2) < x[n-1](2);
    const Vec3 axis = (rev ? x[0] - x[n-1] : x[n-1] - x[0]).normalized();
    const int tipi = rev ? 0 : n-1;
    // length of the segment p-q below the mudline
    auto below = [&](const Vec3 &p, const Vec3 &q)
    {
        const double L = (q - p).norm(), zp = p(2), zq = q(2), zm = s.mudline;
        if(zp <= zm && zq <= zm) return L;
        if(zp >= zm && zq >= zm) return 0.0;
        const double f = (zm - std::min(zp, zq))/std::fabs(zq - zp);
        return f*L;
    };
    int added = 0;
    for(int i=0; i<n; ++i)
    {
        double Lt = 0.0;
        if(i > 0) Lt += below(x[i], 0.5*(x[i] + x[i-1]));
        if(i < n-1) Lt += below(x[i], 0.5*(x[i] + x[i+1]));
        const double X = std::max(0.0, s.mudline - x[i](2));
        if(Lt <= 0.0 && !(qz && i==tipi && X > 0.0))
            continue;
        soil_curve lat, ax, tip;
        if(py) lat = s.py(X, D).scaled(Lt);
        if(tz) ax = s.tz(X, D).scaled(Lt);
        if(qz && i==tipi) tip = s.qz(X, D);
        if(linear)
        {
            if(!lat.empty()) lat = soil_curve::linear(lat.slope(0.0));
            if(!ax.empty()) ax = soil_curve::linear(ax.slope(0.0));
            if(!tip.empty()) {const double k = tip.slope(0.0); tip = soil_curve::linear(k); tip.compression_only = true;}
        }
        if(lat.empty() && ax.empty() && tip.empty())
            continue;
        m.add_load(std::make_unique<soil_spring>(m, nodes[i], axis, lat, ax, tip, name + ":" + std::to_string(i), 0.5*D));
        ++added;
    }
    return added;
}

}
