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

#include"coral/structure/coral_superelement.h"
#include"coral/core/coral_model.h"
#include"coral/core/coral_rotation.h"
#include<algorithm>
#include<cmath>
#include<map>
#include<stdexcept>

namespace coral
{

Eigen::Matrix<double,12,12> beam_stiffness_local(const beam_section &s, double L)
{
    Eigen::Matrix<double,12,12> K = Eigen::Matrix<double,12,12>::Zero();
    auto set = [&K](int i, int j, double v) {K(i,j) = v; K(j,i) = v;};
    set(0, 0, s.EA/L); set(6, 6, s.EA/L); set(0, 6, -s.EA/L);
    set(3, 3, s.GJ/L); set(9, 9, s.GJ/L); set(3, 9, -s.GJ/L);
    // x-y plane: v, theta_z with EIz and shear GAy
    const double pz = s.GAy > 0.0 ? 12.0*s.EIz/(s.GAy*L*L) : 0.0, a = s.EIz/((1.0 + pz)*L*L*L);
    set(1, 1, 12.0*a); set(1, 5, 6.0*L*a); set(1, 7, -12.0*a); set(1, 11, 6.0*L*a);
    set(5, 5, (4.0 + pz)*L*L*a); set(5, 7, -6.0*L*a); set(5, 11, (2.0 - pz)*L*L*a);
    set(7, 7, 12.0*a); set(7, 11, -6.0*L*a); set(11, 11, (4.0 + pz)*L*L*a);
    // x-z plane: w, theta_y with EIy and shear GAz
    const double py = s.GAz > 0.0 ? 12.0*s.EIy/(s.GAz*L*L) : 0.0, b = s.EIy/((1.0 + py)*L*L*L);
    set(2, 2, 12.0*b); set(2, 4, -6.0*L*b); set(2, 8, -12.0*b); set(2, 10, -6.0*L*b);
    set(4, 4, (4.0 + py)*L*L*b); set(4, 8, 6.0*L*b); set(4, 10, (2.0 - py)*L*L*b);
    set(8, 8, 12.0*b); set(8, 10, 6.0*L*b); set(10, 10, (4.0 + py)*L*L*b);
    return K;
}

Eigen::Matrix<double,12,12> beam_mass_local(const beam_section &s, double L)
{
    Eigen::Matrix<double,12,12> M = Eigen::Matrix<double,12,12>::Zero();
    auto add = [&M](int i, int j, double v) {M(i,j) += v; if(i!=j) M(j,i) += v;};
    const double c = s.m*L/420.0;
    add(0, 0, 140.0*c); add(6, 6, 140.0*c); add(0, 6, 70.0*c);
    add(3, 3, s.rIx*L/3.0); add(9, 9, s.rIx*L/3.0); add(3, 9, s.rIx*L/6.0);
    // translation (cubic Hermite)
    add(1, 1, 156.0*c); add(1, 5, 22.0*L*c); add(1, 7, 54.0*c); add(1, 11, -13.0*L*c);
    add(5, 5, 4.0*L*L*c); add(5, 7, 13.0*L*c); add(5, 11, -3.0*L*L*c);
    add(7, 7, 156.0*c); add(7, 11, -22.0*L*c); add(11, 11, 4.0*L*L*c);
    add(2, 2, 156.0*c); add(2, 4, -22.0*L*c); add(2, 8, 54.0*c); add(2, 10, 13.0*L*c);
    add(4, 4, 4.0*L*L*c); add(4, 8, -13.0*L*c); add(4, 10, -3.0*L*L*c);
    add(8, 8, 156.0*c); add(8, 10, 22.0*L*c); add(10, 10, 4.0*L*L*c);
    // rotary inertia of the section
    const double rz = s.rIz/(30.0*L), ry = s.rIy/(30.0*L);
    add(1, 1, 36.0*rz); add(1, 5, 3.0*L*rz); add(1, 7, -36.0*rz); add(1, 11, 3.0*L*rz);
    add(5, 5, 4.0*L*L*rz); add(5, 7, -3.0*L*rz); add(5, 11, -L*L*rz);
    add(7, 7, 36.0*rz); add(7, 11, -3.0*L*rz); add(11, 11, 4.0*L*L*rz);
    add(2, 2, 36.0*ry); add(2, 4, -3.0*L*ry); add(2, 8, -36.0*ry); add(2, 10, -3.0*L*ry);
    add(4, 4, 4.0*L*L*ry); add(4, 8, 3.0*L*ry); add(4, 10, -L*L*ry);
    add(8, 8, 36.0*ry); add(8, 10, 3.0*L*ry); add(10, 10, 4.0*L*L*ry);
    return M;
}

namespace
{
    // local axes of a member (columns x, y, z), as the co-rotational element
    Mat3 member_axes(const Vec3 &p1, const Vec3 &p2, const Vec3 &orient)
    {
        const Vec3 ex = (p2 - p1).normalized();
        Vec3 ey = orient - orient.dot(ex)*ex;
        if(ey.norm() < 1.0e-6)
        {
            const Vec3 alt = std::fabs(ex(2)) < 0.9 ? Vec3::UnitZ() : Vec3::UnitX();
            ey = alt - alt.dot(ex)*ex;
        }
        ey.normalize();
        Mat3 E;
        E.col(0) = ex; E.col(1) = ey; E.col(2) = ex.cross(ey);
        return E;
    }
}

superelement::superelement(model &m, const std::string &name, const std::vector<member_spec> &members, int nmode,
                           double zeta, double beta, const Vec3 &g, double tol)
    : name_(name)
{
    // nodes (merged by position) and their attachment
    std::vector<int> attach;      // -2 free, -1 ground, >= 0 body
    auto node_at = [&](const Vec3 &p)
    {
        for(size_t k=0; k<p0_.size(); ++k)
            if((p0_[k] - p).norm() <= tol*(1.0 + p.norm()))
                return (int)k;
        p0_.push_back(p);
        attach.push_back(-2);
        return (int)p0_.size() - 1;
    };
    auto set_attach = [&](int k, int body, const std::string &mem)
    {
        if(body==-2) return;
        if(attach[k]==-1 || body==-1) {attach[k] = -1; return;}
        if(attach[k] >= 0 && attach[k]!=body)
            throw std::runtime_error("coral::superelement " + name + ": node of member " + mem + " on two bodies");
        attach[k] = body;
    };
    struct elem {int a, b; beam_section s; Mat3 E; double L;};
    std::vector<elem> els;
    for(const auto &ms : members)
    {
        const int n = std::max(1, ms.n);
        std::vector<int> st(n+1);
        for(int i=0; i<=n; ++i)
            st[i] = node_at(ms.p1 + (double)i/n*(ms.p2 - ms.p1));
        mnodes_[ms.name] = st;
        mdiam_[ms.name] = ms.s1.D;
        set_attach(st[0], ms.body1, ms.name);
        set_attach(st[n], ms.body2, ms.name);
        const Mat3 E = member_axes(ms.p1, ms.p2, ms.orient);
        for(int e=0; e<n; ++e)
        {
            elem el;
            el.a = st[e]; el.b = st[e+1];
            el.s = ms.tapered ? beam_section::blend(ms.s1, ms.s2, (e + 0.5)/n) : ms.s1;
            el.E = E;
            el.L = (p0_[el.b] - p0_[el.a]).norm();
            if(el.L <= 0.0)
                throw std::runtime_error("coral::superelement " + name + ": zero length element in member " + ms.name);
            els.push_back(el);
            elem_.emplace_back(el.a, el.b);
        }
    }
    const int nn = (int)p0_.size(), N = 6*nn;
    // global matrices and gravity
    MatX K = MatX::Zero(N, N), M = MatX::Zero(N, N);
    VecX F = VecX::Zero(N);
    for(const auto &el : els)
    {
        Eigen::Matrix<double,12,12> Lm = Eigen::Matrix<double,12,12>::Zero();
        for(int q=0; q<4; ++q)
            Lm.block<3,3>(3*q, 3*q) = el.E.transpose();
        const Eigen::Matrix<double,12,12> Ke = Lm.transpose()*beam_stiffness_local(el.s, el.L)*Lm;
        const Eigen::Matrix<double,12,12> Me = Lm.transpose()*beam_mass_local(el.s, el.L)*Lm;
        // consistent gravity: distributed load w = m g in local axes
        const Vec3 w = el.E.transpose()*(el.s.m*g);
        const double L = el.L;
        Eigen::Matrix<double,12,1> fl = Eigen::Matrix<double,12,1>::Zero();
        fl(0) = fl(6) = 0.5*w(0)*L;
        fl(1) = fl(7) = 0.5*w(1)*L;
        fl(5) = w(1)*L*L/12.0; fl(11) = -w(1)*L*L/12.0;
        fl(2) = fl(8) = 0.5*w(2)*L;
        fl(4) = -w(2)*L*L/12.0; fl(10) = w(2)*L*L/12.0;
        const Eigen::Matrix<double,12,1> fg = Lm.transpose()*fl;
        const int idx[2] = {el.a, el.b};
        for(int i=0; i<2; ++i)
        {
            F.segment<6>(6*idx[i]) += fg.segment<6>(6*i);
            for(int j=0; j<2; ++j)
            {
                K.block<6,6>(6*idx[i], 6*idx[j]) += Ke.block<6,6>(6*i, 6*j);
                M.block<6,6>(6*idx[i], 6*idx[j]) += Me.block<6,6>(6*i, 6*j);
            }
        }
    }
    // interface bodies, interior dofs
    kind_.assign(nn, -2);
    idof_.assign(nn, -1);
    std::map<int,int> bidx;
    for(int k=0; k<nn; ++k)
    {
        if(attach[k]==-1) kind_[k] = -1;
        else if(attach[k] >= 0)
        {
            if(!bidx.count(attach[k]))
            {
                bidx[attach[k]] = (int)bodies_.size();
                bodies_.push_back(attach[k]);
                xb0_.push_back(m.bodies[attach[k]].x);
                Rb0_.push_back(m.bodies[attach[k]].R);
            }
            kind_[k] = bidx[attach[k]];
        }
        else
        {
            idof_[k] = nI_;
            nI_ += 6;
        }
    }
    const int nR = 6*(int)bodies_.size();
    // T: interface body motion -> all node dofs (interface nodes only); S: interior dof selection
    MatX T = MatX::Zero(N, nR), S = MatX::Zero(N, nI_);
    for(int k=0; k<nn; ++k)
    {
        if(kind_[k] >= 0)
        {
            const int b = kind_[k];
            const Vec3 r = p0_[k] - xb0_[b];
            T.block<3,3>(6*k, 6*b) = Mat3::Identity();
            T.block<3,3>(6*k, 6*b+3) = -skew(r);
            T.block<3,3>(6*k+3, 6*b+3) = Mat3::Identity();
        }
        else if(kind_[k]==-2)
            S.block<6,6>(6*k, idof_[k]) = Eigen::Matrix<double,6,6>::Identity();
    }
    const MatX KRR = T.transpose()*K*T, KRI = T.transpose()*K*S, KII = S.transpose()*K*S;
    const MatX MRR = T.transpose()*M*T, MRI = T.transpose()*M*S, MII = S.transpose()*M*S;
    const VecX FRf = T.transpose()*F, FI = S.transpose()*F;
    nm_ = 0;
    if(nI_ > 0)
    {
        const Eigen::LDLT<MatX> kii(KII);
        if(kii.info()!=Eigen::Success)
            throw std::runtime_error("coral::superelement " + name + ": interior stiffness singular (unsupported part)");
        Phi_ = -kii.solve(KRI.transpose());
        nm_ = std::max(0, std::min(nmode, nI_));
        Eigen::GeneralizedSelfAdjointEigenSolver<MatX> es(KII, MII);
        if(es.info()!=Eigen::Success)
            throw std::runtime_error("coral::superelement " + name + ": fixed-interface eigenproblem failed");
        Psi_ = es.eigenvectors().leftCols(nm_);
        om2_ = es.eigenvalues().head(nm_);
    }
    else
    {
        Phi_ = MatX::Zero(0, nR);
        Psi_ = MatX::Zero(0, 0);
        om2_ = VecX::Zero(0);
    }
    KR_ = KRR + (nI_ ? MatX(KRI*Phi_) : MatX::Zero(nR, nR));
    MR_ = MRR;
    if(nI_)
        MR_ += MRI*Phi_ + Phi_.transpose()*MRI.transpose() + Phi_.transpose()*MII*Phi_;
    MRm_ = nI_ ? MatX((MRI + Phi_.transpose()*MII)*Psi_) : MatX::Zero(nR, 0);
    CR_ = beta*KR_;
    FR_ = FRf + (nI_ ? VecX(Phi_.transpose()*FI) : VecX::Zero(nR));
    Fm_ = nI_ ? VecX(Psi_.transpose()*FI) : VecX::Zero(0);
    fmode_.resize(nm_);
    cm_.resize(nm_);
    for(int k=0; k<nm_; ++k)
    {
        const double w = std::sqrt(std::max(om2_(k), 0.0));
        fmode_(k) = w/(2.0*M_PI);
        cm_(k) = 2.0*zeta*w;
    }
    c0_ = m.add_coordinates(nm_, name);
    // displacement per unit modal amplitude: largest translation of the mode shape
    for(int k=0; k<nm_; ++k)
    {
        double mx = 0.0;
        for(int i=0; i<nI_; ++i)
            if(i%6 < 3) mx = std::max(mx, std::fabs(Psi_(i, k)));
        m.cscale(c0_ + k) = mx > 0.0 ? mx : 1.0;
    }
}

void superelement::gather(const model &m, std::vector<bstate> &s) const
{
    s.resize(bodies_.size());
    for(size_t i=0; i<bodies_.size(); ++i)
    {
        const rigid_body &b = m.body(bodies_[i]);
        s[i].x = b.x; s[i].R = b.R; s[i].v = b.v; s[i].w = b.w; s[i].a = b.vd; s[i].wd = b.wd;
    }
}

void superelement::interface_motion(const std::vector<bstate> &s, VecX &u, VecX &v, VecX &a) const
{
    const int nR = 6*(int)bodies_.size();
    u.resize(nR); v.resize(nR); a.resize(nR);
    for(size_t i=0; i<s.size(); ++i)
    {
        u.segment<3>(6*i) = s[i].x - xb0_[i];
        u.segment<3>(6*i+3) = log_so3(s[i].R*Rb0_[i].transpose());
        v.segment<3>(6*i) = s[i].v;
        v.segment<3>(6*i+3) = s[i].R*s[i].w;
        a.segment<3>(6*i) = s[i].a;
        a.segment<3>(6*i+3) = s[i].R*s[i].wd;
    }
}

void superelement::compute(const std::vector<bstate> &s, const VecX &q, const VecX &qd, const VecX &qdd, VecX &fb, VecX &fq) const
{
    VecX u, v, a;
    interface_motion(s, u, v, a);
    VecX fR = FR_ - KR_*u - CR_*v - MR_*a;
    if(nm_)
    {
        fR -= MRm_*qdd;
        fq = Fm_ - om2_.cwiseProduct(q) - cm_.cwiseProduct(qd) - qdd - MRm_.transpose()*a;
    }
    else
        fq.resize(0);
    fb = fR;
    for(size_t i=0; i<s.size(); ++i)
        fb.segment<3>(6*i+3) = s[i].R.transpose()*fR.segment<3>(6*i+3);
}

void superelement::eval(const model &m, double, VecX &f) const
{
    std::vector<bstate> s;
    gather(m, s);
    VecX fb, fq;
    compute(s, m.cq.segment(c0_, nm_), m.cv.segment(c0_, nm_), m.ca.segment(c0_, nm_), fb, fq);
    for(size_t i=0; i<bodies_.size(); ++i)
        f.segment<6>(6*bodies_[i]) += fb.segment<6>(6*i);
    if(nm_)
        f.segment(m.coord_dof(c0_), nm_) += fq;
}

void superelement::tangent(const model &m, double, double cK, double cC, double cM, std::vector<triplet> &Tr) const
{
    const int nb = (int)bodies_.size(), nR = 6*nb;
    std::vector<bstate> s0;
    gather(m, s0);
    const VecX q = m.cq.segment(c0_, nm_), qd = m.cv.segment(c0_, nm_), qdd = m.ca.segment(c0_, nm_);
    // rows: interface bodies (6 each), then the modes; columns likewise
    auto row = [&](int i) {return i < nR ? 6*bodies_[i/6] + i%6 : m.coord_dof(c0_ + i - nR);};
    MatX J = MatX::Zero(nR + nm_, nR + nm_);
    VecX fbp, fqp, fbm, fqm;
    // body columns: central differences through the (nonlinear) map of the body states
    const double e = 1.0e-6;
    for(int j=0; j<nR; ++j)
    {
        const int i = j/6, c = j%6;
        for(int kind=0; kind<3; ++kind)
        {
            const double cf = kind==0 ? cK : kind==1 ? cC : cM;
            if(cf==0.0) continue;
            std::vector<bstate> sp = s0, sm = s0;
            for(int sg=-1; sg<=1; sg+=2)
            {
                bstate &b = (sg < 0 ? sm : sp)[i];
                if(kind==0)
                {
                    if(c < 3) b.x(c) += sg*e; else b.R = b.R*exp_so3(sg*e*Vec3::Unit(c-3));
                }
                else if(kind==1)
                {
                    if(c < 3) b.v(c) += sg*e; else b.w(c-3) += sg*e;
                }
                else
                {
                    if(c < 3) b.a(c) += sg*e; else b.wd(c-3) += sg*e;
                }
            }
            compute(sp, q, qd, qdd, fbp, fqp);
            compute(sm, q, qd, qdd, fbm, fqm);
            J.block(0, j, nR, 1) += cf*(fbp - fbm)/(2.0*e);
            if(nm_) J.block(nR, j, nm_, 1) += cf*(fqp - fqm)/(2.0*e);
        }
    }
    // mode columns: exact (forces linear in q, q', q'')
    for(int k=0; k<nm_; ++k)
    {
        J(nR + k, nR + k) += -cK*om2_(k) - cC*cm_(k) - cM;
        if(cM != 0.0)
        {
            // force on the bodies: -M_Rm q'' (moments to the body frame)
            for(int i=0; i<nb; ++i)
            {
                J.block<3,1>(6*i, nR + k) += -cM*MRm_.block<3,1>(6*i, k);
                J.block<3,1>(6*i+3, nR + k) += -cM*s0[i].R.transpose()*MRm_.block<3,1>(6*i+3, k);
            }
        }
    }
    for(int j=0; j<nR + nm_; ++j)
        for(int i=0; i<nR + nm_; ++i)
            if(J(i, j) != 0.0)
                Tr.emplace_back(row(i), row(j), J(i, j));
}

double superelement::potential(const model &m) const
{
    std::vector<bstate> s;
    gather(m, s);
    VecX u, v, a;
    interface_motion(s, u, v, a);
    const VecX q = m.cq.segment(c0_, nm_);
    double e = 0.5*u.dot(KR_*u) - FR_.dot(u);
    if(nm_) e += 0.5*q.dot(om2_.cwiseProduct(q)) - Fm_.dot(q);
    return e;
}

void superelement::channel_names(std::vector<std::string> &n) const
{
    for(int k=0; k<nm_; ++k)
        n.push_back(name_ + ".q" + std::to_string(k+1));
}

void superelement::channel_values(const model &m, double, std::vector<double> &v) const
{
    for(int k=0; k<nm_; ++k)
        v.push_back(m.cq(c0_ + k));
}

MatX superelement::reduced_mass() const
{
    const int nR = (int)KR_.rows();
    MatX M = MatX::Zero(nR + nm_, nR + nm_);
    M.topLeftCorner(nR, nR) = MR_;
    if(nm_)
    {
        M.topRightCorner(nR, nm_) = MRm_;
        M.bottomLeftCorner(nm_, nR) = MRm_.transpose();
        M.bottomRightCorner(nm_, nm_) = MatX::Identity(nm_, nm_);
    }
    return M;
}

MatX superelement::reduced_stiffness() const
{
    const int nR = (int)KR_.rows();
    MatX K = MatX::Zero(nR + nm_, nR + nm_);
    K.topLeftCorner(nR, nR) = KR_;
    if(nm_)
        K.bottomRightCorner(nm_, nm_) = om2_.asDiagonal();
    return K;
}

VecX superelement::reduced_gravity() const
{
    VecX F(FR_.size() + Fm_.size());
    F << FR_, Fm_;
    return F;
}

Eigen::Matrix<double,6,1> superelement::node_displacement(const model &m, int k) const
{
    Eigen::Matrix<double,6,1> d = Eigen::Matrix<double,6,1>::Zero();
    if(kind_[k]==-1)
        return d;
    std::vector<bstate> s;
    gather(m, s);
    VecX u, v, a;
    interface_motion(s, u, v, a);
    if(kind_[k] >= 0)
    {
        const int b = kind_[k];
        const Vec3 r = p0_[k] - xb0_[b];
        d.head<3>() = u.segment<3>(6*b) + u.segment<3>(6*b+3).cross(r);
        d.tail<3>() = u.segment<3>(6*b+3);
        return d;
    }
    VecX uI = Phi_.middleRows(idof_[k], 6)*u;
    if(nm_)
        uI += Psi_.middleRows(idof_[k], 6)*m.cq.segment(c0_, nm_);
    return uI;
}

std::vector<int> superelement::member_nodes(const std::string &member) const
{
    const auto it = mnodes_.find(member);
    return it==mnodes_.end() ? std::vector<int>() : it->second;
}

double superelement::member_diameter(const std::string &member) const
{
    const auto it = mdiam_.find(member);
    return it==mdiam_.end() ? 0.0 : it->second;
}

void superelement::node_state(const model &m, const VecX &u, const VecX &vr, const VecX &ar, int k, Vec3 &p, Vec3 &v, Vec3 &a) const
{
    p = p0_[k];
    v.setZero();
    a.setZero();
    if(kind_[k]==-1)
        return;
    if(kind_[k] >= 0)
    {
        const int b = kind_[k];
        const Vec3 r = p0_[k] - xb0_[b];
        p += u.segment<3>(6*b) + u.segment<3>(6*b+3).cross(r);
        v = vr.segment<3>(6*b) + vr.segment<3>(6*b+3).cross(r);
        a = ar.segment<3>(6*b) + ar.segment<3>(6*b+3).cross(r);
        return;
    }
    const int i = idof_[k];
    p += Phi_.middleRows(i, 3)*u;
    v = Phi_.middleRows(i, 3)*vr;
    a = Phi_.middleRows(i, 3)*ar;
    if(nm_)
    {
        p += Psi_.middleRows(i, 3)*m.cq.segment(c0_, nm_);
        v += Psi_.middleRows(i, 3)*m.cv.segment(c0_, nm_);
        a += Psi_.middleRows(i, 3)*m.ca.segment(c0_, nm_);
    }
}

void superelement::point_state(const model &m, int k, Vec3 &p, Vec3 &v, Vec3 &a) const
{
    std::vector<bstate> s;
    gather(m, s);
    VecX u, vr, ar;
    interface_motion(s, u, vr, ar);
    node_state(m, u, vr, ar, k, p, v, a);
}

void superelement::point_states(const model &m, const std::vector<int> &k, std::vector<Vec3> &p, std::vector<Vec3> &v,
                                std::vector<Vec3> &a) const
{
    std::vector<bstate> s;
    gather(m, s);
    VecX u, vr, ar;
    interface_motion(s, u, vr, ar);
    p.resize(k.size()); v.resize(k.size()); a.resize(k.size());
    for(size_t i=0; i<k.size(); ++i)
        node_state(m, u, vr, ar, k[i], p[i], v[i], a[i]);
}

void superelement::add_point_force(const model &m, int k, const Vec3 &F, VecX &f) const
{
    if(kind_[k]==-1)
        return;
    if(kind_[k] >= 0)
    {
        Vec3 p, v, a;
        point_state(m, k, p, v, a);
        m.add_force(f, bodies_[kind_[k]], p, F);
        return;
    }
    // interior node: generalized forces through the recovery u_I = Phi u_R + Psi q
    const int i = idof_[k];
    const VecX fR = Phi_.middleRows(i, 3).transpose()*F;
    for(size_t b=0; b<bodies_.size(); ++b)
    {
        const rigid_body &B = m.body(bodies_[b]);
        f.segment<3>(6*bodies_[b]) += fR.segment<3>(6*b);
        f.segment<3>(6*bodies_[b]+3) += B.R.transpose()*fR.segment<3>(6*b+3);
    }
    if(nm_)
        f.segment(m.coord_dof(c0_), nm_) += Psi_.middleRows(i, 3).transpose()*F;
}

void superelement::node_positions(const model &m, std::vector<Vec3> &p) const
{
    p.resize(p0_.size());
    for(size_t k=0; k<p0_.size(); ++k)
        p[k] = p0_[k] + node_displacement(m, (int)k).head<3>();
}

}
