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

#include"coral/linear/coral_modes.h"
#include"coral/core/coral_model.h"
#include"coral/core/coral_rotation.h"
#include"coral/core/coral_fdjac.h"
#include<algorithm>
#include<cmath>
#include<complex>
#include<numeric>
#include<vector>

namespace coral
{

linear_model linearize(model &m, double t)
{
    const int n = m.ndof(), nc = m.ncon(), nb = m.nbody(), nq = m.ncoord();
    linear_model L;
    const model::config c0 = m.save_config();
    std::vector<Vec3> v0, w0, vd0, wd0;
    for(const auto &b : m.bodies)
    {
        v0.push_back(b.v); w0.push_back(b.w); vd0.push_back(b.vd); wd0.push_back(b.wd);
    }
    const VecX cv0 = m.cv, ca0 = m.ca;
    m.zero_rates();
    const bool tang = m.has_tangents();
    VecX f, in, g;
    MatX G;
    m.forces(t, f);
    m.constraints(g, G);
    // multipliers of the static balance in this configuration: G^T lambda = f (least squares)
    VecX lam = VecX::Zero(nc);
    if(nc)
        lam = (G*G.transpose()).ldlt().solve(G*f);
    const VecX zero = VecX::Zero(n);
    auto dense = [&](const std::vector<triplet> &T)
    {
        MatX A = MatX::Zero(n, n);
        for(const auto &x : T)
            if(x.row() < n && x.col() < n)
                A(x.row(), x.col()) += x.value();
        return A;
    };
    // configuration: restoring residual -f + G^T lambda (lambda fixed), constraint rows appended
    std::vector<triplet> T;
    VecX dq = zero;
    fd_jacobian(m, T, VecX(),
        [&](int col, double e) {dq(col) = e;},
        [&](VecX &r)
        {
            m.set_config(c0, dq);
            m.forces(t, f, tang);
            m.constraints(g, G);
            r.resize(n + nc);
            r.head(n) = -f;
            if(nc) {r.head(n) += G.transpose()*lam; r.tail(nc) = g;}
        },
        [&](int col) {return 1.0e-6*(col < 6*nb ? ((col%6) < 3 ? std::max(1.0, c0.x[col/6].cwiseAbs().maxCoeff()) : 1.0)
                                                : std::max(1.0, std::fabs(c0.cq(col - 6*nb))));},
        true, tang);
    m.set_config(c0, zero);
    if(tang) m.tangents(t, -1.0, 0.0, 0.0, T);
    L.K = dense(T);
    // velocities and accelerations
    auto rates = [&](bool acc)
    {
        std::vector<triplet> T2;
        VecX du = zero;
        auto apply = [&]()
        {
            for(int i=0; i<nb; ++i)
            {
                rigid_body &b = m.bodies[i];
                (acc ? b.vd : b.v) = du.segment<3>(6*i);
                (acc ? b.wd : b.w) = du.segment<3>(6*i+3);
            }
            if(nq) (acc ? m.ca : m.cv) = du.tail(nq);
        };
        fd_jacobian(m, T2, VecX(),
            [&](int col, double e) {du(col) = e;},
            [&](VecX &r)
            {
                apply();
                m.forces(t, f, tang);
                r.setZero(n + nc);
                r.head(n) = -f;
                if(acc) {m.inertia_forces(in); r.head(n) += in;}
            },
            [&](int) {return 1.0e-6;},
            true, tang);
        du.setZero();
        apply();
        if(tang) m.tangents(t, 0.0, acc ? 0.0 : -1.0, acc ? -1.0 : 0.0, T2);
        return dense(T2);
    };
    L.C = rates(false);
    L.M = rates(true);
    m.set_config(c0, zero);
    for(int i=0; i<nb; ++i)
    {
        m.bodies[i].v = v0[i]; m.bodies[i].w = w0[i]; m.bodies[i].vd = vd0[i]; m.bodies[i].wd = wd0[i];
    }
    m.cv = cv0;
    m.ca = ca0;
    // null space of G
    if(nc)
    {
        Eigen::JacobiSVD<MatX> svd(G, Eigen::ComputeFullV);
        const VecX sv = svd.singularValues();
        int rank = 0;
        for(int q=0; q<sv.size(); ++q)
            if(sv(q) > 1.0e-10*sv(0)) ++rank;
        L.N = svd.matrixV().rightCols(n - rank);
    }
    else
        L.N = MatX::Identity(n, n);
    L.Mr = L.N.transpose()*L.M*L.N;
    L.Cr = L.N.transpose()*L.C*L.N;
    L.Kr = L.N.transpose()*L.K*L.N;
    return L;
}

mode_result modes(model &m, double t, int nmode, bool damped)
{
    mode_result res;
    const linear_model L = linearize(m, t);
    const int nr = (int)L.Mr.rows();
    if(nr==0)
        return res;
    nmode = std::min(nmode, nr);
    std::vector<double> fr, ze;
    std::vector<VecX> sh;
    const MatX Mi = L.Mr.ldlt().solve(MatX::Identity(nr, nr));
    const double asym = (L.Kr - L.Kr.transpose()).norm()/std::max(L.Kr.norm(), 1.0e-300);
    if(!damped && asym < 1.0e-5)
    {
        // (nearly) symmetric stiffness: symmetric generalized problem, M-orthogonal shapes
        Eigen::GeneralizedSelfAdjointEigenSolver<MatX> es(0.5*(L.Kr + L.Kr.transpose()), 0.5*(L.Mr + L.Mr.transpose()));
        for(int q=0; q<nr; ++q)
        {
            const double w2 = es.eigenvalues()(q);
            fr.push_back(std::sqrt(std::max(w2, 0.0))/(2.0*M_PI));
            ze.push_back(0.0);
            sh.push_back(es.eigenvectors().col(q));
        }
    }
    else if(!damped)
    {
        Eigen::EigenSolver<MatX> es(Mi*L.Kr);
        for(int q=0; q<nr; ++q)
        {
            const double w2 = es.eigenvalues()(q).real();
            fr.push_back(std::sqrt(std::max(w2, 0.0))/(2.0*M_PI));
            ze.push_back(0.0);
            sh.push_back(es.eigenvectors().col(q).real());
        }
    }
    else
    {
        MatX A = MatX::Zero(2*nr, 2*nr);
        A.topRightCorner(nr, nr) = MatX::Identity(nr, nr);
        A.bottomLeftCorner(nr, nr) = -Mi*L.Kr;
        A.bottomRightCorner(nr, nr) = -Mi*L.Cr;
        Eigen::EigenSolver<MatX> es(A);
        for(int q=0; q<2*nr; ++q)
        {
            const std::complex<double> l = es.eigenvalues()(q);
            if(l.imag() < 0.0) continue;    // one of each pair (and the real ones)
            const double a = std::abs(l);
            fr.push_back(a/(2.0*M_PI));
            ze.push_back(a > 0.0 ? -l.real()/a : 1.0);
            sh.push_back(es.eigenvectors().col(q).head(nr).real());
        }
    }
    std::vector<int> idx(fr.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) {return fr[a] < fr[b];});
    nmode = std::min(nmode, (int)idx.size());
    res.frequency.resize(nmode);
    res.damping.resize(nmode);
    res.shape.resize(L.N.rows(), nmode);
    // shapes in the full coordinates; in a cluster of equal frequencies (symmetric structures) the
    // basis is turned so that the modes separate by global direction: generalized eigenvectors of the
    // x-translation kinetic energy against the mass in the cluster's subspace (x first, then the rest)
    const int nall = (int)idx.size();
    MatX S(L.N.rows(), nall);
    for(int q=0; q<nall; ++q)
        S.col(q) = L.N*sh[idx[q]];
    if(!damped)
    {
        MatX Dx = MatX::Zero(L.M.rows(), L.M.cols());
        for(int b=0; b<m.nbody(); ++b)
            Dx(6*b, 6*b) = L.M(6*b, 6*b);
        for(int q=0; q<nmode; )
        {
            int e = q + 1;
            const double f0 = fr[idx[q]];
            while(e < nall && std::fabs(fr[idx[e]] - f0) <= 1.0e-6*std::max(f0, 1.0e-12)) ++e;
            if(e - q > 1)
            {
                const MatX C = S.middleCols(q, e - q);
                const MatX Mc = C.transpose()*L.M*C, Qx = C.transpose()*Dx*C;
                Eigen::GeneralizedSelfAdjointEigenSolver<MatX> ge(0.5*(Qx + Qx.transpose()), 0.5*(Mc + Mc.transpose()));
                if(ge.info()==Eigen::Success)
                {
                    const MatX V = ge.eigenvectors().rowwise().reverse();
                    S.middleCols(q, e - q) = C*V;
                }
            }
            q = e;
        }
    }
    for(int q=0; q<nmode; ++q)
    {
        res.frequency(q) = fr[idx[q]];
        res.damping(q) = ze[idx[q]];
        const VecX s = S.col(q);
        const double mx = s.cwiseAbs().maxCoeff();
        res.shape.col(q) = mx > 0.0 ? VecX(s/mx) : s;
    }
    res.ok = true;
    return res;
}

}
