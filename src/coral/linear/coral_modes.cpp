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
#include<algorithm>
#include<cmath>
#include<complex>
#include<numeric>
#include<vector>

namespace coral
{

linear_model linearize(model &m, double t)
{
    const int n = m.ndof(), nc = m.ncon();
    linear_model L;
    std::vector<Vec3> x0, v0, w0, vd0, wd0;
    std::vector<Mat3> R0;
    for(const auto &b : m.bodies)
    {
        x0.push_back(b.x); R0.push_back(b.R); v0.push_back(b.v); w0.push_back(b.w); vd0.push_back(b.vd); wd0.push_back(b.wd);
    }
    auto restore = [&]()
    {
        for(int i=0; i<m.nbody(); ++i)
        {
            rigid_body &b = m.bodies[i];
            b.x = x0[i]; b.R = R0[i]; b.v.setZero(); b.w.setZero(); b.vd.setZero(); b.wd.setZero();
        }
    };
    restore();
    VecX f, in, g;
    MatX G;
    m.forces(t, f);
    m.constraints(g, G);
    // multipliers of the static balance in this configuration: G^T lambda = f (least squares)
    VecX lam = VecX::Zero(nc);
    if(nc)
        lam = (G*G.transpose()).ldlt().solve(G*f);
    auto rest = [&](VecX &r)    // restoring residual -f + G^T lambda (lambda fixed)
    {
        m.forces(t, f);
        m.constraints(g, G);
        r = -f;
        if(nc) r += G.transpose()*lam;
    };
    VecX rp, rm;
    L.K.setZero(n, n);
    L.C.setZero(n, n);
    L.M.setZero(n, n);
    for(int k=0; k<n; ++k)
    {
        const int i = k/6, c = k%6;
        const double e = c < 3 ? 1.0e-6*std::max(1.0, x0[i].cwiseAbs().maxCoeff()) : 1.0e-6;
        for(int s=-1; s<=1; s+=2)
        {
            restore();
            if(c < 3) m.bodies[i].x(c) += s*e;
            else m.bodies[i].R = R0[i]*exp_so3(s*e*Vec3::Unit(c-3));
            rest(s < 0 ? rm : rp);
        }
        L.K.col(k) = (rp - rm)/(2.0*e);
        // velocities and accelerations
        for(int s=-1; s<=1; s+=2)
        {
            restore();
            if(c < 3) m.bodies[i].v(c) = s*1.0e-6; else m.bodies[i].w(c-3) = s*1.0e-6;
            m.forces(t, f);
            (s < 0 ? rm : rp) = -f;
        }
        L.C.col(k) = (rp - rm)/2.0e-6;
        for(int s=-1; s<=1; s+=2)
        {
            restore();
            if(c < 3) m.bodies[i].vd(c) = s*1.0e-6; else m.bodies[i].wd(c-3) = s*1.0e-6;
            m.forces(t, f);
            m.inertia_forces(in);
            (s < 0 ? rm : rp) = in - f;
        }
        L.M.col(k) = (rp - rm)/2.0e-6;
    }
    restore();
    for(int i=0; i<m.nbody(); ++i)
    {
        m.bodies[i].v = v0[i]; m.bodies[i].w = w0[i]; m.bodies[i].vd = vd0[i]; m.bodies[i].wd = wd0[i];
    }
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
    if(!damped)
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
    for(int q=0; q<nmode; ++q)
    {
        res.frequency(q) = fr[idx[q]];
        res.damping(q) = ze[idx[q]];
        VecX s = L.N*sh[idx[q]];
        const double mx = s.cwiseAbs().maxCoeff();
        res.shape.col(q) = mx > 0.0 ? VecX(s/mx) : s;
    }
    res.ok = true;
    return res;
}

}
