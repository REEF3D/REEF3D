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

#include"coral/hydro/coral_radiation.h"
#include"coral/hydro/coral_hydro_db.h"
#include<algorithm>
#include<cmath>
#include<complex>

namespace coral
{

void retardation::compute(const hydro_db &db, double tmax, double dt)
{
    dt_ = dt;
    tmax_ = tmax;
    const int nd = db.ndof();
    const int ns = (int)std::lround(tmax/dt) + 1;
    K_.assign(ns, MatX::Zero(nd, nd));
    // segments [w0, w1] with B linear; B(0) = 0
    std::vector<double> ws;
    std::vector<MatX> Bs;
    ws.push_back(0.0);
    Bs.push_back(MatX::Zero(nd, nd));
    for(int k=0; k<db.nfreq(); ++k)
    {
        if(db.w[k] <= 0.0) continue;
        ws.push_back(db.w[k]);
        Bs.push_back(db.B[k]);
    }
    for(int s=0; s<ns; ++s)
    {
        const double t = s*dt;
        MatX I = MatX::Zero(nd, nd);
        for(size_t q=0; q+1<ws.size(); ++q)
        {
            const double w0 = ws[q], w1 = ws[q+1];
            const MatX &B0 = Bs[q], &B1 = Bs[q+1];
            if(t*(w1 - w0) < 1.0e-4)
            {
                // small t: trapezoid of B cos(w t) is exact to O((t dw)^2) relative
                I += 0.5*(w1 - w0)*(B0*std::cos(w0*t) + B1*std::cos(w1*t));
            }
            else
            {
                const MatX sl = (B1 - B0)/(w1 - w0);
                // int (B0 + sl (w - w0)) cos(w t) dw = [B(w) sin(w t)/t + sl cos(w t)/t^2]
                I += (B1*std::sin(w1*t) - B0*std::sin(w0*t))/t + sl*(std::cos(w1*t) - std::cos(w0*t))/(t*t);
            }
        }
        K_[s] = (2.0/M_PI)*I;
    }
}

MatX retardation::at(double t) const
{
    if(K_.empty() || t < 0.0 || t > tmax_)
        return K_.empty() ? MatX() : MatX::Zero(K_[0].rows(), K_[0].cols());
    const double x = t/dt_;
    const int k = std::min((int)x, (int)K_.size()-2);
    if(k < 0) return K_[0];
    const double f = x - k;
    return (1.0 - f)*K_[k] + f*K_[k+1];
}

MatX ogilvie_ainf(const hydro_db &db, const retardation &rk)
{
    const int nd = db.ndof(), n = db.nfreq();
    MatX sum = MatX::Zero(nd, nd);
    int cnt = 0;
    for(int k=n/2; k<n; ++k)
    {
        const double w = db.w[k];
        if(w <= 0.0) continue;
        // int_0^T K(t) sin(w t) dt, trapezoid on the samples
        MatX I = MatX::Zero(nd, nd);
        for(int s=1; s<rk.nsample(); ++s)
        {
            const double t0 = (s-1)*rk.dt(), t1 = s*rk.dt();
            I += 0.5*rk.dt()*(rk.sample(s-1)*std::sin(w*t0) + rk.sample(s)*std::sin(w*t1));
        }
        sum += db.A[k];
        sum += I/w;
        ++cnt;
    }
    return cnt ? MatX(sum/cnt) : MatX::Zero(nd, nd);
}

// ---------------------------------------------------------------- state space

void radiation_ss::identify(const retardation &rk, double tol, int rmax, double skip)
{
    pairs.clear();
    ndof_ = rk.ndof();
    const int ns = rk.nsample();
    double kmax = 0.0;
    for(int s=0; s<ns; ++s)
        kmax = std::max(kmax, rk.sample(s).cwiseAbs().maxCoeff());
    if(kmax <= 0.0 || ns < 12)
        return;
    // Hankel matrices m x m of the samples with the stride that spans the whole memory
    const int m = std::min(150, (ns - 2)/2);
    const int stride = std::max(1, (ns - 2)/(2*m));
    const double dts = stride*rk.dt();
    for(int i=0; i<ndof_; ++i)
        for(int j=0; j<ndof_; ++j)
        {
            double kij = 0.0;
            for(int s=0; s<ns; ++s)
                kij = std::max(kij, std::fabs(rk.sample(s)(i,j)));
            if(kij <= skip*kmax)
                continue;
            MatX H0(m, m), H1(m, m);
            for(int a=0; a<m; ++a)
                for(int b=0; b<m; ++b)
                {
                    H0(a,b) = rk.sample(stride*(a+b))(i,j);
                    H1(a,b) = rk.sample(stride*(a+b+1))(i,j);
                }
            Eigen::JacobiSVD<MatX> svd(H0, Eigen::ComputeThinU | Eigen::ComputeThinV);
            const VecX sv = svd.singularValues();
            const double total = sv.sum();
            // smallest order whose discarded singular values are below tol; raised (up to rmax) while
            // the kernel fit is worse than 10 tol
            int r0 = 1;
            double tail = total - sv(0);
            while(r0 < std::min(rmax, (int)sv.size()) && tail > tol*total)
            {
                tail -= sv(r0);
                ++r0;
            }
            pair best;
            best.error = 1.0e30;
            for(int r=r0; r<=std::min(rmax, (int)sv.size()); ++r)
            {
                const MatX U = svd.matrixU().leftCols(r), V = svd.matrixV().leftCols(r);
                const VecX s12 = sv.head(r).cwiseSqrt(), sm12 = s12.cwiseInverse();
                const MatX Ad = sm12.asDiagonal()*U.transpose()*H1*V*sm12.asDiagonal();
                const VecX Bd = s12.asDiagonal()*V.transpose().col(0);
                const VecX Cd = (U.row(0)*s12.asDiagonal()).transpose();
                // continuous poles log(lambda_d)/dts in diagonal form; unstable modes are pulled inside
                Eigen::EigenSolver<MatX> es(Ad);
                Eigen::VectorXcd lam = es.eigenvalues();
                for(int q=0; q<r; ++q)
                {
                    if(std::abs(lam(q)) >= 0.9999)
                        lam(q) *= 0.9999/std::abs(lam(q));
                    lam(q) = std::log(lam(q))/dts;
                }
                const Eigen::MatrixXcd Vc = es.eigenvectors();
                pair p;
                p.i = i;
                p.j = j;
                p.lam = lam;
                p.b = Vc.inverse()*Bd.cast<std::complex<double>>();
                p.c = (Cd.transpose().cast<std::complex<double>>()*Vc).transpose();
                double err = 0.0;
                for(int s=0; s<ns; ++s)
                {
                    std::complex<double> k = 0.0;
                    for(int q=0; q<r; ++q)
                        k += p.c(q)*p.b(q)*std::exp(p.lam(q)*(s*rk.dt()));
                    err = std::max(err, std::fabs(k.real() - rk.sample(s)(i,j)));
                }
                p.error = err/kij;
                if(p.error < best.error)
                    best = p;
                if(best.error < 10.0*tol)
                    break;
            }
            pairs.push_back(best);
        }
}

double radiation_ss::kernel(int p, double t) const
{
    const pair &q = pairs[p];
    std::complex<double> k = 0.0;
    for(int r=0; r<q.lam.size(); ++r)
        k += q.c(r)*q.b(r)*std::exp(q.lam(r)*t);
    return k.real();
}

void radiation_ss::advance(const pair &q, Eigen::VectorXcd &x, double u0, double u1, double h)
{
    for(int r=0; r<q.lam.size(); ++r)
    {
        const std::complex<double> z = q.lam(r)*h;
        const std::complex<double> e = std::exp(z);
        std::complex<double> p1, p2;    // (e^z - 1)/z and (e^z - 1 - z)/z^2
        if(std::abs(z) < 1.0e-4)
        {
            p1 = 1.0 + z/2.0 + z*z/6.0;
            p2 = 0.5 + z/6.0 + z*z/24.0;
        }
        else
        {
            p1 = (e - 1.0)/z;
            p2 = (e - 1.0 - z)/(z*z);
        }
        x(r) = e*x(r) + q.b(r)*h*(u0*p1 + (u1 - u0)*p2);
    }
}

int radiation_ss::nstates() const
{
    int n = 0;
    for(const auto &p : pairs)
        n += (int)p.lam.size();
    return n;
}

}
