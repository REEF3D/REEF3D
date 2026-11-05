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

#include"iqn_ils.h"
#include<algorithm>
#include<cmath>

iqn_ils::iqn_ils(int slots, double omega, int reuse, int maxcols, double eps, bool imvj) : hist_(std::max(1,slots)), slot_(0),
                 omega_(omega), eps_(eps), reuse_(std::max(0,reuse)), maxcols_(std::max(1,maxcols)), ncols_(0), imvj_(imvj),
                 M_(std::max(1,slots)), mset_(std::max(1,slots),0)
{
}

int iqn_ils::qr(const std::vector<Eigen::VectorXd> &Vc, int n, Eigen::MatrixXd &Q, Eigen::MatrixXd &R, std::vector<int> &acc) const
{
    // QR of the columns by modified Gram-Schmidt; a column that is (nearly) a combination of the newer ones is dropped
    Q.resize(n,Vc.size());
    R = Eigen::MatrixXd::Zero(Vc.size(),Vc.size());
    acc.clear();
    
    for(size_t c=0; c<Vc.size(); ++c)
    {
        Eigen::VectorXd w = Vc[c];
        const double n0 = w.norm();
        
        if(!(n0>0.0) || !std::isfinite(n0))
        continue;
        
        const int m = acc.size();
        Eigen::VectorXd rc = Eigen::VectorXd::Zero(m);
        
        for(int i=0; i<m; ++i)
        {
            rc(i) = Q.col(i).dot(w);
            w -= rc(i)*Q.col(i);
        }
        
        const double nw = w.norm();
        
        if(nw < eps_*n0)
        continue;
        
        Q.col(m) = w/nw;
        R.block(0,m,m,1) = rc;
        R(m,m) = nw;
        acc.push_back(c);
    }
    
    return acc.size();
}

void iqn_ils::begin(int slot)
{
    slot_ = std::min(std::max(0,slot), (int)hist_.size()-1);
    xt_.clear();
    r_.clear();
    mr_.clear();
    ncols_ = 0;
}

Eigen::VectorXd iqn_ils::next(const Eigen::VectorXd &x, const Eigen::VectorXd &xt)
{
    if(imvj_)
    return next_imvj(x,xt);
    
    xt_.push_back(xt);
    r_.push_back(xt - x);

    const int k = (int)r_.size()-1;
    const int n = (int)x.size();
    const int mmax = std::min(maxcols_, std::max(1,n/2));

    // columns, newest first: this sequence, then the kept sequences of earlier calls
    std::vector<Eigen::VectorXd> Vc, Wc;

    for(int j=k; j>0 && (int)Vc.size()<mmax; --j)
    {
        Vc.push_back(r_[j] - r_[j-1]);
        Wc.push_back(xt_[j] - xt_[j-1]);
    }

    for(const auto &b : hist_[slot_])
    for(int c=0; c<b.V.cols() && (int)Vc.size()<mmax; ++c)
    if(b.V.rows()==n)
    {
        Vc.push_back(b.V.col(c));
        Wc.push_back(b.W.col(c));
    }

    // QR of V by modified Gram-Schmidt; a column that is (nearly) a combination of the newer ones is dropped
    Eigen::MatrixXd Q(n,Vc.size()), R = Eigen::MatrixXd::Zero(Vc.size(),Vc.size());
    std::vector<int> acc;

    for(size_t c=0; c<Vc.size(); ++c)
    {
        Eigen::VectorXd w = Vc[c];
        const double n0 = w.norm();

        if(!(n0>0.0) || !std::isfinite(n0))
        continue;

        const int m = acc.size();
        Eigen::VectorXd rc = Eigen::VectorXd::Zero(m);

        for(int i=0; i<m; ++i)
        {
            rc(i) = Q.col(i).dot(w);
            w -= rc(i)*Q.col(i);
        }

        const double nw = w.norm();

        if(nw < eps_*n0)
        continue;

        Q.col(m) = w/nw;
        R.block(0,m,m,1) = rc;
        R(m,m) = nw;
        acc.push_back(c);
    }

    ncols_ = acc.size();

    if(ncols_==0)
    return x + omega_*r_[k];

    const int m = ncols_;
    const Eigen::VectorXd rhs = -(Q.leftCols(m).transpose()*r_[k]);
    const Eigen::VectorXd c = R.topLeftCorner(m,m).triangularView<Eigen::Upper>().solve(rhs);

    Eigen::VectorXd xn = xt;
    for(int i=0; i<m; ++i)
    xn += c(i)*Wc[acc[i]];

    return xn;
}

void iqn_ils::restart()
{
    // the model failed (residual growing): forget the columns of earlier sequences of this slot, keep the ones of
    // the current sequence (plain relaxation diverges for a strong added-mass effect)
    hist_[slot_].clear();
    
    if(imvj_ && mset_[slot_])
    {
        M_[slot_].setZero();
        mset_[slot_] = 0;
        
        for(auto &v : mr_)
        v.setZero();
    }
}

void iqn_ils::end()
{
    if(imvj_)
    {
        end_imvj();
        return;
    }
    
    const int k = (int)r_.size()-1;

    if(reuse_>0 && k>0)
    {
        block b;
        b.V.resize(r_[0].size(),k);
        b.W.resize(r_[0].size(),k);

        for(int j=k, c=0; j>0; --j, ++c)
        {
            b.V.col(c) = r_[j] - r_[j-1];
            b.W.col(c) = xt_[j] - xt_[j-1];
        }

        hist_[slot_].push_front(b);

        while((int)hist_[slot_].size()>reuse_)
        hist_[slot_].pop_back();
    }

    xt_.clear();
    r_.clear();
}

Eigen::VectorXd iqn_ils::next_imvj(const Eigen::VectorXd &x, const Eigen::VectorXd &xt)
{
    const int n = (int)x.size();
    
    if(M_[slot_].rows()!=n)
    {
        M_[slot_] = Eigen::MatrixXd::Zero(n,n);
        mset_[slot_] = 0;
    }
    
    xt_.push_back(xt);
    r_.push_back(xt - x);
    mr_.push_back(mset_[slot_] ? Eigen::VectorXd(M_[slot_]*r_.back()) : Eigen::VectorXd(Eigen::VectorXd::Zero(n)));
    
    const int k = (int)r_.size()-1;
    const int mmax = std::min(maxcols_, std::max(1,n/2));
    
    // columns of the current sequence, newest first
    std::vector<Eigen::VectorXd> Vc, Wc, Zc;
    
    for(int j=k; j>0 && (int)Vc.size()<mmax; --j)
    {
        Vc.push_back(r_[j] - r_[j-1]);
        Wc.push_back(xt_[j] - xt_[j-1]);
        Zc.push_back(mr_[j] - mr_[j-1]);
    }
    
    Eigen::MatrixXd Q, R;
    std::vector<int> acc;
    ncols_ = qr(Vc,n,Q,R,acc);
    
    if(ncols_==0 && !mset_[slot_])
    return x + omega_*r_[k];
    
    // M_k r = M_prev r + (W - Z) (V^T V)^-1 V^T r
    Eigen::VectorXd Mr = mr_[k];
    
    if(ncols_>0)
    {
        const int m = ncols_;
        const Eigen::VectorXd c = R.topLeftCorner(m,m).triangularView<Eigen::Upper>().solve(Q.leftCols(m).transpose()*r_[k]);
        
        for(int i=0; i<m; ++i)
        Mr += c(i)*(Wc[acc[i]] - Zc[acc[i]]);
    }
    
    return xt - Mr;
}

void iqn_ils::end_imvj()
{
    // M_prev <- M_prev + (W - Z)(V^T V)^-1 V^T with the columns of the finished sequence
    const int k = (int)r_.size()-1;
    
    if(k>0)
    {
        const int n = (int)r_[0].size();
        std::vector<Eigen::VectorXd> Vc, Wc, Zc;
        
        for(int j=k; j>0 && (int)Vc.size()<std::min(maxcols_, std::max(1,n/2)); --j)
        {
            Vc.push_back(r_[j] - r_[j-1]);
            Wc.push_back(xt_[j] - xt_[j-1]);
            Zc.push_back(mr_[j] - mr_[j-1]);
        }
        
        Eigen::MatrixXd Q, R;
        std::vector<int> acc;
        const int m = qr(Vc,n,Q,R,acc);
        
        if(m>0)
        {
            Eigen::MatrixXd B(n,m);
            for(int i=0; i<m; ++i)
            B.col(i) = Wc[acc[i]] - Zc[acc[i]];
            
            // (V^T V)^-1 V^T = R^-1 Q^T
            const Eigen::MatrixXd C = R.topLeftCorner(m,m).triangularView<Eigen::Upper>().solve(Eigen::MatrixXd(Q.leftCols(m).transpose()));
            
            M_[slot_].noalias() += B*C;
            mset_[slot_] = 1;
        }
    }
    
    xt_.clear();
    r_.clear();
    mr_.clear();
}
