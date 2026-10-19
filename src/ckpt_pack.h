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

#ifndef CKPT_PACK_H_
#define CKPT_PACK_H_

//  The scalars of an object in a checkpoint (P 48, I 45): packed into a vector of doubles in a fixed
//  order (ckpt_put) and taken back in the same order (ckpt_get).  A short vector leaves the rest as it
//  is and clears ok.

#include<vector>
#include<cstddef>
#include<Eigen/Dense>

struct ckpt_put
{
    std::vector<double> &v;

    void d(double x) { v.push_back(x); }
    void i(int x) { v.push_back(double(x)); }
    void b(bool x) { v.push_back(x ? 1.0 : 0.0); }
    template<class T, int R, int C, int O, int MR, int MC>
    void m(const Eigen::Matrix<T,R,C,O,MR,MC> &a) { for(int n=0; n<a.size(); ++n) v.push_back(double(a.data()[n])); }
    void vd(const std::vector<double> &a) { v.push_back(double(a.size())); v.insert(v.end(),a.begin(),a.end()); }
    void ad(const double *a, int n) { for(int q=0; q<n; ++q) v.push_back(a[q]); }
};

struct ckpt_get
{
    const std::vector<double> &v;
    size_t k = 0;
    bool ok = true;

    bool has(size_t n) { if(k+n>v.size()) ok = false; return ok; }
    void d(double &x) { if(has(1)) x = v[k++]; }
    void i(int &x) { if(has(1)) x = int(v[k++]); }
    void b(bool &x) { if(has(1)) x = v[k++]>0.5; }
    template<class T, int R, int C, int O, int MR, int MC>
    void m(Eigen::Matrix<T,R,C,O,MR,MC> &a) { if(has(size_t(a.size()))) for(int n=0; n<a.size(); ++n) a.data()[n] = T(v[k++]); }
    void vd(std::vector<double> &a) { if(!has(1)) return; const size_t n = size_t(v[k++]); if(!has(n)) return; a.assign(v.begin()+k,v.begin()+k+n); k += n; }
    void ad(double *a, int n) { if(has(size_t(n))) for(int q=0; q<n; ++q) a[q] = v[k++]; }
};

#endif
