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

// REEF3D::CORAL: the state of a run as a vector of doubles (restart of REEF3D from a checkpoint, P 48,
// I 45). state_out appends values in a fixed order, state_in takes them back in the same order; the
// sizes of the dynamic arrays are written with them. A vector that is too short leaves the rest as it
// is and clears ok.

#ifndef CORAL_STATE_H_
#define CORAL_STATE_H_

#include"coral/core/coral_types.h"
#include<complex>
#include<deque>
#include<vector>

namespace coral
{
    struct state_out
    {
        std::vector<double> &v;

        void d(double x) {v.push_back(x);}
        void i(long x) {v.push_back(double(x));}
        void b(bool x) {v.push_back(x ? 1.0 : 0.0);}
        template<class M>
        void m(const M &a)
        {
            v.push_back(double(a.rows()));
            v.push_back(double(a.cols()));
            for(index c=0; c<a.cols(); ++c)
                for(index r=0; r<a.rows(); ++r)
                    v.push_back(double(a(r,c)));
        }
        void cx(const Eigen::VectorXcd &a)
        {
            v.push_back(double(a.size()));
            for(index k=0; k<a.size(); ++k)
            {
                v.push_back(a(k).real());
                v.push_back(a(k).imag());
            }
        }
        void vd(const std::vector<double> &a) {v.push_back(double(a.size())); v.insert(v.end(), a.begin(), a.end());}
        template<class T, class F>
        void list(const T &a, F f) {v.push_back(double(a.size())); for(const auto &x : a) f(x);}
    };

    struct state_in
    {
        const std::vector<double> &v;
        size_t k = 0;
        bool ok = true;

        bool has(size_t n) {if(k + n > v.size()) ok = false; return ok;}
        double next() {return has(1) ? v[k++] : 0.0;}
        void d(double &x) {if(has(1)) x = v[k++];}
        void i(int &x) {if(has(1)) x = int(v[k++]);}
        void i(long &x) {if(has(1)) x = long(v[k++]);}
        void i(size_t &x) {if(has(1)) x = size_t(v[k++]);}
        void b(bool &x) {if(has(1)) x = v[k++] > 0.5;}
        template<class M>
        void m(M &a)
        {
            if(!has(2))
                return;
            const index r = index(v[k]), c = index(v[k+1]);
            k += 2;
            if(!has(size_t(r*c)))
                return;
            a.resize(r, c);
            for(index cc=0; cc<c; ++cc)
                for(index rr=0; rr<r; ++rr)
                    a(rr,cc) = typename M::Scalar(v[k++]);
        }
        void cx(Eigen::VectorXcd &a)
        {
            if(!has(1))
                return;
            const index n = index(v[k++]);
            if(!has(size_t(2*n)))
                return;
            a.resize(n);
            for(index q=0; q<n; ++q)
            {
                a(q) = std::complex<double>(v[k], v[k+1]);
                k += 2;
            }
        }
        void vd(std::vector<double> &a)
        {
            if(!has(1))
                return;
            const size_t n = size_t(v[k++]);
            if(!has(n))
                return;
            a.assign(v.begin() + k, v.begin() + k + n);
            k += n;
        }
        // a list: its size, then f(element) for each (the container resized first)
        template<class T, class F>
        void list(T &a, F f)
        {
            if(!has(1))
                return;
            const size_t n = size_t(v[k++]);
            a.resize(n);
            for(auto &x : a)
            {
                if(!ok)
                    return;
                f(x);
            }
        }
    };
}

#endif
