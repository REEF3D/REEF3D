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

#ifndef SFLOW_MOMENTUM_SI_NH_H_
#define SFLOW_MOMENTUM_SI_NH_H_

#include"reefmg2D.h"
#include"matrix2D.h"
#include"vec2D.h"
#include<vector>
#include<memory>
#include<algorithm>
#include<cmath>

// work space of the coupled solve: band LU, FGMRES vectors, the face map of the open boundaries
struct sflow_si_nh_lu
{
    // band LU (row-major, each row holds the columns r-kl .. r+ku+kl)
    int n=0, kl=5, ku=5, W=0;
    std::vector<double> AB;
    std::vector<int> piv;

    double& at(int r, int c) { return AB[(size_t)r*W + (c - r + kl)]; }

    void resize(int nn) { n=nn; W=2*kl+ku+1; AB.assign((size_t)n*W,0.0); piv.assign(n,0); }

    bool factor()
    {
        for(int k=0; k<n; ++k)
        {
            int pr=k; double pm=fabs(at(k,k));
            for(int r=k+1; r<=std::min(n-1,k+kl); ++r)
            if(fabs(at(r,k))>pm) { pm=fabs(at(r,k)); pr=r; }

            if(pm==0.0)
            return false;

            piv[k]=pr;
            const int cmax = std::min(n-1,k+ku+kl);

            if(pr!=k)
            for(int c=k; c<=cmax; ++c)
            std::swap(at(k,c),at(pr,c));

            const double d = at(k,k);

            for(int r=k+1; r<=std::min(n-1,k+kl); ++r)
            {
                const double l = at(r,k)/d;
                at(r,k) = l;

                if(l!=0.0)
                for(int c=k+1; c<=cmax; ++c)
                at(r,c) -= l*at(k,c);
            }
        }
        return true;
    }

    void solve(std::vector<double> &x)
    {
        for(int k=0; k<n; ++k)
        {
            if(piv[k]!=k)
            std::swap(x[k],x[piv[k]]);

            for(int r=k+1; r<=std::min(n-1,k+kl); ++r)
            x[r] -= at(r,k)*x[k];
        }

        for(int k=n-1; k>=0; --k)
        {
            double s = x[k];
            for(int c=k+1; c<=std::min(n-1,k+ku+kl); ++c)
            s -= at(k,c)*x[c];
            x[k] = s/at(k,k);
        }
    }

    // FGMRES
    std::vector<std::vector<double>> V, Z;
    std::vector<double> H, cs, sn, g;

    // open boundary face of x-face / y-face (i,j): index into bcf, -1 none
    std::vector<int> bfx, bfy;
    int nbcf=-1;

    // preconditioner: M5 (the 5-point Schur operator) and A2 = I - div(rho grad) (the A 220 3 term), one V-cycle each
    std::unique_ptr<reefmg2D> mg1, mg2;
    std::unique_ptr<matrix2D> M1, M2;
    std::unique_ptr<vec2D> r1, r2;
};

#endif
