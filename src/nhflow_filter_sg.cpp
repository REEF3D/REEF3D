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

#include"nhflow_filter.h"
#include"lexer.h"
#include"ghostcell.h"
#include"slice.h"
#include<cmath>
#include<cstdio>

// integer power with the convention 0^0 = 1
static inline double ipow(double x, int n)
{
    double r = 1.0;
    for(int q=0; q<n; ++q)
        r *= x;
    return r;
}

// ---------------------------------------------------------------------
// In-place filtering. Two passes so that neighbour reads are never
// contaminated by already-written cells: pass 1 buffers all filtered
// values in loop order, pass 2 writes them back.
// ---------------------------------------------------------------------
void nhflow_filter::filter_sg(lexer *p, ghostcell *pgc, slice &f)
{
    int iterations = 5;
    
    for(int it=0; it<iterations; ++it)
    {
        pgc->gcsl_start4(p, f, 1);

        buffer.clear();

        SLICELOOP4
        {
            double val=0.0;
            int k=0;
            for(int dv=-hw; dv<=hw; ++dv)
            for(int du=-hw; du<=hw; ++du)
            {
                val += kernel[k]*f(i+du, j+dv);
                ++k;
            }
            buffer.push_back(val);
        }

        size_t n=0;
        SLICELOOP4
        {
            f(i,j) = buffer[n];
            ++n;
        }

        pgc->gcsl_start4(p, f, 1);
    }
}

// ---------------------------------------------------------------------
// Filter "in" into "out" (distinct slices). Single pass, no buffer.
// ---------------------------------------------------------------------
void nhflow_filter::filter_sg(lexer *p, ghostcell *pgc, slice &in, slice &out)
{
    pgc->gcsl_start4(p, in, 1);

    SLICELOOP4
    {
        double val=0.0;
        int k=0;
        for(int dv=-hw; dv<=hw; ++dv)
        for(int du=-hw; du<=hw; ++du)
        {
            val += kernel[k]*in(i+du, j+dv);
            ++k;
        }
        out(i,j) = val;
    }

    pgc->gcsl_start4(p, out, 1);
}

// ---------------------------------------------------------------------
// Precompute the smoothing kernel from the 2D least-squares normal
// equations. The fitted value at the window centre equals the constant
// term of the polynomial, so the kernel weight for each window point k is
//
//     w_k = sum_j  (A^T A)^{-1}_{0j} * A_{kj}
//
// where A is the Vandermonde-type design matrix over the monomials
// u^a v^b (a+b <= order) evaluated at the integer offsets (u,v).
// ---------------------------------------------------------------------
void nhflow_filter::build_coeffs(lexer *p)
{
    // enumerate monomials u^a v^b with a+b <= order; index 0 is the
    // constant term (a=b=0), which is what we evaluate at the centre
    std::vector<int> pa, pb;
    for(int t=0; t<=order; ++t)
    for(int a=0; a<=t; ++a)
    {
        pa.push_back(a);
        pb.push_back(t-a);
    }
    int ncoeff = int(pa.size());
    int npts   = nw*nw;

    // guard: need at least as many window points as unknowns
    if(ncoeff > npts)
    {
        if(p->mpirank==0)
            printf("\n nhflow_filter WARNING: order %d too high for %dx%d window; "
                   "results may be singular.\n", order, nw, nw);
    }

    // design matrix A (npts x ncoeff); row order: dv outer, du inner
    std::vector<std::vector<double> > A(npts, std::vector<double>(ncoeff, 0.0));
    int k=0;
    for(int dv=-hw; dv<=hw; ++dv)
    for(int du=-hw; du<=hw; ++du)
    {
        for(int c=0; c<ncoeff; ++c)
            A[k][c] = ipow(double(du), pa[c]) * ipow(double(dv), pb[c]);
        ++k;
    }

    // normal matrix M = A^T A
    std::vector<std::vector<double> > M(ncoeff, std::vector<double>(ncoeff, 0.0));
    for(int r=0; r<ncoeff; ++r)
    for(int c=0; c<ncoeff; ++c)
    {
        double s=0.0;
        for(int q=0; q<npts; ++q)
            s += A[q][r]*A[q][c];
        M[r][c]=s;
    }

    // invert M
    std::vector<std::vector<double> > Minv;
    invert(M, Minv, ncoeff);

    // kernel weights
    kernel.assign(npts, 0.0);
    double sum=0.0;
    for(int q=0; q<npts; ++q)
    {
        double w=0.0;
        for(int c=0; c<ncoeff; ++c)
            w += Minv[0][c]*A[q][c];
        kernel[q]=w;
        sum += w;
    }

    // defensive normalisation to enforce exact constant reproduction
    if(fabs(sum) > 1.0e-30)
        for(int q=0; q<npts; ++q)
            kernel[q] /= sum;
}

// Gauss-Jordan inversion with partial pivoting (n is small)
void nhflow_filter::invert(std::vector<std::vector<double> > M,
                      std::vector<std::vector<double> >& Minv, int n)
{
    Minv.assign(n, std::vector<double>(n, 0.0));
    for(int i=0; i<n; ++i)
        Minv[i][i]=1.0;

    for(int col=0; col<n; ++col)
    {
        // partial pivot
        int piv=col;
        double best=fabs(M[col][col]);
        for(int r=col+1; r<n; ++r)
            if(fabs(M[r][col])>best){ best=fabs(M[r][col]); piv=r; }

        std::swap(M[piv],    M[col]);
        std::swap(Minv[piv], Minv[col]);

        double d = M[col][col];
        if(fabs(d) < 1.0e-300)
            d = 1.0e-300;

        for(int c=0; c<n; ++c)
        {
            M[col][c]    /= d;
            Minv[col][c] /= d;
        }

        for(int r=0; r<n; ++r)
        {
            if(r==col) continue;
            double f=M[r][col];
            for(int c=0; c<n; ++c)
            {
                M[r][c]    -= f*M[col][c];
                Minv[r][c] -= f*Minv[col][c];
            }
        }
    }
}

void nhflow_filter::print_kernel(lexer *p)
{
    if(p->mpirank!=0)
        return;

    printf("\n Savitzky-Golay kernel (order %d, %dx%d window):\n", order, nw, nw);
    int k=0;
    for(int dv=-hw; dv<=hw; ++dv)
    {
        for(int du=-hw; du<=hw; ++du)
            printf(" % .5f", kernel[k++]);
        printf("\n");
    }
}
