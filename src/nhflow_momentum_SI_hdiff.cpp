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

// ---------------------------------------------------------------------------
// NHFLOW semi-implicit step (A 510 11/12/13): implicit horizontal diffusion of the predictor (A 512 2, one grid).
//
// The explicit horizontal diffusion of the predictor (hdiffusion: WL div_h(nu_h grad_h u) per layer, compact Laplacian,
// no flux through walls, to dry cells and through open boundary faces) limits the time step to
// dt <= 0.25/(nu_h (1/dx^2 + 1/dy^2)); with k-omega at dx = 1 cm (abutment flume, nu_t up to 7e-3 m2/s) that was
// 0.002 s, a third of the advective limit. With A 512 2 the same operator acts implicitly (backward Euler) after the
// predictor, per layer k and velocity component u = U, V, W on the sigma grid of t^n:
//
//   u** - dt L_h(u**) = u*,   L_h u = sum_f c_f nu_f (u_nb - u),   UH** = UH* + WL^n (u** - u*)
//
// with u* the velocity of the predictor UH* (desingularised below A 546, as the stage velocities), nu_f the mean of the
// two cells (VISC + EV, EVH with A 561) of t^n; face coefficients a_f = dt nu_f/dn_f between two wet fluid cells,
// 0 at walls, dry cells and open boundary faces (curvilinear: the face lengths over the cell area, as
// nhflow_curv::lapcoef). First order in time as the vertical diffusion of the column operator; the steady state is the
// same. Without diffusion (nu = 0) the predictor is unchanged.
//
// The systems (one per layer and component, 3 K in 3D) are solved together: Jacobi-preconditioned conjugate gradients
// on the symmetric form (rows times the cell area: W u + sum_f a_f l_f (u - u_nb) = W u*, l_f the face length), all
// systems in the same sweeps over the layer arrays, one halo exchange per component and two reductions (of all
// systems) per iteration, starting from u*. The systems are well conditioned (reaction 1: kappa <= 1 + 8 a_f/dx); a
// system stops when |r|_2 <= N 44 |b|_2. Before: 3 K separate REEFMG2D solves per step (30 in the abutment flume with
// 10 layers, each with its own multigrid set-up and messages, and a terminal line each).
//
// Not with wave breaking (A 550: hdiffusion_brk keeps its sub-steps), multi-block grids or mesh refinement (the 2D
// systems would need the composite grid of nhflow_amr_si): there the horizontal diffusion stays explicit.
// ---------------------------------------------------------------------------

#include"nhflow_momentum_SI.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include<cmath>
#include<vector>

bool nhflow_momentum_SI::hdiff_implicit(lexer *p)
{
    return p->A512==2 && p->geo_mb==0 && p->G1==0 && !(p->A550==1 || p->A550>=3);
}

void nhflow_momentum_SI::hdiffusion_implicit(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    const double starttime = pgc->timer();
    const double dt = p->dt;
    const int K = p->knoz;
    const double hv4 = p->A546*p->A546*p->A546*p->A546;
    double *EVh = (p->A561>0) ? d->EVH : d->EV;
    const bool curv = p->geo_curv==1;
    const size_t N = fN;

    // the components: U, V (with j_dir 1), W
    double *comp[3] = {UHp, VHp, WHp};
    int cl[3], nc=0;
    for(int c=0; c<3; ++c)
    if(c!=1 || p->j_dir==1)
    cl[nc++] = c;

    const int ns = nc*K;    // systems: component a (of cl), layer k: a*K + k

    // the arrays are kept from step to step (zero where they are not written: the ghost cells of walls, faces outside
    // SLICELOOP1/2)
    std::vector<double> &fx = hiF[0], &fy = hiF[1], &D = hiF[2];

    if(hiF[0].size()!=N)
    {
    for(int q=0; q<3; ++q)
    hiF[q].assign(N,0.0);

        for(int a=0; a<3; ++a)
        {
        hdX[a].assign(N,0.0);
        hdR[a].assign(N,0.0);
        hdP[a].assign(N,0.0);
        hdA[a].assign(N,0.0);
        }
    }

    // a wet fluid cell (also of a neighbouring rank): the faces between two of them are open
    auto fl = [&](int a, int b)
    {
        const int q = (a-p->imin)*p->jmax + (b-p->jmin);
        return p->flagslice4[q]>0 && p->wet[q]==1;
    };

    // a cell of the diffusion: wet, above the wetting criterion
    auto act = [&]() { return p->wet[IJ]==1 && WL0(i,j)>p->A544; };

    // the velocity of a layer momentum (desingularised below A 546, as si_stage_vel)
    auto vel = [&](double uh)
    {
        const double h = WL0(i,j);

        if(h>=p->A546)
        return uh/h;

        return sqrt(2.0)*h/sqrt(h*h*h*h + hv4)*uh;
    };

    // the cell area (the weight of a row) and the face lengths
    auto area = [&]() { return curv ? p->AREA2D[IJ] : p->DXN[IP]*(p->j_dir==1 ? p->DYN[JP] : 1.0); };

    // face coefficients times the face length: x-faces (i,j)-(i+1,j), y-faces (i,j)-(i,j+1), every layer
    SLICELOOP1
    {
    const bool open = fl(i,j) && fl(i+1,j);
    const double r = open ? dt*(curv ? p->LEN1[IJ] : (p->j_dir==1 ? p->DYN[JP] : 1.0))/fdn(p,1,i,j) : 0.0;

        for(k=0;k<K;++k)
        fx[IJK] = open ? r*0.5*(d->VISC[IJK] + EVh[IJK] + d->VISC[Ip1JK] + EVh[Ip1JK]) : 0.0;
    }

    if(p->j_dir==1)
    SLICELOOP2
    {
    const bool open = fl(i,j) && fl(i,j+1);
    const double r = open ? dt*(curv ? p->LEN2[IJ] : p->DXN[IP])/fdn(p,2,i,j) : 0.0;

        for(k=0;k<K;++k)
        fy[IJK] = open ? r*0.5*(d->VISC[IJK] + EVh[IJK] + d->VISC[IJp1K] + EVh[IJp1K]) : 0.0;
    }

    exchangeF(p,pgc,fx.data());
    exchangeF(p,pgc,fy.data());

    // diagonal, right-hand side b = W u*, x0 = u*, r0 = b - A x0 (for the convergence the norm of b)
    std::vector<double> bb(ns,0.0), rz(ns,0.0), red(2*ns,0.0);
    std::vector<char> on(ns,1);

    SLICELOOP4
    {
    const double W = area();
    const bool ac = act();

        for(k=0;k<K;++k)
        {
        D[IJK] = W + fx[IJK] + fx[Im1JK] + (p->j_dir==1 ? fy[IJK] + fy[IJm1K] : 0.0);

            for(int a=0; a<nc; ++a)
            hdX[a][IJK] = ac ? vel(comp[cl[a]][PIJK]) : 0.0;
        }
    }

    // A x: W x + sum_f c_f (x - x_nb); x with its halo
    auto apply = [&](const std::vector<double> &x, std::vector<double> &y)
    {
        SLICELOOP4
        for(k=0;k<K;++k)
        {
        double s = D[IJK]*x[IJK] - fx[IJK]*x[Ip1JK] - fx[Im1JK]*x[Im1JK];

            if(p->j_dir==1)
            s -= fy[IJK]*x[IJp1K] + fy[IJm1K]*x[IJm1K];

        y[IJK] = s;
        }
    };

    for(int a=0; a<nc; ++a)
    {
    std::vector<double> &Ap = hdA[a];
    exchangeC(p,pgc,hdX[a].data());
    apply(hdX[a],Ap);

        SLICELOOP4
        {
        const double W = area();

            for(k=0;k<K;++k)
            {
            const double b = W*hdX[a][IJK];
            const double r = b - Ap[IJK];
            hdR[a][IJK] = r;
            hdP[a][IJK] = r/D[IJK];
            bb[a*K+k] += b*b;
            rz[a*K+k] += r*r/D[IJK];
            red[a*K+k] += r*r;
            }
        }
    }

    // (b.b, r.r, r.z): one reduction
    {
    std::vector<double> v(3*ns);
    for(int s=0; s<ns; ++s)
    {
    v[s] = bb[s];
    v[ns+s] = red[s];
    v[2*ns+s] = rz[s];
    }
    MPI_Allreduce(MPI_IN_PLACE,v.data(),3*ns,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);
    for(int s=0; s<ns; ++s)
    {
    bb[s] = v[s];
    rz[s] = v[2*ns+s];
    on[s] = (bb[s]>0.0 && v[ns+s]>p->N44*p->N44*bb[s]) ? 1 : 0;
    }
    }

    std::vector<double> pap(ns), alpha(ns), rr(ns), rzn(ns);
    int it=0;

    auto any = [&]() { for(int s=0; s<ns; ++s) if(on[s]) return true; return false; };

    while(any() && it<p->N46)
    {
        ++it;
        std::fill(pap.begin(),pap.end(),0.0);

        // p.Ap of every system
        for(int a=0; a<nc; ++a)
        {
        exchangeC(p,pgc,hdP[a].data());
        apply(hdP[a],hdA[a]);

            SLICELOOP4
            for(k=0;k<K;++k)
            pap[a*K+k] += hdP[a][IJK]*hdA[a][IJK];
        }

        MPI_Allreduce(MPI_IN_PLACE,pap.data(),ns,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);

        for(int s=0; s<ns; ++s)
        alpha[s] = (on[s] && pap[s]>0.0) ? rz[s]/pap[s] : 0.0;

        std::fill(rr.begin(),rr.end(),0.0);
        std::fill(rzn.begin(),rzn.end(),0.0);

        for(int a=0; a<nc; ++a)
        {
        const std::vector<double> &Ap = hdA[a];

            SLICELOOP4
            for(k=0;k<K;++k)
            {
            const double al = alpha[a*K+k];
            hdX[a][IJK] += al*hdP[a][IJK];
            const double r = hdR[a][IJK] - al*Ap[IJK];
            hdR[a][IJK] = r;
            rr[a*K+k] += r*r;
            rzn[a*K+k] += r*r/D[IJK];
            }
        }

        {
        std::vector<double> v(2*ns);
        for(int s=0; s<ns; ++s)
        {
        v[s] = rr[s];
        v[ns+s] = rzn[s];
        }
        MPI_Allreduce(MPI_IN_PLACE,v.data(),2*ns,MPI_DOUBLE,MPI_SUM,pgc->mpi_comm);
        for(int s=0; s<ns; ++s)
        {
        rr[s] = v[s];
        rzn[s] = v[ns+s];
        }
        }

        // p = z + beta p; a converged system stops (its x stays)
        for(int s=0; s<ns; ++s)
        if(on[s] && rr[s]<=p->N44*p->N44*bb[s])
        on[s] = 0;

        for(int a=0; a<nc; ++a)
        SLICELOOP4
        for(k=0;k<K;++k)
        {
        const int s = a*K+k;
            if(on[s])
            {
            const double beta = rz[s]>0.0 ? rzn[s]/rz[s] : 0.0;
            hdP[a][IJK] = hdR[a][IJK]/D[IJK] + beta*hdP[a][IJK];
            }
            else
            hdP[a][IJK] = 0.0;
        }

        for(int s=0; s<ns; ++s)
        rz[s] = rzn[s];
    }

    // UH** = UH* + WL^n (u** - u*)
    for(int a=0; a<nc; ++a)
    {
    double *A = comp[cl[a]];

        SLICELOOP4
        if(act())
        for(k=0;k<K;++k)
        {
        const double u0 = vel(A[PIJK]);
        A[PIJK] += WL0(i,j)*(hdX[a][IJK] - u0);
        }
    }

    // the iterations and time of the step (both passes of tl2; si_begin resets them)
    hd_iter += it;
    hd_sys = ns;
    hd_time += pgc->timer() - starttime;
}
