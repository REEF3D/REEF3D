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
// two cells (VISC + EV, EVH with A 561) of t^n. The system of a layer is the screened 2D Helmholtz problem of the
// eta system, solved with si_eta2D (face coefficients a_f = dt nu_f/dn_f, reaction 1, curvilinear: the face lengths
// over the cell area, as nhflow_curv::lapcoef) and its own 2D solver. First order in time as the vertical diffusion of
// the column operator; the steady state is the same. Without diffusion (nu = 0) the predictor is unchanged.
//
// Not with wave breaking (A 550: hdiffusion_brk keeps its sub-steps), multi-block grids or mesh refinement (the 2D
// systems would need the composite grid of nhflow_amr_si): there the horizontal diffusion stays explicit.
// ---------------------------------------------------------------------------

#include"nhflow_momentum_SI.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"si_eta2D.h"
#include"solver2D.h"
#include<cmath>

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

    if(phd==nullptr)
    {
    phd = new si_eta2D(p);
    psolvHD = new_solver2D(p,pgc);
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

    double *comp[3] = {UHp, VHp, WHp};

    hd_iter = 0;

    for(int kk=0; kk<K; ++kk)
    {
        // face coefficients of layer kk: a_f = dt nu_f/dn_f between two wet fluid cells, 0 at walls, dry cells and open
        // boundary faces (no flux, as hdiffusion)
        k = kk;
        SLICELOOP1
        {
        k = kk;
        hdax(i,j) = (fl(i,j) && fl(i+1,j)) ? dt*0.5*(d->VISC[IJK] + EVh[IJK] + d->VISC[Ip1JK] + EVh[Ip1JK])/fdn(p,1,i,j) : 0.0;
        }

        if(p->j_dir==1)
        SLICELOOP2
        {
        k = kk;
        hday(i,j) = (fl(i,j) && fl(i,j+1)) ? dt*0.5*(d->VISC[IJK] + EVh[IJK] + d->VISC[IJp1K] + EVh[IJp1K])/fdn(p,2,i,j) : 0.0;
        }

        exchange1(p,pgc,hdax);
        exchange2(p,pgc,hday);

        for(int c=0; c<3; ++c)
        {
            if(c==1 && p->j_dir==0)
            continue;

            double *A = comp[c];

            SLICELOOP4
            {
            k = kk;
            hdr(i,j) = act() ? vel(A[PIJK]) : 0.0;
            hdx(i,j) = hdr(i,j);
            }

            phd->solve(p,pgc,psolvHD,hdax,hday,hdr,hdx);
            hd_iter += phd->iter;

            SLICELOOP4
            if(act())
            {
            k = kk;
            A[PIJK] += WL0(i,j)*(hdx(i,j) - hdr(i,j));
            }
        }
    }

    hd_time = pgc->timer() - starttime;
}
