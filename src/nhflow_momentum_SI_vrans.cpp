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

#include"nhflow_momentum_SI.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"vrans_nhflow.h"
#include"vrans_definitions.h"
#include"si_eta2D.h"
#include<cmath>

// ---------------------------------------------------------------------------------------------------------------------
// NHFLOW semi-implicit step (A 510 11/12/13): VRANS, porous structures (B 200 1; B 201, B 202 layers, B 210-230), on one
// Cartesian grid. The conventions of the explicit NHFLOW (vrans_definitions.h): U, V, W, UH, VH, WH and the layer face
// discharges are Darcy (filter) velocities, eta and the pressure are intrinsic,
//
//     du/dt = C ( -div(u u)/n^2 - g grad(eta) - grad(phi) + diffusion + sources ) - C (A + B |u|) u,   C = CPORNH
//     n_s deta/dt + div(Q) = 0,   div(u) = 0
//
// with A, B the Darcy-Forchheimer coefficients of the static structures (vrans_nhflow_f::porous_coeff), n the porosity
// (d->POR), n_s the porosity at the free surface (d->PORS, the storage of the continuity).
//
//  - Drag: implicit in the column operators of the faces and of the cells, dt C (A + B |u|) on the diagonal, |u| of
//    t^n (Picard): unconditionally stable, and gravity, the vertical diffusion and the drag act together, so that the
//    steady Darcy-Forchheimer balance does not depend on dt. |u|: the horizontal speed for UH and VH (on a face its
//    normal velocity from the layer discharge and the mean tangential velocity of the two cells), |w| for WH. The
//    explicit NHFLOW takes the component alone (the same in a flume with one row of cells). B 204 (explicit drag) has
//    no meaning here.
//  - Inertia factor C: the explicit terms of the predictor (vr_predictor) and the rows of the column operator
//    (vertical advection and diffusion, vr_col) times C of the cell (of the face: C of the mean porosity, as CPORNH1);
//    the gravity term with r = A^-1 C instead of A^-1 1 (vr_ones), so the conveyance alpha_f of the eta system
//    contains C and the drag.
//  - Non-hydrostatic pressure: its face terms times C/(1 + dt C (A + B |u|)) of the face layer (vrSx, vrSy: the
//    diagonal part of the column operator, so the coupled system keeps its stencil), its term of w times the same of
//    the cell (vr_sw); in nh_faces, nh_dfaces, nh_residual, nh_assemble_col and nh_cell. The cells take the mean of
//    the face terms as they are (vr_cell_faces: A_c times the mean into the right-hand side of the column solve).
//  - Storage: n_s (WL^{n+1} - WL^n) = -dt div(...) in the eta system (si_eta2D::stor, vr_rhs), in the mass rows of
//    the coupled system and in eta_update. omega: nhflow_sigma::omega_update has the porosity of the layer interfaces
//    (B 221 1).
//  - Advection: the fluxes over the face porosity squared (vr_adv), as the HLL fluxes of the explicit step.
//  - Fourth-order horizontal operators: not in columns with a porous cell nor next to them (ho_begin, vrPc).
// The porosity is cast on the sigma grid of t^n once per step (vrans_nhflow::update, with the force output B 208).
// Not included: the implicit horizontal diffusion (A 512 2) is not scaled by C; the time step does not know the
// advective speed u/n^2 of the pores.
// ---------------------------------------------------------------------------------------------------------------------

void nhflow_momentum_SI::vr_ini(lexer *p)
{
    vr = (p->B200>0);

    if(!vr)
    return;

    vrA.assign(fN,0.0);
    vrB.assign(fN,0.0);
    vrCc.assign(lN,1.0);
    vrDh.assign(lN,0.0);
    vrDw.assign(lN,0.0);
    vrFCx.assign(lN,1.0);
    vrFDx.assign(lN,0.0);
    vrFCy.assign(lN,1.0);
    vrFDy.assign(lN,0.0);
    vrSx.assign(lN,1.0);
    vrSy.assign(lN,1.0);
    vrPc.assign((size_t)p->imax*p->jmax,0);

    if(p->mpirank==0 && !patch)
    cout<<"NHFLOW A 510 11/12/13 with VRANS (B 200): the Darcy-Forchheimer drag implicit in the column operator"<<endl;
}

void nhflow_momentum_SI::vr_begin(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // the porosity on the sigma grid of t^n, the porosity at the free surface; the force on the structures (B 208)
    pvrans->update(p,d,pgc,1.0,1);

    psi->stor = &d->PORS;

    // the coefficients of the cells, with those of the neighbouring ranks
    double a,b;

    LOOP
    {
    pvrans->si_coeff(p,d,i,j,k,a,b);
    vrA[IJK] = a;
    vrB[IJK] = b;
    }

    exchangeC(p,pgc,vrA.data());
    exchangeC(p,pgc,vrB.data());

    // cells: C, the drag of the step with the velocities of t^n
    const double dt = p->dt;

    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    tA(ii,jj) = 0.0;

    LOOP
    {
    const double C = CPORNHval(d->POR[IJK]);
    const double uh = sqrt(d->U[IJK]*d->U[IJK] + d->V[IJK]*d->V[IJK]*p->y_dir);

    vrCc[PIJK] = C;
    vrDh[PIJK] = dt*C*(vrA[IJK] + vrB[IJK]*uh);
    vrDw[PIJK] = dt*C*(vrA[IJK] + vrB[IJK]*fabs(d->W[IJK]));

        if(d->POR[IJK]<1.0-1.0e-12 || vrA[IJK]>0.0 || vrB[IJK]>0.0)
        tA(i,j) = 1.0;
    }

    // columns with a porous cell, with those of the neighbouring ranks
    exchange4(p,pgc,tA);

    for(int ii=p->imin; ii<p->imin+p->imax; ++ii)
    for(int jj=p->jmin; jj<p->jmin+p->jmax; ++jj)
    vrPc[(ii-p->imin)*p->jmax + (jj-p->jmin)] = tA(ii,jj)>0.5 ? 1 : 0;
}

void nhflow_momentum_SI::vr_faces(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // faces: C of the mean porosity, the drag with the mean coefficients of the two cells and the speed of the face
    // (normal: the layer discharge of t^n over the face depth, tangential: the mean of the two cells)
    const double dt = p->dt;

    SLICELOOP1
    for(k=0;k<p->knoz;++k)
    {
    const double C = CPORNHval(0.5*(d->POR[IJK] + d->POR[Ip1JK]));
    const double A = 0.5*(vrA[IJK] + vrA[Ip1JK]);
    const double B = 0.5*(vrB[IJK] + vrB[Ip1JK]);
    const double un = Hx(i,j)>0.0 ? Qx[PIJK]/Hx(i,j) : 0.0;
    const double ut = 0.5*(d->V[IJK] + d->V[Ip1JK])*p->y_dir;

    vrFCx[PIJK] = C;
    vrFDx[PIJK] = dt*C*(A + B*sqrt(un*un + ut*ut));
    vrSx[PIJK] = C/(1.0 + vrFDx[PIJK]);
    }

    if(p->j_dir==1)
    SLICELOOP2
    for(k=0;k<p->knoz;++k)
    {
    const double C = CPORNHval(0.5*(d->POR[IJK] + d->POR[IJp1K]));
    const double A = 0.5*(vrA[IJK] + vrA[IJp1K]);
    const double B = 0.5*(vrB[IJK] + vrB[IJp1K]);
    const double un = Hy(i,j)>0.0 ? Qy[PIJK]/Hy(i,j) : 0.0;
    const double ut = 0.5*(d->U[IJK] + d->U[IJp1K]);

    vrFCy[PIJK] = C;
    vrFDy[PIJK] = dt*C*(A + B*sqrt(un*un + ut*ut));
    vrSy[PIJK] = C/(1.0 + vrFDy[PIJK]);
    }
}

void nhflow_momentum_SI::vr_col(lexer *p, double *cu, double *cl, double *cd)
{
    // the column operator of col_coef in a porous column: its rows (vertical diffusion and advection, also the
    // explicit part of the advection) times C, the drag on the diagonal
    const int K = p->knoz;

    for(int q=0;q<K;++q)
    {
    const double C = vcC[q];

    cu[q] *= C;
    cl[q] *= C;
    cd[q] = C*cd[q] + vcD[q];
    vau[q] *= C;
    val[q] *= C;
    vad[q] *= C;
    }
}

void nhflow_momentum_SI::vr_ones(lexer *p, double *b)
{
    // the right-hand side of r = A^-1 C (the gravity term times the inertia factor)
    for(int q=0;q<p->knoz;++q)
    b[q] = vcC[q];
}

void nhflow_momentum_SI::vr_predictor(lexer *p, fdm_nhf *d, const double *UHs, const double *VHs, const double *WHs,
                                      double *UHo, double *VHo, double *WHo, double a)
{
    // the explicit terms of the predictor stage times C (si_stage_rhs)
    const double alpha = 1.0-a;

    LOOP
    if(p->wet[IJ]==1 && vrCc[PIJK]!=1.0)
    {
    const double C = vrCc[PIJK];

    UHo[PIJK] = a*UH0[PIJK] + alpha*(UHs[PIJK] + p->dt*C*d->F[IJK]);
    VHo[PIJK] = p->j_dir==1 ? a*VH0[PIJK] + alpha*(VHs[PIJK] + p->dt*C*d->G[IJK]) : 0.0;
    WHo[PIJK] = a*WH0[PIJK] + alpha*(WHs[PIJK] + p->dt*C*d->H[IJK]);
    }
}

void nhflow_momentum_SI::vr_adv(lexer *p, fdm_nhf *d)
{
    // the advective fluxes of the faces over the face porosity squared (Darcy velocities)
    SLICELOOP1
    for(k=0;k<p->knoz;++k)
    {
    const double n = 0.5*(d->POR[IJK] + d->POR[Ip1JK]);
    d->Fx[IJK] /= n*n;
    }

    if(p->j_dir==1)
    SLICELOOP2
    for(k=0;k<p->knoz;++k)
    {
    const double n = 0.5*(d->POR[IJK] + d->POR[IJp1K]);
    d->Fy[IJK] /= n*n;
    }
}

void nhflow_momentum_SI::vr_rhs(lexer *p, fdm_nhf *d)
{
    // the storage n_s of the volume equation: n_s WL(eta^{n+1}) - div(a grad eta^{n+1}) = n_s WL^n - dt div(...)
    SLICELOOP4
    rhsW(i,j) += (d->PORS(i,j) - 1.0)*WL0(i,j);
}

void nhflow_momentum_SI::vr_cell_faces(lexer *p, const double *cu, const double *cl, const double *cd, double dtgam,
                                       double *cx, double *cy)
{
    // the phi face terms of a cell of a porous column (nh_cell_faces: the mean of the terms of its faces, which
    // contain the drag of the faces): A_c times them, so that the column solve of the cell returns them as they are
    const int K = p->knoz;
    bool por = false;

    for(int q=0;q<K;++q)
    if(vcD[q]>0.0 || vcC[q]!=1.0)
    por = true;

    if(!por)
    return;

    for(double *c : {cx,cy})
    {
        for(int q=0;q<K;++q)
        {
        double v = (1.0 + cu[q] + cl[q] + cd[q] + (q==0 ? dtgam : 0.0))*c[q];
        if(q<K-1) v -= cu[q]*c[q+1];
        if(q>0) v -= cl[q]*c[q-1];
        wva[q] = v;
        }

        for(int q=0;q<K;++q)
        c[q] = wva[q];
    }
}

double nhflow_momentum_SI::vr_sw(lexer *p)
{
    // the factor of the phi term of w in cell (i,j,k)
    return vrCc[PIJK]/(1.0 + vrDw[PIJK]);
}

void nhflow_momentum_SI::vr_omega(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    // omega of the layer interfaces in a porous column: the water volume a moving interface sweeps is n_f sigma D_t
    // (nhflow_sigma::omega_update has this with B 221 1 only; B 221 0: off). With the storage n_s of the step the
    // flux through the free surface is zero only with it.
    if(p->B221<2)
    return;

    LOOP
    if(p->wet[IJ]==1)
    {
    const double nf = (k==p->knoz-1) ? d->PORS(i,j) : 0.5*(d->POR[IJK] + d->POR[IJKp1]);

    d->omegaF[FIJKp1] += (1.0 - nf)*p->ZN[KP1]*d->detadt(i,j);
    }

    pgc->start7S(p,d->omegaF,17);
}
