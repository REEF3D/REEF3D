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

// ---------------------------------------------------------------------------
// NHFLOW semi-implicit step (A 510 11/12/13): 6DOF bodies with direct forcing (X 10 1 floating body, X 10 2 prescribed
// motion, X 10 4 wavemaker; X 15 0-5; one Cartesian grid).
//
// The bodies use the direct forcing of the fixed solids (df_*, nhflow_momentum_SI.cpp) with the forcing of the body
// object (sixdof_obj_nhflow::update_forcing_nhflow, the same kernel as the explicit NHFLOW):
//
// 1. Body step, once per time step at its start (fb_step): sixdof_nhflow::start_nhflow with the fluid of t^n. The
//    floating body (X 10 1) takes the loads of t^n (pressure and water level at the hull, moorings) and makes the
//    semi-implicit Euler step of the CFD one-step scheme (N 40 5): the momenta with the loads of t^n, the pose with the
//    new momenta (sixdof_obj_nhflow::solve_eqmotion_nhflow). An oscillation keeps its amplitude, the pose advances
//    with the velocity the fluid sees over the step. Prescribed motions (X 10 2) and the wavemaker (X 10 4) move to
//    t^{n+1} with their velocity. The level set FB is cast on the sigma grid of t^n, the cells with their centre in
//    the body get DF = -2 (df_cast).
// 2. Cells: the predictor and the end of the step are forced to the body velocity (df_apply: nhflow_forcing::
//    reforcing with sixdof_nhflow::reforce_nhflow), as the fixed solids.
// 3. Layer faces with a body cell on either side: the discharge of the body motion Q_k = H_f u_b.n inside the solves
//    (fb_blocked): in the column of the face the predictor part QA_k = Q_b and the response to the gravity term R_k =
//    0 (si_assemble), the dissipation (A 538) on the open part of the face, Q_k = Q_b in face_update and no
//    derivatives in the coupled eta-phi system (nh_dfaces). The eta and eta-phi solves see the body at once, its
//    acceleration has the reaction of the water in the same step (with the incremental face force of the fixed solids,
//    df_faces, that reaction built up over several steps: the body fell with too little added mass). u_b is the body
//    velocity of the body cells (fbU, fbV: the forcing of the body on the fluid at rest), the mean of the two cells
//    where both are in the body.
// 4. Columns of a body over the whole water depth (a piston or flap wavemaker, a grounded body) are walls of the
//    implicit step as the full-depth solids (df_columns). The faces between such a column and a wet water column
//    carry the discharge of the moving wall, Q_k = u_b,k.n WL^n of the water side: in the explicit part of the eta
//    system at both time levels (fb_walls_assemble), in the face discharges (fb_walls_faces, face_update) and so in
//    the transport of the step. The water column of the body drains by it as the wall advances (the volume is
//    exact).
// 5. The non-hydrostatic pressure (A 520 > 0): the coupled eta-phi system holds the body velocity in the cells of the
//    body (fb_wcell). The continuity row of the lowest body cell of a column then holds w = w_b on the hull bottom, its
//    phi is the pressure there; the pressure of the step contains the reaction of the water to the motion of the body
//    of the step (the added mass). With the body cells free in the solve and forced after it the body sank with too
//    little added mass (heave decay: the first trough 2x too deep against CFD).
// 6. The body is advanced with the loads of t^n, which contain the reaction of the water to the acceleration of the
//    step before: for a body lighter than its added mass that explicit coupling grows from step to step. The
//    stabilised coupling X 197 of CFD, (M + A) a = F + A a_prev with the added mass A of the equivalent ellipsoid of
//    the hull and the acceleration extrapolated from the last two steps (sixdof_obj::cfd_added_mass), X 197 1 by
//    default for the floating body with A 520 > 0 (check).
// The time step stays within the Courant number N 47 of the body velocity (nhflow_timestep_SI::limit).
//
// Not yet: the ship pressure X 10 3, the pierced columns X 13, X 15 6-11, porous bodies X 16, nets X 320, 6DOF on
// curvilinear, multi-block and refined grids.
// ---------------------------------------------------------------------------

#include"nhflow_momentum_SI.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"6DOF.h"
#include<cmath>
#include<vector>

void nhflow_momentum_SI::fb_step(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    const size_t N = (size_t)p->imax*p->jmax*(p->kmax+2);

    // the water level and the pressure of t^n with the cells of the neighbouring ranks (the semi-implicit step writes
    // d->P without its halo; sixdof_nhflow::reforce_nhflow, called by the forcing, evaluates the loads at the hull as
    // well); the sigma grid of t^n for the cast of the body. The body step itself takes the loads of the end of the
    // step before (df_reforce)
    exchange4(p,pgc,d->WL);
    exchange4(p,pgc,d->eta);
    pgc->start7P(p,d->P,540);

    sigma_update(p,d,pgc,d->WL);

    fbFX.assign(N,0.0);
    fbFY.assign(N,0.0);
    fbFZ.assign(N,0.0);

    LOOP
    d->FHB[IJK] = 0.0;

    SLICELOOP4
    fbfe(i,j) = 0.0;

    // the body to t^{n+1}, its level set on the sigma grid of t^n, its velocity, log and output (finalize)
    p6dof->start_nhflow(p,d,pgc,0,d->U,d->V,d->W,fbFX.data(),fbFY.data(),fbFZ.data(),d->WL,fbfe,true);

    // DF of the step with the solids and the new body
    df_cast(p,d,pgc);

    // the body velocity of the body cells: the forcing of the body on the fluid at rest, H u_b (the full forcing
    // inside the body, the normal part outside it); divided by H (d->FHB of the forcing)
    fbU.assign(N,0.0);
    fbV.assign(N,0.0);
    fbW.assign(N,0.0);

    dfSU.assign(d->U,d->U+N);
    dfSV.assign(d->V,d->V+N);
    dfSW.assign(d->W,d->W+N);

    df_apply(p,d,pgc,fbU.data(),fbV.data(),fbW.data(),d->WL,false);

    std::copy(dfSU.begin(),dfSU.end(),d->U);
    std::copy(dfSV.begin(),dfSV.end(),d->V);
    std::copy(dfSW.begin(),dfSW.end(),d->W);

    LOOP
    {
        const double h = d->WL(i,j)*d->FHB[IJK];

        if(p->DF[IJK]==-2 && h>1.0e-10)
        {
        fbU[IJK] /= h;
        fbV[IJK] = p->j_dir==1 ? fbV[IJK]/h : 0.0;
        fbW[IJK] /= h;
        }

        else
        fbU[IJK] = fbV[IJK] = fbW[IJK] = 0.0;
    }

    pgc->gcparaxV(p,fbU.data(),14);
    pgc->gcparacoxV(p,fbU.data(),14);
    pgc->gcparaxV(p,fbV.data(),14);
    pgc->gcparacoxV(p,fbV.data(),14);
}

// a layer face (dir 1: x, 2: y) of the cell (i,j,k) and its neighbour blocked by a body (no solid cell): its discharge is
// that of the body motion in the solves, the eta and the eta-phi system see the body (si_assemble: QA = Q_b, R = 0,
// face_update, nh_dfaces)
bool nhflow_momentum_SI::fb_blocked(lexer *p, int dir)
{
    const int a = p->DF[IJK];
    const int b = dir==1 ? p->DF[Ip1JK] : p->DF[IJp1K];

    return (a==-2 || b==-2) && a!=-1 && b!=-1;
}

// the coupled eta-phi system (A 520 > 0): the cell (i,j,k) in a body (DF -2) has the vertical velocity of the body,
// w_k = w_b instead of WH*/H - dt (phi_{k+1} - phi_k)/h_k (nh_residual, nh_assemble). The continuity row of the lowest
// body cell of the column then holds w = w_b on the hull bottom, its phi is the pressure there, which the cell below
// takes; a body cell above another one has phi_k = phi_{k-1}. 0: not in a body, 1: the lowest body cell (or one above
// water), 2: a body cell above another one
int nhflow_momentum_SI::fb_wcell(lexer *p)
{
    if(p->DF[IJK]!=-2)
    return 0;

    return (k>0 && p->DF[IJKm1]==-2) ? 2 : 1;
}

// the discharge of the body motion of the layer face (dir 1: x, 2: y) of the cell (i,j,k) and its neighbour, where one of
// the two cells is in a solid or a body: a solid closes the face, a face of the moving wall carries its discharge
// (fb_walls_set), any other face with a body cell H_f u_b.n
double nhflow_momentum_SI::fb_qface(lexer *p, int dir)
{
    const int a = p->DF[IJK];
    const int b = dir==1 ? p->DF[Ip1JK] : p->DF[IJp1K];

    if(a==-1 || b==-1)
    return 0.0;

    if(fbwalls && (dir==1 ? fbwx[IJ] : fbwy[IJ]))
    return dir==1 ? fbQx[PIJK] : fbQy[PIJK];

    const double H = dir==1 ? Hx(i,j) : Hy(i,j);

    if(H<=0.0)
    return 0.0;

    const std::vector<double> &ub = dir==1 ? fbU : fbV;
    const double ua = ub[IJK];
    const double un = dir==1 ? ub[Ip1JK] : ub[IJp1K];

    double u = 0.0;

    if(a==-2 && b==-2)
    u = 0.5*(ua + un);

    else
    u = a==-2 ? ua : un;

    return H*u;
}

// the faces between a column of a body over the whole water depth (a wall of the step, df_columns) and a wet water
// column: the layer discharges of the moving wall, Q_k = u_b,k.n WL^n of the water side
void nhflow_momentum_SI::fb_walls_set(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    fbwalls = false;

    if(!dfcols)
    return;

    const size_t NS = (size_t)p->imax*p->jmax;

    fbwx.assign(NS,0);
    fbwy.assign(NS,0);

    if(fbQx.size()!=lN)
    {
    fbQx.assign(lN,0.0);
    fbQy.assign(lN,0.0);
    }

    auto sidx = [&](int ii, int jj) { return (ii-p->imin)*p->jmax + (jj-p->jmin); };
    auto cidx = [&](int ii, int jj, int kk) { return (ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax + kk-p->kmin; };

    // a column of a body: a wall column with a body cell and no solid cell
    auto bcol = [&](int ii, int jj)
    {
        if(!df_col(p,ii,jj))
        return false;

        bool body=false;

        for(int kk=0; kk<p->knoz; ++kk)
        {
        const int f = p->DF[cidx(ii,jj,kk)];

            if(f==-1)
            return false;

            if(f==-2)
            body = true;
        }

        return body;
    };

    auto water = [&](int ii, int jj)
    {
        const int q = sidx(ii,jj);
        return p->flagslice4[q]>0 && p->wet[q]==1 && !df_col(p,ii,jj);
    };

    int n=0;

    SLICELOOP1
    {
        const bool l = bcol(i,j), r = bcol(i+1,j);

        if(l==r)
        continue;

        const int ib = l ? i : i+1, iw = l ? i+1 : i;

        if(!water(iw,j))
        continue;

        fbwx[IJ] = 1;
        ++n;

        for(k=0;k<p->knoz;++k)
        fbQx[PIJK] = fbU[cidx(ib,j,k)]*WL0(iw,j);
    }

    if(p->j_dir==1)
    SLICELOOP2
    {
        const bool l = bcol(i,j), r = bcol(i,j+1);

        if(l==r)
        continue;

        const int jb = l ? j : j+1, jw = l ? j+1 : j;

        if(!water(i,jw))
        continue;

        fbwy[IJ] = 1;
        ++n;

        for(k=0;k<p->knoz;++k)
        fbQy[PIJK] = fbV[cidx(i,jb,k)]*WL0(i,jw);
    }

    fbwalls = pgc->globalisum(n)>0;
}

// the explicit part of the eta system: the discharge of the moving wall at both time levels (q^E and q^n), and the
// layer discharges of t^n, so that the transport of the step is the discharge of the wall
void nhflow_momentum_SI::fb_walls_assemble(lexer *p, ghostcell *pgc)
{
    if(!fbwalls)
    return;

    SLICELOOP1
    if(fbwx[IJ])
    {
    double q=0.0;

        for(k=0;k<p->knoz;++k)
        {
        Qxn[PIJK] = fbQx[PIJK];
        q += p->DZN[KP]*fbQx[PIJK];
        }

    qxE(i,j) = qxn(i,j) = q;
    }

    if(p->j_dir==1)
    SLICELOOP2
    if(fbwy[IJ])
    {
    double q=0.0;

        for(k=0;k<p->knoz;++k)
        {
        Qyn[PIJK] = fbQy[PIJK];
        q += p->DZN[KP]*fbQy[PIJK];
        }

    qyE(i,j) = qyn(i,j) = q;
    }

    exchange1(p,pgc,qxE);
    exchange1(p,pgc,qxn);
    exchange2(p,pgc,qyE);
    exchange2(p,pgc,qyn);
    pexchF(p,pgc,Qxn);
    pexchF(p,pgc,Qyn);
}

// the layer discharges of t^{n+1} of the faces of the moving wall (face_update, before its exchange and the transport)
void nhflow_momentum_SI::fb_walls_faces(lexer *p)
{
    if(!fbwalls)
    return;

    SLICELOOP1
    if(fbwx[IJ])
    {
    qx(i,j) = 0.0;

        for(k=0;k<p->knoz;++k)
        {
        Qx[PIJK] = fbQx[PIJK];
        qx(i,j) += p->DZN[KP]*Qx[PIJK];
        }
    }

    if(p->j_dir==1)
    SLICELOOP2
    if(fbwy[IJ])
    {
    qy(i,j) = 0.0;

        for(k=0;k<p->knoz;++k)
        {
        Qy[PIJK] = fbQy[PIJK];
        qy(i,j) += p->DZN[KP]*Qy[PIJK];
        }
    }
}
