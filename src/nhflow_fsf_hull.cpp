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

#include"nhflow_fsf_hull.h"
#include"reefckpt.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"
#include"6DOF.h"
#include"6DOF_obj_nhflow.h"
#include"ioflow.h"
#include<cmath>
#include<sys/stat.h>

// Floating body as the top of the water columns (X 13 2, direct forcing X 10 1/2, NHFLOW RK2/RK3).
// Prototype (Ship Hydrodynamics, validation 20), step 1 of the plan of the added-mass summary.
//
// With the direct forcing (X 13 0/1) the sigma grid of a column pierced by the hull runs through the
// body to a level that moves with the hull. The hydrostatic part of the pressure under the keel then
// has the wrong sign, the non-hydrostatic pressure has to make up a jump at the hull edge, and the
// water next to the hull is carried along with it: the heave and pitch added mass of a floating body
// comes out about 1.5 times potential theory (validation 18). Submerged bodies are right.
//
// Here the columns under the hull end at the hull bottom, as in SWASH (Rijnsdorp, Wolgamot, Zijlema
// 2022): the sigma layers run from the bed to the hull bottom S(x,y,t), with the full vertical
// resolution under the keel, no body cells and no forcing.
//
//   footprint   the columns whose centre lies under the hull (vertical ray through the trimesh at
//               the pose of the stage, sixdof_obj_nhflow::keel_nhflow): S is the lowest crossing.
//               A column is covered when the keel is below its level (free surface or piezometric head).
//   level       WL = S - bed, eta = S - bed - depth, detadt = dS/dt = w_s - u_s dS/dx - v_s dS/dy
//               (u_s, v_s, w_s the body velocity at the keel point), so that the kinematic
//               condition at the top of the column (nhflow_fsf_f::kinematic_fsf) and the grid
//               velocity (omega_update) are those of the hull. The water level update of the
//               continuity equation is overwritten; the volume it would have changed is logged.
//   zeta        the level of the hydrostatic pressure terms (momentum flux g (eta^2/2 + eta h),
//               bed slope source, HLL jump of the continuity flux): the free surface in the free
//               columns, the piezometric head under the hull, continuous across the hull edge.
//               d->zeta points to it; eta keeps the physical top of the column.
//               The flux form gives g (zeta + h) dzeta/dx per layer; the footprint cells get the
//               source g (zeta - eta) dzeta/dx so that the hydrostatic acceleration is
//               g WL dzeta/dx with WL = S - bed (upgrad).
//   correct     after the momentum stage: the depth-integrated continuity of a footprint column,
//               dS/dt + div Q = 0, is a constraint on zeta. With Q' = Q - alpha dt g WL grad(dzeta)
//               it is the 2D Poisson equation div(alpha dt g WL grad dzeta) = dS/dt + div Q over the
//               footprint, dzeta = 0 in the free columns (Dirichlet), Neumann at walls, solved by
//               SOR. UH, VH of the footprint cells and their free neighbours get the correction,
//               zeta += dzeta. The non-hydrostatic projection keeps P = 0 at the top of all columns.
//   forces      hull bottom: p = rho g (zeta - eta) on the footprint columns, with the keel slope
//               for the horizontal components; hull sides (faces between footprint and free
//               columns): the hydrostatic and non-hydrostatic pressure of the free column between
//               the keel and its free surface. The loads go to the body (hull_loads_nhflow).
//
// Result (2D box B 1 m, T 0.4 m, forced heave, omega 3.5 rad/s, Frank close-fit A/m 0.975,
// B 63.6 per 0.1 m span): A/m 1.13 at dx 0.1 m and 1.04 at dx 0.05 m (direct forcing 1.50 / 1.35),
// but the damping is 131-132 at both resolutions and the radiated wave 0.78 of the hull amplitude
// (Frank 0.53): the hull radiates twice the energy. Cause: the top of the column jumps at the hull
// edge, and the sigma layers of the footprint column (ending at the keel) exchange mass, momentum
// and pressure with the layers of the free column at the same sigma index, i.e. up to a draft higher
// up - the water under the keel edge is linked to the free surface next to the hull and drives it
// like a piston. A layer-matched exchange across the edge needs links between different layer
// indices, which the 7-point Poisson stencil of the solvers does not hold; cutting the links of the
// upper band (wall) seals the water under the keel edge (A/m 2.6), a deferred correction diverges.
// The rigid sharp membrane bag (X 330, mobility sharp, structure rigid) keeps the sigma grid
// continuous across the hull side and treats the keel as a blocked internal link (staircase): it gives
// A/m 1.13 / 1.06 and B 49 / 55 at dx 0.1 / 0.05 m, converging to the potential-flow values. See the
// README of validation 20 for the way forward (hull as a sharp thin body of the 6DOF geometry).
//
// Log (rank 0, every X 19 time steps): REEF3D_NHFLOW_6DOF/REEF3D_6DOF_hull.dat
//   time, footprint columns, edge columns, SOR sweeps since the last line, SOR factor, max
//   continuity residual of the footprint before and after the correction, volume replaced by the
//   prescribed level since the last line and over the run, Fz and Fx of the hull pressure.
//
// Not supported: the semi-implicit scheme (A 510 4), curvilinear grids (G 50), mesh refinement,
// porous bodies (X 16), nets (X 320), membranes (X 330), ship load models (X 350), bodies at rest
// below the surface. Columns entering or leaving the footprint (waterline moving across columns)
// switch the column top between the free surface and the keel without a ramp.

nhflow_fsf_hull::nhflow_fsf_hull(lexer *p) : cover(p), S(p), St(p), zeta(p), dzeta(p), R(p), ratio(p), usk(p), vsk(p), wsk(p)
{
    mode = (p->X13==2 && (p->X10==1 || p->X10==2) && p->X16==0 && (p->A510==2 || p->A510==3) && p->geo_curv==0) ? 1 : 0;
    told = false;
    pending = false;
    fresh = false;

    gcval_eta = 50;

    if(p->F50==1)
	gcval_eta = 51;

    if(p->F50==2)
	gcval_eta = 52;

    if(p->F50==3)
	gcval_eta = 53;

    if(p->F50==4)
	gcval_eta = 54;

    nfoot = nedge = 0;
    itstep = 0;
    omega = 1.0;
    rmax0 = rmax1 = 0.0;
    dVstep = dVsum = 0.0;
    Fz_p = Fx_p = 0.0;
    tstep = ttot = 0.0;

    SLICELOOP4
    {
    cover(i,j) = 0.0;
    S(i,j) = 0.0;
    St(i,j) = 0.0;
    zeta(i,j) = 0.0;
    dzeta(i,j) = 0.0;
    R(i,j) = 0.0;
    ratio(i,j) = 1.0;
    usk(i,j) = vsk(i,j) = wsk(i,j) = 0.0;
    }
}

nhflow_fsf_hull::~nhflow_fsf_hull()
{
    if(out.is_open())
    out.close();
}

bool nhflow_fsf_hull::wall_x(lexer *p, int ii) const
{
    return (ii<0 && p->nb1<0) || (ii>=p->knox && p->nb4<0);
}

bool nhflow_fsf_hull::wall_y(lexer *p, int jj) const
{
    return p->j_dir==1 && ((jj<0 && p->nb3<0) || (jj>=p->knoy && p->nb2<0));
}

double nhflow_fsf_hull::gsum(ghostcell *pgc, double v)
{
    return pgc->globalsum(v);
}

int nhflow_fsf_hull::gsum_i(ghostcell *pgc, int v)
{
    return pgc->globalisum(v);
}

// ------------------------------------------------------------------ footprint and level of the stage
void nhflow_fsf_hull::level(lexer *p, fdm_nhf *d, ghostcell *pgc, sixdof *p6dof, ioflow *pflow, int iter, slice &WL, double alpha)
{
    if(!told)
    {
        told = true;

        if(p->mpirank==0)
        {
            if(mode==1)
            cout<<"NHFLOW 6DOF X 13 2: water columns end at the hull bottom, hull pressure from the piezometric head"<<endl;

            if(mode==0 && p->X13==2)
            cout<<"NHFLOW 6DOF X 13 2 ignored: needs X 10 1 or 2, A 510 2 or 3, no porous body (X 16), Cartesian grid"<<endl;
        }
    }

    pending = false;

    if(mode==0)
    return;

    starttime = pgc->timer();

    // the level a column is compared with: the free surface of the free columns (the stage value),
    // the piezometric head of the footprint columns (last solve)
    SLICELOOP4
    if(cover(i,j)<0.5)
    zeta(i,j) = d->eta(i,j);

    // hull at the pose of this stage (the body is advanced first): keel height, cover and the body
    // velocity at the keel point of every column
    const bool got = (p6dof!=nullptr) && p6dof->keel_nhflow(p,d,pgc,iter,zeta,S,cover,usk,vsk,wsk);

    if(!got)
    {
        if(p->mpirank==0)
        cout<<"NHFLOW 6DOF X 13 2: no keel from the body, mode switched off"<<endl;

        mode = 0;
        return;
    }

    fresh = true;

    pgc->gcsl_start4(p,cover,1);
    pgc->gcsl_start4(p,S,1);

    footprint(p,d,pgc,WL);

    // dS/dt of the footprint columns: the keel slopes one-sided inside the footprint
    const int nf = int(fi.size());

    for(int n=0; n<nf; ++n)
    {
        i = fi[n];
        j = fj[n];

        double sx=0.0, sy=0.0;

        const bool cE = cover(i+1,j)>0.5 && !wall_x(p,i+1);
        const bool cW = cover(i-1,j)>0.5 && !wall_x(p,i-1);

        if(cE && cW)
        sx = (S(i+1,j) - S(i-1,j))/(p->DXP[IP] + p->DXP[IM1]);
        else if(cE)
        sx = (S(i+1,j) - S(i,j))/p->DXP[IP];
        else if(cW)
        sx = (S(i,j) - S(i-1,j))/p->DXP[IM1];

        if(p->j_dir==1)
        {
        const bool cN = cover(i,j+1)>0.5 && !wall_y(p,j+1);
        const bool cS = cover(i,j-1)>0.5 && !wall_y(p,j-1);

        if(cN && cS)
        sy = (S(i,j+1) - S(i,j-1))/(p->DYP[JP] + p->DYP[JM1]);
        else if(cN)
        sy = (S(i,j+1) - S(i,j))/p->DYP[JP];
        else if(cS)
        sy = (S(i,j) - S(i,j-1))/p->DYP[JM1];
        }

        St(i,j) = wsk(i,j) - usk(i,j)*sx - vsk(i,j)*sy;
    }

    // level, depth and grid velocity of the footprint columns: the top of the column is the keel
    double dV=0.0;

    for(int n=0; n<nf; ++n)
    {
        i = fi[n];
        j = fj[n];

        const double wl0 = WL(i,j);
        const double wl1 = S(i,j) - d->bed(i,j);

        // column entering the footprint: the piezometric head starts from its free-surface level
        if(ratio(i,j)<0.0)
        zeta(i,j) = d->eta(i,j);

        dV += (wl1 - wl0)*p->DXN[IP]*p->DYN[JP];

        d->eta(i,j) = wl1 - d->depth(i,j);
        WL(i,j) = wl1;
        d->detadt(i,j) = St(i,j);

        ratio(i,j) = (wl0>p->A544 && wl1>p->A544) ? wl1/wl0 : 1.0;

        if(ratio(i,j)!=1.0)
        pending = true;
    }

    dVstep += gsum(pgc,dV);

    // the level of the hydrostatic pressure terms
    SLICELOOP4
    if(cover(i,j)<0.5)
    zeta(i,j) = d->eta(i,j);

    d->zeta = &zeta;

    pgc->gcsl_start4(p,d->eta,gcval_eta);
    pgc->gcsl_start4(p,WL,gcval_eta);
    pgc->gcsl_start4(p,d->detadt,1);
    pgc->gcsl_start4(p,zeta,gcval_eta);
    pgc->gcsl_start4(p,St,1);

    tstep += pgc->timer() - starttime;
}

void nhflow_fsf_hull::footprint(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL)
{
    // footprint columns (cover from the ray cast, keel below the water level) and the free columns
    // next to them; ratio < 0 marks a column that enters the footprint in this stage
    fi.clear();
    fj.clear();
    ei.clear();
    ej.clear();

    const size_t nsl = size_t(p->imax)*size_t(p->jmax);

    if(was.size()!=nsl)
    was.assign(nsl,0);

    SLICELOOP4
    {
        ratio(i,j) = 1.0;

        if(p->wet[IJ]==0)
        cover(i,j) = 0.0;

        if(cover(i,j)>0.5)
        {
        fi.push_back(i);
        fj.push_back(j);

        // column entering the footprint (free in the last stage)
        if(was[IJ]==0)
        ratio(i,j) = -1.0;
        }

        was[IJ] = cover(i,j)>0.5 ? 1 : 0;
    }

    pgc->gcsl_start4(p,cover,1);

    SLICELOOP4
    if(cover(i,j)<0.5 && p->wet[IJ]==1)
    {
        bool near = (cover(i+1,j)>0.5 && !wall_x(p,i+1)) || (cover(i-1,j)>0.5 && !wall_x(p,i-1));

        if(p->j_dir==1)
        near = near || (cover(i,j+1)>0.5 && !wall_y(p,j+1)) || (cover(i,j-1)>0.5 && !wall_y(p,j-1));

        if(near)
        {
        ei.push_back(i);
        ej.push_back(j);
        }
    }

    nfoot = gsum_i(pgc,int(fi.size()));
    nedge = gsum_i(pgc,int(ei.size()));
}

// ------------------------------------------------------------------ momentum of the changed columns
void nhflow_fsf_hull::momentum(lexer *p, double *UH, double *VH, double *WH)
{
    // the velocities of the columns whose depth level() changed stay as they are
    if(mode==0 || !pending)
    return;

    const int nf = int(fi.size());

    for(int n=0; n<nf; ++n)
    {
        i = fi[n];
        j = fj[n];

        const double r = ratio(i,j);

        if(r!=1.0 && r>0.0)
        KLOOP
        {
        UH[IJK] *= r;
        VH[IJK] *= r;
        WH[IJK] *= r;
        }
    }

    pending = false;
}

// ------------------------------------------------------------------ hydrostatic source of the footprint
void nhflow_fsf_hull::upgrad(lexer *p, fdm_nhf *d, slice &WL)
{
    // flux form: g (zeta + h) dzeta/dx per layer; the footprint column has the depth S - bed = eta + h
    if(mode==0 || !fresh)
    return;

    const double g = fabs(p->W22);
    const int nf = int(fi.size());

    for(int n=0; n<nf; ++n)
    {
        i = fi[n];
        j = fj[n];

        if(p->wet[IJ]==0)
        continue;

        const double zx = (zeta(i+1,j) - zeta(i-1,j))/(p->DXP[IP] + p->DXP[IM1]);
        const double zy = p->j_dir==1 ? (zeta(i,j+1) - zeta(i,j-1))/(p->DYP[JP] + p->DYP[JM1]) : 0.0;
        const double dz = zeta(i,j) - d->eta(i,j);

        KLOOP
        {
        d->F[IJK] += g*dz*zx;
        d->G[IJK] += g*dz*zy;
        }
    }
}

// ------------------------------------------------------------------ one-sided eta slopes at the footprint edge
bool nhflow_fsf_hull::edge_slope(lexer *p, fdm_nhf *d, int ii, int jj, double &sx, double &sy)
{
    // column (ii,jj) next to the edge of the footprint: the slopes of eta from the neighbours on
    // the same side of the edge only (the top of the column jumps at the edge); false elsewhere
    if(mode==0 || !fresh)
    return false;

    i = ii;
    j = jj;

    const bool c0 = cover(i,j)>0.5;
    const bool cE = cover(i+1,j)>0.5;
    const bool cW = cover(i-1,j)>0.5;
    bool edge = (cE!=c0 || cW!=c0);

    bool cN=c0, cS=c0;

    if(p->j_dir==1)
    {
    cN = cover(i,j+1)>0.5;
    cS = cover(i,j-1)>0.5;
    edge = edge || cN!=c0 || cS!=c0;
    }

    if(!edge)
    return false;

    if(cE==c0 && cW==c0)
    sx = (d->eta(i+1,j) - d->eta(i-1,j))/(p->DXP[IP] + p->DXP[IM1]);
    else if(cE==c0 && !wall_x(p,i+1))
    sx = (d->eta(i+1,j) - d->eta(i,j))/p->DXP[IP];
    else if(cW==c0 && !wall_x(p,i-1))
    sx = (d->eta(i,j) - d->eta(i-1,j))/p->DXP[IM1];
    else
    sx = 0.0;

    sy = 0.0;

    if(p->j_dir==1)
    {
    if(cN==c0 && cS==c0)
    sy = (d->eta(i,j+1) - d->eta(i,j-1))/(p->DYP[JP] + p->DYP[JM1]);
    else if(cN==c0 && !wall_y(p,j+1))
    sy = (d->eta(i,j+1) - d->eta(i,j))/p->DYP[JP];
    else if(cS==c0 && !wall_y(p,j-1))
    sy = (d->eta(i,j) - d->eta(i,j-1))/p->DYP[JM1];
    }

    return true;
}

void nhflow_fsf_hull::slope(lexer *p, fdm_nhf *d)
{
    // the top of the column jumps between a footprint column and a free column: the sigma layers of
    // each side are tilted and curved by their own top only (Ex, Ey; Exx, Eyy within two cells of
    // the edge from the unfiltered eta of the same side)
    if(mode==0 || !fresh)
    return;

    SLICELOOP4
    if(p->wet[IJ]==1)
    {
        double sx=0.0, sy=0.0;

        if(edge_slope(p,d,i,j,sx,sy))
        {
        d->Ex(i,j) = sx;
        d->Ey(i,j) = sy;
        }

        const bool c0 = cover(i,j)>0.5;
        bool near=false;

        for(int q=-2; q<=2 && !near; ++q)
        if(q!=0 && !wall_x(p,i+q) && (cover(i+q,j)>0.5)!=c0)
        near = true;

        if(p->j_dir==1)
        for(int q=-2; q<=2 && !near; ++q)
        if(q!=0 && !wall_y(p,j+q) && (cover(i,j+q)>0.5)!=c0)
        near = true;

        if(near)
        {
            const bool cE = cover(i+1,j)>0.5, cW = cover(i-1,j)>0.5;

            d->Exx(i,j) = (cE==c0 && cW==c0 && !wall_x(p,i+1) && !wall_x(p,i-1))
                        ? ((d->eta(i+1,j) - d->eta(i,j))/p->DXP[IP] - (d->eta(i,j) - d->eta(i-1,j))/p->DXP[IM1])/p->DXN[IP] : 0.0;

            if(p->j_dir==1)
            {
            const bool cN = cover(i,j+1)>0.5, cS = cover(i,j-1)>0.5;

            d->Eyy(i,j) = (cN==c0 && cS==c0 && !wall_y(p,j+1) && !wall_y(p,j-1))
                        ? ((d->eta(i,j+1) - d->eta(i,j))/p->DYP[JP] - (d->eta(i,j) - d->eta(i,j-1))/p->DYP[JM1])/p->DYN[JP] : 0.0;
            }
        }
    }
}

// ------------------------------------------------------------------ implicit piezometric head
void nhflow_fsf_hull::correct(lexer *p, fdm_nhf *d, ghostcell *pgc, int iter, double alpha, double *UH, double *VH, double *WH, slice &WL)
{
    if(mode==0 || !fresh || nfoot==0)
    return;

    starttime = pgc->timer();

    const double a = alpha*p->dt*fabs(p->W22);
    const int nf = int(fi.size());

    // the stage momentum of the halo (the divergence reads the neighbour columns)
    pgc->start4V(p,UH,14);
    pgc->start4V(p,VH,15);

    // depth-integrated flux divergence of the footprint columns with the face flux 0.5 (Q_P + Q_N)
    // (sum over the layers with the sigma thickness), Q = 0 through walls
    auto divQ = [&](int ii, int jj)
    {
        i = ii;
        j = jj;

        double qe=0.0, qw=0.0, qn=0.0, qs=0.0, qp=0.0, qpv=0.0;
        double qN=0.0, qS=0.0, qE=0.0, qW=0.0;

        KLOOP
        {
        qp  += p->DZN[KP]*UH[IJK];
        qpv += p->DZN[KP]*VH[IJK];
        qE  += p->DZN[KP]*UH[Ip1JK];
        qW  += p->DZN[KP]*UH[Im1JK];

        if(p->j_dir==1)
        {
        qN += p->DZN[KP]*VH[IJp1K];
        qS += p->DZN[KP]*VH[IJm1K];
        }
        }

        qe = wall_x(p,i+1) || p->wet[Ip1J]==0 ? 0.0 : 0.5*(qp + qE);
        qw = wall_x(p,i-1) || p->wet[Im1J]==0 ? 0.0 : 0.5*(qp + qW);

        double div = (qe - qw)/p->DXN[IP];

        if(p->j_dir==1)
        {
        qn = wall_y(p,j+1) || p->wet[IJp1]==0 ? 0.0 : 0.5*(qpv + qN);
        qs = wall_y(p,j-1) || p->wet[IJm1]==0 ? 0.0 : 0.5*(qpv + qS);

        div += (qn - qs)/p->DYN[JP];
        }

        return div;
    };

    double r0=0.0;

    SLICELOOP4
    {
    R(i,j) = 0.0;
    dzeta(i,j) = 0.0;
    }

    for(int n=0; n<nf; ++n)
    {
        R(fi[n],fj[n]) = St(fi[n],fj[n]) + divQ(fi[n],fj[n]);
        r0 = MAX(r0,fabs(R(fi[n],fj[n])));
    }

    rmax0 = MAX(rmax0,pgc->globalmax(r0));

    pgc->gcsl_start4(p,dzeta,1);

    solve(p,d,pgc,a,WL);

    // correction of the layer momentum: the mean of the face gradients of the cell
    auto apply = [&](const vector<int> &ci, const vector<int> &cj)
    {
        const int nc = int(ci.size());

        for(int n=0; n<nc; ++n)
        {
            i = ci[n];
            j = cj[n];

            if(p->wet[IJ]==0)
            continue;

            double gx=0.0, gy=0.0;

            if(!wall_x(p,i+1) && p->wet[Ip1J]==1)
            gx += 0.5*(WL(i,j) + WL(i+1,j))*(dzeta(i+1,j) - dzeta(i,j))/p->DXP[IP];

            if(!wall_x(p,i-1) && p->wet[Im1J]==1)
            gx += 0.5*(WL(i,j) + WL(i-1,j))*(dzeta(i,j) - dzeta(i-1,j))/p->DXP[IM1];

            if(p->j_dir==1)
            {
            if(!wall_y(p,j+1) && p->wet[IJp1]==1)
            gy += 0.5*(WL(i,j) + WL(i,j+1))*(dzeta(i,j+1) - dzeta(i,j))/p->DYP[JP];

            if(!wall_y(p,j-1) && p->wet[IJm1]==1)
            gy += 0.5*(WL(i,j) + WL(i,j-1))*(dzeta(i,j) - dzeta(i,j-1))/p->DYP[JM1];
            }

            const double wl = WL(i,j)>p->A544 ? WL(i,j) : 1.0e20;

            KLOOP
            {
            UH[IJK] -= 0.5*a*gx;
            VH[IJK] -= 0.5*a*gy;

            d->U[IJK] = UH[IJK]/wl;
            d->V[IJK] = VH[IJK]/wl;
            }
        }
    };

    apply(fi,fj);
    apply(ei,ej);

    for(int n=0; n<nf; ++n)
    zeta(fi[n],fj[n]) += dzeta(fi[n],fj[n]);

    pgc->gcsl_start4(p,zeta,gcval_eta);
    pgc->start4V(p,UH,14);
    pgc->start4V(p,VH,15);
    pgc->start4V(p,d->U,10);
    pgc->start4V(p,d->V,11);

    // residual after the correction
    double r1=0.0;

    for(int n=0; n<nf; ++n)
    r1 = MAX(r1,fabs(St(fi[n],fj[n]) + divQ(fi[n],fj[n])));

    rmax1 = MAX(rmax1,pgc->globalmax(r1));

    tstep += pgc->timer() - starttime;
}

void nhflow_fsf_hull::solve(lexer *p, fdm_nhf *d, ghostcell *pgc, double a, slice &WL)
{
    // div(a WL grad dzeta) = R over the footprint, dzeta = 0 in the free columns, Neumann at walls:
    // SOR over the footprint columns (finite volume, non-uniform grids). The halo is exchanged every
    // sweep when the footprint spans ranks; the residual is checked every ncheck sweeps.
    const int nf = int(fi.size());
    const int nrank = pgc->globalisum(nf>0 ? 1 : 0);
    const int itmax = 2000;
    const int ncheck = 4;

    double rref=0.0;
    for(int n=0; n<nf; ++n)
    rref = MAX(rref,fabs(R(fi[n],fj[n])));
    rref = pgc->globalmax(rref);

    const double tol = 1.0e-8*MAX(rref,1.0e-12);

    const double N = (p->j_dir==1) ? sqrt(double(nfoot)) : double(nfoot);
    omega = 2.0/(1.0 + sin(PI/(N + 1.0)));
    omega = MIN(MAX(omega,1.0),1.9);

    auto add = [&](int ii, int jj, double w, double &s, double &ws)
    {
        if(p->wet[(ii-p->imin)*p->jmax + (jj-p->jmin)]==0)
        return;

        // free column: Dirichlet 0, the coefficient stays in the diagonal
        s  += w*(cover(ii,jj)>0.5 ? dzeta(ii,jj) : 0.0);
        ws += w;
    };

    int it;

    for(it=0; it<itmax; ++it)
    {
        double rmax = 0.0;

        for(int n=0; n<nf; ++n)
        {
            i = fi[n];
            j = fj[n];

            double s=0.0, ws=0.0;

            if(!wall_x(p,i+1))
            add(i+1,j,a*0.5*(WL(i,j)+WL(i+1,j))/(p->DXN[IP]*p->DXP[IP]),s,ws);

            if(!wall_x(p,i-1))
            add(i-1,j,a*0.5*(WL(i,j)+WL(i-1,j))/(p->DXN[IP]*p->DXP[IM1]),s,ws);

            if(p->j_dir==1)
            {
            if(!wall_y(p,j+1))
            add(i,j+1,a*0.5*(WL(i,j)+WL(i,j+1))/(p->DYN[JP]*p->DYP[JP]),s,ws);

            if(!wall_y(p,j-1))
            add(i,j-1,a*0.5*(WL(i,j)+WL(i,j-1))/(p->DYN[JP]*p->DYP[JM1]),s,ws);
            }

            if(ws>0.0)
            {
            const double res = (s - R(i,j))/ws - dzeta(i,j);
            dzeta(i,j) += omega*res;
            rmax = MAX(rmax,fabs(res*ws));
            }
        }

        if(nrank>1)
        pgc->gcsl_start4(p,dzeta,1);

        if((it+1)%ncheck==0)
        {
            if(nrank>1)
            rmax = pgc->globalmax(rmax);

            if(rmax<tol)
            {
            ++it;
            break;
            }
        }
    }

    if(nrank<=1)
    pgc->gcsl_start4(p,dzeta,1);

    itstep += pgc->globalimax(it);
}

// ------------------------------------------------------------------ hull loads
void nhflow_fsf_hull::forces(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL, sixdof_obj_nhflow *obj, bool finalize)
{
    double Fx=0.0, Fy=0.0, Fz=0.0, Mx=0.0, My=0.0, Mz=0.0;
    double Fzb=0.0, Fxs=0.0;

    if(mode==1 && fresh)
    {
    const double g = fabs(p->W22);
    const double rho = p->W1;
    const double xc = obj->amr_c(0), yc = obj->amr_c(1), zc = obj->amr_c(2);

    auto moment = [&](double x, double y, double z, double fx, double fy, double fz)
    {
        Mx += (y - yc)*fz - (z - zc)*fy;
        My += (z - zc)*fx - (x - xc)*fz;
        Mz += (x - xc)*fy - (y - yc)*fx;
    };

    // hull bottom: p = rho g (zeta - eta) on the footprint columns, force p (-dS/dx, -dS/dy, 1) dA
    const int nf = int(fi.size());

    for(int n=0; n<nf; ++n)
    {
        i = fi[n];
        j = fj[n];

        const double pt = rho*g*(zeta(i,j) - d->eta(i,j));
        const double dA = p->DXN[IP]*p->DYN[JP];

        double sx=0.0, sy=0.0;

        const bool cE = cover(i+1,j)>0.5 && !wall_x(p,i+1);
        const bool cW = cover(i-1,j)>0.5 && !wall_x(p,i-1);

        if(cE && cW)
        sx = (S(i+1,j) - S(i-1,j))/(p->DXP[IP] + p->DXP[IM1]);
        else if(cE)
        sx = (S(i+1,j) - S(i,j))/p->DXP[IP];
        else if(cW)
        sx = (S(i,j) - S(i-1,j))/p->DXP[IM1];

        if(p->j_dir==1)
        {
        const bool cN = cover(i,j+1)>0.5 && !wall_y(p,j+1);
        const bool cS = cover(i,j-1)>0.5 && !wall_y(p,j-1);

        if(cN && cS)
        sy = (S(i,j+1) - S(i,j-1))/(p->DYP[JP] + p->DYP[JM1]);
        else if(cN)
        sy = (S(i,j+1) - S(i,j))/p->DYP[JP];
        else if(cS)
        sy = (S(i,j) - S(i,j-1))/p->DYP[JM1];
        }

        const double fx = -pt*sx*dA;
        const double fy = -pt*sy*dA;
        const double fz = pt*dA;

        Fx += fx;
        Fy += fy;
        Fz += fz;
        Fzb += fz;

        moment(p->XP[IP], p->j_dir==1 ? p->YP[JP] : yc, S(i,j), fx,fy,fz);
    }

    // hull sides: the faces between a footprint column and a free column, from the keel of the
    // footprint column to the free surface of the free column, with the pressure of the free column
    // rho g (top - z) + P(z); the water pushes the hull away from the free column
    auto side = [&](int ii, int jj, int ni, int nj, double nx, double ny, double xf, double yf, double len)
    {
        const double zlo = S(ii,jj);
        i = ni;
        j = nj;
        const double ztop = d->bed(i,j) + WL(i,j);

        if(ztop<=zlo)
        return;

        for(k=0; k<p->knoz; ++k)
        {
            const double z0 = p->ZSN[FIJK];
            const double z1 = p->ZSN[FIJKp1];
            const double za = MAX(z0,zlo);
            const double zb = MIN(z1,ztop);

            if(zb<=za)
            continue;

            const double hs = rho*g*(ztop*(zb - za) - 0.5*(zb*zb - za*za));
            const double pn = 0.5*(d->P[FIJK] + d->P[FIJKp1])*(zb - za);
            const double f = (hs + pn)*len;

            Fx += f*nx;
            Fy += f*ny;
            Fxs += f*nx;

            moment(xf,yf,0.5*(za+zb),f*nx,f*ny,0.0);
        }
    };

    for(int n=0; n<nf; ++n)
    {
        const int ii = fi[n], jj = fj[n];
        i = ii;
        j = jj;

        const double yp = p->j_dir==1 ? p->YP[JP] : yc;

        if(cover(i+1,j)<0.5 && !wall_x(p,i+1) && p->wet[Ip1J]==1)
        side(ii,jj,ii+1,jj,-1.0,0.0,p->XP[IP]+0.5*p->DXN[IP],yp,p->DYN[JP]);

        if(cover(i-1,j)<0.5 && !wall_x(p,i-1) && p->wet[Im1J]==1)
        side(ii,jj,ii-1,jj, 1.0,0.0,p->XP[IP]-0.5*p->DXN[IP],yp,p->DYN[JP]);

        if(p->j_dir==1)
        {
        if(cover(i,j+1)<0.5 && !wall_y(p,j+1) && p->wet[IJp1]==1)
        side(ii,jj,ii,jj+1,0.0,-1.0,p->XP[IP],yp+0.5*p->DYN[JP],p->DXN[IP]);

        if(cover(i,j-1)<0.5 && !wall_y(p,j-1) && p->wet[IJm1]==1)
        side(ii,jj,ii,jj-1,0.0, 1.0,p->XP[IP],yp-0.5*p->DYN[JP],p->DXN[IP]);
        }
    }
    }

    Fx = pgc->globalsum(Fx);
    Fy = pgc->globalsum(Fy);
    Fz = pgc->globalsum(Fz);
    Mx = pgc->globalsum(Mx);
    My = pgc->globalsum(My);
    Mz = pgc->globalsum(Mz);

    Fz_p = pgc->globalsum(Fzb);
    Fx_p = pgc->globalsum(Fxs);

    obj->hull_loads_nhflow(p,pgc,finalize,Fx,Fy,Fz,Mx,My,Mz);
}

// ------------------------------------------------------------------ log
void nhflow_fsf_hull::print(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL)
{
    if(mode==0 || p->count%p->X19!=0)
    return;

    const double dVs = dVstep;

    dVsum += dVs;
    ttot += pgc->globalmax(tstep);

    if(p->mpirank==0)
    {
        if(!out.is_open())
        {
        mkdir("./REEF3D_NHFLOW_6DOF",0777);
        if(reefckpt_open(out,"./REEF3D_NHFLOW_6DOF/REEF3D_6DOF_hull.dat",p->I45))
        out<<"time \t columns \t edge \t sweeps \t omega \t res_before \t res_after \t dV_step \t dV_sum \t Fz_bottom \t Fx_sides"<<endl;
        }

        out<<p->simtime<<" \t "<<nfoot<<" \t "<<nedge<<" \t "<<itstep<<" \t "<<omega<<" \t "<<rmax0<<" \t "<<rmax1
           <<" \t "<<dVs<<" \t "<<dVsum<<" \t "<<Fz_p<<" \t "<<Fx_p<<endl;
    }

    rmax0 = rmax1 = 0.0;
    dVstep = 0.0;
    tstep = 0.0;
    itstep = 0;
}
