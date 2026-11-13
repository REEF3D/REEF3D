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

#include"nhflow_hullbody.h"
#include"reefckpt.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"slice.h"
#include"6DOF_obj_nhflow.h"
#include"nhflow_thinbody.h"
#include<cmath>
#include<cstdlib>
#include<sys/stat.h>

// The 6DOF hull as a sharp thin body (X 15 8).
//
// With the direct forcing (X 15 0-7) the cells inside the hull are forced to the body velocity and
// the sigma column of a pierced cell runs through the hull to a level that moves with it: the
// hydrostatic pressure under the keel has the wrong sign, the projection leaks, and the heave and
// pitch added mass of a floating body comes out about 1.5 times potential theory (Ship
// Hydrodynamics, validations 18 and 20). Water columns that end at the hull bottom (X 13 2) fix the
// hydrostatics but not the exchange across the hull edge in the sigma grid (validation 20).
//
// Here the hull is a sharp internal boundary of the undeformed sigma grid, with the service that
// carries the sharp membranes (X 330, nhflow_thinbody.h): the water inside the hull stays in the grid
// as enclosed water that moves with the hull, so the columns keep their free surface (the inner
// level inside the hull, the outer one outside) and the sigma layers run continuously across the hull
// side at their true height. The hull surface is the staircase of the faces between cells on opposite
// sides of the body (the links that cross the surface):
//   sides       cell centres and nodes (the pressure nodes between the cells of a column) by the sign of
//               the body level set d->FB (the ray cast of the stage pose); a point changes its side only
//               when it lies more than a hysteresis distance beyond the surface (no flicker of cells that
//               sit on a moving wall), a node whose two cells lie on the other side takes their side
//   links       the horizontal links between cell centres and the vertical links between nodes on
//               opposite sides are blocked: wall fluxes of the two sides (nhflow_thinbody::flux_hook),
//               the wall velocity in the projection (projection_rhs) and the link mobility
//               1/(1 + a R_n/l) of the Poisson equation (d->MBX, MBY, MBZ; R_n = 1e5 m/s as the sharp
//               membranes), cut cells and pockets of the staircase at the wall velocity (cut_forcing)
//   wall velocity   the component of the rigid-body velocity at the face centre along the link; the
//               staircase then sweeps the volume of the true surface (sum of the projected areas)
//   under the keel  the cells below a crossing of the column are "lower" cells: their hydrostatic head
//               is the outer free surface carried across the footprint (harmonic extension, etaL), and
//               the lid correction makes their net inflow the keel velocity times the column area
//   loads       the outer pressure, rho g eta_L + P, on the outer side of every blocked link times the
//               face area, towards the body; moments about the centre of gravity. The enclosed inner
//               water takes no part in the body dynamics (its pressure is not on the hull). The side faces
//               of the staircase (horizontal links) are taken at the heights of their layer in still water
//               (face_height): the sigma layers stretch with the free surface, the hull does not.
//               The dynamic part is weighted with the hull normal, (nf . n) n, so that the work of the body
//               against it equals the work of the wall fluxes. In the ship frame (surge, sway and yaw fixed,
//               X 11) the horizontal loads are taken on the face normals themselves (plain_horizontal): the
//               body does no horizontal work, and the solver takes the momentum out of the water on the faces,
//               along their normals. With the weighting, the floors at the keel corners took an x part of the
//               stagnation pressure (a 2D box in a current: -12 N per 0.1 m, a thrust) and the steps of a sloped
//               bow lost most of theirs (Ship Hydrodynamics validation 30).
//   enclosed water  its horizontal velocity is the horizontal rigid-body velocity, every stage
//               (enclosed_water). Left free, it sloshed inside the hull. In a current (a hull at speed in the
//               ship frame, B 531) the inner level of the footprint columns oscillated by 2-3 cm at 5 Hz and
//               moved the sigma layers of the lower cells and the keel pressure: a 2D box fixed in a current
//               of 1 m/s had a heave force noise of 135 N (now 2 N), in forced heave a negative damping, free
//               in heave it oscillated by +-5 cm (now +-2.5 cm at 1 m/s, decaying at 0.5 m/s), and the
//               calm-water force of the KVLCC2-like hull in the ship frame was -67 N with a noise of 35 N
//               (now -22 N, 9 N; Ship Hydrodynamics validation 27). Its vertical velocity stays with the
//               solver, so that the inner level follows the keel through the lid correction.
// The direct forcing is off. Pressure scheme: the incremental projection (A 520 2, default; the old pressure in the
// predictor with the link mobilities, nhflow_pjm_corr) or the full projection (A 520 1, on the mesh refinement
// hierarchy). The full projection damps short waves (Ship Hydrodynamics validation 34).
//
// Result (2D box B 1 m, T 0.4 m, forced heave at 3.5 rad/s; Frank close-fit A/m 0.975, B 63.6 per
// 0.1 m span): see the README of validation 21.
//
// Log (rank 0, every X 19 time steps): REEF3D_NHFLOW_6DOF/REEF3D_6DOF_hullbody.dat
//   time, blocked links (x, y, z), lower cells, inner water volume, Fx, Fy, Fz, Mx, My, Mz of the hull pressure.
//
// Mesh refinement (nhflow_amr, G 1 with the zone around the hull, G 12): every grid holds the hull (a patch its own
// nhflow_hullbody and nhflow_thinbody, amr_patch), the hook of the thin body runs before the flux matching of the
// hierarchy, the composite projection takes the rows, the right-hand side and the corrections of every grid with the
// mobilities; the loads are summed over the blocked faces of the leaf grids (loads with the finest grid at the face
// centre, finish_loads by the hull of level 0). The zone keeps the hull and its band inside the finest level (margin of
// at least 4 level-0 cells); one time step for all levels (G 7 0); on several ranks placed patches (G 40 1, the zone
// whole on one rank). A patch rebuilt by a regrid starts its hull afresh (sides from the level set of the stage).
//
// Not covered: submerged bodies (the cells inside would be lower cells of their column; keep the direct
// forcing), porous bodies (X 16), curvilinear grids (G 50), more than one body, subcycling (G 7 1).
//
// X 15 10: the body as a solid with cut faces (stage_cut, loads_cut). The staircase of X 15 8 changes in steps: a cell
// centre or a node that crosses the surface (the sigma grid moves with the free surface, the body moves) switches a
// whole link at once, the projection gets a pressure kick, the free surface a bump, and the next cell switches (a
// submerged ellipse in a current: a front of node switches moving upstream by one column per step until the run blew
// up; the steps also give a form drag of their own, C_D 0.4 at 5 cells across the ellipse, falling with the grid at
// order 0.8). Here every face near the body has its open fraction (aperture) from the level set, linear between the
// face corners (the corners: the level set averaged from the cell centres), so the wall moves continuously:
//   faces       the x and y faces in a lower and an upper half (node k to the cell centre, the centre to node k+1), the
//               horizontal faces through the cell centres (the vertical links of the projection) and at the nodes
//   cells       the water fraction kap of every cell from its corners; the apertures of the faces of a cell with
//               little water are scaled with min(kap/0.25, 1), so that a cell has no open face left when its water
//               vanishes (a sliver of water beyond a fully open face would otherwise be a whole cell of water that
//               disappears at once). A cell without an open face is solid: it moves with the body (also after the
//               projection, reforce); the switch carries no flux. Corners on the surface within 1e-6 dx count as body.
//   momentum    the momentum of a cut cell is that of its water kap V: its right-hand side and its velocity correction
//               are divided by max(kap, kap0) (the links of the Poisson matrix with it); with the whole cell the water
//               of the cut layer was too heavy (added mass of the hemisphere +20 %, now +9 %)
//   wall part   (1 - a) of a face moves with the body: the wall velocity is the body velocity at the face centre (the
//               open parts then see the swept volume of the true surface, -u_b . sum a n A)
//   service     nhflow_thinbody in cut mode: fluxes a F + (1 - a) F(wall), the projection with the apertures as link
//               mobilities and the blended face velocity in its divergence, the head of the columns the body covers
//               (the open fraction of the column top weights eta and the outer level eta_L), see nhflow_thinbody.h
//   loads       the still-water pressure rho g (wd - z) on the STL below the still water level, exact for planar
//               facets (the buoyancy and the zero horizontal force at rest); the dynamic pressure (rho g head + P,
//               at the surface point nearest to the cell centre) on the wall part of every cell, -p sum_f a_f n_f A_f
//               (the divergence theorem on the water of the cell): the same pressure that the projection and the wall
//               fluxes put on the water. The force file gets the loads of the stages in the weights of the
//               Runge-Kutta scheme (finish_loads), the force that changes the momentum of the water over the step:
//               the loads of the last stage alone made the drag of the ellipse 25 % larger than the momentum balance
//               of a control volume around it (now within 5 %)
// Results (Ship Hydrodynamics validation 31): submerged ellipse in a current C_D 0.07 at 5 cells across (X 15 1: 0.41),
// stable as the grid moves through it; 2D box in forced heave at 3.5 rad/s, dx 0.1 m: A/m 0.94, B 66 (Frank 0.975, 64;
// X 15 8 0.96, 76, with ten times the force noise); the box in waves and in a current of 1 m/s free in heave without the
// noise of the staircase; a symmetric box at rest stays symmetric. In the ship frame (validation 33): the 2D box in a
// current of 1 m/s loses less than half the momentum of X 15 8 (12 N against 28-32 N per 0.1 m at dx 0.1 m, 6 N at
// dx 0.05 m, where X 15 8 is unstable); the KVLCC2-like hull at Fr 0.142, dx 0.3 m, heave and pitch free, runs clean
// over 15 s, the flow around it loses 190 N (control volume; X 15 8 about 300 N, tests 49 N with friction): at that
// grid the excess is in the flow, not in the staircase. On the mesh refinement (validation 40) every grid cuts its
// faces with the body and the loads come from the finest grid at each cell and each STL facet (keep): the semicircle
// in forced heave on dx 0.1 m with one level of 0.05 m around it gives the A/m and B of the uniform dx 0.05 m grid.
// Not covered yet: turbulence walls at the cut faces (the faces keep the diffusion), porous bodies, more than one body.

nhflow_hullbody::nhflow_hullbody(lexer *p, fdm_nhf *d, ghostcell *pgc, sixdof_obj_nhflow *pobj)
{
    obj = pobj;

    const int ncell = p->imax*p->jmax*(p->kmax+2);
    const int nnode = p->imax*p->jmax*p->kmaxF;

    p->Darray(side0C,ncell);
    p->Darray(side0N,nnode);
    p->Darray(zk0,p->imax*p->jmax);
    
    // X 15 10: the body as a solid with cut faces
    cut = (p->X15==10);
    p->Darray(solidC,ncell);
    p->Darray(kapC,ncell);

    for(int q=0; q<p->imax*p->jmax; ++q)
    zk0[q] = 1.0e20;

    for(int q=0; q<ncell; ++q)
    side0C[q] = 0.0;

    for(int q=0; q<nnode; ++q)
    side0N[q] = 0.0;

    // mobility fields of the projection, 1 away from the hull (as the sharp membranes, link mode)
    if(d->MBETA==nullptr)
    p->Darray(d->MBETA,ncell);

    if(d->MCHI==nullptr)
    p->Darray(d->MCHI,ncell);

    if(d->MRCX==nullptr)
    p->Darray(d->MRCX,ncell);

    if(d->MRCY==nullptr)
    p->Darray(d->MRCY,ncell);

    if(d->MBX==nullptr)
    {
    p->Darray(d->MBX,ncell);
    p->Darray(d->MBY,ncell);
    p->Darray(d->MBZ,ncell);
    }

    for(int q=0; q<ncell; ++q)
    {
    d->MBETA[q] = 1.0;
    d->MCHI[q] = 0.0;
    d->MRCX[q] = 0.0;
    d->MRCY[q] = 0.0;
    d->MBX[q] = d->MBY[q] = d->MBZ[q] = 1.0;
    }

    d->MPROJ = 1;

    if(d->thinbody==nullptr)
    d->thinbody = new nhflow_thinbody(p,d,pgc);

    tb = d->thinbody;
    
    if(cut)
    tb->enable_cut(p);

    double dmin = p->DXM;

    for(int q=0; q<p->knox; ++q)
    dmin = MIN(dmin,p->DXN[q+marge]);

    if(p->j_dir==1)
    for(int q=0; q<p->knoy; ++q)
    dmin = MIN(dmin,p->DYN[q+marge]);

    hyst = 0.02*pgc->globalmin(dmin);
    // the horizontal loads on the face normals of the staircase: a hull at rest horizontally (X 11: surge, sway and yaw
    // fixed; the ship frame). The vertical loads keep the hull-normal weighting (heave and pitch damping).
    plain_horizontal = (p->X11_u==0 && p->X11_v==0 && p->X11_r==0);

    // the footprint level: in the ship frame, a hull at rest horizontally in a current (B 514); relaxation time [s]
    bool current = false;
    for(int n=0; n<p->B514; ++n)
    if(fabs(p->B514_U[n]) + fabs(p->B514_V[n])>0.0)
    current = true;
    fp_on = plain_horizontal && current;
    fp_tau = 0.5;

    // the lid correction of the lower region drains through its outflow edges (ship frame with a current)
    tb->lid_outflow = fp_on;
    Rn = 1.0e5;

    nbx = nby = nbz = nlink = 0;
    told = marked = false;

    for(int q=0; q<6; ++q)
    Fp[q] = 0.0;

    tstage = tsum = 0.0;
}

nhflow_hullbody::~nhflow_hullbody()
{
    delete [] side0C;
    delete [] side0N;
    delete [] zk0;
    delete [] solidC;
    delete [] kapC;

    if(out.is_open())
    out.close();
}

void nhflow_hullbody::body_velocity(double x, double y, double z, double &u, double &v, double &w) const
{
    // rigid-body velocity at the point (x,y,z)
    const double cx = obj->amr_c(0), cy = obj->amr_c(1), cz = obj->amr_c(2);

    u = obj->amr_u(0) + obj->amr_u(4)*(z - cz) - obj->amr_u(5)*(y - cy);
    v = obj->amr_u(1) + obj->amr_u(5)*(x - cx) - obj->amr_u(3)*(z - cz);
    w = obj->amr_u(2) + obj->amr_u(3)*(y - cy) - obj->amr_u(4)*(x - cx);
}

bool nhflow_hullbody::hull_normal(lexer *p, fdm_nhf *d, int q0, int q1, double &nx, double &ny, double &nz) const
{
    // unit normal of the hull at the face between the cells q0 and q1 (q0 = q1: at the cell), pointing into
    // the body: -grad FB (FB > 0 in the water), the gradient averaged over the two cells; false where the
    // level set has no gradient
    const int qs[2] = {q0,q1};
    nx = ny = nz = 0.0;

    for(int c=0; c<2; ++c)
    {
        const int q = qs[c];
        const int ii = q/(p->jmax*p->kmax) + p->imin;
        const int rem = q - (ii-p->imin)*p->jmax*p->kmax;
        const int jj = rem/p->kmax + p->jmin;
        const int kk = rem - (jj-p->jmin)*p->kmax + p->kmin;
        const int qe = q + p->jmax*p->kmax, qw = q - p->jmax*p->kmax;
        const int qn = q + p->kmax, qs2 = q - p->kmax;
        const double wl = p->WL[(ii-p->imin)*p->jmax + (jj-p->jmin)];

        nx += -(d->FB[qe] - d->FB[qw])/(p->DXP[ii+marge] + p->DXP[ii-1+marge]);

        if(p->j_dir==1)
        ny += -(d->FB[qn] - d->FB[qs2])/(p->DYP[jj+marge] + p->DYP[jj-1+marge]);

        nz += -(d->FB[q+1] - d->FB[q-1])/((p->DZP[kk+marge] + p->DZP[kk-1+marge])*MAX(wl,1.0e-20));
    }

    const double nn = sqrt(nx*nx + ny*ny + nz*nz);

    if(nn<1.0e-20)
    return false;

    nx/=nn;
    ny/=nn;
    nz/=nn;

    return true;
}

// ------------------------------------------------------------------ sides, blocked links, wall velocities
void nhflow_hullbody::stage(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, double *UH, double *VH, double *WH, slice &WL)
{
    const double starttime = pgc->timer();

    if(!told)
    {
        told = true;

        if(p->mpirank==0 && !amr_patch)
        {
        if(cut)
        cout<<"NHFLOW 6DOF X 15 10: the body as a solid with cut faces (open fractions of the faces from the level set, wall part with the body velocity, loads from the pressure on the wall parts)"<<endl;
        else
        cout<<"NHFLOW 6DOF X 15 8: the hull as a sharp thin body (blocked links, wall velocity, loads from the outer pressure)"<<endl;
        }
    }
    
    if(cut)
    {
    stage_cut(p,d,pgc,alpha,UH,VH,WH,WL);
    tstage += pgc->timer() - starttime;
    return;
    }

    const int ncell = p->imax*p->jmax*(p->kmax+2);
    const double a = alpha*p->dt;

    for(int q=0; q<ncell; ++q)
    {
    d->MBETA[q] = 1.0;
    d->MCHI[q] = 0.0;
    d->MBX[q] = d->MBY[q] = d->MBZ[q] = 1.0;
    }

    // the free surface of the footprint columns in the ship frame with a current: relaxed towards the outer level carried
    // under the hull plus the keel displacement
    if(fp_on)
    footprint_level(p,d,pgc,WL,a);

    tb->begin(p);

    // sides of the cell centres: the sign of the level set, with hysteresis
    auto keep = [&](double fb, double old)
    {
        const double s = fb>=0.0 ? 1.0 : -1.0;
        return (old!=0.0 && s!=old && fabs(fb)<hyst) ? old : s;
    };

    LOOP
    {
        tb->sideC[IJK] = keep(d->FB[IJK],side0C[IJK]);
        side0C[IJK] = tb->sideC[IJK];
    }

    pgc->start4V(p,tb->sideC,1);

    // sides of the nodes (node k at the bottom face of cell k), level set interpolated between the cells
    for(i=0; i<p->knox; ++i)
    for(j=0; j<p->knoy; ++j)
    for(k=0; k<=p->knoz; ++k)
    {
        double fb;

        if(k==0)
        fb = d->FB[IJK];
        else if(k==p->knoz)
        fb = d->FB[IJKm1];
        else
        fb = 0.5*(d->FB[IJKm1] + d->FB[IJK]);

        tb->sideN[FIJK] = keep(fb,side0N[FIJK]);

        // a node whose two cells lie on the other side (a wall thinner than a cell) takes their side
        if(k>0 && k<p->knoz && tb->sideC[IJKm1]*tb->sideC[IJK]>0.0 && tb->sideN[FIJK]*tb->sideC[IJK]<0.0)
        tb->sideN[FIJK] = tb->sideC[IJK];

        side0N[FIJK] = tb->sideN[FIJK];
    }

    // blocked links with the wall velocity of the face centre and the link mobility of the projection.
    // Wall velocity: the normal part (u_b . n) n of the rigid-body velocity at the face, n the hull normal
    // from the level set (as the sharp membranes). The staircase faces then sweep the volume of the true
    // surface, and a sloped surface in heave also pushes the water sideways through the vertical faces of
    // its steps; with the plain velocity components the steps act as horizontal pistons between fixed walls
    // (hemisphere in forced heave at KR 1: the radiated wave gave a damping of 0.10 against Hulme's 0.25,
    // with the normal part 0.23)
    auto beta = [&](double l) {return 1.0/(1.0 + a*Rn/MAX(l,1.0e-20));};

    auto normal_vel = [&](int q0, int q1, double x, double y, double z, double &un, double &vn, double &wn)
    {
        double nx,ny,nz;
        double u,v,w;
        body_velocity(x,y,z,u,v,w);

        if(!hull_normal(p,d,q0,q1,nx,ny,nz))
        {
        un=u; vn=v; wn=w;
        return;
        }

        const double ub = u*nx + v*ny + w*nz;

        un = ub*nx;
        vn = ub*ny;
        wn = ub*nz;
    };

    nbx = nby = nbz = 0;

    LOOP
    if(p->wet[IJ]==1)
    {
        const double s = tb->sideC[IJK];
        double u,v,w;

        // +x link to (i+1,j,k)
        if(tb->sideC[Ip1JK]*s<0.0 && p->flag4[Ip1JK]>0 && p->wet[Ip1J]==1)
        {
            normal_vel(IJK,Ip1JK, p->XP[IP]+0.5*p->DXN[IP], p->YP[JP], 0.5*(p->ZSP[IJK]+p->ZSP[Ip1JK]), u,v,w);

            const double b = beta(p->DXP[IP]);
            d->MBX[IJK] = MIN(d->MBX[IJK],b);
            d->MBETA[IJK] = MIN(d->MBETA[IJK],b);
            d->MBETA[Ip1JK] = MIN(d->MBETA[Ip1JK],b);
            tb->bx[IJK] = 1.0;
            tb->ux[IJK] = u;
            ++nbx;
        }

        // +y link to (i,j+1,k)
        if(p->j_dir==1 && tb->sideC[IJp1K]*s<0.0 && p->flag4[IJp1K]>0 && p->wet[IJp1]==1)
        {
            normal_vel(IJK,IJp1K, p->XP[IP], p->YP[JP]+0.5*p->DYN[JP], 0.5*(p->ZSP[IJK]+p->ZSP[IJp1K]), u,v,w);

            const double b = beta(p->DYP[JP]);
            d->MBY[IJK] = MIN(d->MBY[IJK],b);
            d->MBETA[IJK] = MIN(d->MBETA[IJK],b);
            d->MBETA[IJp1K] = MIN(d->MBETA[IJp1K],b);
            tb->by[IJK] = 1.0;
            tb->uy[IJK] = v;
            ++nby;
        }

        // vertical link through the cell: nodes k and k+1
        if(tb->sideN[FIJK]*tb->sideN[FIJKp1]<0.0)
        {
            normal_vel(IJK,IJK, p->XP[IP], p->YP[JP], p->ZSP[IJK], u,v,w);

            const double b = beta(p->DZN[KP]*WL(i,j));
            d->MBZ[IJK] = MIN(d->MBZ[IJK],b);
            d->MBETA[IJK] = MIN(d->MBETA[IJK],b);
            tb->bz[IJK] = 1.0;
            tb->uz[IJK] = w;
            ++nbz;
        }
    }

    pgc->start4V(p,d->MBETA,1);
    pgc->start4V(p,d->MBX,1);
    pgc->start4V(p,d->MBY,1);
    pgc->start4V(p,d->MBZ,1);

    tb->finish(p,d,pgc);

    nlink = pgc->globalisum(nbx+nby+nbz);

    // cells that crossed the hull with the moving sigma grid, cut cells and pockets at the wall velocity
    tb->side_change(p,d,UH,VH,WH,WL);
    tb->cut_forcing(p,d,UH,VH,WH,WL);
    
    // the enclosed water moves horizontally with the hull
    enclosed_water(p,d,pgc,UH,VH,WL);

    marked = true;

    tstage += pgc->timer() - starttime;
}

// ------------------------------------------------------------------ enclosed water
void nhflow_hullbody::enclosed_water(lexer *p, fdm_nhf *d, ghostcell *pgc, double *UH, double *VH, slice &WL)
{
    // cells inside the hull (side -1) that are not cut by a floor (no blocked vertical link, their W is the wall
    // velocity of cut_forcing): U, V of the rigid-body motion at the cell centre
    LOOP
    if(p->wet[IJ]==1 && tb->sideC[IJK]<0.0 && tb->bz[IJK]<0.5)
    {
        double u,v,w;
        body_velocity(p->XP[IP],p->YP[JP],p->ZSP[IJK],u,v,w);
        
        d->U[IJK] = u;
        d->V[IJK] = v;
        UH[IJK] = u*WL(i,j);
        VH[IJK] = v*WL(i,j);
    }
    
    pgc->start4V(p,d->U,10);
    pgc->start4V(p,d->V,11);
}

// ------------------------------------------------------------------ footprint level (ship frame)
// The footprint columns (the free surface of the column inside the hull, nhflow_thinbody::top_inside) carry the level of
// the enclosed water, set by the column's own continuity and the lid correction. Their sigma layers then end at that
// inner level, while the outer columns at the hull edge end at the outer free surface. In a current the bow pile-up
// and the stern trough put the layers below the keel at different heights on the two sides of the hull edge: a 2D box
// free in heave in a current of 1 m/s kept a limit cycle of +-2.5 cm (validation 27), and the heave force of the
// KVLCC2-like hull in the ship frame a noise of 70-200 N.
// In the ship frame with a current the level of a footprint column relaxes, with the time fp_tau, towards the outer
// level carried under the hull (eta_L, the harmonic extension of the outer levels, the head of the lower cells) plus
// the keel displacement since the column entered the footprint:
//     d eta/dt = (eta_L + z_keel - z_keel,0 - eta)/fp_tau
// so that the layers continue across the hull edge, and the lower region still follows the keel. Set to that level at
// every stage instead, the column fought its own continuity: on the 3D hull a column near the keel took a flux of
// several m/s and the time step fell to a third. Only with a current (B 514): at zero speed the level of the enclosed
// water gives the radiation of the heaving hull (2D box, forced heave at 3.5 rad/s: A/m 0.922, B 63.0 per span, Frank
// 0.975, 63.6; with the relaxation 0.986, 78.1). Ship Hydrodynamics validation 30.
void nhflow_hullbody::footprint_level(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL, double a)
{
    const double r = MIN(a/fp_tau, 1.0);

    SLICELOOP4
    {
        const int q = (i-p->imin)*p->jmax + (j-p->jmin);

        if(p->wet[IJ]==1 && tb->top_inside(i,j))
        {
            const double zk = keel_height(p,d,i,j);

            if(zk0[q]>1.0e19)
            zk0[q] = zk;

            const double dz = (zk<1.0e19 && zk0[q]<1.0e19) ? zk - zk0[q] : 0.0;

            d->eta(i,j) += (tb->etaL_at(i,j) + dz - d->eta(i,j))*r;
            WL(i,j) = d->eta(i,j) + d->depth(i,j);
        }
        else
        zk0[q] = 1.0e20;
    }

    pgc->gcsl_start4(p,d->eta,50);
    pgc->gcsl_start4(p,WL,50);
}

// lowest crossing of the hull in the column, from the level set at the cell centres (1e20: none)
double nhflow_hullbody::keel_height(lexer *p, fdm_nhf *d, int ii, int jj) const
{
    int i=ii, j=jj, k;

    for(k=1; k<p->knoz; ++k)
    if(d->FB[IJKm1]>=0.0 && d->FB[IJK]<0.0)
    {
        const double fa = d->FB[IJKm1], fb = d->FB[IJK];
        const double za = p->ZSP[IJKm1], zb = p->ZSP[IJK];
        return za + (zb - za)*fa/(fa - fb);
    }

    return 1.0e20;
}

// ------------------------------------------------------------------ loads: the outer pressure on the blocked links
void nhflow_hullbody::forces(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL, bool finalize)
{
    double F[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    
    loads(p,d,WL,nullptr,F,nullptr);
    finish_loads(p,d,pgc,F,finalize,nullptr);
}

// the loads of the blocked links of this grid and rank, added to F (Fx, Fy, Fz, Mx, My, Mz); with keep only the
// faces whose centre passes keep(x,y) (mesh refinement: the leaf part of the grid), S[0] += the number of those
// faces. No collective.
void nhflow_hullbody::loads(lexer *p, fdm_nhf *d, slice &WL, const std::function<bool(double,double)> *keep, double *F, double *S)
{
    if(cut)
    {
    loads_cut(p,d,WL,keep,F,S);
    return;
    }
    
    double Fx=0.0, Fy=0.0, Fz=0.0, Mx=0.0, My=0.0, Mz=0.0;
    int nf=0;

    if(marked)
    {
    const double rg = p->W1*fabs(p->W22);
    const double *P = d->P;
    const double xc = obj->amr_c(0), yc = obj->amr_c(1), zc = obj->amr_c(2);

    auto add = [&](double x, double y, double z, double fx, double fy, double fz)
    {
        if(keep!=nullptr && !(*keep)(x,y))
        return;
        
        ++nf;
        Fx += fx;
        Fy += fy;
        Fz += fz;
        Mx += (y - yc)*fz - (z - zc)*fy;
        My += (z - zc)*fx - (x - xc)*fz;
        Mz += (x - xc)*fy - (y - yc)*fx;
    };

    // the pressure of the outer cell at the face: the still-water part -rho g z' (z' above the still water
    // level) and the dynamic part rho g head + P (the head of the column above the still water level, the
    // non-hydrostatic pressure)
    auto pstill = [&](int ii, int jj, double z)
    {
        return -rg*(z - d->bed(ii,jj) - d->depth(ii,jj));
    };
    
    auto pdyn = [&](int ii, int jj, int kk, double pn)
    {
        return rg*tb->head(p,d,ii,jj,kk) + pn;
    };
    
    // the force on the blocked face of area A with the face normal nf (unit, into the body): the still-water
    // pressure on the face itself, ps A nf (the buoyancy of the staircase), the dynamic pressure as
    // pd A (nf . n) n with n the hull normal of the wall velocity, so that the work of the body against the
    // dynamic loads equals the work of the wall fluxes on the water (the energy the waves carry away); the
    // face normal itself where the level set has no gradient. The sum over the staircase of A (nf . n) n is
    // the surface vector of the true hull, as the sum of A nf is. (The still-water part is conservative and
    // large: weighted with n it would carry the errors of the normals at the waterline into the heave force.)
    auto face = [&](int q0, int q1, double x, double y, double z, double ps, double pd, double A, double nfx, double nfy, double nfz)
    {
        double nx,ny,nz;
        
        if(hull_normal(p,d,q0,q1,nx,ny,nz))
        {
            const double f = pd*A*(nfx*nx + nfy*ny + nfz*nz);
            
            // ship frame: the horizontal part on the face normal (the momentum the solver takes from the water at the
            // face), the vertical part weighted
            if(plain_horizontal)
            add(x,y,z,(ps+pd)*A*nfx, (ps+pd)*A*nfy, ps*A*nfz + f*nz);
            else
            add(x,y,z,ps*A*nfx + f*nx, ps*A*nfy + f*ny, ps*A*nfz + f*nz);
        }
        else
        add(x,y,z,(ps+pd)*A*nfx,(ps+pd)*A*nfy,(ps+pd)*A*nfz);
    };
    
    // the vertical extent of a side face (horizontal link, layer k of the outer column): the layer at its height in
    // still water, the top layer from there to the free surface. With the sigma heights (layer k * WL) the walls of
    // the staircase stretched with the free surface: the bottom of a hull side, the layer of the keel, rose and fell
    // with the outer level, and the still-water part of the face loads (rho g depth times the face area) changed at
    // first order in the wave elevation where the true hull has none. Bow and stern side differently: KVLCC2-like hull
    // in the ship frame at Fr 0.142 (Ship Hydrodynamics validation 28), forced heave and pitch at zero speed in calm
    // water: mean surge force +19.6 N (+1.9 N on still-water heights), the still-water part of the surge force
    // oscillating by +-90 N (+-1 N); in head waves the added resistance -20 N (+25.5 N, tests 15-17 N); the bow
    // pile-up and the stern trough of a hull in a current did not act on the walls (calm-water force -22 N, now +118
    // N, near the 125-135 N of the direct forcing X 15 1). Height zf of the face centre, width times extent in A.
    auto face_height = [&](int ii, int jj, int kk, double w, double &zf, double &A)
    {
        const double D = d->depth(ii,jj), bd = d->bed(ii,jj);
        
        if(kk<p->knoz-1)
        {
            zf = bd + p->ZP[kk+marge]*D;
            A = w*p->DZN[kk+marge]*D;
        }
        else
        {
            const double z0 = bd + (1.0 - p->DZN[kk+marge])*D, z1 = bd + D + d->eta(ii,jj);
            zf = 0.5*(z0 + z1);
            A = w*MAX(z1 - z0,0.0);
        }
    };
    
    LOOP
    if(p->wet[IJ]==1)
    {
        // x link: the outer cell pushes towards the body, pressure at the face centre
        if(tb->bx[IJK]>0.5)
        {
            const bool outer_i = tb->sideC[IJK]>0.0;
            const int ii = outer_i ? i : i+1;
            double zf, A;
            face_height(ii,j,k,p->DYN[JP],zf,A);

            face(IJK,Ip1JK, p->XP[IP]+0.5*p->DXN[IP], p->YP[JP], zf, pstill(ii,j,zf), pdyn(ii,j,k,tb->pcell(p,P,ii,j,k)), A, (outer_i ? 1.0 : -1.0), 0.0, 0.0);
        }

        if(p->j_dir==1 && tb->by[IJK]>0.5)
        {
            const bool outer_j = tb->sideC[IJK]>0.0;
            const int jj = outer_j ? j : j+1;
            double zf, A;
            face_height(i,jj,k,p->DXN[IP],zf,A);

            face(IJK,IJp1K, p->XP[IP], p->YP[JP]+0.5*p->DYN[JP], zf, pstill(i,jj,zf), pdyn(i,jj,k,tb->pcell(p,P,i,jj,k)), A, 0.0, (outer_j ? 1.0 : -1.0), 0.0);
        }
        
        // vertical link through cell k (the surface between the nodes k and k+1): the cell centres below (kb)
        // and above (ka) the crossing as in nhflow_thinbody::link_force; the pressure at the surface height
        // (level set interpolated between the two cell centres), the non-hydrostatic pressure of the node on
        // the outer side
        if(tb->bz[IJK]>0.5)
        {
            int kb=k-1, ka=k+1;
            int wf[6];
            tb->wall_faces(p,i,j,k,wf);
            
            if(wf[4])
            {
            kb=k-1;
            ka=k;
            }
            else
            if(wf[5])
            {
            kb=k;
            ka=k+1;
            }
            
            kb = MAX(kb,0);
            ka = MIN(ka,p->knoz-1);
            
            const double A = p->DXN[IP]*p->DYN[JP];
            const int qb = (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + kb-p->kmin;
            const int qa = (i-p->imin)*p->jmax*p->kmax + (j-p->jmin)*p->kmax + ka-p->kmin;
            
            const double fb = d->FB[qb], fa = d->FB[qa];
            const double t = (fabs(fb-fa)>1.0e-20) ? MIN(MAX(fb/(fb-fa),0.0),1.0) : 0.5;
            const double zs = p->ZSP[qb] + t*(p->ZSP[qa] - p->ZSP[qb]);
            
            if(tb->sideC[qb]>0.0)
            face(IJK,IJK, p->XP[IP], p->YP[JP], zs, pstill(i,j,zs), pdyn(i,j,kb,P[FIJK]), A, 0.0, 0.0, 1.0);

            if(tb->sideC[qa]>0.0)
            face(IJK,IJK, p->XP[IP], p->YP[JP], zs, pstill(i,j,zs), pdyn(i,j,ka,P[FIJKp1]), A, 0.0, 0.0, -1.0);
        }
    }
    }
    
    F[0] += Fx;
    F[1] += Fy;
    F[2] += Fz;
    F[3] += Mx;
    F[4] += My;
    F[5] += Mz;
    
    if(S!=nullptr)
    S[0] += double(nf);
}

// the loads of all grids and ranks (F: the local sums of loads) into the body; S (mesh refinement): the local
// sums of the blocked faces, lower cells and inner water volume of the leaf grids for the log, nullptr: this grid
void nhflow_hullbody::finish_loads(lexer *p, fdm_nhf *d, ghostcell *pgc, const double *F, bool finalize, const double *S)
{
    double Fx,Fy,Fz,Mx,My,Mz;
    
    Fp[0] = Fx = pgc->globalsum(F[0]);
    Fp[1] = Fy = pgc->globalsum(F[1]);
    Fp[2] = Fz = pgc->globalsum(F[2]);
    Fp[3] = Mx = pgc->globalsum(F[3]);
    Fp[4] = My = pgc->globalsum(F[4]);
    Fp[5] = Mz = pgc->globalsum(F[5]);

    // X 15 10: the printed loads are the stage loads in the weights of the Runge-Kutta scheme (the loads of a stage
    // follow its pressure projection): the force by which the time step changes the momentum of the water. The
    // pressure of the single stages differs by more than the mean load on a body in a current; only the weighted
    // sum enters the flow. The body motion takes the loads of the stage as before.
    if(cut)
    {
        for(int s=0; s<2; ++s)
        for(int q=0; q<6; ++q)
        Fst[s][q] = Fst[s+1][q];
        
        for(int q=0; q<6; ++q)
        Fst[2][q] = Fp[q];
        
        nst = MIN(nst+1,3);
        
        const int ns = p->A534==1 ? 1 : (p->A510==3 ? 3 : (p->A510==2 ? 2 : 1));
        const double w3[3] = {1.0/6.0, 1.0/6.0, 2.0/3.0}, w2[3] = {0.0, 0.5, 0.5}, w1[3] = {0.0, 0.0, 1.0};
        const double *w = ns==3 ? w3 : (ns==2 ? w2 : w1);
        double Fw[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
        
        // after a start or restart with fewer stages stored: their mean
        for(int s=3-MIN(ns,nst); s<3; ++s)
        for(int q=0; q<6; ++q)
        Fw[q] += (nst>=ns ? w[s] : 1.0/double(MIN(ns,nst)))*Fst[s][q];
        
        obj->hull_loads_nhflow(p,pgc,finalize,Fx,Fy,Fz,Mx,My,Mz,Fw);
    }
    
    else
    obj->hull_loads_nhflow(p,pgc,finalize,Fx,Fy,Fz,Mx,My,Mz);

    if(finalize)
    print(p,d,pgc,S);
}

// ------------------------------------------------------------------ regrid (mesh refinement)
int nhflow_hullbody::amr_nstate(lexer *p) const
{
    return p->knoz + (p->knoz+1) + tb->amr_nstate(p) + 2*p->knoz + 1;
}

void nhflow_hullbody::amr_get(lexer *p, fdm_nhf *d, int ii, int jj, double *v)
{
    const int K = p->knoz;
    const int qc = (ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax - p->kmin;
    const int qf = (ii-p->imin)*p->jmax*p->kmaxF + (jj-p->jmin)*p->kmaxF - p->kmin;
    int n=0;
    
    for(int kk=0; kk<K; ++kk)
    v[n++] = side0C[qc+kk];
    
    for(int kk=0; kk<=K; ++kk)
    v[n++] = side0N[qf+kk];
    
    tb->amr_get(p,ii,jj,v+n);
    n += tb->amr_nstate(p);
    
    for(int kk=0; kk<K; ++kk)
    {
    v[n++] = d->MRCX[qc+kk];
    v[n++] = d->MRCY[qc+kk];
    }

    v[n++] = zk0[(ii-p->imin)*p->jmax + (jj-p->jmin)];
}

void nhflow_hullbody::amr_set(lexer *p, fdm_nhf *d, int ii, int jj, const double *v)
{
    const int K = p->knoz;
    const int qc = (ii-p->imin)*p->jmax*p->kmax + (jj-p->jmin)*p->kmax - p->kmin;
    const int qf = (ii-p->imin)*p->jmax*p->kmaxF + (jj-p->jmin)*p->kmaxF - p->kmin;
    int n=0;
    
    for(int kk=0; kk<K; ++kk)
    side0C[qc+kk] = v[n++];
    
    for(int kk=0; kk<=K; ++kk)
    side0N[qf+kk] = v[n++];
    
    tb->amr_set(p,ii,jj,v+n);
    n += tb->amr_nstate(p);
    
    for(int kk=0; kk<K; ++kk)
    {
    d->MRCX[qc+kk] = v[n++];
    d->MRCY[qc+kk] = v[n++];
    }

    zk0[(ii-p->imin)*p->jmax + (jj-p->jmin)] = v[n++];
}

// ------------------------------------------------------------------ log
void nhflow_hullbody::print(lexer *p, fdm_nhf *d, ghostcell *pgc, const double *S)
{
    if(p->count%p->X19!=0)
    return;

    // mesh refinement (S): the blocked faces of the leaf grids in the column of links_z (the split by direction is
    // not kept), links_x and links_y 0
    int nx,ny,nz,nl;
    double vin;
    
    if(S==nullptr)
    {
    nx = pgc->globalisum(nbx);
    ny = pgc->globalisum(nby);
    nz = pgc->globalisum(nbz);
    nl = tb->nlower(p,pgc);
    vin = tb->inner_volume(p,d);
    }
    else
    {
    nx = ny = 0;
    nz = int(pgc->globalsum(S[0])+0.5);
    nl = int(pgc->globalsum(S[1])+0.5);
    vin = pgc->globalsum(S[2]);
    }
    tsum += pgc->globalmax(tstage);

    if(p->mpirank==0)
    {
        if(!out.is_open())
        {
        mkdir("./REEF3D_NHFLOW_6DOF",0777);
        if(reefckpt_open(out,"./REEF3D_NHFLOW_6DOF/REEF3D_6DOF_hullbody.dat",p->I45))
        out<<"time \t links_x \t links_y \t links_z \t lower_cells \t V_inner \t Fx \t Fy \t Fz \t Mx \t My \t Mz"<<endl;
        }

        out<<p->simtime<<" \t "<<nx<<" \t "<<ny<<" \t "<<nz<<" \t "<<nl<<" \t "<<vin
           <<" \t "<<Fp[0]<<" \t "<<Fp[1]<<" \t "<<Fp[2]<<" \t "<<Fp[3]<<" \t "<<Fp[4]<<" \t "<<Fp[5]<<endl;
    }

    tstage = 0.0;
}

// X 15 10: the cells without an open face at the rigid-body velocity after the projection (nhflow_forcing::reforcing)
void nhflow_hullbody::reforce(lexer *p, fdm_nhf *d, ghostcell *pgc, double *UH, double *VH, double *WH, slice &WL)
{
    if(cut && marked)
    solid_velocity(p,d,pgc,UH,VH,WH,WL);
}

// the cells of the solid at the rigid-body velocity of their centre
void nhflow_hullbody::solid_velocity(lexer *p, fdm_nhf *d, ghostcell *pgc, double *UH, double *VH, double *WH, slice &WL)
{
    LOOP
    if(p->wet[IJ]==1 && solidC[IJK]>0.5)
    {
        double u,v,w;
        body_velocity(p->XP[IP],p->YP[JP],p->ZSP[IJK],u,v,w);
        
        d->U[IJK] = u;
        d->V[IJK] = v;
        d->W[IJK] = w;
        UH[IJK] = u*WL(i,j);
        VH[IJK] = v*WL(i,j);
        WH[IJK] = w*WL(i,j);
    }
    
    pgc->start4V(p,d->U,10);
    pgc->start4V(p,d->V,11);
    pgc->start4V(p,d->W,12);
}

// ------------------------------------------------------------------ X 15 10: cut faces
// fluid fraction of a segment, a triangle and a quadrilateral (corners in cyclic order) with the level set linear
// between the corner values (> 0 water); continuous in the corner values
static double cut_seg(double a, double b)
{
    if(a<=0.0 && b<=0.0)
    return 0.0;
    
    if(a>=0.0 && b>=0.0)
    return 1.0;
    
    return (a>0.0 ? a : b)/(fabs(a)+fabs(b));
}

static double cut_tri(double a, double b, double c)
{
    const int np = (a>0.0) + (b>0.0) + (c>0.0);
    
    if(np==3)
    return 1.0;
    
    if(np==0)
    return 0.0;
    
    // the corner alone on its side: v, the two others w1, w2; fraction of its side v^2/((v-w1)(v-w2))
    double v,w1,w2;
    
    if(np==1)
    {
        if(a>0.0)      {v=a; w1=b; w2=c;}
        else if(b>0.0) {v=b; w1=a; w2=c;}
        else           {v=c; w1=a; w2=b;}
        
        return v*v/((v-w1)*(v-w2));
    }
    
    if(a<=0.0)      {v=a; w1=b; w2=c;}
    else if(b<=0.0) {v=b; w1=a; w2=c;}
    else            {v=c; w1=a; w2=b;}
    
    return 1.0 - v*v/((v-w1)*(v-w2));
}

static double cut_quad(double a, double b, double c, double d)
{
    return 0.25*(cut_tri(a,b,c) + cut_tri(a,c,d) + cut_tri(a,b,d) + cut_tri(b,c,d));
}

// the height of the horizontal wall part of cut cell (i,j,k) between its two nodes, as a fraction of the cell: the water
// fraction below a wall above the water (more open at the bottom), one minus it above a wall below the water
double nhflow_hullbody::cut_lam(lexer *p, int i, int j, int k) const
{
    const double a0 = tb->azn[IJK], a1 = k<p->knoz-1 ? tb->azn[IJKp1] : tb->atop_(i,j);
    const double kap = MAX(MIN(kapC[IJK],1.0),0.0);
    
    return a0>a1+1.0e-6 ? kap : (a1>a0+1.0e-6 ? 1.0-kap : 0.5);
}

// the vertical link mobility of cut cell (i,j,k): (1 - lam) a_k+1 + lam a_k
double nhflow_hullbody::cut_zlink(lexer *p, int i, int j, int k) const
{
    const double a0 = tb->azn[IJK], a1 = k<p->knoz-1 ? tb->azn[IJKp1] : tb->atop_(i,j);
    
    if(fabs(a0-a1)<=1.0e-6)
    return 0.5*(a0+a1);
    
    const double lam = cut_lam(p,i,j,k);
    
    return (1.0-lam)*a1 + lam*a0;
}

void nhflow_hullbody::stage_cut(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, double *UH, double *VH, double *WH, slice &WL)
{
    const int ncell = p->imax*p->jmax*(p->kmax+2);
    const double a = alpha*p->dt;
    
    for(int q=0; q<ncell; ++q)
    {
    d->MBETA[q] = 1.0;
    d->MCHI[q] = 0.0;
    d->MBX[q] = d->MBY[q] = d->MBZ[q] = 1.0;
    }
    
    tb->begin(p);
    
    // level set at the cell corners: (x_ii+1/2, y_jj+1/2) at the level of the cell centres kk (vC) or of node kk (vN)
    // a corner on the surface within rounding (a body face that coincides with a grid face) counts as body, on every
    // face alike
    const double ftol = 1.0e-6*p->DXM;
    auto snap = [&](double v) {return fabs(v)<ftol ? 0.0 : v;};
    auto FBc = [&](int ii, int jj, int kk) {return d->FB[tb->cidx(p,ii,jj,kk)];};
    
    auto vC = [&](int ii, int jj, int kk)
    {
        if(p->j_dir==0)
        return snap(0.5*(FBc(ii,jj,kk) + FBc(ii+1,jj,kk)));
        
        return snap(0.25*(FBc(ii,jj,kk) + FBc(ii+1,jj,kk) + FBc(ii,jj+1,kk) + FBc(ii+1,jj+1,kk)));
    };
    
    auto vN = [&](int ii, int jj, int kk)
    {
        if(kk<=0)
        return vC(ii,jj,0);
        
        if(kk>=p->knoz)
        return vC(ii,jj,p->knoz-1);
        
        return snap(0.5*(vC(ii,jj,kk-1) + vC(ii,jj,kk)));
    };
    
    const int jo = p->j_dir==1 ? 1 : 0;   // corner offset in y (2D: the y extent is one cell)
    
    // x faces i+1/2, i = -1 .. knox-1
    for(i=-1; i<p->knox; ++i)
    for(j=0; j<p->knoy; ++j)
    for(k=0; k<p->knoz; ++k)
    {
        if(p->j_dir==0)
        {
        tb->axl[IJK] = cut_seg(vN(i,j,k),vC(i,j,k));
        tb->axh[IJK] = cut_seg(vC(i,j,k),vN(i,j,k+1));
        }
        else
        {
        tb->axl[IJK] = cut_quad(vN(i,j-jo,k),vN(i,j,k),vC(i,j,k),vC(i,j-jo,k));
        tb->axh[IJK] = cut_quad(vC(i,j-jo,k),vC(i,j,k),vN(i,j,k+1),vN(i,j-jo,k+1));
        }
    }
    
    // y faces j+1/2, j = -1 .. knoy-1
    if(p->j_dir==1)
    for(i=0; i<p->knox; ++i)
    for(j=-1; j<p->knoy; ++j)
    for(k=0; k<p->knoz; ++k)
    {
        tb->ayl[IJK] = cut_quad(vN(i-1,j,k),vN(i,j,k),vC(i,j,k),vC(i-1,j,k));
        tb->ayh[IJK] = cut_quad(vC(i-1,j,k),vC(i,j,k),vN(i,j,k+1),vN(i-1,j,k+1));
    }
    
    // horizontal faces through the cell centre (vertical link of the projection) and at node k (bottom of cell k)
    LOOP
    {
        if(p->j_dir==0)
        {
        tb->azc[IJK] = cut_seg(vC(i-1,j,k),vC(i,j,k));
        tb->azn[IJK] = cut_seg(vN(i-1,j,k),vN(i,j,k));
        }
        else
        {
        tb->azc[IJK] = cut_quad(vC(i-1,j-1,k),vC(i,j-1,k),vC(i,j,k),vC(i-1,j,k));
        tb->azn[IJK] = cut_quad(vN(i-1,j-1,k),vN(i,j-1,k),vN(i,j,k),vN(i-1,j,k));
        }
    }
    
    // water fraction of the cells (2D: the cell between its four corners; 3D: the mean of the three mid sections) and
    // the effective apertures: a face next to a cell with little water closes with it, a (kap/kap0), so that a cell
    // whose water vanishes has no open face left when it does (a sliver of water beyond a fully open face would
    // otherwise be a whole cell of water that disappears at once)
    const double kap0 = 0.25;
    
    for(i=-1; i<=p->knox; ++i)
    for(j=(p->j_dir==1?-1:0); j<=(p->j_dir==1?p->knoy:p->knoy-1); ++j)
    for(k=0; k<p->knoz; ++k)
    {
        double kap;
        
        if(p->j_dir==0)
        kap = cut_quad(vN(i-1,j,k),vN(i,j,k),vN(i,j,k+1),vN(i-1,j,k+1));
        else
        {
        auto mx = [&](int jj, int kk) {return 0.5*(vN(i-1,jj,kk)+vN(i,jj,kk));};
        auto my = [&](int ii, int kk) {return 0.5*(vN(ii,j-1,kk)+vN(ii,j,kk));};
        const double sx = cut_quad(mx(j-1,k),mx(j,k),mx(j,k+1),mx(j-1,k+1));
        const double sy = cut_quad(my(i-1,k),my(i,k),my(i,k+1),my(i-1,k+1));
        const double sz = cut_quad(vC(i-1,j-1,k),vC(i,j-1,k),vC(i,j,k),vC(i-1,j,k));
        kap = (sx+sy+sz)/3.0;
        }
        
        kapC[IJK] = kap;
    }
    
    auto gk = [&](int q) {return kap0>0.0 ? MIN(kapC[q]/kap0,1.0) : 1.0;};
    
    for(i=-1; i<p->knox; ++i)
    for(j=0; j<p->knoy; ++j)
    for(k=0; k<p->knoz; ++k)
    {
        const double g = gk(IJK)*gk(Ip1JK);
        tb->axl[IJK] *= g;
        tb->axh[IJK] *= g;
    }
    
    if(p->j_dir==1)
    for(i=0; i<p->knox; ++i)
    for(j=-1; j<p->knoy; ++j)
    for(k=0; k<p->knoz; ++k)
    {
        const double g = gk(IJK)*gk(IJp1K);
        tb->ayl[IJK] *= g;
        tb->ayh[IJK] *= g;
    }
    
    LOOP
    {
        tb->azc[IJK] *= gk(IJK);
        tb->azn[IJK] *= gk(IJK)*(k>0 ? gk(IJKm1) : 1.0);
    }
    
    pgc->start4V(p,tb->azc,1);
    pgc->start4V(p,tb->azn,1);
    
    // open fraction of the free surface of the columns (the face at node knoz, with the top cell's water)
    SLICELOOP4
    {
        k = p->knoz-1;
        double a;
        
        if(p->j_dir==0)
        a = cut_seg(vN(i-1,j,p->knoz),vN(i,j,p->knoz));
        else
        a = cut_quad(vN(i-1,j-1,p->knoz),vN(i,j-1,p->knoz),vN(i,j,p->knoz),vN(i-1,j,p->knoz));
        
        tb->atop(i,j) = a*gk(IJK);
    }
    
    pgc->gcsl_start4(p,tb->atop,1);
    
    auto AX = [&](int q) {return 0.5*(tb->axl[q]+tb->axh[q]);};
    auto AY = [&](int q) {return 0.5*(tb->ayl[q]+tb->ayh[q]);};
    
    // cells without an open face: no water, they move with the body
    
    LOOP
    {
        // (a tolerance: the two sides of a symmetric body differ in the last digits of the level set)
        const double atol = 1.0e-10;
        bool s = AX(IJK)<=atol && AX(Im1JK)<=atol && tb->azc[IJK]<=atol && tb->azn[IJK]<=atol && (k==p->knoz-1 || tb->azn[IJKp1]<=atol);
        
        if(p->j_dir==1)
        s = s && AY(IJK)<=atol && AY(IJm1K)<=atol;
        
        solidC[IJK] = s ? 1.0 : 0.0;
        tb->sideC[IJK] = d->FB[IJK]>=0.0 ? 1.0 : -1.0;
    }
    
    pgc->start4V(p,solidC,1);
    pgc->start4V(p,tb->sideC,1);
    tb->solidC = solidC;
    
    // the momentum of a cut cell is that of its water (nhflow_thinbody::ikap): 1/max(kap, kap0). Below kap0 the faces
    // of a cell close with its water, so that its fluxes, and with them its acceleration, stay bounded
    for(int q=0; q<p->imax*p->jmax*(p->kmax+2); ++q)
    tb->ikap[q] = 1.0;
    
    LOOP
    if(solidC[IJK]<0.5)
    tb->ikap[IJK] = 1.0/MAX(MIN(kapC[IJK],1.0),kap0);
    
    pgc->start4V(p,tb->ikap,1);
    
    // body velocity at the faces (the wall part of a face moves with the body: its flux is the body velocity), the
    // link mobilities of the projection: the apertures, at least the mobility of a closed link
    auto beta = [&](double l) {return 1.0/(1.0 + a*Rn/MAX(l,1.0e-20));};
    double u,v,w;
    
    nbx = nby = nbz = 0;
    
    for(i=-1; i<p->knox; ++i)
    for(j=0; j<p->knoy; ++j)
    for(k=0; k<p->knoz; ++k)
    {
        const double A = AX(IJK);
        
        if(A<1.0)
        {
        body_velocity(p->XP[IP]+0.5*p->DXN[IP], p->YP[JP], p->ZSP[IJK], u,v,w);
        tb->ux[IJK] = u;
        ++nbx;
        }
        
        d->MBX[IJK] = MAX(A, beta(p->DXP[IP]));
    }
    
    if(p->j_dir==1)
    for(i=0; i<p->knox; ++i)
    for(j=-1; j<p->knoy; ++j)
    for(k=0; k<p->knoz; ++k)
    {
        const double A = AY(IJK);
        
        if(A<1.0)
        {
        body_velocity(p->XP[IP], p->YP[JP]+0.5*p->DYN[JP], p->ZSP[IJK], u,v,w);
        tb->uy[IJK] = v;
        ++nby;
        }
        
        d->MBY[IJK] = MAX(A, beta(p->DYP[JP]));
    }
    
    LOOP
    {
        if(tb->azc[IJK]<1.0 || solidC[IJK]>0.5)
        {
        body_velocity(p->XP[IP], p->YP[JP], p->ZSP[IJK], u,v,w);
        tb->uz[IJK] = w;
        ++nbz;
        }
        
        // the vertical link of the cell: its vertical pressure force as a finite volume, with the pressure on the
        // horizontal wall part at the height of the wall between the two nodes, p_w = P_k + lam (P_k+1 - P_k) (lam = kap
        // below a wall above the water, 1 - kap above a wall below it): ((1 - lam) a_k+1 + lam a_k)(P_k+1 - P_k), a link
        // with that mobility in the correction (wcorr) and the matrix alike; the loads take the same wall pressure
        // (loads_cut)
        d->MBZ[IJK] = MAX(cut_zlink(p,i,j,k), beta(p->DZN[KP]*WL(i,j)));
    }
    
    pgc->start4V(p,d->MBX,1);
    pgc->start4V(p,d->MBY,1);
    pgc->start4V(p,d->MBZ,1);
    
    LOOP
    {
        double b = MIN(MIN(d->MBX[IJK],d->MBX[Im1JK]),d->MBZ[IJK]);
        
        if(p->j_dir==1)
        b = MIN(b,MIN(d->MBY[IJK],d->MBY[IJm1K]));
        
        d->MBETA[IJK] = MIN(b,1.0);
    }
    
    pgc->start4V(p,d->MBETA,1);
    
    tb->finish_cut(p,d,pgc);
    
    nlink = pgc->globalisum(nbx+nby+nbz);
    
    solid_velocity(p,d,pgc,UH,VH,WH,WL);
    
    marked = true;
}

// loads of the cut faces: the body takes -p sum_f a_f n_f A_f of every cell with water (the wall part of the cell
// surface, divergence theorem), p the pressure of the cell centre (still water, head, non-hydrostatic)
void nhflow_hullbody::loads_cut(lexer *p, fdm_nhf *d, slice &WL, const std::function<bool(double,double)> *keep, double *F, double *S)
{
    double Fx=0.0, Fy=0.0, Fz=0.0, Mx=0.0, My=0.0, Mz=0.0;
    int nf=0;
    
    if(marked)
    {
    const double rg = p->W1*fabs(p->W22);
    const double *P = d->P;
    const double xc = obj->amr_c(0), yc = obj->amr_c(1), zc = obj->amr_c(2);
    
    auto AX = [&](int q) {return 0.5*(tb->axl[q]+tb->axh[q]);};
    auto AY = [&](int q) {return 0.5*(tb->ayl[q]+tb->ayh[q]);};
    
    LOOP
    if(p->wet[IJ]==1 && solidC[IJK]<0.5)
    {
        const double ax = AX(IJK) - AX(Im1JK);
        const double ay = p->j_dir==1 ? AY(IJK) - AY(IJm1K) : 0.0;
        const double az = (k<p->knoz-1 ? tb->azn[IJKp1] : tb->atop(i,j)) - tb->azn[IJK];
        
        // the horizontal part of the open surface of the cell, sum a n A, as the finite volume of the pressure force of
        // the solver (nhflow_cut_gx): the x / y faces with their heights, the sigma faces with their slope
        const double He = 0.5*(WL(i,j)+WL(i+1,j))*p->DZN[KP], Hw = 0.5*(WL(i-1,j)+WL(i,j))*p->DZN[KP];
        const double a0 = tb->azn[IJK], a1 = k<p->knoz-1 ? tb->azn[IJKp1] : tb->atop(i,j);
        const double ox = (AX(IJK)*He - AX(Im1JK)*Hw)*p->DYN[JP] - p->DXN[IP]*p->DYN[JP]*WL(i,j)*(p->sigx[FIJK]*a0 - p->sigx[FIJKp1]*a1);
        double oy = 0.0;
        
        if(p->j_dir==1)
        {
        const double Hn = 0.5*(WL(i,j)+WL(i,j+1))*p->DZN[KP], Hs = 0.5*(WL(i,j-1)+WL(i,j))*p->DZN[KP];
        oy = (AY(IJK)*Hn - AY(IJm1K)*Hs)*p->DXN[IP] - p->DXN[IP]*p->DYN[JP]*WL(i,j)*(p->sigy[FIJK]*a0 - p->sigy[FIJKp1]*a1);
        }
        
        if(ax==0.0 && ay==0.0 && az==0.0 && fabs(ox)<1.0e-14 && fabs(oy)<1.0e-14)
        continue;
        
        // the pressure at the surface point nearest to the cell centre, x_c - FB grad FB (the still-water part there,
        // so that the load does not jump by rho g dz when the wall passes from one cell layer to the next)
        double gx = (d->FB[Ip1JK] - d->FB[Im1JK])/(p->DXP[IP] + p->DXP[IM1]);
        double gy = p->j_dir==1 ? (d->FB[IJp1K] - d->FB[IJm1K])/(p->DYP[JP] + p->DYP[JM1]) : 0.0;
        double gz = (d->FB[IJKp1] - d->FB[IJKm1])/((p->DZP[KP] + p->DZP[KM1])*MAX(WL(i,j),1.0e-20));
        const double gn = sqrt(gx*gx + gy*gy + gz*gz);
        
        if(gn>1.0e-20)
        {
        gx/=gn;
        gy/=gn;
        gz/=gn;
        }
        
        const double fb = d->FB[IJK];
        const double x = p->XP[IP] - fb*gx, y = p->YP[JP] - fb*gy, z = p->ZSP[IJK] - fb*gz;
        
        if(keep!=nullptr && !(*keep)(p->XP[IP],p->YP[JP]))
        continue;
        
        // the dynamic part (head above the still water level, non-hydrostatic pressure); the still-water part below.
        // Horizontally both as the solver puts them on the water: the head as the wall flux of the hydrostatic part,
        // g (e^2/2 + e h) per layer fraction (nhflow_thinbody::cut_faces), the non-hydrostatic pressure on the open surface
        // (ox, oy); vertically the pressure on the open fractions of the bottom and the top face
        const double eh = tb->head(p,d,i,j,k), pn = tb->pcell(p,P,i,j,k);
        const double hh = rg*(0.5*eh*eh + eh*d->depth(i,j));
        
        const double fx = -hh*ax*p->DYN[JP]*p->DZN[KP] - pn*ox;
        const double fy = -hh*ay*p->DXN[IP]*p->DZN[KP] - pn*oy;
        const double fz = -(rg*eh + P[FIJK] + cut_lam(p,i,j,k)*(P[FIJKp1]-P[FIJK]))*az*p->DXN[IP]*p->DYN[JP];
        
        ++nf;
        Fx += fx;
        Fy += fy;
        Fz += fz;
        Mx += (y - yc)*fz - (z - zc)*fy;
        My += (z - zc)*fx - (x - xc)*fz;
        Mz += (x - xc)*fy - (y - yc)*fx;
    }
    }
    
    // still-water pressure rho g (wd - z) on the STL below the still water level, exact for planar facets (edge-midpoint
    // rule on the clipped facets): the buoyancy and the zero horizontal force of a body at rest come out exactly, not
    // with the discretisation error of the cut faces; triangles owned by the rank of their centroid, on the mesh
    // refinement by the finest grid there
    if(marked)
    {
    const double rg = p->W1*fabs(p->W22);
    const double xc = obj->amr_c(0), yc = obj->amr_c(1), zc = obj->amr_c(2);
    double **tx = obj->amr_tri(0), **ty = obj->amr_tri(1), **tz = obj->amr_tri(2);
    const double wd = p->wd;
    
    // orientation of the STL: the normals of the facets point out of the body when the enclosed volume is positive
    double vol=0.0;
    
    for(int n=0; n<obj->amr_tricount(); ++n)
    vol += (tx[n][0]*(ty[n][1]*tz[n][2]-tz[n][1]*ty[n][2]) - ty[n][0]*(tx[n][1]*tz[n][2]-tz[n][1]*tx[n][2])
          + tz[n][0]*(tx[n][1]*ty[n][2]-ty[n][1]*tx[n][2]))/6.0;
    
    const double sgn = vol<0.0 ? -1.0 : 1.0;
    
    for(int n=0; n<obj->amr_tricount(); ++n)
    {
        const double cxx = (tx[n][0]+tx[n][1]+tx[n][2])/3.0, cyy = (ty[n][0]+ty[n][1]+ty[n][2])/3.0;
        
        // mesh refinement: the facet on the finest grid at its centroid (and on the rank of that grid) only; the
        // subdomain of the rank otherwise
        if(keep!=nullptr)
        {
            if(!(*keep)(cxx,p->j_dir==1 ? cyy : 0.5*(p->originy+p->endy)))
            continue;
        }
        else
        {
            if(!(cxx>=p->originx && cxx<p->endx))
            continue;
            
            if(p->j_dir==1 && !(cyy>=p->originy && cyy<p->endy))
            continue;
        }

        
        double nx = (ty[n][1]-ty[n][0])*(tz[n][2]-tz[n][0]) - (ty[n][2]-ty[n][0])*(tz[n][1]-tz[n][0]);
        double ny = (tx[n][2]-tx[n][0])*(tz[n][1]-tz[n][0]) - (tx[n][1]-tx[n][0])*(tz[n][2]-tz[n][0]);
        double nz = (tx[n][1]-tx[n][0])*(ty[n][2]-ty[n][0]) - (tx[n][2]-tx[n][0])*(ty[n][1]-ty[n][0]);
        const double nn = sqrt(nx*nx+ny*ny+nz*nz);
        
        if(nn<1.0e-20)
        continue;
        
        nx*=sgn/nn; ny*=sgn/nn; nz*=sgn/nn;
        
        // clip at z = wd (keep z <= wd); in 2D the end caps too (no x or z part when they are planar; with them the
        // surface stays closed, so that the horizontal force at rest is zero to round-off)
        double px[4],py[4],pz[4];
        int np=0;
        
        for(int q=0; q<3; ++q)
        {
            const int r = (q+1)%3;
            const double za = tz[n][q], zb = tz[n][r];
            
            if(za<=wd)
            {
            px[np]=tx[n][q]; py[np]=ty[n][q]; pz[np]=za; ++np;
            }
            
            if((za<wd && zb>wd) || (za>wd && zb<wd))
            {
            const double t = (wd-za)/(zb-za);
            px[np]=tx[n][q]+t*(tx[n][r]-tx[n][q]); py[np]=ty[n][q]+t*(ty[n][r]-ty[n][q]); pz[np]=wd; ++np;
            }
        }
        
        for(int q=1; q+1<np; ++q)
        {
            const int a=0, b=q, c=q+1;
            const double ux_=px[b]-px[a], uy_=py[b]-py[a], uz_=pz[b]-pz[a];
            const double vx_=px[c]-px[a], vy_=py[c]-py[a], vz_=pz[c]-pz[a];
            const double cx_=uy_*vz_-uz_*vy_, cy_=uz_*vx_-ux_*vz_, cz_=ux_*vy_-uy_*vx_;
            const double As = 0.5*sqrt(cx_*cx_+cy_*cy_+cz_*cz_);
            
            const int e0[3]={a,b,c}, e1[3]={b,c,a};
            
            for(int m=0; m<3; ++m)
            {
                const double mx = 0.5*(px[e0[m]]+px[e1[m]]), my = 0.5*(py[e0[m]]+py[e1[m]]), mz = 0.5*(pz[e0[m]]+pz[e1[m]]);
                const double ps = rg*(wd - mz);
                const double fx = -ps*nx*As/3.0, fy = -ps*ny*As/3.0, fz = -ps*nz*As/3.0;
                
                Fx += fx;
                Fy += fy;
                Fz += fz;
                Mx += (my - yc)*fz - (mz - zc)*fy;
                My += (mz - zc)*fx - (mx - xc)*fz;
                Mz += (mx - xc)*fy - (my - yc)*fx;
            }
        }
    }
    }
    
    F[0] += Fx;
    F[1] += Fy;
    F[2] += Fz;
    F[3] += Mx;
    F[4] += My;
    F[5] += Mz;
    
    if(S!=nullptr)
    S[0] += double(nf);
}
