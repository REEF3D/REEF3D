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
//               water takes no part in the body dynamics (its pressure is not on the hull).
// The direct forcing is off. The pressure scheme is the full projection (A 520 1).
//
// Result (2D box B 1 m, T 0.4 m, forced heave at 3.5 rad/s; Frank close-fit A/m 0.975, B 63.6 per
// 0.1 m span): see the README of validation 21.
//
// Log (rank 0, every X 19 time steps): REEF3D_NHFLOW_6DOF/REEF3D_6DOF_hullbody.dat
//   time, blocked links (x, y, z), lower cells, inner water volume, Fx, Fy, Fz, Mx, My, Mz of the hull pressure.
//
// Not covered: submerged bodies (the cells inside would be lower cells of their column; keep the direct
// forcing), porous bodies (X 16), curvilinear grids (G 50), mesh refinement, more than one body.

nhflow_hullbody::nhflow_hullbody(lexer *p, fdm_nhf *d, ghostcell *pgc, sixdof_obj_nhflow *pobj)
{
    obj = pobj;

    const int ncell = p->imax*p->jmax*(p->kmax+2);
    const int nnode = p->imax*p->jmax*p->kmaxF;

    p->Darray(side0C,ncell);
    p->Darray(side0N,nnode);

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

    double dmin = p->DXM;

    for(int q=0; q<p->knox; ++q)
    dmin = MIN(dmin,p->DXN[q+marge]);

    if(p->j_dir==1)
    for(int q=0; q<p->knoy; ++q)
    dmin = MIN(dmin,p->DYN[q+marge]);

    hyst = 0.02*pgc->globalmin(dmin);
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

        if(p->mpirank==0)
        cout<<"NHFLOW 6DOF X 15 8: the hull as a sharp thin body (blocked links, wall velocity, loads from the outer pressure)"<<endl;
    }

    const int ncell = p->imax*p->jmax*(p->kmax+2);
    const double a = alpha*p->dt;

    for(int q=0; q<ncell; ++q)
    {
    d->MBETA[q] = 1.0;
    d->MCHI[q] = 0.0;
    d->MBX[q] = d->MBY[q] = d->MBZ[q] = 1.0;
    }

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

    marked = true;

    tstage += pgc->timer() - starttime;
}

// ------------------------------------------------------------------ loads: the outer pressure on the blocked links
void nhflow_hullbody::forces(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL, bool finalize)
{
    double Fx=0.0, Fy=0.0, Fz=0.0, Mx=0.0, My=0.0, Mz=0.0;

    if(marked)
    {
    const double rg = p->W1*fabs(p->W22);
    const double *P = d->P;
    const double xc = obj->amr_c(0), yc = obj->amr_c(1), zc = obj->amr_c(2);

    auto add = [&](double x, double y, double z, double fx, double fy, double fz)
    {
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
            add(x,y,z,ps*A*nfx + f*nx, ps*A*nfy + f*ny, ps*A*nfz + f*nz);
        }
        else
        add(x,y,z,(ps+pd)*A*nfx,(ps+pd)*A*nfy,(ps+pd)*A*nfz);
    };
    
    LOOP
    if(p->wet[IJ]==1)
    {
        // x link: the outer cell pushes towards the body, pressure at the face centre
        if(tb->bx[IJK]>0.5)
        {
            const bool outer_i = tb->sideC[IJK]>0.0;
            const int ii = outer_i ? i : i+1;
            const int q = outer_i ? IJK : Ip1JK;
            const double A = p->DYN[JP]*p->DZN[KP]*WL(ii,j);

            face(IJK,Ip1JK, p->XP[IP]+0.5*p->DXN[IP], p->YP[JP], p->ZSP[q], pstill(ii,j,p->ZSP[q]), pdyn(ii,j,k,tb->pcell(p,P,ii,j,k)), A, (outer_i ? 1.0 : -1.0), 0.0, 0.0);
        }

        if(p->j_dir==1 && tb->by[IJK]>0.5)
        {
            const bool outer_j = tb->sideC[IJK]>0.0;
            const int jj = outer_j ? j : j+1;
            const int q = outer_j ? IJK : IJp1K;
            const double A = p->DXN[IP]*p->DZN[KP]*WL(i,jj);

            face(IJK,IJp1K, p->XP[IP], p->YP[JP]+0.5*p->DYN[JP], p->ZSP[q], pstill(i,jj,p->ZSP[q]), pdyn(i,jj,k,tb->pcell(p,P,i,jj,k)), A, 0.0, (outer_j ? 1.0 : -1.0), 0.0);
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
    
    Fp[0] = Fx = pgc->globalsum(Fx);
    Fp[1] = Fy = pgc->globalsum(Fy);
    Fp[2] = Fz = pgc->globalsum(Fz);
    Fp[3] = Mx = pgc->globalsum(Mx);
    Fp[4] = My = pgc->globalsum(My);
    Fp[5] = Mz = pgc->globalsum(Mz);

    obj->hull_loads_nhflow(p,pgc,finalize,Fx,Fy,Fz,Mx,My,Mz);

    if(finalize)
    print(p,d,pgc);
}

// ------------------------------------------------------------------ log
void nhflow_hullbody::print(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    if(p->count%p->X19!=0)
    return;

    const int nx = pgc->globalisum(nbx), ny = pgc->globalisum(nby), nz = pgc->globalisum(nbz);
    const int nl = tb->nlower(p,pgc);
    const double vin = tb->inner_volume(p,d);
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
