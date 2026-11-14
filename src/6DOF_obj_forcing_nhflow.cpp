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

#include"6DOF_obj_nhflow.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"vrans_definitions.h"
#include"nhflow_bodyflux.h"

void sixdof_obj_nhflow::update_forcing_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, 
                             double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, slice &fe, int iter)
{
    // porous floating body: Darcy-Forchheimer resistance instead of rigid direct forcing
    if(p->X16==1)
    {
    update_forcing_nhflow_porous(p,d,pgc,U,V,W,FX,FY,FZ,WL,iter);
    return;
    }
    
    forcing_kernel_nhflow(p,d,pgc,U,V,W,FX,FY,FZ,WL,iter,false);
}

void sixdof_obj_nhflow::update_forcing_nhflow_wavemaker(lexer *p, fdm_nhf *d, ghostcell *pgc, 
                             double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, slice &fe, int iter)
{
    // wavemaker (X 10 4): the same direct forcing (X 15 0/1/2) as the floating body, with the
    // paddle velocity of the stage as the target
    forcing_kernel_nhflow(p,d,pgc,U,V,W,FX,FY,FZ,WL,iter,true);
}

void sixdof_obj_nhflow::target_velocity(lexer *p, bool wavemaker, double &uf, double &vf, double &wf)
{
    // target velocity of the direct forcing at the cell (i,j,k)
    if(wavemaker)
    {
        uf = wm_u(p,p->pos_x(),p->pos_z());
        vf = 0.0;
        wf = 0.0;
        return;
    }
    
    uf = u_fb(0) + u_fb(4)*(p->pos_z() - c_(2)) - u_fb(5)*(p->pos_y() - c_(1));
    vf = u_fb(1) + u_fb(5)*(p->pos_x() - c_(0)) - u_fb(3)*(p->pos_z() - c_(2));
    wf = u_fb(2) + u_fb(3)*(p->pos_y() - c_(1)) - u_fb(4)*(p->pos_x() - c_(0));
}

void sixdof_obj_nhflow::forcing_kernel_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, 
                             double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, int iter, bool wavemaker)
{
    // X 13 2: no body cells in the water columns, no direct forcing
    // X 15 8: the hull is a sharp thin body (blocked links), no direct forcing
    if(d->hull!=nullptr || d->hullbody!=nullptr)
    return;
    
    // X 15 11: direct forcing with cell and face fractions and the body in the pressure projection
    // (the wavemaker, X 10 4, with the forcing of X 15 4)
    if(p->X15==11 && !wavemaker)
    {
    forcing_kernel_fraction(p,d,pgc,U,V,W,FX,FY,FZ,WL,iter,wavemaker);
    return;
    }
    
    // Calculate forcing fields
    double H, uf, vf, wf;
    
    double du,dv,dw,udotn;
    double beta;   // 0.0 = free slip (tangential preserved), 1.0 = no slip
    
    if(p->X15==0)
    beta=1.0;
    
    if(p->X15>=1)
    beta=0.0;
    
    
    if(p->X15==0)
    LOOP
    {
        H = Hsolidface_nhflow(p,d,0,0,0);
        
        target_velocity(p,wavemaker,uf,vf,wf);
         
        d->FHB[IJK] = MIN(d->FHB[IJK] + H, 1.0); 
        
        FX[IJK] += H*(uf - U[IJK])/(alpha[iter]*p->dt);
        FY[IJK] += H*(vf - V[IJK])/(alpha[iter]*p->dt);
        FZ[IJK] += H*(wf - W[IJK])/(alpha[iter]*p->dt);
    }
    
    // X 15 3: as X 15 1, but only the cells inside the body are forced. A body moving through the
    // grid at speed (towed or self-propelled ship) otherwise pushes the water in the forcing band ahead
    // of the bow up to the body velocity; that water passes under the hull and is carried along with it.
    // X 15 5: as X 15 1 with the sharp switch at FB = 0 (full forcing for FB <= 0, normal-weighted outside),
    // the forcing before the shell ramp was introduced (9ce616dab)
    // X 15 6 / 7: the forcing of X 15 1 / 3; the body faces are reflective walls in the HLL fluxes (nhflow_HLL.cpp)
    if(p->X15==1 || p->X15==3 || p->X15==5 || p->X15==6 || p->X15==7)
    LOOP
    {
        H = Hsolidface_nhflow(p,d,0,0,0);
        
        target_velocity(p,wavemaker,uf,vf,wf);
         
        d->FHB[IJK] = MIN(d->FHB[IJK] + H, 1.0); 
        
        if((p->X15==3 || p->X15==7) && d->FB[IJK]>=0.0)
        continue;
        
    // Normal vectors calculation 
        if(p->geo_curv==1)
        fbnormal_curv(p,d,nx,ny);
        else
        {
		nx = -(d->FB[Ip1JK] - d->FB[Im1JK])/(p->DXP[IP] + p->DXP[IM1])
            - 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
        
        if(p->j_dir==0)
        ny = 0.0;
        
        if(p->j_dir==1)
		ny = -(d->FB[IJp1K] - d->FB[IJm1K])/(p->DYP[JP] + p->DYP[JM1])
            - 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
        }
            
		nz = -(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP]*WL(i,j) + p->DZP[KM1]*WL(i,j));

		norm = sqrt(nx*nx + ny*ny + nz*nz);
                
		nx /= norm > 1.0e-20 ? norm : 1.0e20;
		ny /= norm > 1.0e-20 ? norm : 1.0e20;
		nz /= norm > 1.0e-20 ? norm : 1.0e20;
        
        
        // full forcing inside the body, normal-weighted forcing (tangential velocity kept) outside,
        // ramped over a thin shell |FB| < dsh around the surface. A sharp switch at FB = 0 jumps
        // for cell centres on or next to the surface (a flat side through a row of cell centres):
        // the smallest rotation moves the row partly inside and partly outside, and the jump in
        // the forcing kicks a symmetric hull into sway and yaw.
        // X 15 5: the sharp switch at FB = 0 (2D cases and grid-aligned hulls, where the kick does not
        // arise; reproduces the floating-body motions of the code before the ramp)
        double sb;
        
        if(p->X15==5)
        sb = d->FB[IJK]<=0.0 ? 1.0 : 0.0;
        
        else
        {
        const double dsh = 0.25*Hpsi_nhflow(p,d);
        const double tb = MIN(MAX(0.5*(1.0 - d->FB[IJK]/dsh), 0.0), 1.0);
        sb = tb*tb*(3.0 - 2.0*tb);
        }
        
        FX[IJK] += (sb + (1.0-sb)*fabs(nx))*H*(uf - U[IJK])/(alpha[iter]*p->dt);
        FY[IJK] += (sb + (1.0-sb)*fabs(ny))*H*(vf - V[IJK])/(alpha[iter]*p->dt);
        FZ[IJK] += (sb + (1.0-sb)*fabs(nz))*H*(wf - W[IJK])/(alpha[iter]*p->dt);
    
    }
    
    
    if(p->X15==2)
    LOOP
    {
        H = Hsolidface_nhflow(p,d,0,0,0);
        
        d->FHB[IJK] = MIN(d->FHB[IJK] + H, 1.0); 
        
        target_velocity(p,wavemaker,uf,vf,wf);
        
        du = uf - U[IJK];
        dv = vf - V[IJK];
        dw = wf - W[IJK];
        
        if(d->FB[IJK] < 0.0)
        {
            // interior of the body: enforce the full rigid-body velocity
            FX[IJK] += H*du/(alpha[iter]*p->dt);
            FY[IJK] += H*dv/(alpha[iter]*p->dt);
            FZ[IJK] += H*dw/(alpha[iter]*p->dt);
        }
         
        else
        {
            // Normal vectors calculation 
            if(p->geo_curv==1)
            fbnormal_curv(p,d,nx,ny);
            else
            {
            nx = -(d->FB[Ip1JK] - d->FB[Im1JK])/(p->DXP[IP] + p->DXP[IM1])
                - 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
                
            if(p->j_dir==0)
            ny = 0.0;
            
            if(p->j_dir==1)
            ny = -(d->FB[IJp1K] - d->FB[IJm1K])/(p->DYP[JP] + p->DYP[JM1])
                - 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
            }
                
            nz = -(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP]*WL(i,j) + p->DZP[KM1]*WL(i,j));

            norm = sqrt(nx*nx + ny*ny + nz*nz);
                    
            if(norm > 1.0e-10)
            {
            nx /= norm;
            ny /= norm;
            nz /= norm;
            
            udotn = nx*du + ny*dv + nz*dw;
            
            FX[IJK] += H*(beta*du + (1.0-beta)*udotn*nx)/(alpha[iter]*p->dt);
            FY[IJK] += H*(beta*dv + (1.0-beta)*udotn*ny)/(alpha[iter]*p->dt);
            FZ[IJK] += H*(beta*dw + (1.0-beta)*udotn*nz)/(alpha[iter]*p->dt);
            }
            
            else
            {
            // degenerate gradient (flat level set): fall back to full forcing
            FX[IJK] += H*du/(alpha[iter]*p->dt);
            FY[IJK] += H*dv/(alpha[iter]*p->dt);
            FZ[IJK] += H*dw/(alpha[iter]*p->dt);
            }
        }
    
    }
    
    // X 15 4: normal projection outside the body (as X 15 2), ramped over the thin shell |FB| < dsh of X 15 1:
    //   F = H [sb du + (1 - sb)(du.n) n] / (alpha dt CPORNH)
    // sb = 1 inside the body (full rigid-body velocity), 0 outside (only the normal velocity is imposed,
    // the tangential velocity is kept). Unlike the component weighting |n_i| of X 15 1, the projection does
    // not depend on the orientation of the hull to the grid (rolled hull, bilges): the tangential velocity
    // outside is not damped. Unlike the sharp switch of X 15 2 at FB = 0, the forcing stays continuous when
    // the hull moves across cell centres. On grid-aligned faces X 15 1 and 4 are identical.
    // The forcing is divided by CPORNH as in nhflow_forcing::solid_forcing, so that the body velocity is
    // fully imposed also inside a porous structure (without porous media CPORNH = 1).
    if(p->X15==4 || (p->X15==11 && wavemaker))
    LOOP
    {
        H = Hsolidface_nhflow(p,d,0,0,0);
        
        d->FHB[IJK] = MIN(d->FHB[IJK] + H, 1.0); 
        
        if(H<=0.0)
        continue;
        
        target_velocity(p,wavemaker,uf,vf,wf);
        
        du = uf - U[IJK];
        dv = vf - V[IJK];
        dw = wf - W[IJK];
        
        const double dsh = 0.25*Hpsi_nhflow(p,d);
        const double tb = MIN(MAX(0.5*(1.0 - d->FB[IJK]/dsh), 0.0), 1.0);
        double sb = tb*tb*(3.0 - 2.0*tb);
        
        nx=ny=nz=0.0;
        udotn=0.0;
        
        if(sb<1.0)
        {
            // Normal vectors calculation 
            if(p->geo_curv==1)
            fbnormal_curv(p,d,nx,ny);
            else
            {
            nx = -(d->FB[Ip1JK] - d->FB[Im1JK])/(p->DXP[IP] + p->DXP[IM1])
                - 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
                
            if(p->j_dir==0)
            ny = 0.0;
            
            if(p->j_dir==1)
            ny = -(d->FB[IJp1K] - d->FB[IJm1K])/(p->DYP[JP] + p->DYP[JM1])
                - 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
            }
                
            nz = -(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP]*WL(i,j) + p->DZP[KM1]*WL(i,j));

            norm = sqrt(nx*nx + ny*ny + nz*nz);
            
            if(norm > 1.0e-10)
            {
            nx /= norm;
            ny /= norm;
            nz /= norm;
            
            udotn = nx*du + ny*dv + nz*dw;
            }
            
            // degenerate gradient (flat level set): fall back to full forcing
            else
            sb = 1.0;
        }
        
        const double fac = H/(alpha[iter]*p->dt*CPORNH);
        
        FX[IJK] += fac*(sb*du + (1.0-sb)*udotn*nx);
        FY[IJK] += fac*(sb*dv + (1.0-sb)*udotn*ny);
        FZ[IJK] += fac*(sb*dw + (1.0-sb)*udotn*nz);
    }
   
    pgc->start5V(p,d->FHB,50);
}
    
void sixdof_obj_nhflow::fbnormal_curv(lexer *p, fdm_nhf *d, double &gx, double &gy)
{
    // curvilinear grid (G 50): the horizontal part of -grad(FB) at the cell (i,j,k), the gradient along the sigma
    // level by the chain rule on the cell centres, then the sigma correction (sigx, sigy are physical)
    double fx,fy;
    p->geo_grad(i,j,d->FB[Ip1JK],d->FB[Im1JK],d->FB[IJp1K],d->FB[IJm1K],fx,fy);
    const double fs = (d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);

    gx = -fx - 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*fs;
    gy = -fy - 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*fs;
}

double sixdof_obj_nhflow::Hpsi_nhflow(lexer *p, fdm_nhf *d)
{
    // half width of the smoothed Heaviside of the direct forcing (as in Hsolidface_nhflow)
    if(p->geo_curv==1)
    return p->A526*p->geo_width(i,j);

    if(p->j_dir==0)
    return p->A526*p->DXN[IP];
    
    return p->A526*0.5*(p->DXN[IP] + p->DYN[JP]);
}

double sixdof_obj_nhflow::Hsolidface_nhflow(lexer *p, fdm_nhf *d, int aa, int bb, int cc)
{
    double psi, H, phival_fb;
    
    if(p->j_dir==0)
    psi = p->A526*(1.0/1.0)*(p->DXN[IP] + 0.0*p->DZN[KP]*p->WL[IJ]);
    
    if(p->j_dir==1)
    psi = p->A526*(1.0/2.0)*(p->DXN[IP] + p->DYN[JP] + 0.0*p->DZN[KP]*p->WL[IJ]);

    // curvilinear grid (G 50): the physical mean width of the cell
    if(p->geo_curv==1)
    psi = p->A526*p->geo_width(i,j);


    // Construct solid heaviside function
    phival_fb = d->FB[IJK];
    
    if(-phival_fb > psi)
    H = 1.0;

    if(-phival_fb < -psi)
    H = 0.0;

    if(fabs(phival_fb)<=psi)
    H = 0.5*(1.0 + (-phival_fb)/psi + (1.0/PI)*sin((PI*(-phival_fb))/psi));
    
    return H;
}

// X 15 11: direct forcing with the body in the pressure projection
//
//   fraction    kap of every cell: the part of the cell inside the body, from the signed distance FB of the
//               cell centre and its gradient, kap = 1/2 - FB/|(nx dx, ny dy, nz dz)| limited to 0..1 (exact
//               for the faces of the hull along the grid). It changes linearly while the surface moves through
//               the cell: no switching of cells, the volume of the body on the grid is that of the hull.
//   forcing     u = u* + tau (u_b - u*), tau = max(2 kap - 1, 0): only where the cell centre is inside the
//               hull, the rigid-body velocity in the cells of the body; continuous in the position of the hull.
//               The cells with their centre in the water are not forced: the body acts on them through the
//               fluxes of its faces and the pressure. A forcing of their normal velocity with the fraction
//               works against the projection, which has made the water flow around the hull: the pressure then
//               holds a part that only balances the forcing, and the pressure loads a mean force (a box free in
//               waves of 2 s drifted downwave with 4 times the drift force of its reflected and transmitted
//               waves, the mean force on the fixed box of 2 m in waves of 1.4 s was 63 % of that of its waves,
//               now 93 %). The reaction of the forcing as a load on the body did not help: the heave decay and
//               the mean force on the fixed box got worse.
//               Open: a box free in waves of 2 s (H 0.04 m) still gains about 0.002 m/s2 downwave, 3 times the
//               drift force of its waves (NHFLOW Floating validation 04).
//   projection  open fractions of the links of the Poisson equation: 1 - kap of the cell for the vertical
//               link through it (W of the cell), 1 - the fraction of the face from the mean FB of its two
//               cells for the horizontal links (d->MBX, MBY, MBZ; nhflow_membrane_beta.h, link mode). The
//               pressure correction then does not move water through the body, the pressure below the hull
//               follows from the normal velocity of the hull (its added mass is in the pressure), and the
//               Dirichlet condition at the top of a column covered by the body does not reach the water below.
//               The nodes inside the body hold no pressure (nhflow_membrane_row).
//   fluxes      the faces of the body carry the body through their closed part, in the divergence of the
//               projection, the continuity flux and the momentum fluxes (nhflow_bodyflux.h, nhflow_HLL.cpp)
//   level       the water level of the columns covered by the body is that of the surrounding free surface
//               (nhflow_fsf_body, as X 13 1)
//   The velocity is not forced again after the projection (the projected field keeps its continuity).
//   Loads: 6DOF_obj_nhflow_forces_stl.cpp, coupling: 6DOF_obj_nhflow_addedmass.cpp
//
//   Not yet: mesh refinement (G 1), curvilinear grids (G 50), the semi-implicit time step (A 510 11-13), several
//   bodies; checkpoints restart with the added mass and the acceleration of the first step after the restart
double sixdof_obj_nhflow::cellfraction_nhflow(lexer *p, fdm_nhf *d, slice &WL, double &ex, double &ey, double &ez)
{
    double gx,gy,gz;
    
    if(p->geo_curv==1)
    fbnormal_curv(p,d,gx,gy);
    else
    {
    gx = -(d->FB[Ip1JK] - d->FB[Im1JK])/(p->DXP[IP] + p->DXP[IM1])
        - 0.5*(p->sigx[FIJK]+p->sigx[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
    
    gy = 0.0;
    
    if(p->j_dir==1)
    gy = -(d->FB[IJp1K] - d->FB[IJm1K])/(p->DYP[JP] + p->DYP[JM1])
        - 0.5*(p->sigy[FIJK]+p->sigy[FIJKp1])*(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP] + p->DZP[KM1]);
    }
    
    const double wl = MAX(WL(i,j),1.0e-6);
    
    gz = -(d->FB[IJKp1] - d->FB[IJKm1])/(p->DZP[KP]*wl + p->DZP[KM1]*wl);
    
    const double g = sqrt(gx*gx + gy*gy + gz*gz);
    
    const double dx = p->geo_curv==1 ? p->geo_width(i,j) : p->DXN[IP];
    const double dy = p->geo_curv==1 ? p->geo_width(i,j) : p->DYN[JP];
    const double dz = p->DZN[KP]*wl;
    
    // away from the surface (FB at the limit of its band) or a degenerate gradient: inside or outside
    if(g<0.1)
    {
    ex = ey = ez = 0.0;
    return d->FB[IJK]<0.0 ? 1.0 : 0.0;
    }
    
    ex = gx/g;
    ey = gy/g;
    ez = gz/g;
    
    // width of the cell in the direction of the normal. With the sum of the sides, |nx| dx + |ny| dy + |nz| dz (a
    // plane through a box cell), the cells at an edge of the hull, where the normal of the distance points along the
    // diagonal, stay open by a part although they lie inside: nodes with nearly closed links, whose pressure is
    // not determined by the water. With the length of (nx dx, ny dy, nz dz) the fraction is exact for the faces
    // along the grid and 1 in the cells at the edges
    const double w = sqrt(ex*ex*dx*dx + (p->j_dir==1 ? ey*ey*dy*dy : 0.0) + ez*ez*dz*dz);
    
    const double kap = MIN(MAX(0.5 - d->FB[IJK]/w, 0.0), 1.0);
    
    // no cells and links that are open or closed by less than 1 %: their pressure nodes are nearly without water
    return kap>0.99 ? 1.0 : (kap<0.01 ? 0.0 : kap);
}

void sixdof_obj_nhflow::fraction_ini(lexer *p, fdm_nhf *d)
{
    const int n = p->imax*p->jmax*(p->kmax+2);
    
    double **mb[] = {&d->MBETA,&d->MCHI,&d->MRCX,&d->MRCY,&d->MBX,&d->MBY,&d->MBZ};
    
    if(FNX==nullptr)
    {
    p->Darray(FNX,n);
    p->Darray(FNY,n);
    p->Darray(FNZ,n);
    }
    
    if(d->FBUX==nullptr)
    {
    p->Darray(d->FBUX,n);
    p->Darray(d->FBUY,n);
    p->Darray(d->FBWZ,n);
    p->Iarray(d->FBNEAR,p->imax*p->jmax);
    }
    
    // open fractions of the links (the mobilities of nhflow_membrane_beta.h, link mode) and the face corrections
    // of the last projection, 1 and 0 away from the body
    for(int q=0; q<7; ++q)
    if(*mb[q]==nullptr)
    {
        p->Darray(*mb[q],n);
        
        for(int m=0; m<n; ++m)
        (*mb[q])[m] = (q==1 || q==2 || q==3) ? 0.0 : 1.0;
    }
    
    d->MPROJ = 1;
}

void sixdof_obj_nhflow::forcing_kernel_fraction(lexer *p, fdm_nhf *d, ghostcell *pgc, 
                             double *U, double *V, double *W, double *FX, double *FY, double *FZ, slice &WL, int iter, bool wavemaker)
{
    // smallest open fraction of a link: the Poisson equation stays regular inside the body
    const double bmin = 1.0e-4;
    
    double uf,vf,wf,ex,ey,ez;
    
    fraction_ini(p,d);
    
    // acceleration of the time step (the loads of prescribed degrees of freedom, force_calc_stl)
    if(iter==0 && !fraction_reforce && p->simtime!=tfd_)
    {
        if(tfd_>=0.0 && p->simtime>tfd_)
        afd_ = (u_fb - ufd_)/(p->simtime - tfd_);
        
        ufd_ = u_fb;
        tfd_ = p->simtime;
    }
    
    LOOP
    {
        ex = ey = ez = 0.0;
        
        const double kap = (p->wet[IJ]==1) ? cellfraction_nhflow(p,d,WL,ex,ey,ez) : 0.0;
        
        FNX[IJK] = ex;
        FNY[IJK] = ey;
        FNZ[IJK] = ez;
        
        d->FHB[IJK] = MIN(d->FHB[IJK] + kap, 1.0);
        d->MBETA[IJK] = MAX(1.0 - kap, bmin);
        d->MBZ[IJK] = d->MBETA[IJK];
        
        // forcing only where the cell centre is inside the hull
        if(kap<=0.5 || fraction_reforce)
        continue;
        
        target_velocity(p,wavemaker,uf,vf,wf);
        
        const double du = uf - U[IJK];
        const double dv = vf - V[IJK];
        const double dw = wf - W[IJK];
        
        const double tau = 2.0*kap - 1.0;
        const double fac = 1.0/(alpha[iter]*p->dt*CPORNH);
        
        FX[IJK] += fac*tau*du;
        FY[IJK] += fac*tau*dv;
        FZ[IJK] += fac*tau*dw;
    }
    
    pgc->start5V(p,d->FHB,50);
    pgc->start4V(p,d->MBETA,1);
    pgc->start4V(p,d->MBZ,1);
    pgc->start4V(p,FNX,1);
    pgc->start4V(p,FNY,1);
    pgc->start4V(p,FNZ,1);
    
    // open fractions of the horizontal links
    LOOP
    {
        d->MBX[IJK] = 1.0;
        d->MBY[IJK] = 1.0;
        d->FBUX[IJK] = 0.0;
        d->FBUY[IJK] = 0.0;
        d->FBWZ[IJK] = 0.0;
        
        if(p->wet[IJ]==0)
        continue;
        
        const double dz = p->DZN[KP]*MAX(WL(i,j),1.0e-6);
        
        // normal velocity of the body at the centres of the x and y faces, vertical velocity at the cell centre
        {
        const double ry = p->pos_y() - c_(1);
        const double rx = p->pos_x() - c_(0);
        
        d->FBWZ[IJK] = u_fb(2) + u_fb(3)*ry - u_fb(4)*rx;
        
        if(p->wet[Ip1J]==1)
        d->FBUX[IJK] = u_fb(0) + u_fb(4)*(0.5*(p->ZSP[IJK] + p->ZSP[Ip1JK]) - c_(2)) - u_fb(5)*ry;
        
        if(p->j_dir==1 && p->wet[IJp1]==1)
        d->FBUY[IJK] = u_fb(1) + u_fb(5)*rx - u_fb(3)*(0.5*(p->ZSP[IJK] + p->ZSP[IJp1K]) - c_(2));
        }
        
        if(p->wet[Ip1J]==1 && p->flag4[Ip1JK]>0)
        d->MBX[IJK] = linkfraction_nhflow(d->FB[IJK],d->FB[Ip1JK],d->MBETA[IJK],d->MBETA[Ip1JK],
                                          (p->j_dir==1 ? 0.5*fabs(FNY[IJK]+FNY[Ip1JK])*p->DYN[JP] : 0.0) + 0.5*fabs(FNZ[IJK]+FNZ[Ip1JK])*dz, bmin);
        
        if(p->j_dir==1 && p->wet[IJp1]==1 && p->flag4[IJp1K]>0)
        d->MBY[IJK] = linkfraction_nhflow(d->FB[IJK],d->FB[IJp1K],d->MBETA[IJK],d->MBETA[IJp1K],
                                          0.5*fabs(FNX[IJK]+FNX[IJp1K])*p->DXN[IP] + 0.5*fabs(FNZ[IJK]+FNZ[IJp1K])*dz, bmin);
    }
    
    // columns near the body
    for(int q=0; q<p->imax*p->jmax; ++q)
    d->FBNEAR[q] = 0;
    
    LOOP
    if(d->MBETA[IJK]<1.0)
    d->FBNEAR[IJ] = 1;
    
    pgc->gcsl_start4Vint(p,d->FBNEAR,1);
    
    for(int qn=0; qn<nhflow_bodyflux::nnear; ++qn)
    {
        SLICELOOP4
        if(d->FBNEAR[IJ]==0 && (d->FBNEAR[Ip1J]==1 || d->FBNEAR[Im1J]==1 
           || (p->j_dir==1 && (d->FBNEAR[IJp1]==1 || d->FBNEAR[IJm1]==1))))
        d->FBNEAR[IJ] = 2;
        
        SLICELOOP4
        if(d->FBNEAR[IJ]==2)
        d->FBNEAR[IJ] = 1;
        
        pgc->gcsl_start4Vint(p,d->FBNEAR,1);
    }
    
    // no pressure at the nodes inside the body (between two cells of the body)
    LOOP
    if(k>0 && d->MBZ[IJK]<=bmin && d->MBZ[IJKm1]<=bmin)
    d->P[FIJK] = 0.0;
    
    pgc->start4V(p,d->MBX,1);
    pgc->start4V(p,d->MBY,1);
    pgc->start4V(p,d->FBUX,1);
    pgc->start4V(p,d->FBUY,1);
    pgc->start4V(p,d->FBWZ,1);
}

double sixdof_obj_nhflow::linkfraction_nhflow(double fa, double fb, double ba, double bb, double w, double bmin)
{
    // the face between two cells: FB at its centre from the two cell centres, w the change of FB over the face
    // (the components of the normal in the plane of the face times its sides). The open fraction of the face is the
    // part with FB > 0. A face parallel to the surface (w = 0) would switch when the surface passes it: the face
    // is open by at most twice the open fraction of its cells, so that it closes linearly while the surface
    // moves through the last half cell in front of it.
    if(ba>=1.0 && bb>=1.0)
    return 1.0;
    
    const double fm = 0.5*(fa + fb);
    const double af = w>1.0e-12 ? MIN(MAX(0.5 + fm/w, 0.0), 1.0) : (fm>0.0 ? 1.0 : 0.0);
    
    const double b = MIN(af, MIN(2.0*MIN(ba,bb), 1.0));
    
    return b>0.99 ? 1.0 : (b<0.01 ? bmin : b);
}
