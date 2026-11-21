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
#include"gradient.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"nhflow_fsf_hull.h"
#include"nhflow_hullbody.h"
#include"nhflow_thinbody.h"

void sixdof_obj_nhflow::hydrodynamic_forces_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL, bool finalize)
{
	// forcecalc
    // X 13 2: the hull is the top of the water columns, the loads are the column integrals
    if(d->hull!=nullptr)
    {
    d->hull->forces(p,d,pgc,WL,this,finalize);
    return;
    }
    
    // X 15 8: the hull is a sharp thin body, the loads are the outer pressure on the blocked links
    // mesh refinement: every grid holds the hull, every blocked face counts on the finest grid at its centre
    // (and with placed patches, G 40 1, on the rank of that grid)
    if(d->hullbody!=nullptr && amr_grids_nhflow)
    {
    double F[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    double S[3] = {0.0,0.0,0.0};
    
    amr_grids_nhflow([&](const nhflow_grid &g)
    {
        if(g.d->hullbody==nullptr)
        return;
        
        const std::function<bool(double,double)> keep = [&](double x, double y)
        {
            if(amr_owner_nhflow && !amr_owner_nhflow(x,y))
            return false;
            return amr_grid_nhflow(x,y).d==g.d;
        };
        
        g.d->hullbody->loads(g.p,g.d,*g.WL,&keep,F,S);
        S[1] += double(g.d->hullbody->thinbody()->nlower_local(g.p,keep));
        S[2] += g.d->hullbody->thinbody()->inner_volume_local(g.p,g.d,keep);
    });
    
    d->hullbody->finish_loads(p,d,pgc,F,finalize,S);
    return;
    }
    
    if(d->hullbody!=nullptr)
    {
    d->hullbody->forces(p,d,pgc,WL,finalize);
    return;
    }
    
    // porous body (X 16) always uses the STL integration: the (1-n) scaled envelope pressure
    // plus the volume drag from update_forcing_nhflow_porous()
    if(p->X60==1 || p->X16==1)
    force_calc_stl(p,d,pgc,WL,finalize);
    
    
    if(p->X60==2 && p->X16==0)
    {
    triangulation(p,d,pgc);
	reconstruct(p,d);
    force_calc_lsm(p,d,pgc,WL);
        
    deallocate(p,d,pgc);
    }
} 


void sixdof_obj_nhflow::force_calc_stl(lexer* p, fdm_nhf *d, ghostcell *pgc, slice &WL, bool finalize)
{
    // Pressure + viscous load on the wetted part of the STL surface.
    //
    // Each triangle is clipped against the local free surface z = wd + eta (Sutherland-Hodgman,
    // so vertices lying exactly on the surface are handled), the wetted polygon is fan-triangulated
    // and every sub-triangle is integrated with the 3-point edge-midpoint rule. That rule is exact
    // for quadratic integrands, so for a planar facet both the force (linear hydrostatic pressure)
    // and the moment (lever arm x pressure) are integrated exactly, independent of STL resolution.
    // Moment arms are the quadrature points, not the centroid of the unclipped triangle.
    
    double Fv_x,Fv_y,Fv_z;
    double Xe_p,Ye_p,Ze_p,Xe_v,Ye_v,Ze_v;
    double A;
    double vx[3],vy[3],vz[3];
    double px[4],py[4],pz[4];   // wetted polygon: a triangle clipped by a plane has at most 4 vertices
    
    A=0.0;
    Xe=Ye=Ze=Ke=Me=Ne=0.0;
    
    // distributed loads (fluid_member): the hull forces binned along the member, this rank's part
    external_strips_zero();
    
    // porous body (X 16): the pressure acts on the solid fraction only.
    // -int_V (1-n) grad p dV = -(1-n_fb) oint p n dA  (divergence theorem, n_fb const inside)
    // -> buoyancy and Froude-Krylov load on the skeleton
    const double pfac = (p->X16==1) ? (1.0 - p->X16_n) : 1.0;
    Xe_p=Ye_p=Ze_p=Xe_v=Ye_v=Ze_v=0.0;
    
    // X 15 11: the pressure sampled in the water, one and two cells off the hull (X 42 cells, if given), and brought
    // to the hull with a parabola through the two samples with the normal pressure gradient of the wall,
    // dp/dn = -rho a_b.n:  p_wall = (4 p1 - p2)/3 + (2/3) rho doff a_b.n. The last term is the load -G a_b, G the
    // mass of 2/3 of the water layer between hull and first sample points in the normal direction (kept in Gdf_).
    // The cells cut by the hull and the cells of the body hold no pressure of the water. With one sample and the
    // gradient of the wall the heave added mass of a box grew with the sample distance (the acceleration of the
    // water falls off from the hull): 1.10, 1.13, 1.17, 1.22 of the displaced mass for 0.5, 1, 1.5, 2 cells.
    // load_P: the loads of the pressure field load_P alone (added mass, 6DOF_obj_nhflow_addedmass.cpp)
    const bool frac = (p->X15==11);
    double Gl[6][6];
    
    for(int qa=0; qa<6; ++qa)
    for(int qb=0; qb<6; ++qb)
    Gl[qa][qb] = 0.0;
    
    // Set new time
    // the loads of the last stage belong to the end of the step (simtime is advanced after it)
    curr_time = finalize ? p->simtime + p->dt : p->simtime;
    
    for(int n=0; n<tricount; ++n)
    {     
        for(int q=0; q<3; ++q)
        {
            vx[q] = tri_x[n][q];
            vy[q] = tri_y[n][q];
            vz[q] = tri_z[n][q];
        }
        
		// Center of triangle
		const double xc = (vx[0] + vx[1] + vx[2])/3.0;
		const double yc = (vy[0] + vy[1] + vy[2])/3.0;
        
        // Ownership by the triangle centroid.
        // NHFLOW is decomposed in x and y only, and originz/endz are the bounds of the (sigma)
        // mesh, not physical heights -> no z test. In 2D the STL may be wider than the one-cell
        // strip -> no y test either (same convention as ray_cast_z and band_distance).
        if(amr_owner_nhflow)
        {
            if(!amr_owner_nhflow(xc,yc))
            continue;
        }
        else if(p->geo_curv==1)
        {
            // curvilinear grid (G 50): the rank of the cell that holds the centroid (grid_locate.cpp)
            int ci,cj;
            if(!p->geo_cell(xc,yc,ci,cj))
            continue;
        }
        else
        {
        if(!(xc >= p->originx && xc < p->endx))
        continue;
        
        if(p->j_dir==1 && !(yc >= p->originy && yc < p->endy))
        continue;
        }
        
        // mesh refinement (nhflow_amr): the triangle is sampled on the finest grid at its centroid
        lexer *pl = p;
        fdm_nhf *dl = d;
        slice *WLl = &WL;
        double dsm = DSM;
        if(amr_grid_nhflow)
        {
            nhflow_grid g = amr_grid_nhflow(xc,yc);
            pl = g.p;
            dl = g.d;
            WLl = g.WL;
            dsm = (g.hr>0.0) ? DSM*g.hr : DSM*pl->DXM/p->DXM;
        }
        lexer *plc = pl;
        fdm_nhf *dlc = dl;
        slice *WLc = WLl;
        
        // Normal vector (pointing outwards)
        double nx = (vy[1] - vy[0])*(vz[2] - vz[0]) - (vy[2] - vy[0])*(vz[1] - vz[0]);
        double ny = (vx[2] - vx[0])*(vz[1] - vz[0]) - (vx[1] - vx[0])*(vz[2] - vz[0]); 
        double nz = (vx[1] - vx[0])*(vy[2] - vy[0]) - (vx[2] - vx[0])*(vy[1] - vy[0]);
        const double norm = sqrt(nx*nx + ny*ny + nz*nz);
        
        if(norm<1.0e-20)
        continue;
        
        nx /= norm;
        ny /= norm;
        nz /= norm;
        
        // 2D: the y-normal end caps carry no load
        if(p->j_dir==0)
        {
            if(fabs(ny)>0.9)
            continue;
            
            ny = 0.0;
        }
        
        // X 15 11: the water level at the hull from the water side, as the pressure: at one and two cells off the hull
        // along the horizontal part of the normal, (4 eta1 - eta2)/3. At the hull itself the interpolation takes in
        // the level of the columns covered by the body, the surface extended from all around the body, which at a
        // side of the hull in waves is not the level of the water at that side. Below the keel (vertical normal)
        // it is the level of the column.
        auto eta_hull = [&](double x, double y, double z)
        {
            if(!frac)
            return pl->ccslipol4(dl->eta,x,y);
            
            const double off = sample_dist_nhflow(pl,*WLl,x,y,z,nx,ny,nz);
            
            return (4.0*pl->ccslipol4(dl->eta,x + nx*off,y + ny*off) - pl->ccslipol4(dl->eta,x + 2.0*nx*off,y + 2.0*ny*off))/3.0;
        };
        
        const double fsf_z = p->wd + eta_hull(xc,yc,(vz[0] + vz[1] + vz[2])/3.0);
        
        // Clip the triangle to the wetted side z <= fsf_z
        int np=0;
        for(int q=0; q<3; ++q)
        {
            const int r = (q+1)%3;
            const bool in_q = (vz[q] <= fsf_z);
            const bool in_r = (vz[r] <= fsf_z);
            
            if(in_q)
            {
                px[np] = vx[q];
                py[np] = vy[q];
                pz[np] = vz[q];
                ++np;
            }
            
            if(in_q != in_r)   // edge crosses the surface; vz[r]!=vz[q] is guaranteed here
            {
                const double t = (fsf_z - vz[q])/(vz[r] - vz[q]);
                px[np] = vx[q] + t*(vx[r] - vx[q]);
                py[np] = vy[q] + t*(vy[r] - vy[q]);
                pz[np] = vz[q] + t*(vz[r] - vz[q]);
                ++np;
            }
        }
        
        if(np<3)
        continue;
        
        // Fan triangulation of the wetted polygon
        for(int s=1; s<np-1; ++s)
        {
            const double ax = px[0],   ay = py[0],   az = pz[0];
            const double bx = px[s],   by = py[s],   bz = pz[s];
            const double cx = px[s+1], cy = py[s+1], cz = pz[s+1];
            
            const double crx = (by - ay)*(cz - az) - (cy - ay)*(bz - az);
            const double cry = (cx - ax)*(bz - az) - (bx - ax)*(cz - az);
            const double crz = (bx - ax)*(cy - ay) - (cx - ax)*(by - ay);
            const double A_sub = 0.5*sqrt(crx*crx + cry*cry + crz*crz);
            
            if(A_sub<1.0e-30)
            continue;
            
            // Pressure: 3-point edge-midpoint rule
            const double mx[3] = {0.5*(ax + bx), 0.5*(bx + cx), 0.5*(cx + ax)};
            const double my[3] = {0.5*(ay + by), 0.5*(by + cy), 0.5*(cy + ay)};
            const double mz[3] = {0.5*(az + bz), 0.5*(bz + cz), 0.5*(cz + az)};
            double pq[3];       // pfac p at the quadrature points, for the strips
            
            // multi-block grid: every quadrature point on the grid of this rank that holds it (amr_point_nhflow)
            if(amr_point_grid)
            {
            for(int q=0; q<3; ++q)
            {
                const double sx = mx[q] + p->X42*nx*dsm, sy = my[q] + p->X42*ny*dsm;
                amr_point_nhflow(p,sx,sy,plc,dlc,WLc,pl,dl,WLl);
                
                const double pval = pl->ccipol7V(dl->P, *WLl, dl->bed, sx, sy, mz[q] + p->X42*nz*dsm);
                const double hsp = MAX(0.0, (p->wd + pl->ccslipol4(dl->eta,mx[q],my[q]) - mz[q])*p->W1*fabs(p->W22));
                
                const double w  = A_sub/3.0;
                pq[q] = pfac*(pval + hsp);
                const double fx = -pfac*(pval + hsp)*w*nx;
                const double fy = -pfac*(pval + hsp)*w*ny;
                const double fz = -pfac*(pval + hsp)*w*nz;
                
                const double rx = mx[q] - c_(0);
                const double ry = my[q] - c_(1);
                const double rz = mz[q] - c_(2);
                
                Xe += fx;
                Ye += fy;
                Ze += fz;
                Ke += ry*fz - rz*fy;
                Me += rz*fx - rx*fz;
                Ne += rx*fy - ry*fx;
                
                Xe_p += fx;
                Ye_p += fy;
                Ze_p += fz;
            }
            pl = plc;
            dl = dlc;
            WLl = WLc;
            }
            
            else
            for(int q=0; q<3; ++q)
            {
                // non-hydrostatic pressure, optionally sampled X42 mean cell sizes off the wall
                const double doff = frac ? sample_dist_nhflow(pl,*WLl,mx[q],my[q],mz[q],nx,ny,nz) : p->X42*dsm;
                
                double pval = pl->ccipol7V(load_P!=nullptr ? load_P : dl->P, *WLl, dl->bed, mx[q] + nx*doff, 
                                                                    my[q] + ny*doff, 
                                                                    mz[q] + nz*doff);
                
                // X 15 11: second sample at twice the distance, p_wall = (4 p1 - p2)/3 - (2/3) doff dp/dn
                if(frac)
                pval = (4.0*pval - pl->ccipol7V(load_P!=nullptr ? load_P : dl->P, *WLl, dl->bed, mx[q] + 2.0*nx*doff, 
                                                                    my[q] + 2.0*ny*doff, 
                                                                    mz[q] + 2.0*nz*doff))/3.0;
                // hydrostatic pressure at the quadrature point
                const double hsp = load_P!=nullptr ? 0.0 : MAX(0.0, (p->wd + eta_hull(mx[q],my[q],mz[q]) - mz[q])*p->W1*fabs(p->W22));
                
                if(frac)
                {
                    const double rgx = mx[q] - c_(0), rgy = my[q] - c_(1), rgz = mz[q] - c_(2);
                    const double N[6] = {nx, ny, nz, rgy*nz - rgz*ny, rgz*nx - rgx*nz, rgx*ny - rgy*nx};
                    
                    for(int qa=0; qa<6; ++qa)
                    for(int qb=0; qb<6; ++qb)
                    Gl[qa][qb] += p->W1*(2.0/3.0)*doff*(A_sub/3.0)*N[qa]*N[qb];
                }
                
                const double w  = A_sub/3.0;
                pq[q] = pfac*(pval + hsp);
                const double fx = -pfac*(pval + hsp)*w*nx;
                const double fy = -pfac*(pval + hsp)*w*ny;
                const double fz = -pfac*(pval + hsp)*w*nz;
                
                const double rx = mx[q] - c_(0);
                const double ry = my[q] - c_(1);
                const double rz = mz[q] - c_(2);
                
                Xe += fx;
                Ye += fy;
                Ze += fz;
                Ke += ry*fz - rz*fy;
                Me += rz*fx - rx*fz;
                Ne += rx*fy - ry*fx;
                
                Xe_p += fx;
                Ye_p += fy;
                Ze_p += fz;
            }
            
            if(external_strips_on())
            {
                const double ta[3] = {ax,ay,az}, tb[3] = {bx,by,bz}, tc[3] = {cx,cy,cz};
                external_strips_tri(ta,tb,tc,pq,nx,ny,nz);
            }
            
            A += A_sub;
            
            // Viscous forces at the sub-triangle centroid
            // (not for a porous body: the envelope is not a wall, the volume drag replaces it)
            if(p->X16==1 || load_P!=nullptr)
            continue;
            
            const double gx = (ax + bx + cx)/3.0;
            const double gy = (ay + by + cy)/3.0;
            const double gz = (az + bz + cz)/3.0;
            
            hydrodynamic_viscous_forces_nhflow(pl, dl, pgc, *WLl, Fv_x, Fv_y, Fv_z, A_sub, gx, gy, gz, nx, ny, nz);
            
            Xe += Fv_x;
            Ye += Fv_y;
            Ze += Fv_z;
            external_strips_add(gx,gy,gz,Fv_x,Fv_y,Fv_z);
            Ke += (gy - c_(1))*Fv_z - (gz - c_(2))*Fv_y;
            Me += (gz - c_(2))*Fv_x - (gx - c_(0))*Fv_z;
            Ne += (gx - c_(0))*Fv_y - (gy - c_(1))*Fv_x;
            
            Xe_v += Fv_x;
            Ye_v += Fv_y;
            Ze_v += Fv_z;
        }
	}
    
	// Communication with other processors
    A = pgc->globalsum(A);
    
	Xe = pgc->globalsum(Xe);
	Ye = pgc->globalsum(Ye);
	Ze = pgc->globalsum(Ze);
	Ke = pgc->globalsum(Ke);
	Me = pgc->globalsum(Me);
	Ne = pgc->globalsum(Ne);
    
	Xe_p = pgc->globalsum(Xe_p);
	Ye_p = pgc->globalsum(Ye_p);
	Ze_p = pgc->globalsum(Ze_p);
	Xe_v = pgc->globalsum(Xe_v);
	Ye_v = pgc->globalsum(Ye_v);
	Ze_v = pgc->globalsum(Ze_v);
    
    // X 15 11: the layer of water between the hull and the sample points
    if(frac)
    {
        for(int qa=0; qa<6; ++qa)
        for(int qb=0; qb<6; ++qb)
        Gdf_(qa,qb) = pgc->globalsum(Gl[qa][qb]);
        
        // the loads of a pressure field alone: no layer, no weight, no output
        if(load_P!=nullptr)
        {
            if(p->j_dir==0)
            Ye = Ke = Ne = 0.0;
            
            return;
        }
        
        // the acceleration of the body in the projection of this stage (the one of its last RK stage)
        // (prescribed degrees of freedom: from the velocities of the last two time steps)
        Eigen::Matrix<double,6,1> ab = am_prev_;
        
        for(int qa=0; qa<6; ++qa)
        if(rb.fixed(qa))
        ab(qa) = afd_(qa);
        
        const Eigen::Matrix<double,6,1> Fg = -Gdf_*ab;
        
        Xe += Fg(0);
        Ye += Fg(1);
        Ze += Fg(2);
        Ke += Fg(3);
        Me += Fg(4);
        Ne += Fg(5);
        
        Xe_p += Fg(0);
        Ye_p += Fg(1);
        Ze_p += Fg(2);
    }
    
    // porous body: add the Darcy-Forchheimer drag reaction (already global sums).
    // It is written to the viscous columns of the force file.
    if(p->X16==1)
    {
        Xe += Xd;
        Ye += Yd;
        Ze += Zd;
        Ke += Kd;
        Me += Md;
        Ne += Nd;
        
        Xe_v = Xd;
        Ye_v = Yd;
        Ze_v = Zd;
    }
    
    // 2D: only surge, heave and pitch exist. Out-of-plane moments from an asymmetric STL
    // triangulation would otherwise enter h_ and tilt the trimesh out of the x-z plane,
    // although u_fb(3) and u_fb(5) are zeroed in update_fbvel.
    if(p->j_dir==0)
    {
        Ye = Ke = Ne = 0.0;
        Ye_p = Ye_v = 0.0;
    }
    
	// Add gravity force
	Xe += p->W20*Mass_fb;
	Ye += p->W21*Mass_fb;
	Ze += p->W22*Mass_fb;
    
    
    // Print results	
    if (p->mpirank==0 && finalize==1) 
    {
        printforce<<curr_time<<" \t "<<Xe<<" \t "<<Ye<<" \t "<<Ze<<" \t "<<Ke
        <<" \t "<<Me<<" \t "<<Ne<<" \t "<<Xe_p<<" \t "<<Ye_p<<" \t "<<Ze_p<<" \t "<<Xe_v<<" \t "<<Ye_v<<" \t "<<Ze_v<<endl;   
    }
}

// multi-block grid (amr_point_grid): the grid of this rank that holds the quadrature point (x,y): a patch or block piece,
// level 0 (the native grid with its ring) where it holds the point, else the grid of the centroid (plc, dlc, WLc)
void sixdof_obj_nhflow::amr_point_nhflow(lexer *p, double x, double y, lexer *plc, fdm_nhf *dlc, slice *WLc,
                                         lexer *&pl, fdm_nhf *&dl, slice *&WLl)
{
    nhflow_grid g = amr_grid_nhflow(x,y);
    int a,b;
    
    if(g.p!=p || plc==p || p->geo_cell(x,y,a,b))
    {
        pl = g.p;
        dl = g.d;
        WLl = g.WL;
        return;
    }
    
    pl = plc;
    dl = dlc;
    WLl = WLc;
}

// X 15 11: distance of the pressure sample point from the hull at (x,y,z) with the outward normal n: X 42 cells
// (1 if not given) in the direction of the normal, the cell of the point with its height in the water column
double sixdof_obj_nhflow::sample_dist_nhflow(lexer *pl, slice &WLl, double x, double y, double z, double nx, double ny, double nz)
{
    const int marge = increment::marge;
    const double c = pl->X42>0.0 ? pl->X42 : 1.0;
    
    if(pl->geo_curv==1)
    return c*DSM;
    
    const int ii = MAX(0, MIN(pl->knox-1, pl->posc_i(x)));
    const int jj = (pl->j_dir==1) ? MAX(0, MIN(pl->knoy-1, pl->posc_j(y))) : 0;
    const int kk = MAX(0, MIN(pl->knoz-1, pl->posc_sig(ii,jj,z)));
    
    return c*(fabs(nx)*pl->DXN[ii+marge] + (pl->j_dir==1 ? fabs(ny)*pl->DYN[jj+marge] : 0.0) + fabs(nz)*pl->DZN[kk+marge]*MAX(WLl(ii,jj),1.0e-6));
}
