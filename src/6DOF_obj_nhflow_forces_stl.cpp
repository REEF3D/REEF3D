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

#include"6DOF_obj.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void sixdof_obj::hydrodynamic_forces_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &WL, bool finalize)
{
	// forcecalc
    if(p->X60==1)
    force_calc_stl(p,d,pgc,WL,finalize);
    
    
    if(p->X60==2)
    {
    triangulation(p,d,pgc);
	reconstruct(p,d);
    force_calc_lsm(p,d,pgc,WL);
        
    deallocate(p,d,pgc);
    }
} 


void sixdof_obj::force_calc_stl(lexer* p, fdm_nhf *d, ghostcell *pgc, slice &WL, bool finalize)
{
    // Pressure + viscous load on the wetted part of the STL surface.
    //
    // Each triangle is clipped against the local free surface z = wd + eta (Sutherland-Hodgman,
    // so vertices lying exactly on the surface are handled), the wetted polygon is fan-triangulated
    // and every sub-triangle is integrated with the 3-point edge-midpoint rule. That rule is exact
    // for quadratic integrands, so for a planar facet both the force (linear hydrostatic pressure)
    // and the moment (lever arm x pressure) are integrated exactly, independent of STL resolution.
    // Moment arms are the quadrature points, not the centroid of the unclipped triangle.
    
    double Fx,Fy,Fz,Fv_x,Fv_y,Fv_z;
    double Xe_p,Ye_p,Ze_p,Xe_v,Ye_v,Ze_v;
    double A;
    double vx[3],vy[3],vz[3];
    double px[4],py[4],pz[4];   // wetted polygon: a triangle clipped by a plane has at most 4 vertices
    
    A=0.0;
    Xe=Ye=Ze=Ke=Me=Ne=0.0;
    Xe_p=Ye_p=Ze_p=Xe_v=Ye_v=Ze_v=0.0;
    
    // Set new time
    curr_time = p->simtime;
    
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
        if(!(xc >= p->originx && xc < p->endx))
        continue;
        
        if(p->j_dir==1 && !(yc >= p->originy && yc < p->endy))
        continue;
        
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
        
        const double fsf_z = p->wd + p->ccslipol4(d->eta,xc,yc);
        
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
            
            for(int q=0; q<3; ++q)
            {
                // non-hydrostatic pressure, optionally sampled X42 mean cell sizes off the wall
                const double pval = p->ccipol7V(d->P, WL, d->bed, mx[q] + p->X42*nx*DSM, 
                                                                  my[q] + p->X42*ny*DSM, 
                                                                  mz[q] + p->X42*nz*DSM);
                // hydrostatic pressure at the quadrature point
                const double hsp = MAX(0.0, (p->wd + p->ccslipol4(d->eta,mx[q],my[q]) - mz[q])*p->W1*fabs(p->W22));
                
                const double w  = A_sub/3.0;
                const double fx = -(pval + hsp)*w*nx;
                const double fy = -(pval + hsp)*w*ny;
                const double fz = -(pval + hsp)*w*nz;
                
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
            
            // Viscous forces at the sub-triangle centroid
            const double gx = (ax + bx + cx)/3.0;
            const double gy = (ay + by + cy)/3.0;
            const double gz = (az + bz + cz)/3.0;
            
            hydrodynamic_viscous_forces_nhflow(p, d, pgc, WL, Fv_x, Fv_y, Fv_z, A_sub, gx, gy, gz, nx, ny, nz);
            
            Xe += Fv_x;
            Ye += Fv_y;
            Ze += Fv_z;
            Ke += (gy - c_(1))*Fv_z - (gz - c_(2))*Fv_y;
            Me += (gz - c_(2))*Fv_x - (gx - c_(0))*Fv_z;
            Ne += (gx - c_(0))*Fv_y - (gy - c_(1))*Fv_x;
            
            Xe_v += Fv_x;
            Ye_v += Fv_y;
            Ze_v += Fv_z;
            
            A += A_sub;
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
    
    // 2D: only surge, heave and pitch exist. Out-of-plane moments from an asymmetric STL
    // triangulation would otherwise enter h_ and tilt the trimesh out of the x-z plane,
    // although u_fb(3) and u_fb(5) are zeroed in update_fbvel.
    if(p->j_dir==0)
    {
        Ye = Ke = Ne = 0.0;
        Ye_p = Ye_v = 0.0;
    }
    
    Fx = Xe;
    Fy = Ye;
    Fz = Ze;
    
	// Add gravity force
	Xe += p->W20*Mass_fb;
	Ye += p->W21*Mass_fb;
	Ze += p->W22*Mass_fb;
    
    if(p->mpirank==0)
    {
    cout<<"Mass_fb: "<<Mass_fb<<" G_fb: "<<p->W22*Mass_fb<<endl;
    cout<<"Fx: "<<Fx<<" Fy: "<<Fy<<" Fz: "<<Fz<<endl;
    cout<<"A_tot: "<<A<<endl;
    cout<<"Xe: "<<Xe<<" Ye: "<<Ye<<" Ze: "<<Ze<<" Ke: "<<Ke<<" Me: "<<Me<<" Ne: "<<Ne<<endl;
    }
    
    // Print results	
    if (p->mpirank==0 && finalize==1) 
    {
        printforce<<curr_time<<" \t "<<Xe<<" \t "<<Ye<<" \t "<<Ze<<" \t "<<Ke
        <<" \t "<<Me<<" \t "<<Ne<<" \t "<<Xe_p<<" \t "<<Ye_p<<" \t "<<Ze_p<<" \t "<<Xe_v<<" \t "<<Ye_v<<" \t "<<Ze_v<<endl;   
    }
}
