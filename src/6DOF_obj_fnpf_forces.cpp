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
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<algorithm>

void sixdof_obj::ray_cast_fnpf(lexer *p, fdm_fnpf *c, ghostcell *pgc, double *FBF, slice &foot)
{
    // Vertical ray per sigma column: the z-crossings with the trimesh give an exact
    // inside/outside test (parity) for every sigma node of the column.
    // Marks FBF = n6DOF+1 inside the body and foot = n6DOF+1 where the top node
    // (the free-surface node k=knoz) lies inside, i.e. the body footprint.
    // The caller zeroes FBF/foot beforehand and exchanges them afterwards.
    
    const double id = double(n6DOF+1);
    const int ny = p->knoy;
    
    vector<vector<double> > hits(p->knox*p->knoy);
    
    // 2D: the STL is extruded in y, cast in its mid-plane
    double ymid=0.0;
    
    if(p->j_dir==0 && tricount>0)
    {
        double ylo=1.0e20, yhi=-1.0e20;
        
        for(n=0; n<tricount; ++n)
        for(int q=0; q<3; ++q)
        {
        ylo = MIN(ylo,tri_y[n][q]);
        yhi = MAX(yhi,tri_y[n][q]);
        }
        
        ymid = 0.5*(ylo+yhi);
    }
    
    // small irrational offset against rays through edges and vertices
    const double ex = 0.1234567891e-6*p->DXM;
    const double ey = 0.3141592653e-6*p->DXM;
    
    for(n=0; n<tricount; ++n)
    {
        const double x0 = tri_x[n][0], x1 = tri_x[n][1], x2 = tri_x[n][2];
        const double y0 = tri_y[n][0], y1 = tri_y[n][1], y2 = tri_y[n][2];
        const double z0 = tri_z[n][0], z1 = tri_z[n][1], z2 = tri_z[n][2];
        
        const double det = (y1-y2)*(x0-x2) + (x2-x1)*(y0-y2);
        
        // vertical facets are never crossed by a vertical ray
        if(fabs(det)<1.0e-30)
        continue;
        
        const double txmin = MIN(x0,MIN(x1,x2));
        const double txmax = MAX(x0,MAX(x1,x2));
        const double tymin = MIN(y0,MIN(y1,y2));
        const double tymax = MAX(y0,MAX(y1,y2));
        
        ILOOP
        {
            const double xr = p->XP[IP] + ex;
            
            if(xr<txmin || xr>txmax)
            continue;
            
            JLOOP
            {
                const double yr = ((p->j_dir==0) ? ymid : p->YP[JP]) + ey;
                
                if(yr<tymin || yr>tymax)
                continue;
                
                const double l0 = ((y1-y2)*(xr-x2) + (x2-x1)*(yr-y2))/det;
                const double l1 = ((y2-y0)*(xr-x2) + (x0-x2)*(yr-y2))/det;
                const double l2 = 1.0 - l0 - l1;
                
                if(l0<0.0 || l1<0.0 || l2<0.0)
                continue;
                
                hits[i*ny+j].push_back(l0*z0 + l1*z1 + l2*z2);
            }
        }
    }
    
    ILOOP
    JLOOP
    {
        vector<double> &h = hits[i*ny+j];
        
        if(h.empty())
        continue;
        
        // a ray through a shared edge or vertex hits every adjacent facet:
        // count coincident crossings once
        sort(h.begin(),h.end());
        
        const double tol = 1.0e-9*p->DXM;
        size_t nu=1;
        
        for(size_t m=1; m<h.size(); ++m)
        if(h[m]-h[nu-1]>tol)
        h[nu++] = h[m];
        
        h.resize(nu);
        
        FKLOOP
        if(p->flag7[FIJK]>0)
        {
            const double z = p->ZSN[FIJK];
            
            int cnt=0;
            for(size_t m=0; m<h.size(); ++m)
            if(h[m]>z)
            ++cnt;
            
            if(cnt%2==1)
            FBF[FIJK] = id;
        }
        
        k = p->knoz;
        if(fabs(FBF[FIJK]-id)<0.5)
        foot(i,j) = id;
    }
}

void sixdof_obj::face_data_fnpf(lexer *p, fdm_fnpf *c, ghostcell *pgc, int mode, double *FBF, double *U, double *V, double *W)
{
    // Neumann data of the staircase body faces, stored at the body nodes:
    // the Laplace assembly imposes d(f)/dx_face = (U,V,W).e_face on every
    // fluid/body face.
    //  mode -1: rigid-body velocity           u_c + w x r     (Laplace for phi)
    //  mode -2: centripetal acceleration      w x (w x r)     (psi_0 = phi_t, m-terms neglected)
    //  mode 0-2: unit translation             e_j             (added-mass modes)
    //  mode 3-5: unit rotation                e_j x r
    
    const double id = double(n6DOF+1);
    const Eigen::Vector3d w(u_fb(3),u_fb(4),u_fb(5));
    const Eigen::Vector3d uc(u_fb(0),u_fb(1),u_fb(2));
    
    ILOOP
    JLOOP
    FKLOOP
    if(fabs(FBF[FIJK]-id)<0.5)
    {
        const Eigen::Vector3d r(p->XP[IP]-c_(0), p->YP[JP]-c_(1), p->ZSN[FIJK]-c_(2));
        Eigen::Vector3d v;
        
        if(mode==-1)
        v = uc + w.cross(r);
        
        else if(mode==-2)
        v = w.cross(w.cross(r));
        
        else if(mode<3)
        {
        v.setZero();
        v(mode) = 1.0;
        }
        
        else
        {
        Eigen::Vector3d e = Eigen::Vector3d::Zero();
        e(mode-3) = 1.0;
        v = e.cross(r);
        }
        
        U[FIJK] = v(0);
        V[FIJK] = v(1);
        W[FIJK] = v(2);
    }
}

void sixdof_obj::forces_fnpf(lexer *p, fdm_fnpf *c, ghostcell *pgc, double *psi0, double **psi, bool computeA)
{
    // Pressure from Bernoulli with psi = phi_t split as psi = psi_0 + sum_j a_j psi_j:
    //   p = -rho*(psi_0 + 0.5|grad phi|^2) + rho*g*(wd - z)  - rho*sum_j a_j psi_j
    // F_0 = -int p_0 N dS goes to Xe..Ne, the last term is -A a with
    //   A_ij = -rho int psi_j N_i dS,   N = (n, r x n),  n pointing into the fluid.
    // Integration over the trimesh clipped at the local free surface (as NHFLOW),
    // fluid quantities sampled half a mean cell off the hull.
    
    const double rho = p->W1;
    const double grav = fabs(p->W22);
    const double del = 0.5*DSM;
    
    double F[3]  = {0.0,0.0,0.0};
    double Mo[3] = {0.0,0.0,0.0};
    double Am[36];
    double Atot=0.0;
    
    for(int m=0; m<36; ++m)
    Am[m]=0.0;
    
    double vx[3],vy[3],vz[3];
    double px[4],py[4],pz[4];
    
    curr_time = p->simtime;
    
    for(n=0; n<tricount; ++n)
    {     
        for(int q=0; q<3; ++q)
        {
            vx[q] = tri_x[n][q];
            vy[q] = tri_y[n][q];
            vz[q] = tri_z[n][q];
        }
        
		const double xc = (vx[0] + vx[1] + vx[2])/3.0;
		const double yc = (vy[0] + vy[1] + vy[2])/3.0;
        
        // ownership by the triangle centroid (decomposition in x and y only)
        if(!(xc >= p->originx && xc < p->endx))
        continue;
        
        if(p->j_dir==1 && !(yc >= p->originy && yc < p->endy))
        continue;
        
        double nx = (vy[1] - vy[0])*(vz[2] - vz[0]) - (vy[2] - vy[0])*(vz[1] - vz[0]);
        double ny = (vx[2] - vx[0])*(vz[1] - vz[0]) - (vx[1] - vx[0])*(vz[2] - vz[0]); 
        double nz = (vx[1] - vx[0])*(vy[2] - vy[0]) - (vx[2] - vx[0])*(vy[1] - vy[0]);
        const double norm = sqrt(nx*nx + ny*ny + nz*nz);
        
        if(norm<1.0e-20)
        continue;
        
        nx /= norm;
        ny /= norm;
        nz /= norm;
        
        if(p->j_dir==0)
        {
            if(fabs(ny)>0.9)
            continue;
            
            ny = 0.0;
        }
        
        const double fsf_z = p->wd + p->ccslipol4(c->eta,xc,yc);
        
        // wetted part z <= fsf_z
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
            
            if(in_q != in_r)
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
            
            const double mx[3] = {0.5*(ax + bx), 0.5*(bx + cx), 0.5*(cx + ax)};
            const double my[3] = {0.5*(ay + by), 0.5*(by + cy), 0.5*(cy + ay)};
            const double mz[3] = {0.5*(az + bz), 0.5*(bz + cz), 0.5*(cz + az)};
            
            for(int q=0; q<3; ++q)
            {
                const double xp = mx[q] + del*nx;
                const double yp = my[q] + del*ny;
                const double zp = mz[q] + del*nz;
                
                const double uval = p->ccipol7V(c->U, c->WL, c->bed, xp, yp, zp);
                const double vval = p->ccipol7V(c->V, c->WL, c->bed, xp, yp, zp);
                const double wval = p->ccipol7V(c->W, c->WL, c->bed, xp, yp, zp);
                const double ps0  = p->ccipol7V(psi0,  c->WL, c->bed, xp, yp, zp);
                
                const double pres = -rho*(ps0 + 0.5*(uval*uval + vval*vval + wval*wval)) 
                                  + rho*grav*(p->wd - mz[q]);
                
                const double wgt = A_sub/3.0;
                
                const double fx = -pres*wgt*nx;
                const double fy = -pres*wgt*ny;
                const double fz = -pres*wgt*nz;
                
                const double rx = mx[q] - c_(0);
                const double ry = my[q] - c_(1);
                const double rz = mz[q] - c_(2);
                
                F[0] += fx;
                F[1] += fy;
                F[2] += fz;
                Mo[0] += ry*fz - rz*fy;
                Mo[1] += rz*fx - rx*fz;
                Mo[2] += rx*fy - ry*fx;
                
                if(computeA)
                {
                    const double N[6] = {nx, ny, nz, ry*nz - rz*ny, rz*nx - rx*nz, rx*ny - ry*nx};
                    
                    for(int jj=0; jj<6; ++jj)
                    if(psi[jj]!=nullptr)
                    {
                        const double pj = p->ccipol7V(psi[jj], c->WL, c->bed, xp, yp, zp);
                        
                        for(int ii=0; ii<6; ++ii)
                        Am[ii*6+jj] += -rho*pj*N[ii]*wgt;
                    }
                }
            }
            
            Atot += A_sub;
        }
	}
    
    Atot = pgc->globalsum(Atot);
    
    for(int m=0; m<3; ++m)
    {
    F[m]  = pgc->globalsum(F[m]);
    Mo[m] = pgc->globalsum(Mo[m]);
    }
    
    if(p->j_dir==0)
    F[1] = Mo[0] = Mo[2] = 0.0;
    
    Fx = F[0];
    Fy = F[1];
    Fz = F[2];
    
    Xe = F[0] + p->W20*Mass_fb;
    Ye = F[1] + p->W21*Mass_fb;
    Ze = F[2] + p->W22*Mass_fb;
    Ke = Mo[0];
    Me = Mo[1];
    Ne = Mo[2];
    
    if(computeA)
    {
        for(int m=0; m<36; ++m)
        Am[m] = pgc->globalsum(Am[m]);
        
        for(int ii=0; ii<6; ++ii)
        for(int jj=0; jj<6; ++jj)
        Aadd_(ii,jj) = 0.5*(Am[ii*6+jj] + Am[jj*6+ii]);
        
        am_on_ = true;
    }
    
    if(p->mpirank==0 && (p->count%p->P12==0))
    {
    cout<<"Fx: "<<Fx<<" Fy: "<<Fy<<" Fz: "<<Fz<<"  A_tot: "<<Atot<<endl;
    cout<<"Xe: "<<Xe<<" Ye: "<<Ye<<" Ze: "<<Ze<<" Ke: "<<Ke<<" Me: "<<Me<<" Ne: "<<Ne<<endl;
    
    if(computeA)
    cout<<"A11: "<<Aadd_(0,0)<<" A33: "<<Aadd_(2,2)<<" A55: "<<Aadd_(4,4)<<endl;
    }
}
