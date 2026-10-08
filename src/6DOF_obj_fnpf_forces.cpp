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

#include"6DOF_obj_fnpf.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include<algorithm>

void sixdof_obj_fnpf::ray_cast_fnpf(lexer *p, fdm_fnpf *c, ghostcell *pgc, double *FBF, slice &foot)
{
    // Vertical ray per sigma column: the z-crossings with the trimesh give an exact
    // inside/outside test (parity) for every sigma node of the column.
    // Marks FBF = n6DOF+1 inside the body and foot = n6DOF+1 where the top node
    // (the free-surface node k=knoz) lies inside, i.e. the body footprint.
    // The caller zeroes FBF/foot beforehand and exchanges them afterwards.
    
    const double id = double(n6DOF+1);
    const int ny = p->knoy;
    
    // sorted, unique crossings of the column rays (geometry core)
    vector<vector<double> > hits;
    
    georay.column_hits(p,tri_x,tri_y,tri_z,0,tricount,p->DXM,hits);
    
    ILOOP
    JLOOP
    {
        vector<double> &h = hits[i*ny+j];
        
        if(h.empty())
        continue;
        
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

void sixdof_obj_fnpf::face_data_fnpf(lexer *p, fdm_fnpf *c, ghostcell *pgc, int mode, double *FBF, double *U, double *V, double *W)
{
    // Neumann data of the staircase body faces, stored at the body nodes:
    // the Laplace assembly imposes d(f)/dx_face = (U,V,W).e_face on every
    // fluid/body face.
    //  mode -1: rigid-body velocity           u_c + w x r     (Laplace for phi)
    //  mode -2: psi_0, the velocity-dependent part of the time derivative of phi:
    //           chi_on (one body): psi_0 = chi = phi_t + V.grad(phi), V = u_c + w x r, the time
    //           derivative following the body points. Its Neumann data on the rigid body are
    //           exactly  w x (w x r) + V x w  (from d/dt (grad(phi).n - V.n) = 0 along the body,
    //           grad(V.grad(phi)) = (V.grad)grad(phi) - w x grad(phi) for a rigid V): no second
    //           derivatives of phi, no m-terms to neglect
    //           otherwise psi_0 = phi_t with  w x (w x r)  (m-terms neglected)
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
        {
        v = w.cross(w.cross(r));
        
        if(chi_on)
        v += (uc + w.cross(r)).cross(w);
        }
        
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

bool sixdof_obj_fnpf::chi_on = false;

void sixdof_obj_fnpf::chi_fsf(lexer *p, fdm_fnpf *c, slice &D, slice &foot)
{
    if(!chi_on)
    return;
    
    const Eigen::Vector3d w(u_fb(3),u_fb(4),u_fb(5));
    const Eigen::Vector3d uc(u_fb(0),u_fb(1),u_fb(2));
    
    FFILOOP4
    if(foot(i,j)<0.5)
    {
        const Eigen::Vector3d r(p->XP[IP]-c_(0), p->YP[JP]-c_(1), p->ZSN[FIJK]-c_(2));
        const Eigen::Vector3d V = uc + w.cross(r);
        
        D(i,j) += V(0)*c->U[FIJK] + V(1)*c->V[FIJK] + V(2)*c->W[FIJK];
    }
}

void sixdof_obj_fnpf::chi_mean(lexer *p, ghostcell *pgc, slice &D, slice &foot)
{
    // X 18 tau: in a statistically steady state (steady manoeuvre, regular waves) the time mean
    // of the exact chi = d phi/dt following a body point is zero. Next to a moving
    // surface-piercing hull the free-surface data psi_D + V.grad(phi) are not: columns entering
    // and leaving the footprint leave a mean residual near the waterline, which biases the hull
    // pressure (in a steady turn it gives a drift-yaw interaction of about a third of the Munk
    // moment, validation 11). The running mean of chi is kept on a body-frame grid around the
    // body (horizontal extent of the body plus 5 cells) and removed from the free-surface data
    // there. The mean is a second-order low-pass (two exponential stages with time constant
    // tau each): the oscillating part (waves, heave and pitch) is kept, at a frequency w the mean
    // takes up about 1/(w tau)^2 of it with a phase near 180 deg, i.e. no spurious damping. (A
    // first-order mean takes up 1/(w tau) at 90 deg, which acts as negative damping in heave and
    // pitch and made a free hull unstable.)
    if(!chi_on || p->X18<=0.0)
    return;
    
    const double psi = atan2(R_(1,0),R_(0,0));
    const double cp = cos(psi), sp = sin(psi);
    
    if(chim_nx==0)
    {
        // horizontal extent of the body in the body frame
        double xmin=1.0e20, xmax=-1.0e20, ymin=1.0e20, ymax=-1.0e20;
        for(int n=0; n<tricount; ++n)
        for(int q=0; q<3; ++q)
        {
            const double rx = tri_x[n][q]-c_(0), ry = tri_y[n][q]-c_(1);
            const double xb =  cp*rx + sp*ry;
            const double yb = -sp*rx + cp*ry;
            xmin = MIN(xmin,xb); xmax = MAX(xmax,xb);
            ymin = MIN(ymin,yb); ymax = MAX(ymax,yb);
        }
        xmin = pgc->globalmin(xmin); ymin = pgc->globalmin(ymin);
        xmax = pgc->globalmax(xmax); ymax = pgc->globalmax(ymax);
        
        chim_h = p->DXM;
        const double margin = 5.0*chim_h;
        chim_x0 = xmin - margin;
        chim_y0 = ymin - margin;
        chim_nx = int((xmax - xmin + 2.0*margin)/chim_h) + 2;
        chim_ny = int((ymax - ymin + 2.0*margin)/chim_h) + 2;
        chim_.assign(chim_nx*chim_ny,0.0);
        chim2_.assign(chim_nx*chim_ny,0.0);
    }
    
    // update once per time step (forces() is called in every stage)
    const bool update = (p->count != chim_count);
    chim_count = p->count;
    const double alpha = MIN(1.0, p->dt/p->X18);
    
    SLICELOOP4
    if(foot(i,j)<0.5)
    {
        const double rx = p->XP[IP]-c_(0), ry = p->YP[JP]-c_(1);
        const double xb =  cp*rx + sp*ry;
        const double yb = -sp*rx + cp*ry;
        const int bi = int(floor((xb - chim_x0)/chim_h + 0.5));
        const int bj = int(floor((yb - chim_y0)/chim_h + 0.5));
        
        if(bi<0 || bj<0 || bi>=chim_nx || bj>=chim_ny)
        continue;
        
        double &m1 = chim_[bj*chim_nx + bi];
        double &m2 = chim2_[bj*chim_nx + bi];
        
        if(update)
        {
        m1 += alpha*(D(i,j) - m1);
        m2 += alpha*(m1 - m2);
        }
        
        D(i,j) -= m2;
    }
    
    pgc->gcsl_start4(p,D,50);
}

void sixdof_obj_fnpf::forces_fnpf(lexer *p, fdm_fnpf *c, ghostcell *pgc, double *psi0, double **psi, bool computeA)
{
    // Pressure from Bernoulli with psi = phi_t split as psi = psi_0 + sum_j a_j psi_j:
    //   p = -rho*(psi_0 + 0.5|grad phi|^2) + rho*g*(wd - z)  - rho*sum_j a_j psi_j
    // (chi_on: psi_0 - V.grad(phi) instead of psi_0, see face_data_fnpf mode -2)
    // F_0 = -int p_0 N dS goes to Xe..Ne, the last term is -A a with
    //   A_ij = -rho int psi_j N_i dS,   N = (n, r x n),  n pointing into the fluid.
    // Integration over the trimesh clipped at the local free surface (as NHFLOW),
    // fluid quantities sampled half a mean cell off the hull.
    
    fnpf_force_sum S;
    forces_fnpf_zero(p,S);
    
    forces_fnpf_sum(p,c,psi0,psi,computeA,0.5*DSM,nullptr,S);
    forces_fnpf_set(p,pgc,S,computeA);
}

void sixdof_obj_fnpf::forces_fnpf_zero(lexer *p, fnpf_force_sum &S)
{
    curr_time = p->simtime;
    
    for(int m=0; m<3; ++m)
    S.F[m] = S.Mo[m] = 0.0;
    
    for(int m=0; m<36; ++m)
    S.Am[m]=0.0;
    
    S.Atot=0.0;
}

void sixdof_obj_fnpf::forces_fnpf_sum(lexer *p, fdm_fnpf *c, double *psi0, double **psi, bool computeA, double del,
                                 const std::function<bool(double,double)> *own, fnpf_force_sum &S)
{
    const double rho = p->W1;
    const double grav = fabs(p->W22);
    
    double *F = S.F;
    double *Mo = S.Mo;
    double *Am = S.Am;
    double &Atot = S.Atot;
    
    double vx[3],vy[3],vz[3];
    double px[4],py[4],pz[4];
    
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
        
        // ownership by the triangle centroid (decomposition in x and y only); with several
        // grids own() decides (rank box and the finest grid that holds the centroid)
        if(own!=nullptr)
        {
            if(!(*own)(xc,yc))
            continue;
        }
        else
        {
            if(!(xc >= p->originx && xc < p->endx))
            continue;
            
            if(p->j_dir==1 && !(yc >= p->originy && yc < p->endy))
            continue;
        }
        
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
                
                // chi_on: phi_t = chi - V.grad(phi) with V the rigid-body velocity at the sample point
                double vgp = 0.0;
                
                if(chi_on)
                {
                    const Eigen::Vector3d Vs = Eigen::Vector3d(u_fb(0),u_fb(1),u_fb(2))
                                             + Eigen::Vector3d(u_fb(3),u_fb(4),u_fb(5)).cross(Eigen::Vector3d(xp-c_(0),yp-c_(1),zp-c_(2)));
                    vgp = Vs(0)*uval + Vs(1)*vval + Vs(2)*wval;
                }
                
                const double pres = -rho*(ps0 - vgp + 0.5*(uval*uval + vval*vval + wval*wval)) 
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
}

void sixdof_obj_fnpf::forces_fnpf_set(lexer *p, ghostcell *pgc, fnpf_force_sum &S, bool computeA)
{
    double *F = S.F;
    double *Mo = S.Mo;
    double *Am = S.Am;
    double Atot = S.Atot;
    
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
