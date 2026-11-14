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


#include"6DOF_obj_cfd.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"solver.h"
#include"field1.h"
#include"field2.h"
#include"field3.h"
#include"field4.h"
#include"forcing_share.h"
#include"solid_closed.h"
#include<sys/stat.h>

// X 14 11: the added mass of the hull on the grid (CFD, as NHFLOW X 15 11, 6DOF_obj_nhflow_addedmass.cpp)
//
// With the closed projection (X 49) the pressure on the hull contains the reaction of the fluid to the acceleration
// of the body, -A a, A the added mass. The body is advanced with the loads of the last projection, i.e. with -A a
// of its last acceleration: for a light body that exchange grows from step to step. With the added mass on both
// sides (sixdof_obj::lagged_added_mass),
//
//      (M + A) a = F + A a_prev,
//
// a is the acceleration that satisfies the equation of motion with the reaction of the fluid to the new
// acceleration. X 197 takes a multiple of the displaced fluid (or of the equivalent ellipsoid) for A; here A is the
// one of the grid:
//
// for every free degree of freedom j the forcing sets the share w of a face (forcing_share.h) to the velocity of the
// body, so that the unit acceleration e_j changes the predicted velocity by w e_j dt. Its divergence is the
// right-hand side of the Poisson equation of the projection that was just solved (the matrix a->M with the open
// shares 1-w, the density and the boundaries of the step),
//
//      M P_j = -div(w e_j),
//
// and the loads of P_j on the hull, taken as the loads of the pressure (X 60 1: at the probe points X 42 cells off
// the hull, to the hull with the gradient rho e_j.n of X 34 2), are -A(:,j). One Poisson solve per free degree of
// freedom after the projection of every X 28 th time step (default every step); the solves start from the pressure
// of the last one.
//
// One body, the one-step scheme N 40 5, no mesh refinement (driver_logic_cfd). Until the first added mass is there
// (the first step, a start from a checkpoint) and after a failed solve the coupling is that of X 197.
//
// Log: REEF3D_CFD_6DOF/REEF3D_6DOF_addedmass_<n>.dat, the time, the diagonal of A, A15, A24 and the mass

void sixdof_obj_cfd::added_mass_cfd(lexer *p, fdm *a, ghostcell *pgc, solver *psolv)
{
    if(p->X28>1 && acfd_on_ && p->count%p->X28!=0)
    return;
    
    bool free_dof[6];
    int nfree=0;
    
    for(int qn=0; qn<6; ++qn)
    {
    free_dof[qn] = !p_fixed_dof(p,qn);
    
    if(p->j_dir==0 && (qn==1 || qn==3 || qn==5))
    free_dof[qn] = false;
    
    if(free_dof[qn])
    ++nfree;
    }
    
    if(nfree==0)
    return;
    
    const double starttime = pgc->timer();
    
    if(am_u_==nullptr)
    {
    am_u_ = new field1(p);
    am_v_ = new field2(p);
    am_w_ = new field3(p);
    }
    
    field &ua = *am_u_;
    field &va = *am_v_;
    field &wa = *am_w_;
    
    const int itsave = p->solveriter;
    int itsum = 0;
    
    Eigen::Matrix<double,6,6> A = Eigen::Matrix<double,6,6>::Zero();
    
    for(int qn=0; qn<6; ++qn)
    if(free_dof[qn])
    {
        // velocity of the body at (x,y,z) for the unit velocity of degree of freedom qn (rotations about the
        // centre of gravity)
        auto unit = [&](double x, double y, double z, int c)
        {
            const double r[3] = {x - c_(0), y - c_(1), z - c_(2)};
            
            if(qn<3)
            return (c==qn) ? 1.0 : 0.0;
            
            const int m = qn-3;
            
            // e_m x r
            if(c==(m+1)%3)
            return -r[(m+2)%3];
            
            if(c==(m+2)%3)
            return r[(m+1)%3];
            
            return 0.0;
        };
        
        // the share of the faces taken by the forcing of this body, times the unit velocity
        ULOOP
        ua(i,j,k) = forcing_share_face(p,a,i,j,k,0,false)*unit(p->pos1_x(),p->pos1_y(),p->pos1_z(),0);
        
        VLOOP
        va(i,j,k) = (p->j_dir==1) ? forcing_share_face(p,a,i,j,k,1,false)*unit(p->pos2_x(),p->pos2_y(),p->pos2_z(),1) : 0.0;
        
        WLOOP
        wa(i,j,k) = forcing_share_face(p,a,i,j,k,2,false)*unit(p->pos3_x(),p->pos3_y(),p->pos3_z(),2);
        
        pgc->start1(p,ua,10);
        pgc->start2(p,va,11);
        pgc->start3(p,wa,12);
        
        // right-hand side: the divergence, as pjm_corr::rhs with alpha dt = 1 (D 1 2: through the open fractions
        // of the faces at the solids)
        const bool ap = (p->D1==2 && solid_walls(p));
        int n=0;
        
        LOOP
        {
            const double Ae = (ap && p->flag4[Ip1JK]>0) ? a->apx(i,j,k) : 1.0;
            const double Aw = (ap && p->flag4[Im1JK]>0) ? a->apx(i-1,j,k) : 1.0;
            const double An = (ap && p->j_dir==1 && p->flag4[IJp1K]>0) ? a->apy(i,j,k) : 1.0;
            const double As = (ap && p->j_dir==1 && p->flag4[IJm1K]>0) ? a->apy(i,j-1,k) : 1.0;
            const double At = (ap && p->flag4[IJKp1]>0) ? a->apz(i,j,k) : 1.0;
            const double Ab = (ap && p->flag4[IJKm1]>0) ? a->apz(i,j,k-1) : 1.0;
            
            a->rhsvec.V[n] = -(Ae*ua(i,j,k) - Aw*ua(i-1,j,k))/p->DXN[IP]
                             -(An*va(i,j,k) - As*va(i,j-1,k))/p->DYN[JP]*p->y_dir
                             -(At*wa(i,j,k) - Ab*wa(i,j,k-1))/p->DZN[KP];
            
            // a cell without an open face (inside the body or a solid): the row 1 x = 0 (poisson_pcorr)
            if(a->M.n[n]==0.0 && a->M.s[n]==0.0 && a->M.e[n]==0.0 && a->M.w[n]==0.0 && a->M.t[n]==0.0 && a->M.b[n]==0.0)
            a->rhsvec.V[n] = 0.0;
            
        ++n;
        }
        
        // the pressure of the unit acceleration, from that of the last solve
        if(am_p_[qn]==nullptr)
        {
        am_p_[qn] = new field4(p);
        
        LOOP
        (*am_p_[qn])(i,j,k) = 0.0;
        }
        
        field &PA = *am_p_[qn];
        
        pgc->start4(p,PA,1);
        psolv->start(p,a,pgc,PA,a->rhsvec,5);
        pgc->start4(p,PA,40);
        
        itsum += p->solveriter;
        
        // its loads on the hull
        double Fu[6];
        added_mass_loads(p,a,pgc,PA,qn,Fu);
        
        for(int qm=0; qm<6; ++qm)
        if(free_dof[qm])
        A(qm,qn) = -Fu[qm];
    }
    
    p->solveriter = itsave;
    
    // symmetric, with a positive diagonal
    A = (0.5*(A + A.transpose())).eval();
    
    bool ok = true;
    
    for(int qn=0; qn<6; ++qn)
    if(free_dof[qn] && !(A(qn,qn)>0.0))
    ok = false;
    
    for(int qn=0; qn<6; ++qn)
    for(int qm=0; qm<6; ++qm)
    if(A(qn,qm)!=A(qn,qm))
    ok = false;
    
    // a failed solve: the last added mass is kept, the pressures start from 0 again
    if(ok)
    {
    Acfd_ = A;
    acfd_on_ = true;
    }
    
    else
    for(int qn=0; qn<6; ++qn)
    if(am_p_[qn]!=nullptr)
    LOOP
    (*am_p_[qn])(i,j,k) = 0.0;
    
    if(p->mpirank==0)
    {
        if(p->count%p->P12==0)
        cout<<"6DOF added mass: "<<nfree<<" solves  iter: "<<itsum<<"  time: "<<setprecision(3)<<pgc->timer()-starttime<<(ok ? "" : "  !!! failed, the last one kept")<<endl;
        
        if(!amout_.is_open())
        {
            char str[200];
            mkdir("./REEF3D_CFD_6DOF",0777);
            snprintf(str,200,"./REEF3D_CFD_6DOF/REEF3D_6DOF_addedmass_%i.dat",n6DOF);
            amout_.open(str);
            amout_<<"time \t A11 \t A22 \t A33 \t A44 \t A55 \t A66 \t A15 \t A24 \t mass"<<endl;
        }
        
        amout_<<setprecision(9)<<p->simtime<<" \t "<<A(0,0)<<" \t "<<A(1,1)<<" \t "<<A(2,2)<<" \t "<<A(3,3)<<" \t "<<A(4,4)<<" \t "<<A(5,5)
              <<" \t "<<A(0,4)<<" \t "<<A(1,3)<<" \t "<<Mass_fb<<endl;
    }
}

// the loads (forces, moments about the centre of gravity) of the pressure P of the unit acceleration of degree of
// freedom qn on the hull: the pressure part of forces_stl (X 60 1) without the weight of the fluid, the acceleration
// of the hull at the triangle that of the unit acceleration
void sixdof_obj_cfd::added_mass_loads(lexer *p, fdm *a, ghostcell *pgc, field &P, int qn, double *F)
{
    double Fl[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    
    const double x42 = (p->X49>0 && p->X42<=0.0) ? 1.5 : p->X42;
    const int x34 = (p->X49>0 && p->X34==0) ? 2 : p->X34;
    
    for(int n=0; n<tricount; ++n)
    {
        const double x0 = tri_x[n][0], y0 = tri_y[n][0], z0 = tri_z[n][0];
        const double x1 = tri_x[n][1], y1 = tri_y[n][1], z1 = tri_z[n][1];
        const double x2 = tri_x[n][2], y2 = tri_y[n][2], z2 = tri_z[n][2];
        
        const double xc = (x0 + x1 + x2)/3.0;
        const double yc = (y0 + y1 + y2)/3.0;
        const double zc = (z0 + z1 + z2)/3.0;
        
        if(!(xc >= p->originx && xc < p->endx && yc >= p->originy && yc < p->endy && zc >= p->originz && zc < p->endz))
        continue;
        
        i = p->posc_i(xc);
        j = p->posc_j(yc);
        k = p->posc_k(zc);
        
        // area (Heron) and outward normal
        const double at = sqrt(pow(x1-x0,2.0) + pow(y1-y0,2.0) + pow(z1-z0,2.0));
        const double bt = sqrt(pow(x1-x2,2.0) + pow(y1-y2,2.0) + pow(z1-z2,2.0));
        const double ct = sqrt(pow(x2-x0,2.0) + pow(y2-y0,2.0) + pow(z2-z0,2.0));
        const double st = 0.5*(at+bt+ct);
        const double A_triang = sqrt(MAX(0.0,st*(st-at)*(st-bt)*(st-ct)));
        
        double nx = (y1 - y0) * (z2 - z0) - (y2 - y0) * (z1 - z0);
        double ny = (x2 - x0) * (z1 - z0) - (x1 - x0) * (z2 - z0); 
        double nz = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
        
        const double norm = sqrt(nx*nx + ny*ny + nz*nz);
        
        nx /= norm > 1.0e-20 ? norm : 1.0e20;
        ny /= norm > 1.0e-20 ? norm : 1.0e20;
        nz /= norm > 1.0e-20 ? norm : 1.0e20;
        
        if(p->j_dir==0)
        ny = 0.0;
        
        const double ddx = x42*nx*p->DXP[IP], ddy = x42*ny*p->DYP[JP], ddz = x42*nz*p->DZP[KP];
        
        double p_int = p->ccipol4a(P,xc+ddx,yc+ddy,zc+ddz);
        
        // to the hull with the pressure gradient at the wall, dp/dn = -rho n.a
        if(x34==2 && x42>0.0)
        {
            const double rho_m = p->ccipol4a(a->ro,xc+0.5*ddx,yc+0.5*ddy,zc+0.5*ddz);
            const double rx = xc - c_(0), ry = yc - c_(1), rz = zc - c_(2);
            double ax = 0.0, ay = 0.0, az = 0.0;
            
            if(qn==0) ax = 1.0;
            if(qn==1) ay = 1.0;
            if(qn==2) az = 1.0;
            if(qn==3) {ay = -rz; az =  ry;}
            if(qn==4) {ax =  rz; az = -rx;}
            if(qn==5) {ax = -ry; ay =  rx;}
            
            p_int += rho_m*(ddx*ax + ddy*ay + ddz*az);
        }
        
        const double Fx = -nx*p_int*A_triang;
        const double Fy = (p->j_dir==0) ? 0.0 : -ny*p_int*A_triang;
        const double Fz = -nz*p_int*A_triang;
        
        Fl[0] += Fx;
        Fl[1] += Fy;
        Fl[2] += Fz;
        Fl[3] += (yc - c_(1))*Fz - (zc - c_(2))*Fy;
        Fl[4] += (zc - c_(2))*Fx - (xc - c_(0))*Fz;
        Fl[5] += (xc - c_(0))*Fy - (yc - c_(1))*Fx;
    }
    
    for(int q=0; q<6; ++q)
    F[q] = pgc->globalsum(Fl[q]);
}
