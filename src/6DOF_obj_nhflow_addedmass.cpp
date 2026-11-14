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

#include"6DOF_obj_nhflow.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"solver.h"
#include"nhflow_curv.h"
#include"nhflow_bodyflux.h"
#include<sys/stat.h>

// X 15 11: coupling of the body and the water in the projection
//
// With the body in the projection (open fractions of the links, 6DOF_obj_forcing_nhflow.cpp) the pressure on the
// hull contains the reaction of the water to the acceleration of the body, -A a, A the added mass. The body is
// advanced with the loads of the last projection, i.e. with -A a of its last acceleration: for A > M, a floating
// body, that exchange grows from stage to stage. With the added mass on both sides,
//
//      (M + A) a = F + A a_prev,
//
// a is the acceleration of the body that satisfies F_rest - A a = M a, the equation of motion with the reaction of
// the water to the new acceleration: the implicit coupling, without iterations, for A the added mass of the grid.
// An error dA of A leaves a lag (dA/(M + A)) (a - a_prev), stable for |dA| < M + A.
//
// A on the grid: for every free degree of freedom j the velocity field of the direct forcing for the unit
// acceleration e_j of the body (the weights of forcing_kernel_fraction), its divergence as the right-hand side of
// the Poisson equation of the stage (the matrix of the projection that was just solved, with the open fractions
// and the water level of the stage), and the loads of that pressure on the hull (force_calc_stl with load_P),
// A(:,j) = -F(P_j) + G(:,j), G the water layer of the pressure samples. One Poisson solve per free degree of
// freedom, at the first stage of every X 28 th time step (default every step).
//
// Log: REEF3D_NHFLOW_6DOF/REEF3D_6DOF_addedmass_<n>.dat, the time and the diagonal of A

void sixdof_obj_nhflow::added_mass_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, solver *psolv, slice &WL)
{
    const int nc = p->imax*p->jmax*(p->kmax+2);
    
    if(UA==nullptr)
    {
        p->Darray(UA,nc);
        p->Darray(VA,nc);
        p->Darray(WA,nc);
        p->Darray(PA,nc);
        p->Darray(UBA,nc);
        p->Darray(VBA,nc);
        p->Darray(WBA,nc);
    }
    
    const double Fsave[6] = {Xe,Ye,Ze,Ke,Me,Ne};
    const Eigen::Matrix<double,6,6> Gsave = Gdf_;
    const int itsave = p->solveriter;
    
    Eigen::Matrix<double,6,6> A = Eigen::Matrix<double,6,6>::Zero();
    bool free_dof[6];
    
    for(int qn=0; qn<6; ++qn)
    free_dof[qn] = !rb.fixed(qn);
    
    for(int qn=0; qn<6; ++qn)
    if(free_dof[qn])
    {
        // velocity of the direct forcing for the unit velocity of degree of freedom qn
        auto unit = [&](double rx, double ry, double rz, double *v)
        {
            v[0] = v[1] = v[2] = 0.0;
            
            if(qn<3)
            v[qn] = 1.0;
            
            if(qn==3)
            {
            v[1] = -rz;
            v[2] =  ry;
            }
            
            if(qn==4)
            {
            v[0] =  rz;
            v[2] = -rx;
            }
            
            if(qn==5)
            {
            v[0] = -ry;
            v[1] =  rx;
            }
        };
        
        LOOP
        {
            UA[IJK] = VA[IJK] = WA[IJK] = 0.0;
            UBA[IJK] = VBA[IJK] = WBA[IJK] = 0.0;
            
            if(p->wet[IJ]==0)
            continue;
            
            const double rx = p->pos_x() - c_(0);
            const double ry = p->pos_y() - c_(1);
            const double rz = p->pos_z() - c_(2);
            
            double v[3];
            
            // velocity of the body at the cell centre and at the centres of the x and y faces
            unit(rx,ry,rz,v);
            WBA[IJK] = v[2];
            
            if(p->wet[Ip1J]==1)
            {
            unit(rx + 0.5*p->DXP[IP], ry, 0.5*(p->ZSP[IJK] + p->ZSP[Ip1JK]) - c_(2), v);
            UBA[IJK] = v[0];
            }
            
            if(p->j_dir==1 && p->wet[IJp1]==1)
            {
            unit(rx, ry + 0.5*p->DYP[JP], 0.5*(p->ZSP[IJK] + p->ZSP[IJp1K]) - c_(2), v);
            VBA[IJK] = v[1];
            }
            
            const double kap = 1.0 - d->MBZ[IJK];
            
            if(kap<=1.0e-3)
            continue;
            
            unit(rx,ry,rz,v);
            
            const double tau = MAX(2.0*kap - 1.0, 0.0);
            
            UA[IJK] = tau*v[0];
            VA[IJK] = tau*v[1];
            WA[IJK] = tau*v[2];
        }
        
        pgc->start4V(p,UA,10);
        pgc->start4V(p,VA,11);
        pgc->start4V(p,WA,12);
        pgc->start4V(p,UBA,1);
        pgc->start4V(p,VBA,1);
        pgc->start4V(p,WBA,1);
        
        // divergence at the pressure nodes, as nhflow_pjm::rhs with alpha dt = 1
        int n=0;
        
        FBASELOOP
        {
        d->rhsvec.V[n] = 0.0;
        ++n;
        }
        
        n=0;
        LOOP
        {
            if(p->wet[IJ]==1 && p->deep[IJ]==1 && p->flag7[FIJK]>0)
            {
                const bool bednode = (k==0);
                const double fac = p->DZN[KM1]/(p->DZN[KP]+p->DZN[KM1]);
                
                const double Ub1 = bednode ? 0.0 : UA[Im1JKm1], Ub2 = bednode ? 0.0 : UA[Ip1JKm1];
                const double Vb1 = bednode ? 0.0 : VA[IJm1Km1], Vb2 = bednode ? 0.0 : VA[IJp1Km1];
                const double Ubc = bednode ? 0.0 : UA[IJKm1],   Vbc = bednode ? 0.0 : VA[IJKm1];
                
                const double U1 = (1.0-fac)*UA[Im1JK] + fac*Ub1;
                const double U2 = (1.0-fac)*UA[Ip1JK] + fac*Ub2;
                const double V1 = (1.0-fac)*VA[IJm1K] + fac*Vb1;
                const double V2 = (1.0-fac)*VA[IJp1K] + fac*Vb2;
                
                const double dUdz = (UA[IJK] - Ubc)/p->DZP[KM1];
                const double dVdz = (VA[IJK] - Vbc)/p->DZP[KM1];
                const double dWdz = p->sigz[IJ]*(WA[IJK] - WA[IJKm1])/p->DZP[KM1];
                
                if(p->geo_curv==1)
                d->rhsvec.V[n] = -(nhflow_curv::divh_level(p,d,UA,VA,i,j,k) + p->sigx[FIJK]*dUdz + p->sigy[FIJK]*dVdz + dWdz);
                else
                d->rhsvec.V[n] = -((U2-U1)/(p->DXP[IP] + p->DXP[IM1]) + p->sigx[FIJK]*dUdz
                                 + (V2-V1)/(p->DYP[JP] + p->DYP[JM1]) + p->sigy[FIJK]*dVdz + dWdz);
            }
        ++n;
        }
        
        // the flux of the body through the closed part of its faces
        nhflow_bodyflux::rhs(p,d,UA,VA,WA,UBA,VBA,WBA,1.0,false);
        
        // no pressure at the nodes inside the body (rows without neighbours)
        n=0;
        LOOP
        {
        if(d->M.n[n]==0.0 && d->M.s[n]==0.0 && d->M.t[n]==0.0 && d->M.b[n]==0.0 && d->M.e[n]==0.0 && d->M.w[n]==0.0)
        d->rhsvec.V[n] = 0.0;
        ++n;
        }
        
        for(int q=0; q<nc; ++q)
        PA[q] = 0.0;
        
        psolv->startF(p,pgc,PA,d->rhsvec,d->M,8);
        pgc->start7P(p,PA,540);
        
        // loads of the pressure of the unit acceleration
        load_P = PA;
        force_calc_stl(p,d,pgc,WL,false);
        load_P = nullptr;
        
        const double Fu[6] = {Xe,Ye,Ze,Ke,Me,Ne};
        
        for(int qm=0; qm<6; ++qm)
        if(free_dof[qm])
        A(qm,qn) = -Fu[qm] + Gdf_(qm,qn);
    }
    
    Xe = Fsave[0];
    Ye = Fsave[1];
    Ze = Fsave[2];
    Ke = Fsave[3];
    Me = Fsave[4];
    Ne = Fsave[5];
    Gdf_ = Gsave;
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
    
    // a failed solve: the last added mass is kept
    if(ok)
    {
    Adf_ = A;
    adf_on_ = true;
    }
    
    if(p->mpirank==0 && p->count%p->X19==0)
    {
        if(!amout.is_open())
        {
            char str[200];
            mkdir("./REEF3D_NHFLOW_6DOF",0777);
            snprintf(str,200,"./REEF3D_NHFLOW_6DOF/REEF3D_6DOF_addedmass_%i.dat",n6DOF);
            amout.open(str);
            amout<<"time \t A11 \t A22 \t A33 \t A44 \t A55 \t A66 \t A15 \t A24 \t mass"<<endl;
        }
        
        amout<<setprecision(9)<<p->simtime<<" \t "<<A(0,0)<<" \t "<<A(1,1)<<" \t "<<A(2,2)<<" \t "<<A(3,3)<<" \t "<<A(4,4)<<" \t "<<A(5,5)
             <<" \t "<<A(0,4)<<" \t "<<A(1,3)<<" \t "<<Mass_fb<<endl;
    }
}

void sixdof_obj_nhflow::fraction_stabilisation(lexer *p)
{
    if(!adf_on_)
    return;
    
    lagged_added_mass(p,Adf_);
}
