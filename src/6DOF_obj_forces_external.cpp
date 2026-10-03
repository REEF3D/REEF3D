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
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"mooring.h"
#include"net_interface.h"

void sixdof_obj::externalForces_cfd(lexer *p, fdm* a, ghostcell *pgc, double alpha, bool finalize)
{
    Xext = Yext = Zext = Kext = Mext = Next = 0.0;
    
    // Mooring forces
	if (p->X310>0)
	mooringForces(p,pgc,alpha);

    // Net forces
	if (p->X320>0)
	netForces_cfd(p,a,pgc,alpha,finalize);
    
    // VRANS forces
}

void sixdof_obj::externalForces_nhflow(lexer *p, fdm_nhf* d, ghostcell *pgc, double alpha, bool finalize)
{
    Xext = Yext = Zext = Kext = Mext = Next = 0.0;

    // Mooring forces
	if (p->X310>0)
	mooringForces(p,pgc,alpha);

    // Net forces
	if (p->X320>0)
	netForces_nhflow(p,d,pgc,alpha,finalize);
    
    // Membrane forces (X 330): load of rigid or flexible membranes attached to the body
    if (p->X330>0)
    {
        double X,Y,Z,K,M,N;
        pnetinter->membraneForces_nhflow(p,c_,quatRotMat,X,Y,Z,K,M,N);
        
        Xext += X;
        Yext += Y;
        Zext += Z;
        Kext += K;
        Mext += M;
        Next += N;
    }
    
    // VRANS forces
}

void sixdof_obj::membrane_forcing_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, 
                                         double *UH, double *VH, double *WH, slice &WL)
{
    // current body kinematics for the attached membrane nodes (u_fb: also prescribed motions, X 11 = 2)
    pnetinter->membrane_body_nhflow(c_, quatRotMat, Eigen::Vector3d(u_fb(0),u_fb(1),u_fb(2)), Eigen::Vector3d(u_fb(3),u_fb(4),u_fb(5)));
    
    pnetinter->membrane_forcing_nhflow(p,d,pgc,alpha,UH,VH,WH,WL);
}

void sixdof_obj::membrane_stabilisation(lexer *p, int iter)
{
    // A membrane bag carries much more water than a floating collar weighs. The loads on the body (membrane
    // load, and the pressure on the hull, computed after the projection of the last stage) contain the inertia
    // of that water as -M_a,true a, lagged by a stage, and the explicit exchange diverges for M_a,true >> M
    // (added-mass instability). In addition the top edge of a flexible membrane acts on the collar as a stiff
    // spring, F = F^n + J q (net_membrane::attach_response), which has to be integrated implicitly together with
    // the added mass. The translation of the stage (h = alpha dt) is therefore integrated as
    //
    //   (M + M_a - h^2 J) dv = h (F + M_a a^n + h J v)
    //
    // and passed to the RK update as the equivalent force M dv/h. M_a a^n (acceleration of the last time step)
    // leaves the equilibrium unchanged; stable for M_a > (M_a,true - M)/2.
    // M_a: membrane.dat 'bodyaddedmass', default 2 rho V_bag (water enclosed below the still water level).
    // Rotations of a rigid bag in the same way with the added inertia I_a (net_membrane::body_addedinertia, 2 rho x
    // inertia of the bag water about the centre of gravity), world frame, small rotations (no gyroscopic coupling):
    //
    //   (I_w + I_a) dw = h (M + I_a alpha^n),   I_w = R I R^T,
    //
    // passed on as the equivalent moment I_w dw/h. A flexible bag needs none (its top edge is coupled implicitly).
    const double Ma = pnetinter->membrane_addedmass_nhflow(p);
    Eigen::Matrix3d J = pnetinter->membrane_stiffness_nhflow(p);
    const Eigen::Matrix3d Ia = pnetinter->membrane_addedinertia_nhflow(p);
    
    if(Ma<=0.0 && J.norm()<=0.0 && Ia.norm()<=0.0)
    return;
    
    // acceleration of the last time step, from the velocity at the start of each step
    if(iter==0 && p->simtime!=tmem_n_)
    {
        const Eigen::Vector3d u(u_fb(0),u_fb(1),u_fb(2));
        const Eigen::Vector3d w(u_fb(3),u_fb(4),u_fb(5));
        
        if(tmem_n_>=0.0 && p->simtime>tmem_n_)
        {
        amem_n_ = (u - umem_n_)/(p->simtime - tmem_n_);
        almem_n_ = (w - wmem_n_)/(p->simtime - tmem_n_);
        }
        
        umem_n_ = u;
        wmem_n_ = w;
        tmem_n_ = p->simtime;
    }
    
    if(Ia.norm()>0.0)
    {
        const int fw[3] = {p->X11_p==1 && p->j_dir==1, p->X11_q==1, p->X11_r==1 && p->j_dir==1};
        const Eigen::Matrix3d Iw = quatRotMat*I_*quatRotMat.transpose();
        
        Eigen::Matrix3d A = Iw + Ia;
        Eigen::Vector3d b = Mfb_ + Ia*almem_n_;
        
        for(int r=0; r<3; ++r)
        if(!fw[r])
        {
            A.row(r).setZero();
            A.col(r).setZero();
            A(r,r) = 1.0;
            b(r) = 0.0;
        }
        
        const Eigen::Vector3d M = Iw*A.lu().solve(b);
        
        for(int r=0; r<3; ++r)
        if(fw[r])
        Mfb_(r) = M(r);
    }
    
    if(Ma<=0.0 && J.norm()<=0.0)
    return;
    
    // free translations only
    const int fr[3] = {p->X11_u==1, p->X11_v==1 && p->j_dir==1, p->X11_w==1};
    
    for(int r=0; r<3; ++r)
    for(int c=0; c<3; ++c)
    if(!fr[r] || !fr[c])
    J(r,c) = 0.0;
    
    const double h = alpha[iter]*p->dt;
    const Eigen::Vector3d v(u_fb(0),u_fb(1),u_fb(2));
    
    Eigen::Matrix3d A = (Mass_fb + Ma)*Eigen::Matrix3d::Identity() - h*h*J;
    Eigen::Vector3d b = Ffb_ + Ma*amem_n_ + h*J*v;
    
    for(int r=0; r<3; ++r)
    if(!fr[r])
    {
        A.row(r).setZero();
        A.col(r).setZero();
        A(r,r) = 1.0;
        b(r) = 0.0;
    }
    
    const Eigen::Vector3d F = Mass_fb*A.lu().solve(b);
    
    for(int r=0; r<3; ++r)
    if(fr[r])
    Ffb_(r) = F(r);
}

bool sixdof_obj::membrane_iterated()
{
    return pnetinter->membrane_iterated();
}

void sixdof_obj::membrane_reforce_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, double *UH, double *VH, double *WH, slice &WL)
{
    pnetinter->membrane_reforce_nhflow(p,d,pgc,alpha,UH,VH,WH,WL);
}

bool sixdof_obj::membrane_couple_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, int iter, double alpha, slice &WL, int it)
{
    return pnetinter->membrane_couple_nhflow(p,d,pgc,iter,alpha,WL,it);
}

void sixdof_obj::membrane_reaction_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, slice &WL, bool finalize)
{
    pnetinter->membrane_reaction_nhflow(p,d,pgc,alpha,WL,finalize);
}

void sixdof_obj::mooringForces(lexer *p, ghostcell *pgc, double alpha)
{
	for (int ii=0; ii<p->mooring_count; ii++)
	{
		// Update coordinates of end point
        Eigen::Vector3d point(X311_xen[ii], X311_yen[ii], X311_zen[ii]);
					
        point = R_*point;
					
        p->X311_xe[ii] = point(0) + c_(0);
        p->X311_ye[ii] = point(1) + c_(1);
        p->X311_ze[ii] = point(2) + c_(2);

        // Advance in time
        pmooring[ii]->start(p, pgc);
                
        // Get forces
        pmooring[ii]->mooringForces(Xme[ii],Yme[ii],Zme[ii]);
                
        // Calculate moments
        Kme[ii] = (p->X311_ye[ii] - c_(1))*Zme[ii] - (p->X311_ze[ii] - c_(2))*Yme[ii];
        Mme[ii] = (p->X311_ze[ii] - c_(2))*Xme[ii] - (p->X311_xe[ii] - c_(0))*Zme[ii];
        Nme[ii] = (p->X311_xe[ii] - c_(0))*Yme[ii] - (p->X311_ye[ii] - c_(1))*Xme[ii];
            
        // Distribute forces and moments to all processors
        pgc->bcast_double(&Xme[ii],1);
        pgc->bcast_double(&Yme[ii],1);
        pgc->bcast_double(&Zme[ii],1);
        pgc->bcast_double(&Kme[ii],1);
        pgc->bcast_double(&Mme[ii],1);
        pgc->bcast_double(&Nme[ii],1);	
        
        // Add to external forces
        Xext += Xme[ii];
        Yext += Yme[ii];
        Zext += Zme[ii];
        
        Kext += Kme[ii];
        Mext += Mme[ii];
        Next += Nme[ii];
    }
}

void sixdof_obj::netForces_cfd(lexer *p, fdm* a, ghostcell *pgc, double alpha, bool finalize)
{    
    pnetinter->netForces_cfd(p,a,pgc,alpha,quatRotMat,Xne,Yne,Zne,Kne,Mne,Nne,finalize);
    
    NETLOOP
    {
    // Add to external forces
        Xext += Xne[n];
        Yext += Yne[n];
        Zext += Zne[n];
        Kext += Kne[n];
        Mext += Mne[n];
        Next += Nne[n];
    }
}

void sixdof_obj::netForces_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double alpha, bool finalize)
{
    pnetinter->netForces_nhflow(p,d,pgc,alpha,quatRotMat,Xne,Yne,Zne,Kne,Mne,Nne,finalize);
    
    NETLOOP
    {
    // Add to external forces
        Xext += Xne[n];
        Yext += Yne[n];
        Zext += Zne[n];
        Kext += Kne[n];
        Mext += Mne[n];
        Next += Nne[n];
    }
}	

