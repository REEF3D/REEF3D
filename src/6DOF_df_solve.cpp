/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Tobias Martin
--------------------------------------------------------------------*/

#include"6DOF_df_object.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_df_object::solve_eqmotion(lexer *p, fdm *a, ghostcell *pgc, double alpha, double gamma, double zeta, vrans *pvrans, vector<net*>& pnet)
{
    externalForces(p, a, pgc, alpha, pvrans, pnet);
    
    updateForces(a);
   
    rk3(p,a,pgc,alpha,gamma,zeta);
}

void sixdof_df_object::rk3(lexer *p, fdm *a, ghostcell *pgc, double alpha, double gamma, double zeta)
{
    get_trans(p,a,pgc, dp_, dc_, p_, c_);    
    get_rot(dh_, de_, h_, e_);

    p_ = p_ + gamma*p->dt*dp_ + zeta*p->dt*pk_;
    c_ = c_ + gamma*p->dt*dc_ + zeta*p->dt*ck_;
    h_ = h_ + gamma*p->dt*dh_ + zeta*p->dt*hk_;
    e_ = e_ + gamma*p->dt*de_ + zeta*p->dt*ek_;
    
    pk_ = dp_;
    ck_ = dc_;
    hk_ = dh_;
    ek_ = de_;
}

void sixdof_df_object::get_trans(lexer *p, fdm *a, ghostcell *pgc, Eigen::Vector3d& dp, Eigen::Vector3d& dc, const Eigen::Vector3d& pp, const Eigen::Vector3d& c)
{
    dp = Ffb_; 
    dc = pp/Mass_fb;

	// Prescribed motions
	prescribedMotion(p,a,pgc,dp,dc);
} 

void sixdof_df_object::get_rot(Eigen::Vector3d& dh, Eigen::Vector4d& de, const Eigen::Vector3d& h, const Eigen::Vector4d& e)
{
    // Update Euler parameter matrices
    quat_matrices(e);

    // RHS of e
    de = 0.5*G_.transpose()*I_.inverse()*h;
    
    // RHS of h
    // Transforming torsion into body fixed system (Shivarama and Schwab)
    Gdot_ << -de(1), de(0), de(3),-de(2),
             -de(2),-de(3), de(0), de(1),
             -de(3), de(2),-de(1), de(0); 
   
    dh = 2.0*Gdot_*G_.transpose()*h + Rinv_*Mfb_;
} 

void sixdof_df_object::prescribedMotion(lexer *p, fdm *a, ghostcell *pgc, Eigen::Vector3d& dp, Eigen::Vector3d& dc)
{
    
    if (p->X11_u == 2)
    {
        dp(0) = 0.0; 
        dc(0) = Uext;
    }
    
    if (p->X11_v == 2)
    {
        dp(1) = 0.0; 
        dc(1) = Vext;
    }

    if (p->X11_w == 2)
    {
        dp(2) = 0.0; 
        dc(2) = Wext;
    }
    
    if(p->X11_p==2)
    {
        //Pext;
        cout<<"not implemented yet"<<endl;
    }
                
    if(p->X11_q==2)
    {
        //Qext;
        cout<<"not implemented yet"<<endl;
    }
    
    if(p->X11_r==2)
    {
        //Rext;
        cout<<"not implemented yet"<<endl;
    }
    
}





