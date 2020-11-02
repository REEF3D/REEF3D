/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
--------------------------------------------------------------------*/

#include"6DOF_df.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"


void sixdof_df::predictor(lexer *p, fdm *a, ghostcell *pgc, double alpha)
{
}
    

bool sixdof_df::corrector(lexer *p, fdm *a, ghostcell *pgc, double alpha, vrans *pvrans, vector<net*>& pnet)
{
    double err_norm = 0.0;
    
    double Fxold = Ffb_(0); 
    double Fyold = Ffb_(1); 
    double Fzold = Ffb_(2); 
    double Mxold = Mfb_(0); 
    double Myold = Mfb_(1); 
    double Mzold = Mfb_(2); 

    p_ = pn1_;
    c_ = cn1_;
    h_ = hn1_;
    e_ = en1_;


    if (alpha == 1.0)
    {
        externalForces(p, a, pgc, alpha, pvrans, pnet);
    }
    
    updateForces(a);
    
    rk3(p,a,pgc,alpha);

    err_norm = 
          fabs(Fxold - Ffb_(0)) + fabs(Fyold - Ffb_(1)) + fabs(Fzold - Ffb_(2)) 
        + fabs(Mxold - Mfb_(0)) + fabs(Myold - Mfb_(1)) + fabs(Mzold - Mfb_(2)); 

    if (err_norm < 1e-7) 
    {
        if (p->mpirank == 0)
        {
            //cout<<setprecision(10)<<"Converged FSI after "<<nCorr<<" iterations with error: "<<err_norm<<endl;
        }

        return true;
    }
    else
    {
        if (p->mpirank == 0)
        {
            //cout<<setprecision(10)<<"Iteration loop "<<nCorr<<" with error: "<<err_norm<<endl;
        }
        
        nCorr++;
        
        return false;
    }
}




void sixdof_df::get_trans
(
    lexer *p,
    fdm *a,
    ghostcell *pgc,
    Eigen::Vector3d& dp,
    Eigen::Vector3d& dc, 
    const Eigen::Vector3d& pp,
    const Eigen::Vector3d& c
)
{
    dp = Ffb_; 
    dc = pp/Mass_fb;

	// Prescribed motions
	prescribedMotion(p,a,pgc,dp,dc);
} 


void sixdof_df::get_rot
(
    Eigen::Vector3d& dh, 
    Eigen::Vector4d& de, 
    const Eigen::Vector3d& h, 
    const Eigen::Vector4d& e
)
{
    // Update Euler parameter matrices
    quat_matrices(e);

    // RHS of e
    de = 0.5*G_.transpose()*I_.inverse()*h;
    
    // RHS of h
    // Transforming torsion into body fixed system (Shivarama and Schwab)
    Gdot_ << -de(1), de(0), de(3), -de(2),
            -de(2), -de(3), de(0), de(1),
            -de(3), de(2), -de(1), de(0); 
   
    dh = 2.0*Gdot_*G_.transpose()*h + Rinv_*Mfb_;
} 


void sixdof_df::prescribedMotion(lexer *p, fdm *a, ghostcell *pgc, Eigen::Vector3d& dp, Eigen::Vector3d& dc)
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





