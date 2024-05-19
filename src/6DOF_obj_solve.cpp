/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
Authors: Tobias Martin, Hans Bihs
--------------------------------------------------------------------*/

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_obj::solve_eqmotion(lexer *p, fdm *a, ghostcell *pgc, int iter, vrans *pvrans, vector<net*>& pnet)
{
    externalForces(p, a, pgc, alpha[0], pvrans, pnet);
    
    update_forces(p);
    
    if(p->N40==2 || p->N40==22)
    rk2(p,pgc,iter);
    
    if(p->N40==3 || p->N40==23 || p->N40==33)
    rk3(p,pgc,iter);
   
    if(p->N40==4)
    rkls3(p,pgc,iter);
}

void sixdof_obj::rkls3(lexer *p, ghostcell *pgc, int iter)
{
    get_trans(p,pgc, dp_, dc_, p_, c_);    
    get_rot(p,dh_, de_, h_, e_);

    p_ = p_ + gamma[iter]*p->dt*dp_ + zeta[iter]*p->dt*pk_;
    c_ = c_ + gamma[iter]*p->dt*dc_ + zeta[iter]*p->dt*ck_;
    h_ = h_ + gamma[iter]*p->dt*dh_ + zeta[iter]*p->dt*hk_;
    e_ = e_ + gamma[iter]*p->dt*de_ + zeta[iter]*p->dt*ek_;
    
    pk_ = dp_;
    ck_ = dc_;
    hk_ = dh_;
    ek_ = de_;
}
    
void sixdof_obj::rk3(lexer *p, ghostcell *pgc, int iter)
{   
    if(iter==0)
    {
        pk_ = p_;
        ck_ = c_;
        hk_ = h_;
        ek_ = e_;
        
        get_trans(p,pgc, dp_, dc_, p_, c_);    
        get_rot(p,dh_, de_, h_, e_);
        
        p_ = p_ + p->dt*dpn1_;
        c_ = c_ + p->dt*dcn1_;
        h_ = h_ + p->dt*dhn1_;
        e_ = e_ + p->dt*den1_;
    }
    
    if(iter==1)
    {
        get_trans(p,pgc, dp_, dc_, p_, c_);    
        get_rot(p,dh_, de_, h_, e_);
    
        p_ = 0.75*pk_ + 0.25*p_ + 0.25*p->dt*dpn1_;
        c_ = 0.75*ck_ + 0.25*c_ + 0.25*p->dt*dcn1_;
        h_ = 0.75*hk_ + 0.25*h_ + 0.25*p->dt*dhn1_;
        e_ = 0.75*ek_ + 0.25*e_ + 0.25*p->dt*den1_;    
    }  
    
    else
    {
        get_trans(p,pgc, dp_, dc_, p_, c_);    
        get_rot(p,dh_, de_, h_, e_);
        
        p_ = (1.0/3.0)*pk_ + (2.0/3.0)*p_ + (2.0/3.0)*p->dt*dpn1_;
        c_ = (1.0/3.0)*ck_ + (2.0/3.0)*c_ + (2.0/3.0)*p->dt*dcn1_;
        h_ = (1.0/3.0)*hk_ + (2.0/3.0)*h_ + (2.0/3.0)*p->dt*dhn1_;
        e_ = (1.0/3.0)*ek_ + (2.0/3.0)*e_ + (2.0/3.0)*p->dt*den1_;         
    }
}

void sixdof_obj::rk2(lexer *p, ghostcell *pgc, int iter)
{   
    get_trans(p,pgc, dp_, dc_, p_, c_);    
    get_rot(p,dh_, de_, h_, e_);
        
    if(iter==0)
    {
        pk_ = p_;
        ck_ = c_;
        hk_ = h_;
        ek_ = e_;

        p_ = p_ + p->dt*dpn1_;
        c_ = c_ + p->dt*dcn1_;
        h_ = h_ + p->dt*dhn1_;
        e_ = e_ + p->dt*den1_;
    }
    
    else
    {  
        p_ = 0.5*pk_ + 0.5*p_ + 0.5*p->dt*dpn1_;
        c_ = 0.5*ck_ + 0.5*c_ + 0.5*p->dt*dcn1_;
        h_ = 0.5*hk_ + 0.5*h_ + 0.5*p->dt*dhn1_;
        e_ = 0.5*ek_ + 0.5*e_ + 0.5*p->dt*den1_;         
    }
}

void sixdof_obj::solve_eqmotion_oneway(lexer *p, ghostcell *pgc)
{
    get_trans(p,pgc, dp_, dc_, p_, c_);    
    get_rot(p,dh_, de_, h_, e_);
        
        pk_ = p_;
        ck_ = c_;
        hk_ = h_;
        ek_ = e_;

        p_ = p_ + p->dt*dp_;
        c_ = c_ + p->dt*dc_;
        h_ = h_ + p->dt*dh_;
        e_ = e_ + p->dt*de_;

        p_ = 0.5*pk_ + 0.5*p_ + 0.5*p->dt*dp_;
        c_ = 0.5*ck_ + 0.5*c_ + 0.5*p->dt*dc_;
        h_ = 0.5*hk_ + 0.5*h_ + 0.5*p->dt*dh_;
        e_ = 0.5*ek_ + 0.5*e_ + 0.5*p->dt*de_;         
}


