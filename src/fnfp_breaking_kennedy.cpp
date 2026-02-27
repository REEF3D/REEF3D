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

#include"fnpf_breaking.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"

void fnpf_breaking::breaking_kennedy(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &eta_n, slice &Fifsf, double alpha)
{
    // ini
    
    // detect breaking
    SLICELOOP4
    eta_t(i, j) = (c->eta(i, j) - eta_n(i, j)) /(alpha* p->dt);

    
     double threshold_I, threshold_F;

    SLICELOOP4
    {
        
        //---------------
        // Shallow Water
        //---------------
        if((p->A351==1 || p->A351==3) && p->count>1)
        {
        // Onset threshold: eta_t > eta_I * sqrt(g * h)
        threshold_I = eta_I * sqrt(9.81 * c->WL(i,j));

        // Cessation threshold: eta_t < eta_F * sqrt(g * h)
        threshold_F = eta_F * sqrt(9.81 * c->WL(i,j));

        // Check for breaking onset
        if (c->breaking(i,j) == 0 && eta_t(i,j) > threshold_I)
        {
            c->breaking(i,j) = 1;
            t_break(i,j)  = p->simtime;
        }
        
        // Check for breaking cessation
        if(c->breaking(i,j) == 1 && eta_t(i,j) < threshold_F)
        {
            c->breaking(i,j) = 0;
            B_coeff(i,j)  = 0.0;
            c->vb(i,j)     = 0.0;
        }
        }
        
        //---------------
        // Deep Water
        //---------------
        // steepness detection
        if((p->A351==2 || p->A351==3) && p->count>1)
        {
        if (c->breaking(i,j) == 0 && (fabs(c->Ex(i,j)) > p->A355 || fabs(c->Ey(i,j)) > p->A355))
        {
            c->breaking(i,j) = 2;
            t_break(i,j)  = p->simtime;
        }
        
        // Check for breaking cessation
        if(c->breaking(i,j) == 2 && (fabs(c->Ex(i,j)) < p->A356*p->A355 || fabs(c->Ey(i,j)) < p->A356*p->A355))
        {
            c->breaking(i,j) = 0;
            B_coeff(i,j)  = 0.0;
            c->vb(i,j)     = 0.0;
        }
        
        
        }
    }
    
    
    // calculate viscosity
    double T_ramp, dt_break;

    SLICELOOP4
    {
        if(c->breaking(i,j)==1 || c->breaking(i,j)==2)
        {
   
            // Ramp-up timescale: T* = T_star_coeff * sqrt(h/g)
            T_ramp = T_star * sqrt(c->WL(i,j)/ 9.81);

            // Time since breaking onset
            dt_break = p->simtime + p->dt - t_break(i,j);

            // Temporal ramp coefficient B(t)
            if (dt_break < T_ramp && T_ramp > 1.0e-10)
            B_coeff(i, j) = dt_break / T_ramp;
            
            else
            B_coeff(i, j) = 1.0;

            B_coeff(i, j) = MIN(B_coeff(i, j), 1.0);
            B_coeff(i, j) = MAX(B_coeff(i, j), 0.0);

            // Eddy viscosity: c->vb = B * delta_b^2 * h * |eta_t|
            // Using |eta_t| to ensure positive viscosity
            //c->vb(i,j) = B_coeff(i,j) * delta_b * delta_b * c->WL(i,j) * fabs(eta_t(i, j));
            c->vb(i,j) = B_coeff(i,j)*p->A365;
            
            c->vb(i-1,j) = B_coeff(i,j)*p->A365;
            c->vb(i+1,j) = B_coeff(i,j)*p->A365;
        
            if(p->j_dir==1)
            {
            c->vb(i,j-1) = B_coeff(i,j)*p->A365;
            c->vb(i,j+1) = B_coeff(i,j)*p->A365;
            }
            
            //cout<<"c->vb(i,j): "<<c->vb(i,j)<<"  B_coeff(i, j): "<<B_coeff(i, j)<<" dt_break: "<<dt_break<<" t_break(i,j): "<<t_break(i,j)<<endl;
        }
        
        /*else
        if(c->breaking(i,j)==2)
        {
        c->vb(i,j) = p->A365;
        c->vb(i-1,j) = p->A365;
        c->vb(i+1,j) = p->A365;
        
        if(p->j_dir==1)
        {
        c->vb(i,j-1) = p->A365;
        c->vb(i,j+1) = p->A365;
        }
        }*/
        
        else
        {
            c->vb(i, j)    = 0.0;
            B_coeff(i, j) = 0.0;
        }
    }

    pgc->gcsl_start4(p, c->vb, 1);
    
    
    SLICELOOP4
    c->breaklog(i,j)=0;
    
    // breaklog
    int count=0; 
    
    SLICELOOP4
    if(c->breaking(i,j)>0)
    {
    c->breaklog(i,j)=1;
    ++count;
    }
    
    /*
    LOOP
    {
    if(c->breaking(i,j)>0)
    c->test[IJK] = 1.0;
    
    else
    c->test[IJK] = 0.0;
    }*/
    
    LOOP
    c->test[IJK] = c->vb(i,j);
    
    count=pgc->globalisum(count);
    
    if(p->mpirank==0 && (p->count%p->P12==0))
    cout<<"breaking: "<<count<<endl;
    
}