/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include"ptf_fsfbc.h"
#include"lexer.h"
#include"fdm_ptf.h"
#include"ghostcell.h"
#include"solver2D.h"

void ptf_fsfbc::damping(lexer *p, fdm_ptf *a, ghostcell *pgc, slice &f, int gcval, double alpha)
{
    double starttime=pgc->timer();
    
    int check=0;
    
    // Fifsf damping
    if((gcval==60 || gcval==160) && (p->A357==1 || p->A357==2))
    check=1;
    
    // eta damping
    if((gcval==55 || gcval==155) && (p->A357==1 || p->A357==3))
    check=1;
    
    if(p->A350==1 && check==1)
    {
        n=0;
        SLICELOOP4
        {
         visc = a->vb(i,j);
            
         a->N.p[n] =   visc/(p->DXP[IM1]*p->DXN[IP])*p->x_dir
                     + visc/(p->DXP[IP]*p->DXN[IP])*p->x_dir
                     + visc/(p->DYP[JM1]*p->DYN[JP])*p->y_dir
                     + visc/(p->DYP[JP]*p->DYN[JP])*p->y_dir
                       
                       + 1.0/(alpha*p->dt);
        
         a->rvec.V[n] =   f(i,j)/(alpha*p->dt);
         
         a->N.s[n] = -visc/(p->DXP[IM1]*p->DXN[IP])*p->x_dir;
         a->N.n[n] = -visc/(p->DXP[IP]*p->DXN[IP])*p->x_dir;
         
         a->N.e[n] = -visc/(p->DYP[JM1]*p->DYN[JP])*p->y_dir;
         a->N.w[n] = -visc/(p->DYP[JP]*p->DYN[JP])*p->y_dir;
     
         ++n;
        }
        
        
        n=0;
        SLICELOOP4
        {
            if(p->flagslice4[Im1J]<0)
            {
            a->rvec.V[n] -= a->N.s[n]*f(i-1,j);
            a->N.s[n] = 0.0;
            }
            
            if(p->flagslice4[Ip1J]<0)
            {
            a->rvec.V[n] -= a->N.n[n]*f(i+1,j);
            a->N.n[n] = 0.0;
            }
            
            if(p->flagslice4[IJm1]<0)
            {
            a->rvec.V[n] -= a->N.e[n]*f(i,j-1);
            a->N.e[n] = 0.0;
            }
            
            if(p->flagslice4[IJp1]<0)
            {
            a->rvec.V[n] -= a->N.w[n]*f(i,j+1);
            a->N.w[n] = 0.0;
            }
     
        ++n;
        }
        
        pgc->gcsl_start4(p,f,gcval);
        
        psolv->start(p,pgc,f,a->N,a->xvec,a->rvec,4);
        

        double time=pgc->timer()-starttime;
        
        if(p->mpirank==0 && p->count%p->P12==0 && p->D21==1)
        cout<<"fsfbc_damping: "<<p->solveriter<<"  fsfbc_damping_time: "<<setprecision(3)<<time<<endl;
    }
}

void ptf_fsfbc::damping_wd(lexer *p, fdm_ptf *a, ghostcell *pgc, slice &f, int gcval, double alpha)
{
    double starttime=pgc->timer();
    
    int check=0;
    
    // Fifsf damping
    if((gcval==60 || gcval==160) && (p->A357==1 || p->A357==2))
    check=1;
    
    // eta damping
    if((gcval==55 || gcval==155) && (p->A357==1 || p->A357==3))
    check=1;
    
    if(p->A350==1)
    {
        n=0;
        SLICELOOP4
        {
            if(p->wet[IJ]==1 || p->A343==2)
            {
             visc = a->vb(i,j);
                
             a->N.p[n] =   visc/(p->DXP[IM1]*p->DXN[IP])*p->x_dir
                         + visc/(p->DXP[IP]*p->DXN[IP])*p->x_dir
                         + visc/(p->DYP[JM1]*p->DYN[JP])*p->y_dir
                         + visc/(p->DYP[JP]*p->DYN[JP])*p->y_dir
                           
                           + 1.0/(alpha*p->dt);
            
             a->rvec.V[n] =   f(i,j)/(alpha*p->dt);
             
             a->N.s[n] = -visc/(p->DXP[IM1]*p->DXN[IP])*p->x_dir;
             a->N.n[n] = -visc/(p->DXP[IP]*p->DXN[IP])*p->x_dir;
             
             a->N.e[n] = -visc/(p->DYP[JM1]*p->DYN[JP])*p->y_dir;
             a->N.w[n] = -visc/(p->DYP[JP]*p->DYN[JP])*p->y_dir;
            }

            if(p->wet[IJ]==0 && p->A343==1)
            {
             a->N.p[n] =  1.0;
            
             a->rvec.V[n] =   0.0;
             
             a->N.s[n] = 0.0;
             a->N.n[n] = 0.0;
             
             a->N.e[n] = 0.0;
             a->N.w[n] = 0.0;
            }
     
         ++n;
        }
        
        
        n=0;
        SLICELOOP4
        {
            if(p->wet[IJ]==1 || p->A343==2)
            {
                if(p->flagslice4[Im1J]<0 || (p->wet[Im1J]==0 && p->A343==1))
                {
                a->rvec.V[n] -= a->N.s[n]*f(i,j);
                a->N.s[n] = 0.0;
                }
                
                if(p->flagslice4[Ip1J]<0 || (p->wet[Ip1J]==0 && p->A343==1))
                {
                a->rvec.V[n] -= a->N.n[n]*f(i,j);
                a->N.n[n] = 0.0;
                }
                
                if(p->flagslice4[IJm1]<0 || (p->wet[IJm1]==0 && p->A343==1))
                {
                a->rvec.V[n] -= a->N.e[n]*f(i,j);
                a->N.e[n] = 0.0;
                }
                
                if(p->flagslice4[IJp1]<0 || (p->wet[IJp1]==0 && p->A343==1))
                {
                a->rvec.V[n] -= a->N.w[n]*f(i,j);
                a->N.w[n] = 0.0;
                }
            }
     
        ++n;
        }
        
        pgc->gcsl_start4(p,f,gcval);
        
        psolv->start(p,pgc,f,a->N,a->xvec,a->rvec,4);
        
        double time=pgc->timer()-starttime;
        
        if(p->mpirank==0 && p->count%p->P12==0 && p->D21==1)
        cout<<"fsfbc_damping: "<<p->solveriter<<"  fsfbc_damping_time: "<<setprecision(3)<<time<<endl;
    }
}
