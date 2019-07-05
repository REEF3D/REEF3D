/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"fnpf_sg_fsfbc_wd.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"solver2D.h"

void fnpf_sg_fsfbc_wd::damping(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &f, int gcval, double alpha)
{
    double starttime=pgc->timer();
    double visc; 
    
    if(p->A350==1)
    {
        n=0;
        SLICELOOP4
        {
        visc = c->vb(i,j);

            
         c->N.p[n] =   0.5*(c->vb(i,j)+c->vb(i-1,j))/(p->DXP[IM1]*p->DXN[IP])
                     + 0.5*(c->vb(i,j)+c->vb(i+1,j))/(p->DXP[IP]*p->DXN[IP])
                     + 0.5*(c->vb(i,j)+c->vb(i,j-1))/(p->DYP[JM1]*p->DYN[JP])
                     + 0.5*(c->vb(i,j)+c->vb(i,j+1))/(p->DYP[JP]*p->DYN[JP])
                       
                       + 1.0/(alpha*p->dt);
        
         c->rvec.V[n] =    c->N.p[n]*f(i,j)*(1.0/p->N54-1.0)
                             
                             + (f(i,j))/(alpha*p->dt);

         c->N.p[n] /= p->N54;
         
         c->N.s[n] = -0.5*(c->vb(i,j)+c->vb(i-1,j))/(p->DXP[IM1]*p->DXN[IP]);
         c->N.n[n] = -0.5*(c->vb(i,j)+c->vb(i+1,j))/(p->DXP[IP]*p->DXN[IP]);
         
         c->N.e[n] = -0.5*(c->vb(i,j)+c->vb(i,j-1))/(p->DYP[JM1]*p->DYN[JP]);
         c->N.w[n] = -0.5*(c->vb(i,j)+c->vb(i,j+1))/(p->DYP[JP]*p->DYN[JP]);
     
         ++n;
        }
        
        
        n=0;
        SLICELOOP4
        {
            if(p->flagslice4[Im1J]<0)
            {
            c->rvec.V[n] -= c->N.s[n]*f(i-1,j);
            c->N.s[n] = 0.0;
            }
            
            if(p->flagslice4[Ip1J]<0)
            {
            c->rvec.V[n] -= c->N.n[n]*f(i+1,j);
            c->N.n[n] = 0.0;
            }
            
            if(p->flagslice4[IJm1]<0)
            {
            c->rvec.V[n] -= c->N.e[n]*f(i,j-1);
            c->N.e[n] = 0.0;
            }
            
            if(p->flagslice4[IJp1]<0)
            {
            c->rvec.V[n] -= c->N.w[n]*f(i,j+1);
            c->N.w[n] = 0.0;
            }
     
        ++n;
        }
        
        psolv->start(p,pgc,f,c->N,c->xvec,c->rvec,4,50,p->D29,c->C4);
        
        pgc->gcsl_start4(p,f,gcval);
        
        double time=pgc->timer()-starttime;
        if(p->mpirank==0 && innercounter==p->N50-1 && p->D21==1 && (p->count%p->P12==0))
        cout<<"fsfbc_damping: "<<p->solveriter<<"  fsfbc_damping_time: "<<setprecision(3)<<time<<endl;
    }
}