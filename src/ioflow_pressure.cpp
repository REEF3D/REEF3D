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

#include"ioflow_f.h"
#include"lexer.h"
#include"fdm.h"
#include"patchBC_interface.h"
#include"heaviside.h"

void ioflow_f::pressure_io(lexer *p, fdm* a, ghostcell *pgc)
{
    pressure_inlet(p,a,pgc);
    pressure_outlet(p,a,pgc);
    
    
    pBC->patchBC_pressure(p,a,pgc,a->press);
}

void ioflow_f::pressure_inlet(lexer *p, fdm *a, ghostcell *pgc)
{
    double pval=0.0;

    if(p->B76==0)
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];
		
		if(a->phi(i,j,k)>=0.0)
        pval=(p->phimean - p->pos_z())*a->ro(i,j,k)*fabs(p->W22);
		
		if(a->phi(i,j,k)<0.0)
        pval = a->press(i,j,k);

        a->press(i-1,j,k)=pval;
        a->press(i-2,j,k)=pval;
        a->press(i-3,j,k)=pval;
    }
    
    if(p->B76==3)
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];
    
		
		if(a->phi(i,j,k)>=0.0)
        pval=a->press(i,j,k) + p->Ui*p->DXP[IM1]; 
		
		if(a->phi(i,j,k)<0.0)
        pval = a->press(i,j,k);
    
        a->press(i-1,j,k)=pval;
        a->press(i-2,j,k)=pval;
        a->press(i-3,j,k)=pval;
    }
}


void ioflow_f::pressure_outlet(lexer *p, fdm *a, ghostcell *pgc)
{
    double pval=0.0;
    double diff;
    double eps,H,roval;
    
    /*
    if(p->count!=iter0)
    {
    diff = p->phiout-p->fsfout;
    
    p->fsfoutval -= 0.1*diff;
    
    iter0=p->count;
    
    // cout<<p->mpirank<<" fsfout: "<<p->fsfout<<" diff: "<<diff<<" fsfoutval: "<<p->fsfoutval<<" phiout: "<<p->phiout<<endl;
    }*/
    
    
        for(n=0;n<p->gcout_count;++n)
        {
        i=p->gcout[n][0];
        j=p->gcout[n][1];
        k=p->gcout[n][2];
        pval=0.0;
        
        
        if(p->B77==0)
        {
        pval = a->press(i,j,k); 
        a->press(i+1,j,k)=pval;
        a->press(i+2,j,k)=pval;
        a->press(i+3,j,k)=pval;
        }
		
        
			if(p->B77==1)
			{
                
                
            /*
            eps = 2.1*(1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);
        
            H = heaviside(a->phi(i,j,k),eps);
            
            //pval=H*pval + (1.0-H)*a->press(i,j,k);
            
            roval = p->W1*H +   p->W3*(1.0-H);*/
            
            
                if(p->F50==2 || p->F50==3)
                pval=(p->fsfout - p->pos_z())*a->ro(i,j,k)*fabs(p->W22);
                
                if(p->F50==1 || p->F50==4)
                pval=a->press(i,j,k);
                
                pval=a->press(i,j,k);
            
			a->press(i+1,j,k)=pval;
			a->press(i+2,j,k)=pval;
			a->press(i+3,j,k)=pval;
			}
            
            if(p->B77==2)
			{
                pval=a->press(i,j,k);
            
			a->press(i+1,j,k)=pval;
			a->press(i+2,j,k)=pval;
			a->press(i+3,j,k)=pval;
			}
		
        
        
			if(p->B77==10)
			{
            eps = 0.6*(1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);
        
            H = heaviside(a->phi(i,j,k),eps);
        
            pval=(1.0-H)*a->press(i,j,k);
            
			a->press(i+1,j,k)=pval;
			a->press(i+2,j,k)=pval;
			a->press(i+3,j,k)=pval;
			}
        }
}
