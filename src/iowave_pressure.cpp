/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"patchBC_interface.h"

void iowave::pressure_io(lexer *p, fdm* a, ghostcell *pgc)
{
    pressure_inlet(p,a,pgc);
    pressure_outlet(p,a,pgc);
    
    pBC->patchBC_pressure(p,a,pgc,a->press);
}

void iowave::pressure_outlet(lexer *p, fdm *a, ghostcell *pgc)
{
    double pval=0.0;

        for(n=0;n<p->gcout_count;++n)
        {
        i=p->gcout[n][0];
        j=p->gcout[n][1];
        k=p->gcout[n][2];
		pval=0.0;
		
			if(p->B77==1 && p->B99==0)
			{
                if(p->F50==2 || p->F50==3)
                pval=(p->fsfout - p->pos_z())*a->ro(i,j,k)*fabs(p->W22);
                
                if(p->F50==1 || p->F50==4)
                pval=a->press(i,j,k);
            
			a->press(i+1,j,k)=pval;
			a->press(i+2,j,k)=pval;
			a->press(i+3,j,k)=pval;
			}
            
            if(p->B77==1 && (p->B99==1 || p->B99==2))
			{
                pval=a->press(i,j,k);
            
			a->press(i+1,j,k)=pval;
			a->press(i+2,j,k)=pval;
			a->press(i+3,j,k)=pval;
			}
            
            if(p->B77==0 && p->A10==5)
			{
			pval=0.0;
			
			a->press(i+1,j,k)=pval;
			a->press(i+2,j,k)=pval;
			a->press(i+3,j,k)=pval;
			}
		
			if(p->B77==2)
			{
			double eps,H;
                
            eps = 0.6*(1.0/3.0)*(p->DXN[IP] + p->DYN[JP] + p->DZN[KP]);
        
            if(a->phi(i,j,k)>eps)
            H=1.0;

            if(a->phi(i,j,k)<-eps)
            H=0.0;

            if(fabs(a->phi(i,j,k))<=eps)
            H=0.5*(1.0 + a->phi(i,j,k)/eps + (1.0/PI)*sin((PI*a->phi(i,j,k))/eps));
        
        
    
            pval=(1.0-H)*a->press(i,j,k);
			

			a->press(i+1,j,k)=pval;
			a->press(i+2,j,k)=pval;
			a->press(i+3,j,k)=pval;
			}
			
        }
}

void iowave::pressure_inlet(lexer *p, fdm *a, ghostcell *pgc)
{
    double pval=0.0;
    
    if(p->B76==0 && p->A10 != 5)
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
    
    if(p->B76==0 && p->A10==5)
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];
		
        pval = 0.0;

        a->press(i-1,j,k)=pval;
        a->press(i-2,j,k)=pval;
        a->press(i-3,j,k)=pval;
    }
}

void iowave::pressure_wall(lexer *p, fdm *a, ghostcell *pgc)
{
    double pval=0.0;

    GC4LOOP
    if(p->gcb4[n][3] != 5 && p->gcb4[n][3] != 6 && (p->gcb4[n][4] ==3 && p->gcb4[n][4] ==21 && p->gcb4[n][4] ==22))
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];


        if(a->phi(i,j,k)>0.0 || p->I56==0)
        {
        pval=(p->phimean - p->pos_z())*a->ro(i,j,k)*fabs(p->W22);

            if(p->gcb4[n][3]==1)
            {
            a->press(i-1,j,k)=pval;
            a->press(i-2,j,k)=pval;
            a->press(i-3,j,k)=pval;
            }

            if(p->gcb4[n][3]==2)
            {
            a->press(i,j+1,k)=pval;
            a->press(i,j+2,k)=pval;
            a->press(i,j+3,k)=pval;
            }

            if(p->gcb4[n][3]==3)
            {
            a->press(i,j-1,k)=pval;
            a->press(i,j-2,k)=pval;
            a->press(i,j-3,k)=pval;
            }

            if(p->gcb4[n][3]==4)
            {
            a->press(i+1,j,k)=pval;
            a->press(i+2,j,k)=pval;
            a->press(i+3,j,k)=pval;
            }
        }
    }
}

void iowave::pressure_bed(lexer *p, fdm *a, ghostcell *pgc)
{
    double pval=0.0;

    GC4LOOP
    if(p->gcb4[n][3] == 5 )
    {
    i=p->gcb4[n][0];
    j=p->gcb4[n][1];
    k=p->gcb4[n][2];


        if(a->phi(i,j,k)>0.0 || p->I56==0)
        {
        pval=(local_fsf(p,a,pgc) - p->pos_z())*a->ro(i,j,k)*fabs(p->W22);

        a->press(i,j,k-1)=pval;
        a->press(i,j,k-2)=pval;
        a->press(i,j,k-3)=pval;
        }
    }
}

double iowave::local_fsf(lexer *p, fdm *a, ghostcell *pgc)
{
    double wsf=-1.0e20;
    int count=0;

        KLOOP
        PCHECK
        {
            if(a->phi(i,j,k)>=0.0 && a->phi(i,j,k+1)<0.0)
            wsf=MAX(wsf,-(a->phi(i,j,k)*p->DXM)/(a->phi(i,j,k+1)-a->phi(i,j,k)) + p->pos_z());

            if(a->phi(i,j,k)>0.0)
            ++count;
        }

    //wsf=pgc->globalmax(wsf);

    if(count==0)
    wsf=0.0;

    return wsf;
}

