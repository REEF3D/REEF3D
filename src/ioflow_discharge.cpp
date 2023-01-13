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

#include"ioflow_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"patchBC_interface.h"

void ioflow_f::discharge(lexer *p, fdm* a, ghostcell* pgc)
{
    Qin(p,a,pgc);
    Qout(p,a,pgc);

	if(p->count==0)
    if(p->mpirank==0 && (p->count%p->P12==0))
    {
    cout<<"Inflow_0:  "<<setprecision(5)<<p->W10<<" Ui: "<<p->Ui<<" Hi: "<<p->phimean<<endl;
    cout<<"Outflow_0: "<<setprecision(5)<<p->W10<<" Uo: "<<p->Uo<<" Ho: "<<p->phiout<<endl;
    }
	
	if(p->count>0)
	if(p->mpirank==0 && (p->count%p->P12==0))
    {
    cout<<"Inflow:  "<<setprecision(5)<<p->Qi<<" Ui: "<<p->Ua<<" Hi: "<<p->phimean<<endl;
    cout<<"Outflow: "<<setprecision(5)<<p->Qo<<" Uo: "<<p->Uo<<" Ho: "<<p->phiout<<endl;
    }
    
    // patchBC
    pBC->patchBC_discharge(p,a,pgc);
}

void ioflow_f::Qin(lexer *p, fdm* a, ghostcell* pgc)
{
    area=0.0;
    Ai=0.0;
    p->Qi=0.0;
    p->Ui=0.0;

    // in
    for(n=0;n<p->gcin_count;n++)
    if(p->gcin[n][3]>0)
    {
    area=0.0;
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];

        if(a->phi(i,j,k)>-0.5*p->DZN[KP]-1.0e-20 && a->topo(i,j,k)>0.0)
        {
            if(a->phi(i,j,k)>=0.5*p->DZN[KP])
            area=p->DYN[JP]*p->DZN[KP];

            if(a->phi(i,j,k)<0.5*p->DZN[KP] && a->phi(i,j,k)>0.0)
            area=p->DYN[JP]*(p->DZN[KP]*0.5 + a->phi(i,j,k));
			
			if(a->phi(i,j,k)>=-0.5*p->DZN[KP] -1.0e-20 && a->phi(i,j,k)<=0.0)
            area=p->DYN[JP]*(p->DZN[KP]*0.5 - fabs(a->phi(i,j,k)));


            Ai+=area;
            p->Qi+=area*a->u(i-1,j,k);
        }
    }
    Ai=pgc->globalsum(Ai);
    p->Qi=pgc->globalsum(p->Qi);
    
    if(p->B60==1)
    p->Ui=p->W10/(Ai>1.0e-20?Ai:1.0e20); 
    
    if(p->B60==2 || p->B60==4)
    p->Ui=hydrograph_ipol(p,pgc,hydro_in,hydro_in_count)/(Ai>1.0e-20?Ai:1.0e20);    
    
    if(p->mpirank==0 && (p->B60==2 || p->B60==4))
    cout<<"Qi_ipol: "<<hydrograph_ipol(p,pgc,hydro_in,hydro_in_count)<<endl;
        
    p->Ua=p->Qi/Ai;
	
}

void ioflow_f::Qout(lexer *p, fdm* a, ghostcell* pgc)
{
    area=0.0;
    Ao=0.0;
    p->Qo=0.0;
    p->Uo=0.0;

    // out
    for(n=0;n<p->gcout_count;n++)
    if(p->gcout[n][3]>0)
    {
    area=0.0;
    i=p->gcout[n][0];
    j=p->gcout[n][1];
    k=p->gcout[n][2];

        if(a->phi(i+1,j,k)>-0.5*p->DZN[KP]-1.0e-20 && a->topo(i,j,k)>0.0)
        {

            if(a->phi(i+1,j,k)>=0.5*p->DZN[KP])
            area=p->DYN[JP]*p->DZN[KP];

            if(a->phi(i+1,j,k)<0.5*p->DZN[KP] && a->phi(i+1,j,k)>0.0)
            area=p->DYN[JP]*(p->DZN[KP]*0.5 + a->phi(i+1,j,k));
			
			if(a->phi(i+1,j,k)>=-0.5*p->DZN[KP]-1.0e-20 && a->phi(i+1,j,k)<=0.0)
            area=p->DYN[JP]*(p->DZN[KP]*0.5 - fabs(a->phi(i+1,j,k)));

            Ao+=area;
            p->Qo+=area*a->u(i,j,k);
        }
    }
    
    Ao=pgc->globalsum(Ao);
    p->Qo=pgc->globalsum(p->Qo);
	
	if(p->B60==1)
    {
	p->Uo=p->Qo/(Ao>1.0e-20?Ao:1.0e20);
    
    if(p->count==0 && p->I11==1)
    p->Uo=p->W10/(Ao>1.0e-20?Ao:1.0e20);
    }
	
	if(p->B60==2)
	p->Uo=p->Qo/(Ao>1.0e-20?Ao:1.0e20);
	
	if(p->B60==3 || p->B60==4)
	p->Uo=hydrograph_ipol(p,pgc,hydro_out,hydro_out_count)/(Ao>1.0e-20?Ao:1.0e20); 
	
	if(p->mpirank==0 && (p->B60==3 || p->B60==4))
    cout<<"Qo_ipol: "<<hydrograph_ipol(p,pgc,hydro_out,hydro_out_count)<<endl;
}


