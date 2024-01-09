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

#include"ioflow_f.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"patchBC_interface.h"

void ioflow_f::discharge2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
    Qin2D(p,b,pgc);
    Qout2D(p,b,pgc);

	if(p->count==0)
    if(p->mpirank==0)
    {
    cout<<"Inflow_0:  "<<setprecision(5)<<p->W10<<" Ui: "<<p->Ui<<" Ai: "<<Ai<<" Hi: "<<Hi<<endl;
    cout<<"Outflow_0: "<<setprecision(5)<<p->W10<<" Uo: "<<p->Uo<<" Ao: "<<Ao<<" Ho: "<<Ho<<endl;
    }
	
	if(p->count>0)
	if(p->mpirank==0)
    {
    cout<<"Inflow:  "<<setprecision(5)<<p->Qi<<" Ui: "<<p->Ui<<" Ai: "<<Ai<<" Hi: "<<Hi<<endl;
    cout<<"Outflow: "<<setprecision(5)<<p->Qo<<" Uo: "<<p->Uo<<" Ao: "<<Ao<<" Ho: "<<Ho<<endl;
    }
    
    // patchBC
    pBC->patchBC_discharge2D(p,b,pgc,b->P,b->Q,b->eta,b->bed);
}

void ioflow_f::Qin2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
    area=0.0;
    hval=0.0;
    hcount=0;
    Ai=0.0;
    p->Qi=0.0;
    p->Ui=0.0;

    // in
    for(n=0;n<p->gcslin_count;n++)
    {
    area=0.0;
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];
    
        if(p->wet[IJ]==1)
        {
    
        area = p->DYN[JP]*b->hp(i-1,j);
        
        Ai+=area;
        p->Qi+=area*b->P(i,j);
        
        hval += b->hp(i,j);
        ++hcount;
        }

    }
    
    Ai=pgc->globalsum(Ai);
    p->Qi=pgc->globalsum(p->Qi);
    Hi=pgc->globalsum(hval);
    hcount=pgc->globalisum(hcount);
    
    Hi = Hi/(hcount>1.0e-20?hcount:1.0e20); 
    
    if(p->B60==1)
    p->Ui=p->W10/(Ai>1.0e-20?Ai:1.0e20); 
    
    if(p->B60==2 || p->B60==4)
    p->Ui=hydrograph_ipol(p,pgc,hydro_in,hydro_in_count)/(Ai>1.0e-20?Ai:1.0e20);    
    
    if(p->mpirank==0 && (p->B60==2 || p->B60==4))
    cout<<"Qi_ipol: "<<hydrograph_ipol(p,pgc,hydro_in,hydro_in_count)<<endl;
        
    p->Ua=p->Qi/Ai;
}

void ioflow_f::Qout2D(lexer *p, fdm2D* b, ghostcell* pgc)
{
    area=0.0;
    Ao=0.0;
    hval=0.0;
    hcount=0;
    p->Qo=0.0;
    p->Uo=0.0;

    // out
    for(n=0;n<p->gcslout_count;n++)
    {
    area=0.0;
    i=p->gcslout[n][0];
    j=p->gcslout[n][1];
    
        area = p->DYN[JP]*b->hp(i,j);
        
        Ao+=area;
        p->Qo+=area*b->P(i,j);
        
        hval += b->hp(i,j);
        ++hcount;
    }
    
    Ao=pgc->globalsum(Ao);
    p->Qo=pgc->globalsum(p->Qo);
    Ho=pgc->globalsum(hval);
    hcount=pgc->globalisum(hcount);
    
    Ho = Ho/(hcount>1.0e-20?hcount:1.0e20); 
	
	if(p->B60==1)
	p->Uo=p->W10/(Ao>1.0e-20?Ao:1.0e20);
	
	if(p->B60==2)
	p->Uo=p->Qo/(Ao>1.0e-20?Ao:1.0e20);
	
	if(p->B60==3 || p->B60==4)
	p->Uo=hydrograph_ipol(p,pgc,hydro_out,hydro_out_count)/(Ao>1.0e-20?Ao:1.0e20); 
	
	if(p->mpirank==0 && (p->B60==3 || p->B60==4))
    cout<<"Qo_ipol: "<<hydrograph_ipol(p,pgc,hydro_out,hydro_out_count)<<endl;
}


