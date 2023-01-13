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

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include"patchBC_interface.h"

void iowave::inflow(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{
    if(p->I230==0)
    {
    if(p->B98==0 || (p->W10>0.0 && p->B98==2))
    inflow_plain(p,a,pgc,u,v,w);

    
	if(p->B98==3)
	dirichlet_wavegen(p,a,pgc,u,v,w);
	
	if(p->B98==4)
	active_wavegen(p,a,pgc,u,v,w);
	}
    
	if(p->B99==3||p->B99==4||p->B99==5)
	active_beach(p,a,pgc,u,v,w);
    
    if(p->I230>0)
    ff_inflow(p,a,pgc,u,v,w);
    
    pBC->patchBC_ioflow(p,a,pgc,u,v,w);
}

void iowave::rkinflow(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];
    
    u(i-1,j,k) = u(i-2,j,k) = u(i-3,j,k) = a->u(i-1,j,k);
    v(i-1,j,k) = v(i-2,j,k) = v(i-3,j,k) = a->v(i-1,j,k);
    w(i-1,j,k) = w(i-2,j,k) = w(i-3,j,k) = a->w(i-1,j,k);
    }
    
    pBC->patchBC_rkioflow(p,a,pgc,u,v,w);
}

void iowave::inflow_plain(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{
    for(n=0;n<p->gcin_count;n++)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];

        u(i-1,j,k)=p->Ui;
        u(i-2,j,k)=p->Ui;
        u(i-3,j,k)=p->Ui;
		
		v(i-1,j,k)=0.0;
        v(i-2,j,k)=0.0;
        v(i-3,j,k)=0.0;
		
		w(i-1,j,k)=0.0;
        w(i-2,j,k)=0.0;
        w(i-3,j,k)=0.0;
        
        // Air inflow
        if(p->W50_air==1 && a->phi(i,j,k)<-0.6*p->DXM)
        {
        u(i-1,j,k)+=p->W50;
        u(i-2,j,k)+=p->W50;
        u(i-3,j,k)+=p->W50;
        }
    
    }
}

void iowave::flowfile(lexer *p, fdm* a, ghostcell* pgc, turbulence *pturb)
{
    if(p->gcin_count>0 && p->I230>0)
    flowfile_start(p,a,pgc,pturb);
}
