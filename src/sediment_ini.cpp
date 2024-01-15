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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"sediment_f.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm2D.h"
#include"ghostcell.h"
#include"sediment_fdm.h"
#include"patchBC_interface.h"

void sediment_f::ini_cfd(lexer *p, fdm *a,ghostcell *pgc)
{
	double h=0;

	ILOOP
    JLOOP
	{
		KLOOP
		PBASECHECK
		{
        if(a->topo(i,j,k-1)<0.0 && a->topo(i,j,k)>0.0)
        h = -(a->topo(i,j,k-1)*p->DZP[KP])/(a->topo(i,j,k)-a->topo(i,j,k-1)) + p->pos_z()-p->DZP[KP];
		}
		
		s->bedzh(i,j)=h;
        s->bedzh0(i,j)=h;
	}
	
	pgc->gcsl_start4(p,s->bedzh,1);
	
    active_ini_cfd(p,a,pgc);
    
    topo_zh_update(p,a,pgc,s);
    
    ini_parameters(p,pgc);
    ini_guard(p,pgc);
    log_ini(p);
}

void sediment_f::ini_sflow(lexer *p, fdm2D *b, ghostcell *pgc)
{
    //relax(p,b,pgc);
    SLICELOOP4
    {
    s->ks(i,j) = p->S20;
    
    s->bedzh(i,j)=b->topobed(i,j);
    s->bedzh0(i,j)=b->topobed(i,j);
    }
    
    SLICELOOP4
    b->bed(i,j) = MAX(b->topobed(i,j),b->solidbed(i,j));
    
    active_ini_sflow(p,b,pgc);
    
    ini_parameters(p,pgc);
    ini_guard(p,pgc);
    log_ini(p);
}

void sediment_f::ini_parameters(lexer *p, ghostcell *pgc)
{
    // Fredsøe, p.199
    double rhosed=p->S22;
    double rhowat=p->W1;
    double g=9.81;
    double d50=p->S20;
    double visc=p->W2;
    double Ls = p->S20;
    double cd = 1.5;
    
    //s->ws=1.1*(rhosed/rhowat-1.0)*g*d50*d50;
    
    if(p->S23==0)
    s->ws = sqrt(4.0*(rhosed/rhowat-1.0)*g*d50/(3.0*cd));
    
    if(p->S23==1)
    s->ws = p->S23_val;
    
    if(p->mpirank==0)
    cout<<"ws: "<<s->ws<<endl;
    
}

void sediment_f::ini_guard(lexer *p, ghostcell *pgc)
{
    SLICELOOP4
    s->guard(i,j)=1.0;
    
    
    
    if(p->S78==1)
    {
        for(n=0;n<p->gcin_count;++n)
        {
        i=p->gcin[n][0];
        j=p->gcin[n][1];

        s->guard(i,j)=0.0;
        }
            
        for(int qq=0;qq<pBC->obj_count;++qq)
        for(n=0;n<pBC->patch[qq]->gcb_count;++n)
        {
        i=pBC->patch[qq]->gcb[n][0];
        j=pBC->patch[qq]->gcb[n][1];
            
            
        s->guard(i,j)=0.0;
        }
    }
    
    
    pgc->gcsl_start4(p,s->guard,1);
}