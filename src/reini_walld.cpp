/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"reini_walld.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"reinidisc_f.h"

reini_walld::reini_walld(lexer* p, fdm *a):gradient(p),dab(p)
{
	prdisc = new reinidisc_f(p);
}

reini_walld::~reini_walld()
{
}

void reini_walld::start(fdm* a,lexer* p, field &f, ghostcell* pgc,ioflow* pflow)
{
	starttime=pgc->timer();

	pgc->start4(p,f,50);
    
    int qq;
	QQGC4LOOP
    if(p->gcb4[qq][4]==5|| p->gcb4[qq][4]==21|| p->gcb4[qq][4]==22)
    {
        i=p->gcb4[qq][0];
        j=p->gcb4[qq][1];
        k=p->gcb4[qq][2];
        
        if(p->gcb4[qq][3]==1)
        {
        f.V[Im1JK]=-0.5*p->DXN[IP] - 0.0*p->DXN[IP];  
        f.V[Im2JK]=-0.5*p->DXN[IP] - 1.0*p->DXN[IP];  
        f.V[Im3JK]=-0.5*p->DXN[IP] - 2.0*p->DXN[IP];  
        }

        if(p->gcb4[qq][3]==4)
        { 
        f.V[Ip1JK]=-0.5*p->DXN[IP] - 0.0*p->DXN[IP];  
        f.V[Ip2JK]=-0.5*p->DXN[IP] - 1.0*p->DXN[IP];  
        f.V[Ip3JK]=-0.5*p->DXN[IP] - 2.0*p->DXN[IP];  
        }

        if(p->gcb4[qq][3]==3)
        {
        f.V[IJm1K]=-0.5*p->DYN[JP] - 0.0*p->DYN[JP];
        f.V[IJm2K]=-0.5*p->DYN[JP] - 1.0*p->DYN[JP];
        f.V[IJm3K]=-0.5*p->DYN[JP] - 2.0*p->DYN[JP];
        }
        
        if(p->gcb4[qq][3]==2)
        {
        f.V[IJp1K]=-0.5*p->DYN[JP] - 0.0*p->DYN[JP];
        f.V[IJp2K]=-0.5*p->DYN[JP] - 1.0*p->DYN[JP];
        f.V[IJp3K]=-0.5*p->DYN[JP] - 2.0*p->DYN[JP];
        }
        
        if(p->gcb4[qq][3]==5)
        {
        f.V[IJKm1]=-0.5*p->DZN[KP] - 0.0*p->DZN[KP];
        f.V[IJKm2]=-0.5*p->DZN[KP] - 1.0*p->DZN[KP];
        f.V[IJKm3]=-0.5*p->DZN[KP] - 2.0*p->DZN[KP];
        }
                
        if(p->gcb4[qq][3]==6)
        {
        f.V[IJKp1]=-0.5*p->DZN[KP] - 0.0*p->DZN[KP];
        f.V[IJKp2]=-0.5*p->DZN[KP] - 1.0*p->DZN[KP];
        f.V[IJKp3]=-0.5*p->DZN[KP] - 2.0*p->DZN[KP];
        }
    }
    
    QQGC4LOOP
    if(p->gcb4[qq][4]==1|| p->gcb4[qq][4]==2|| p->gcb4[qq][4]==3)
    {
        i=p->gcb4[qq][0];
        j=p->gcb4[qq][1];
        k=p->gcb4[qq][2];
        n=p->gcb4[qq][5];
        
        if(p->gcb4[qq][3]==1)
        {
        f.V[Im1JK] = f.V[IJK];  
        f.V[Im2JK] = f.V[IJK];   
        f.V[Im3JK] = f.V[IJK];   
        }

        if(p->gcb4[qq][3]==4)
        { 
        f.V[Ip1JK] = f.V[IJK];   
        f.V[Ip2JK] = f.V[IJK];   
        f.V[Ip3JK] = f.V[IJK];   
        }

        if(p->gcb4[qq][3]==3)
        {
        f.V[IJm1K] = f.V[IJK]; 
        f.V[IJm2K] = f.V[IJK]; 
        f.V[IJm3K] = f.V[IJK]; 
        }
        
        if(p->gcb4[qq][3]==2)
        {
        f.V[IJp1K] = f.V[IJK]; 
        f.V[IJp2K] = f.V[IJK]; 
        f.V[IJp3K] = f.V[IJK]; 
        }
        
        if(p->gcb4[qq][3]==5)
        {
        f.V[IJKm1] = f.V[IJK]; 
        f.V[IJKm2] = f.V[IJK]; 
        f.V[IJKm3] = f.V[IJK]; 
        }
                
        if(p->gcb4[qq][3]==6)
        {
        f.V[IJKp1] = f.V[IJK]; 
        f.V[IJKp2] = f.V[IJK]; 
        f.V[IJKp3] = f.V[IJK]; 
        }
    }
    
    double dx;
    dt=1e9;
    if(p->N50==1)
    LOOP
    {
    dx = MIN3(p->DXN[IP],p->DYN[JP],p->DZN[KP]);

    dt = MIN(dt,0.5*dx);
    }
    
    reiniter=2*int(p->maxlength/(dt));
    
    reiniter = pgc->globalimax(reiniter);
  
	pgc->gcparax(p,f,4);
   
	for(int q=0;q<reiniter;++q)
	{

		prdisc->start(p,a,pgc,f,a->L,4);

		if(q==0)
		LOOP
		dab.V[IJK]=a->L.V[IJK];


		LOOP
		{
		f.V[IJK] += dt*0.5*(3.0*a->L.V[IJK] - dab.V[IJK]);

		dab.V[IJK]=a->L.V[IJK];
		}
        
	QQGC4LOOP
    if(p->gcb4[qq][4]==5|| p->gcb4[qq][4]==21|| p->gcb4[qq][4]==22)
    {
    i=p->gcb4[qq][0];
    j=p->gcb4[qq][1];
    k=p->gcb4[qq][2];
    n=p->gcb4[qq][5];
    
        if(p->gcb4[qq][3]==1 || p->gcb4[qq][3]==4)
        f.V[IJK] = 0.5*p->DXN[IP];  
        
        if(p->gcb4[qq][3]==3 || p->gcb4[qq][3]==2)
        f.V[IJK] = 0.5*p->DYN[JP];  
        
        if(p->gcb4[qq][3]==5 || p->gcb4[qq][3]==6)
        f.V[IJK] = 0.5*p->DZN[KP];  
	}
	
	pgc->gcparax(p,f,4);
	}

}


void reini_walld::step(fdm* a, lexer *p)
{
	

}




