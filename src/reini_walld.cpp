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

#include"reini_walld.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"reinidisc_f.h"


reini_walld::reini_walld(lexer* p, fdm *a):gradient(p),f(p),dab(p),L(p)
{
	prdisc = new reinidisc_f(p);
}

reini_walld::~reini_walld()
{
}

void reini_walld::start(fdm* a,lexer* p,field& b, ghostcell* pgc,ioflow* pflow)
{
	starttime=pgc->timer();
	
	sizeM=p->sizeM4;
	
	n=0;
	LOOP
	{
	f.V[n]=b(i,j,k);
	++n;
	}

	pgc->start4vec(p,f,50);
    
    int qq;
	QQGC4LOOP
    if(p->gcb4[qq][4]==5|| p->gcb4[qq][4]==21|| p->gcb4[qq][4]==22)
    {
        i=p->gcb4[qq][0];
        j=p->gcb4[qq][1];
        k=p->gcb4[qq][2];
        n=p->gcb4[qq][5];
        
        if(p->gcb4[qq][3]==1)
        {
        f.V[Im1_J_K_4]=-0.5*p->DXM - 0.0*p->DXM;  
        f.V[Im2_J_K_4]=-0.5*p->DXM - 1.0*p->DXM;  
        f.V[Im3_J_K_4]=-0.5*p->DXM - 2.0*p->DXM;  
        }

        if(p->gcb4[qq][3]==4)
        { 
        f.V[Ip1_J_K_4]=-0.5*p->DXM - 0.0*p->DXM;  
        f.V[Ip2_J_K_4]=-0.5*p->DXM - 1.0*p->DXM;  
        f.V[Ip3_J_K_4]=-0.5*p->DXM - 2.0*p->DXM;  
        }

        if(p->gcb4[qq][3]==3)
        {
        f.V[I_Jm1_K_4]=-0.5*p->DXM - 0.0*p->DXM;
        f.V[I_Jm2_K_4]=-0.5*p->DXM - 1.0*p->DXM;
        f.V[I_Jm3_K_4]=-0.5*p->DXM - 2.0*p->DXM;
        }
        
        if(p->gcb4[qq][3]==2)
        {
        f.V[I_Jp1_K_4]=-0.5*p->DXM - 0.0*p->DXM;
        f.V[I_Jp2_K_4]=-0.5*p->DXM - 1.0*p->DXM;
        f.V[I_Jp3_K_4]=-0.5*p->DXM - 2.0*p->DXM;
        }
        
        if(p->gcb4[qq][3]==5)
        {
        f.V[I_J_Km1_4]=-0.5*p->DXM - 0.0*p->DXM;
        f.V[I_J_Km2_4]=-0.5*p->DXM - 1.0*p->DXM;
        f.V[I_J_Km3_4]=-0.5*p->DXM - 2.0*p->DXM;
        }
                
        if(p->gcb4[qq][3]==6)
        {
        f.V[I_J_Kp1_4]=-0.5*p->DXM - 0.0*p->DXM;
        f.V[I_J_Kp2_4]=-0.5*p->DXM - 1.0*p->DXM;
        f.V[I_J_Kp3_4]=-0.5*p->DXM - 2.0*p->DXM;
        }
    }

	pgc->gcparaxvec(p,f,4);
	dt=0.5*p->DXM;
	reiniter=2*int(p->maxlength/(dt));

	for(int q=0;q<reiniter;++q)
	{

		prdisc->start(p,a,pgc,f,L,4);

		if(q==0)
		NLOOP
		dab.V[n]=L.V[n];


		NLOOP
		{
		f.V[n] += dt*0.5*(3.0*L.V[n] - dab.V[n]);

		dab.V[n]=L.V[n];
		}
	
	QQGC4LOOP
    if(p->gcb4[qq][4]==5|| p->gcb4[qq][4]==21|| p->gcb4[qq][4]==22)
    {
    i=p->gcb4[qq][0];
    j=p->gcb4[qq][1];
    k=p->gcb4[qq][2];
    n=p->gcb4[qq][5];
    
	f.V[I_J_K_4]=0.5*p->DXM;
	}
	
	pgc->gcparaxvec(p,f,4);
	}

    // backfill
	n=0;
	LOOP
	{
	b(i,j,k)=f.V[n];
	++n;
	}
}

void reini_walld::startV(fdm* a,lexer* p,vec &f, ghostcell* pgc,ioflow* pflow)
{ 
    
}

void reini_walld::step(fdm* a, lexer *p)
{
	

}




