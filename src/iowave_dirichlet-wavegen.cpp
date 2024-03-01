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

void iowave::dirichlet_wavegen(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{
        count=0;
		for(n=0;n<p->gcin_count;n++)
		{
		i=p->gcin[n][0];
		j=p->gcin[n][1];
		k=p->gcin[n][2];		

        uvel=uval[count]*ramp(p);
        vvel=vval[count]*ramp(p);
        wvel=wval[count]*ramp(p);
            
			if(a->phi(i-1,j,k)>=0.0)
			{
			u(i-1,j,k)=uvel+p->Ui;
			u(i-2,j,k)=uvel+p->Ui;
			u(i-3,j,k)=uvel+p->Ui;
            
            v(i-1,j,k)=vvel;
			v(i-2,j,k)=vvel;
			v(i-3,j,k)=vvel;
			
			w(i-1,j,k)=wvel;
			w(i-2,j,k)=wvel;
			w(i-3,j,k)=wvel;
			}

			if(a->phi(i-1,j,k)<0.0 && a->phi(i-1,j,k)>=-p->F45*p->DZP[KP])
			{
			fac= p->B122*(1.0 - fabs(a->phi(i-1,j,k))/(p->F45*p->DZP[KP]));
            
			u(i-1,j,k)=uvel*fac + p->Ui;
			u(i-2,j,k)=uvel*fac + p->Ui;
			u(i-3,j,k)=uvel*fac + p->Ui;
            
            v(i-1,j,k)=vvel*fac;
			v(i-2,j,k)=vvel*fac;
			v(i-3,j,k)=vvel*fac;
            
			w(i-1,j,k)=wvel*fac;
			w(i-2,j,k)=wvel*fac;
			w(i-3,j,k)=wvel*fac;
			}

			if(a->phi(i-1,j,k)<-p->F45*p->DXM)
			{
			u(i-1,j,k)=0.0 + p->Ui;
			u(i-2,j,k)=0.0 + p->Ui;
			u(i-3,j,k)=0.0 + p->Ui;
            
            if(p->W50_air==1)
            {
            u(i-1,j,k)+=p->W50;
            u(i-2,j,k)+=p->W50;
            u(i-3,j,k)+=p->W50;
            }
            
            v(i-1,j,k)=0.0;
			v(i-2,j,k)=0.0;
			v(i-3,j,k)=0.0;

			w(i-1,j,k)=0.0;
			w(i-2,j,k)=0.0;
			w(i-3,j,k)=0.0;
			}
        ++count;
		}
        
        
        if(p->B98==3||p->B98==4||p->B99==3||p->B99==4||p->B99==5)
		{
		for(int q=0;q<4;++q)
		for(n=0;n<p->gcin_count;++n)
		{
		i=p->gcin[n][0]+q;
		j=p->gcin[n][1];
		k=p->gcin[n][2];

		if(a->phi(i,j,k)<0.0)
		a->eddyv(i,j,k)=MIN(a->eddyv(i,j,k),1.0e-4);
		}
		pgc->start4(p,a->eddyv,24);
		}
        
    // PTF
    /*if(p->A10==4)
    {
        for(n=0;n<p->gcslin_count;n++)
        {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        k=a->etaloc(i,j);
        {
        a->Fifsf(i-1,j) = a->Fifsf(i,j) - u(i-1,j,k)*1.0*p->DXP[IM1];
        a->Fifsf(i-2,j) = a->Fifsf(i,j) - u(i-1,j,k)*2.0*p->DXP[IM1];
        a->Fifsf(i-3,j) = a->Fifsf(i,j) - u(i-1,j,k)*3.0*p->DXP[IM1];
        }
         }
    }*/
    
    
}
