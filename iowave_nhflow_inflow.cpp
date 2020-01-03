/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"iowave.h"
#include"lexer.h"
#include"ghostcell.h"
 
void iowave::nhflow_inflow(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
{
    if(p->I230==0)
    {
    if(p->B98==0)
    nhflow_inflow_plain(p,a,pgc,u,v,w);
    
	if(p->B98==3)
	nhflow_dirichlet_wavegen(p,a,pgc,u,v,w);
	
	if(p->B98==5)
	nhflow_active_wavegen(p,a,pgc,u,v,w);
	}
    
	if(p->B99==3||p->B99==4||p->B99==5)
	nhflow_active_beach(p,a,pgc,u,v,w);
    
    if(p->I230>0)
    ff_inflow(p,a,pgc,u,v,w);
}

void iowave::nhflow_inflow_plain(lexer *p, fdm* a, ghostcell* pgc, field& u, field& v, field& w)
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
    }
}
