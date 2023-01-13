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
#include"ghostcell.h"
#include"fdm_fnpf.h"
 
void iowave::inflow_fnpf(lexer *p, fdm_fnpf *c, ghostcell *pgc, double *Fi, double *Uin, slice &Fifsf, slice &eta)
{
    if(p->B98==3 || p->B98==4)
	dirichlet_wavegen_fnpf(p,c,pgc,Fi,Uin,Fifsf,eta);
    
    if(p->B99==3||p->B99==4||p->B99==5)
	active_beach_fnpf(p,c,pgc,Fi,Uin,Fifsf,eta);
}

void iowave::rkinflow_fnpf(lexer *p, fdm_fnpf*, ghostcell *pgc, slice &frk, slice &f)
{
    for(n=0;n<p->gcslin_count;n++)
    {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        frk(i-1,j) = f(i-1,j);
        frk(i-2,j) = f(i-2,j);
        frk(i-3,j) = f(i-3,j);
    }
    
}