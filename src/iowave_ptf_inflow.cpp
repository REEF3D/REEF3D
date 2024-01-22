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
#include"fdm_ptf.h"
#include"ghostcell.h"
void iowave::inflow_ptf(lexer *p, fdm_ptf *a, ghostcell *pgc, double *Fi_, double *Uin_, slice &Fifsf, slice &eta)
{
    if(p->B98==3 || p->B98==4)
	dirichlet_wavegen_ptf(p,a,pgc,Fi_,Uin_,Fifsf,eta);
    
    if(p->B99==3||p->B99==4||p->B99==5)
	active_beach_ptf(p,a,pgc,Fi_,Uin_,Fifsf,eta);
}

void iowave::rkinflow_ptf(lexer *p, fdm_ptf *a, ghostcell *pgc, slice &frk, slice &f)
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