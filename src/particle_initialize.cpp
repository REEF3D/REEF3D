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

#include"particle_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"

void particle_f::ini(lexer* p, fdm* a, ghostcell* pgc, ioflow *pflow)
{
    LOOP
    active(i,j,k) = 0.0;
    
    // Box
    cellcount=0;
    for(qn=0;qn<p->Q110;++qn)
    LOOP
	if(p->XN[IP]>=p->Q110_xs[qn] && p->XN[IP]<p->Q110_xe[qn]
	&& p->YN[JP]>=p->Q110_ys[qn] && p->YN[JP]<p->Q110_ye[qn]
	&& p->ZN[KP]>=p->Q110_zs[qn] && p->ZN[KP]<p->Q110_ze[qn])
	{
	active(i,j,k) = 1.0;
    ++cellcount;
	}
    
    /*if(p->B139>0)
    srand(p->B139);

    if(p->B139==0)
    srand((unsigned)time(0));

    // make phases
	for(int n=0;n<p->wN;++n)
	ei[n]  = double(rand() % 628)/100.0;*/
    
    
    // guess particle demand
    
    
    // allocated
    
    // seed: distribute





} 
