/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"sediment_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"topo_relax.h"
#include"ioflow.h"


void sediment_f::update(lexer *p, fdm *a,ghostcell *pgc, ioflow *pflow)
{
    for(int qn=0;qn<3;++qn)
    prelax->start(p,a,pgc);
    
    p->sedtime+=p->dtsed;
    pgc->topo_update(p,a);
    pflow->gcio_update(p,a,pgc);
}

void sediment_f::bedlevel(lexer *p, fdm *a, ghostcell *pgc)
{
    p->bedmin=1.0e15;
    p->bedmax=-1.0e15;

    SLICELOOP4
    {
        p->bedmin = MIN(p->bedmin, a->bedzh(i,j));
        p->bedmax = MAX(p->bedmax, a->bedzh(i,j));
    }
	
    p->bedmin=pgc->globalmin(p->bedmin);
    p->bedmax=pgc->globalmax(p->bedmax);

    if(p->mpirank==0)
    {
    cout<<"bedmin: "<<setprecision(4)<<p->bedmin<<endl;
    cout<<"bedmax: "<<setprecision(4)<<p->bedmax<<endl<<endl;
    }
}

void sediment_f::relax(lexer *p, fdm *a,ghostcell *pgc)
{
    prelax->start(p,a,pgc);
}

void sediment_f::topo_zh_update(lexer *p, fdm *a,ghostcell *pgc)
{
	pgc->gcsl_start4(p,a->bedzh,1);
	
    ALOOP
    {
    if(p->pos_x()>p->S77_xs && p->pos_x()<p->S77_xe)
    a->topo(i,j,k)=-a->bedzh(i,j)+p->pos_z();
    }
	
	pgc->start4a(p,a->topo,150);
}
