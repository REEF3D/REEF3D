/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the B117, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/liceonephases/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"onephase_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"ioflow.h"

onephase_f::onephase_f(lexer *p, fdm *a, ghostcell *pgc) : ddweno_f_nug(p),uf(p),vf(p),wf(p),urk1(p),vrk1(p),wrk1(p),xphi(p),yphi(p),zphi(p)
{
    gcval_u=10;
	gcval_v=11;
	gcval_w=12;
    
    
    dt=1.0e8;
	FLUIDLOOP
	{
	dt = MIN(dt,0.5*MIN3(p->DXP[IP],p->DYP[JP],p->DZP[KP]));
	}
    
    dt=pgc->timesync(dt);
}

onephase_f::~onephase_f()
{
}

void onephase_f::update(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow)
{   
    FLUIDLOOP
    {
    if(a->phi(i,j,k)<0.0)
    p->flag4[IJK]=AIR;

    if(a->phi(i,j,k)>=0.0)
    p->flag4[IJK]=WATER;
    }
    
    pgc->flagx(p,p->flag4);
    
    activenum=0;
    LOOP
	++activenum;
    
    activenum=pgc->globalisum(activenum);
    
    if(p->mpirank==0)
    cout<<"active number of cells: "<<activenum<<endl;
    
    FLUIDLOOP
	{
	p->flag1[IJK]=p->flag4[IJK];
	p->flag2[IJK]=p->flag4[IJK];
	p->flag3[IJK]=p->flag4[IJK];
	}
    
    FLUIDLOOP
	{
    if(p->flag4[IJK]==WATER && p->flag4[Ip1JK]==AIR)
    p->flag1[IJK]=AIR;
    
    if(p->flag4[IJK]==WATER && p->flag4[IJp1K]==AIR)
    p->flag2[IJK]=AIR;
    
    if(p->flag4[IJK]==WATER && p->flag4[IJKp1]==AIR)
    p->flag3[IJK]=AIR;
	}
    
    pgc->flagx(p,p->flag1);
    pgc->flagx(p,p->flag2);
    pgc->flagx(p,p->flag3);
}

void onephase_f::fsf_update(lexer *p, fdm *a, ghostcell *pgc)
{
}

void onephase_f::ini(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow)
{
}

