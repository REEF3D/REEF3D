/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"ioflow_gcio.h"
#include"ioflow_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"vrans.h"
#include"patchBC_interface.h"

void ioflow_f::gcio_update(lexer *p, fdm *a, ghostcell *pgc)
{
    ioflow_gcio_lists(p,p->flagsf4,nullptr);

    if(p->I10==1 && p->count==0)
    velini(p,a,pgc);

    ioflow_gcio_marks(p,p->flagsf4,pBC);
}


void ioflow_f::gcio_update_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc)
{
    ioflow_gcio_lists(p,p->DF,p->wet);

    ioflow_gcio_marks(p,p->DF,pBC);
}

void ioflow_f::inflow_walldist(lexer *p, fdm *a, ghostcell *pgc, convection *pconvec, reini *preini, ioflow *pflow)
{

	p->del_Darray(walldin, walldin_size);
	p->del_Darray(walldout, walldout_size);
	
	walldin_size=p->gcin_count;
	walldout_size=p->gcout_count;
	
	p->Darray(walldin, walldin_size);
    p->Darray(walldout, walldout_size);


    for(n=0;n<p->gcin_count;++n)
    {
    i=p->gcin[n][0];
    j=p->gcin[n][1];
    k=p->gcin[n][2];
    walldin[n] = a->walld(i,j,k);
    }

    for(n=0;n<p->gcout_count;++n)
    {
    i=p->gcout[n][0];
    j=p->gcout[n][1];
    k=p->gcout[n][2];

    walldout[n] = a->walld(i,j,k);
    }
}

void ioflow_f::veltimesave(lexer *p, fdm *a, ghostcell *pgc, vrans *pvrans)
{
    pvrans->veltimesave(p,a,pgc);
}

void ioflow_f::vrans_sed_update(lexer *p,fdm *a,ghostcell *pgc, vrans *pvrans)
{
    pvrans->sed_update(p,a,pgc);
}


