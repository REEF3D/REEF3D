/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"sediment_part.h"
#include"sediment_fdm.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"ioflow.h"
#include"topo.h"
#include"bedshear.h"
#include"patchBC_interface.h"
#include"bedslope.h"
#include <sys/stat.h>

sediment_part::sediment_part(lexer *p, fdm *a, ghostcell *pgc, turbulence *ppturb, patchBC_interface *ppBC) : por(p), d50(p)
{
    pBC = ppBC;
    pturb = ppturb;
    
    sediment_logic(p,a,pgc,pturb);

	p->gcin4a_count=p->gcin_count;
	p->gcout4a_count=p->gcout_count;
	
    
    volume_token=0;
    
    if(p->F50==1)
	gcval_eta = 51;
    
    if(p->F50==2)
	gcval_eta = 52;
    
    if(p->F50==3)
	gcval_eta = 53;
    
    if(p->F50==4)
	gcval_eta = 54;
    
    // Create Folder
	if(p->mpirank==0 && p->Q180>0 && (p->Q181>0||p->Q182>0))
    mkdir("./REEF3D_CFD_SedPart",0777);
}

sediment_part::~sediment_part()
{
}

void sediment_part::start_cfd(lexer *p, fdm *a, ghostcell *pgc, ioflow *pflow, reinitopo *preto, solver *psolv)
{
    sedcalc=0;
    
	if((p->S41==1 && p->count>=p->S43) || (p->S41==2 && p->simtime>=p->S45) || (p->S41==3 && p->simtime/p->wT>=p->S47))
	{
		sediment_algorithm_cfd(p,a,pgc,pflow,preto,pturb);
		
    
    sedcalc=1;
	}
    
    if(sedcalc==0)
    {
    fill_bedk(p,a,pgc);
    waterlevel(p,a,pgc);
    pbedshear->taubed(p,a,pgc,s);
    }
}