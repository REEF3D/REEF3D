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

#include"driver.h"
#include"ghostcell.h"
#include"lexer.h"
#include"fdm.h"
#include"fdm2D.h"
#include"freesurface_header.h"
#include"turbulence_header.h"
#include"momentum_header.h"
#include"pressure_header.h"
#include"fdm_header.h"
#include"sediment_header.h"
#include"convection_header.h"
#include"solver_header.h"
#include"field_header.h"
#include"heat_header.h"
#include"concentration_header.h"
#include"benchmark_header.h"
#include"6DOF_header.h"
#include"waves_header.h"
#include"sflow_header.h"
#include<sys/stat.h>
#include<sys/types.h>

void driver::driver_ini_sflow()
{
p->count=0;

p->cellnumtot=pgc->globalisum(p->cellnum);
p->pointnumtot=pgc->globalisum(p->pointnum);


if(p->mpirank==0)
cout<<"number of cells: "<<p->cellnumtot<<endl;

	log_ini();


if(p->mpirank==0)
cout<<"starting driver_ini_PFLOW"<<endl;
    
    pflow->ini2D(p,b,pgc);
	
    if(p->toporead>0 ||p->solidread==1)
    {
    geotopo gtopo(p,a,pgc);
    gtopo.start(p,a,pgc,pflow,preto,pvrans);
    }
	
    ptstep->ini(a,p,pgc);
    preini->start(a,p,a->phi, pgc, pflow);
    pflow->gcio_update(p,a,pgc);
    //ppf->ini(p,a,pgc,pflow,preini,pfsfdisc);
    pflow->fi_relax(p,pgc,a->Fi,a->phi);
    pgc->start4(p,a->Fi,250);
    
	
	pflow->inflow(p,a,pgc,a->u,a->v,a->w);
    

    pprint->start(a,p,pgc,pturb,pheat,pflow,psolv,pdata,pconc,pmp,psed);

	
	p->gctime=0.0;
    p->xtime=0.0;
	p->wavetime=0.0;
	p->field4time=0.0;

if(p->mpirank==0)
cout<<"starting mainloop.PFLOW"<<endl;

}




