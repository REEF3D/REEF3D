/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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

#include"driver.h"
#include"ghostcell.h"
#include"lexer.h"
#include"fdm_fnpf.h"
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
#include"lexer.h"
#include"cart1.h"
#include"cart2.h"
#include"cart3.h"
#include"cart4.h"
#include"cart4a.h"
#include<sys/stat.h>
#include<sys/types.h>

#define WLVL (fabs(c->WL(i,j))>1.0e-20?c->WL(i,j):1.0-20)

void driver::driver_ini_fnpf_sg()
{
    // count cells
    int count=0;
	p->pointnum=0;
	p->cellnum=0;
	p->tpcellnum=0;

	TPLOOP
	{
	++count;
	++p->pointnum;
    c->nodeval(i,j,k)=count;
	}

	FLOOP
	++p->cellnum;
    
    LOOP
    ++p->tpcellnum;
    
    //cout<<p->mpirank<<" POINTNUM: "<<p->pointnum<<" CELLNUM: "<<p->cellnum<<endl;
    p->count=0;
    
// --
    p->cellnumtot=pgc->globalisum(p->cellnum);
    p->pointnumtot=pgc->globalisum(p->pointnum);


    if(p->mpirank==0)
    cout<<"number of cells: "<<p->cellnumtot<<endl;

      //  log_ini();

    if(p->mpirank==0)
    cout<<"starting driver_ini_FNPF"<<endl;

    
    /*
    // Solid
    if(p->G39==1)
    {
    solid solid_object(p,a,pgc);
    solid_object.start(p,a,pgc,pflow,pconvec,preto);
    }
    
    // Geotopo
    if((p->G50>0 && p->G51>0) || p->G60>0 || p->G61>0)
    {
    geotopo gtopo(p,a,pgc);
    gtopo.start(p,a,pgc,pflow,pconvec,preto);
    }*/
    
    // bed ini
    SLICELOOP4
	c->bed(i,j) = p->bed[IJ];
    
    pgc->gcsl_start4(p,c->bed,50);
    
    // eta ini
	SLICELOOP4
	c->eta(i,j) = 0.0;

    pgc->gcsl_start4(p,c->eta,50);
    
     SLICELOOP4
    c->WL(i,j) = MAX(0.0,c->eta(i,j) + p->wd - c->bed(i,j));
    
    p->Darray(p->sigz,p->imax*p->jmax);
    
    SLICELOOP4
    p->sigz[IJ] = 1.0/WLVL;
    
    pflow->ini_fnpf(p,c,pgc);
    pftstep->ini(c,p,pgc);

    ppfsg->ini(p,c,pgc,pflow,preini,poneph);  // --- 
    pflow->eta_relax(p,pgc,c->eta);
    pflow->fivec_relax(p,pgc,c->Fi);
    pflow->fifsf_relax(p,pgc,c->Fifsf);
    pgc->gcsl_start4(p,c->eta,50);
    pgc->gcsl_start4(p,c->Fifsf,50);
    pgc->start7V(p,c->Fi,250);
    
    
    ppfsg->inidisc(p,c,pgc);

	pfprint->start(p,c,pgc,pflow);

    
	p->gctime=0.0;
    p->xtime=0.0;
	p->wavetime=0.0;
	p->field4time=0.0;

if(p->mpirank==0)
cout<<"starting mainloop.FNPF"<<endl;

}