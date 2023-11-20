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

#define WLVL (fabs(c->WL(i,j))>1.0e-20?c->WL(i,j):1.0e20)

void driver::driver_ini_fnpf()
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

	LOOP
	++p->cellnum;
    
    LOOP
    ++p->tpcellnum;
    
    p->count=0;
    
// --
    p->cellnumtot=pgc->globalisum(p->cellnum);
    p->pointnumtot=pgc->globalisum(p->pointnum);


    if(p->mpirank==0)
    cout<<"number of cells: "<<p->cellnumtot<<endl;
    
    
    
    // maxcoor

    p->maxlength=-1.0e9;
    p->xcoormax=-1.0e9;
    p->xcoormin=1.0e9;
    p->ycoormax=-1.0e9;
    p->ycoormin=1.0e9;
    p->zcoormax=-1.0e9;
    p->zcoormin=1.0e9;

    LOOP
    {
        p->xcoormax = MAX(p->xcoormax,p->XN[IP1]);
        p->xcoormin = MIN(p->xcoormin,p->XN[IP]);
        p->ycoormax = MAX(p->ycoormax,p->YN[JP1]);
        p->ycoormin = MIN(p->ycoormin,p->YN[JP]);
        p->zcoormax = MAX(p->zcoormax,p->ZN[KP1]);
        p->zcoormin = MIN(p->zcoormin,p->ZN[KP]);
     }

     p->maxlength=MAX(p->maxlength,p->xcoormax-p->xcoormin);
     p->maxlength=MAX(p->maxlength,p->ycoormax-p->ycoormin);
     p->maxlength=MAX(p->maxlength,p->zcoormax-p->zcoormin);

     p->maxlength=pgc->globalmax(p->maxlength);
	 
	 p->xcoormax=pgc->globalmax(p->xcoormax);
	 p->ycoormax=pgc->globalmax(p->ycoormax);
	 p->zcoormax=pgc->globalmax(p->zcoormax);
	 
	 p->xcoormin=pgc->globalmin(p->xcoormin);
	 p->ycoormin=pgc->globalmin(p->ycoormin);
	 p->zcoormin=pgc->globalmin(p->zcoormin);
     
    if(p->mpirank==0)
    cout<<"starting driver_ini_FNPF"<<endl;

	
    // bed ini
    SLICELOOP4
	c->bed(i,j) = p->bed[IJ];
    
    pgc->gcsl_start4(p,c->bed,50);
    
    // bc ini
    SLICELOOP4
	c->bc(i,j) = 0;
    
    pgc->gcsl_start4int(p,c->bc,50);
   
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];

    c->bc(i-1,j) = 1;
    }
    
    for(n=0;n<p->gcslout_count;n++)
    {
    i=p->gcslout[n][0];
    j=p->gcslout[n][1];
    
    c->bc(i+1,j) = 2;
    }
   
    // 2D mesh
    count=0;
	p->pointnum2D=0;
	p->cellnum2D=0;
	p->polygon_sum=0;
    
   
    TPSLICELOOP  
	{
	++count;
	++p->pointnum2D;
	c->nodeval2D(i,j)=count;
    }
	
	SLICEBASELOOP
	++p->polygon_sum;
	
	p->polygon_sum *=2;

	SLICELOOP4
	++p->cellnum2D;
    
    SLICELOOP4
	++p->cellnum2D;

    p->cellnumtot2D=pgc->globalisum(p->cellnum2D);
    
    
    // eta ini
	SLICELOOP4
	c->eta(i,j) = 0.0;

    pgc->gcsl_start4(p,c->eta,50);
    
     SLICELOOP4
    c->WL(i,j) = MAX(0.0,c->eta(i,j) + p->wd - c->bed(i,j));
    
    
    
    SLICELOOP4
    p->sigz[IJ] = 1.0/WLVL;
    
    ppfsg->ini_wetdry(p,c,pgc);    // ini wetdry and coastline
    
    pflow->gcio_update(p,a,pgc);
    pflow->ini_fnpf(p,c,pgc);  // including fullini

    ppfsg->ini(p,c,pgc,pflow,preini,poneph);  // --- 
    
    pgc->start7V(p,c->Fi,c->bc,250);
    ppfsg->inidisc(p,c,pgc,pflow,psolv);    // ini wetdry and coastline
    
    pflow->eta_relax(p,pgc,c->eta);
    pflow->fivec_relax(p,pgc,c->Fi);
    pflow->fifsf_relax(p,pgc,c->Fifsf);
    pgc->gcsl_start4(p,c->eta,50);
    pgc->gcsl_start4(p,c->Fifsf,50);
    pgc->start7V(p,c->Fi,c->bc,250);
    
    pftstep->ini(c,p,pgc);
	
	pfprint->start(p,c,pgc,pflow);
	
    
	p->gctime=0.0;
    p->xtime=0.0;
	p->wavetime=0.0;
	p->field4time=0.0;

    
     if(p->mpirank==0)
    cout<<"starting mainloop.FNPF"<<endl;

}




