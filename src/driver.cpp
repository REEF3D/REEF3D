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

#include"driver.h"
#include"regression_dump.h"
#include"ghostcell.h"
#include"fdm.h"
#include"fdm2D.h"
#include"fdm_fnpf.h"
#include"fdm_nhf.h"
#include"lexer.h"
#include"waves_header.h"
#include"patchBC.h"
#include"runlog.h"
#include"seastate_f.h"

driver::driver(int& argc, char **argv)
{
	p = new lexer;
	pgc = new ghostcell(argc,argv,p);

	if(p->mpirank==0)
    {
    cout<<endl<<"REEF3D (c) 2008-2026 Hans Bihs"<<endl;
    sprintf(version,"v_261009");
    cout<<endl<<":: Open-Source Hydrodynamics" <<endl;
    cout<<endl<<version<<endl;
    cout<<endl<<"github branch: "<<BRANCH<<endl;
    cout<<endl<<"github version: "<<VERSION<<endl;
    }

    pgc->mpi_check(p);
    
	p->lexer_read(pgc);
    p->vellast();
    
	pgc->gc_ini(p);
    pgc->gcx_ini(p);
    
    p->gridini(pgc);
    patchBC_logic();

    // run log: REEF3D_Case/REEF3D_<SOLVER>_run.jsonl (rank 0 writes)
    p->plog = new runlog(p,VERSION);


    if(p->mpirank==0)
    {
    if(p->A10==2)
    cout<<endl<<"REEF3D::SFLOW" <<endl<<endl;

    if(p->A10==3)
    cout<<endl<<"REEF3D::FNPF" <<endl<<endl;

    if(p->A10==5)
    cout<<endl<<"REEF3D::NHFLOW"<<endl<<endl;

    if(p->A10==6)
    cout<<endl<<"REEF3D::CFD" <<endl<<endl;

    if(p->A10==7)
    cout<<endl<<"REEF3D::SEASTATE" <<endl<<endl;
    }

    // PTF (A 10 4) was removed
    if(p->A10==4)
    {
        if(p->mpirank==0)
        cout<<endl<<"A 10 4: REEF3D::PTF has been removed, use REEF3D::FNPF (A 10 3) or REEF3D::CFD (A 10 6)."<<endl<<endl;

        pgc->final(true);
    }

// 2D Framework - SFLOW
    if(p->A10==2)
    {
        p->flagini2D();
        p->gridini2D();
        makegrid2D(p,pgc);
        pBC->patchBC_ini(p,pgc);
        sflow_driver();
    }

// 2D Framework - SEASTATE
    if(p->A10==7)
    {
        p->flagini2D();
        p->gridini2D();
        makegrid2D(p,pgc);
        seastate_driver();
    }

// 3D Framework
    // sigma grid - FNPF
    if(p->A10==3)
    {
        p->flagini();
        pgc->flagfield(p);
        makegrid_sigma(p,pgc);
        makegrid2D_basic(p,pgc);

        fnpf_driver();
    }

    // sigma grid - NHFLOW
    if(p->A10==5)
    {
        BASELOOP
        if(p->flagslice4[IJ]<0)
        p->flag4[IJK]=-10;

        p->flagini();
        pgc->flagfield(p);
        makegrid_sigma(p,pgc);
        makegrid2D(p,pgc);

        nhflow_driver();
    }

    // fixed grid - CFD
    if(p->A10==6)
    {
        p->flagini();
        pgc->flagfield(p);
        makegrid(p,pgc);
        makegrid2D(p,pgc);

        cfd_driver();
    }
}

void driver::sflow_driver()
{
    if(p->mpirank==0)
	cout<<"initialize fdm"<<endl;

    b=new fdm2D(p);
    bb=b;
    
    pgc->fdm2D_update(b);

    psflow = new sflow_f(p,b,pgc,pBC);

    makegrid2D_cds(p,pgc,b);

    // Start SFLOW
	psflow->start(p,b,pgc);
}

void driver::seastate_driver()
{
    // 2D grid set-up of SFLOW without the SFLOW fdm (makegrid2D_cds)
    p->flagini2D();
    p->gridini2D();
    pgc->sizeS_update(p);

    pseastate = new seastate_f(p,pgc);

    // Start SEASTATE
    pseastate->start(p,pgc);
}

void driver::fnpf_driver()
{
    if(p->mpirank==0)
	cout<<"initialize fdm"<<endl;

    p->grid2Dsize();

    c=new fdm_fnpf(p);

    pgc->fdm_fnpf_update(c);

    makegrid_sigma_cds(p,pgc);

    logic_fnpf();

    driver_ini_fnpf();

    preg = new regression_dump(p);
    preg->fnpf_ini(p,c,pgc);

    // Start MAINLOOP
    loop_fnpf();
}

void driver::nhflow_driver()
{
    if(p->mpirank==0)
	cout<<"initialize fdm"<<endl;

	d=new fdm_nhf(p);

    pgc->fdm_nhf_update(d);

    makegrid_sigma_cds(p,pgc);

    logic_nhflow();

    driver_ini_nhflow();

    preg = new regression_dump(p);
    preg->nhflow_ini(p,d,pgc);

    // Start MAINLOOP
    loop_nhflow();
}

void driver::cfd_driver()
{
    if(p->mpirank==0)
	cout<<"initialize fdm "<<endl;

    a=new fdm(p);

	aa=a;
    pgc->fdm_update(a);

    logic_cfd();

    driver_ini_cfd();

    preg = new regression_dump(p);
    preg->cfd_ini(p,a,pgc,pturb,pconc);

    // Start MAINLOOP
    if(p->X10==0 && p->Z10==0 && p->N40==14)
    loop_cfd_sf(a);

    else
    if((p->X10==1  || p->Z10!=0) && (p->N40==14))
    loop_cfd_df(a);

    else
    loop_cfd(a);
}

driver::~driver()
{
}
