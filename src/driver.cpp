/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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
#include"fdm.h"
#include"fdm2D.h"
#include"fdm_fnpf.h"
#include"fdm_nhf.h"
#include"lexer.h"
#include"waves_header.h"
#include"patchBC.h"

driver::driver(int& argc, char **argv)
{
	p = new lexer;
	pgc = new ghostcell(argc,argv,p);

	if(p->mpirank==0)
    {
    cout<<endl<<"REEF3D (c) 2008-2024 Hans Bihs"<<endl;
    sprintf(version,"v_240218");
    cout<<endl<<":: Open-Source Hydrodynamics" <<endl;
    cout<<endl<<version<<endl<<endl;
    }

	p->lexer_read(pgc);
    p->vellast();
	pgc->gcini(p);
    p->gridini(pgc);
    patchBC_logic();


    if(p->mpirank==0)
    {
    if(p->A10==2)
    cout<<endl<<"REEF3D::SFLOW" <<endl<<endl;

    if(p->A10==3)
    cout<<endl<<"REEF3D::FNPF" <<endl<<endl;

    if(p->A10==4)
    cout<<endl<<"REEF3D::PTF" <<endl<<endl;

    if(p->A10==5)
    cout<<endl<<"REEF3D::NHFLOW"<<endl<<endl;

    if(p->A10==6)
    cout<<endl<<"REEF3D::CFD" <<endl<<endl;
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

// 3D Framework
    // sigma grid - FNPF
    if(p->A10==3)
    {
        p->flagini();
        p->gridini_patchBC();
        pgc->flagfield(p);
        pgc->tpflagfield(p);
        makegrid_sigma(p,pgc);
        makegrid2D_basic(p,pgc);

        pgc->ndflag_update(p);

        fnpf_driver();
    }

    // sigma grid - NHFLOW
    if(p->A10==5)
    {
        BASELOOP
        if(p->flagslice4[IJ]<0)
        p->flag4[IJK]=-10;

        p->flagini();
        p->gridini_patchBC();
        pgc->flagfield(p);
        pgc->tpflagfield(p);
        makegrid_sigma(p,pgc);
        makegrid2D(p,pgc);

        pgc->ndflag_update(p);

        nhflow_driver();
    }

    // fixed grid - PTF & NSEWAVE & CFD
    if(p->A10==4 || p->A10==6)
    {
        p->flagini();
        p->gridini_patchBC();
        pgc->flagfield(p);
        pgc->tpflagfield(p);
        makegrid(p,pgc);
        makegrid2D(p,pgc);

        pgc->ndflag_update(p);

        if(p->A10==4)
        ptf_driver();

        if(p->A10==6)
        cfd_driver();
    }
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

    // Start MAINLOOP
    if(p->X10==0 && p->Z10==0 && p->G3==1 && p->N40==4)
    loop_cfd_sf(a);

    else
    if((p->X10==1  || p->Z10!=0) && (p->N40==4))
    loop_cfd_df(a);

    else
    loop_cfd(a);
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

    // Start MAINLOOP
    loop_nhflow();
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

    // Start MAINLOOP
    loop_fnpf();
}

void driver::ptf_driver()
{
    if(p->mpirank==0)
	cout<<"initialize fdm"<<endl;

    a=new fdm(p);

    aa=a;
    pgc->fdm_update(a);

    logic_ptf();

    driver_ini_ptf();

    // Start MAINLOOP
    loop_ptf(a);
}

void driver::sflow_driver()
{
    if(p->mpirank==0)
	cout<<"initialize fdm"<<endl;

    b=new fdm2D(p);
    bb=b;

    psflow = new sflow_f(p,b,pgc,pBC);

    makegrid2D_cds(p,pgc,b);

    // Start SFLOW
	psflow->start(p,b,pgc);
}

driver::~driver()
{
}
