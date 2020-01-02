/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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
#include"fdm.h"
#include"fdm2D.h"
#include"fdm_fnpf.h"
#include"lexer.h"
#include"waves_header.h"

driver::driver(int& argc, char **argv)
{
	p = new lexer;
	pgc = new ghostcell(argc,argv,p);

	if(p->mpirank==0)
    {
    cout<<endl<<"REEF3D (c) 2008-2019 Hans Bihs"<<endl;
    cout<<endl<<":: Open-Source Hydrodynamics" <<endl;
    cout<<endl<<"v_200102" <<endl<<endl;
    }

	p->lexer_read(pgc);
	pgc->gcini(p);
    p->gridini(pgc);

    if(p->mpirank==0)
    {
    if(p->A10==2)
    cout<<endl<<"REEF3D::SFLOW" <<endl<<endl;

    if(p->A10==3)
    cout<<endl<<"REEF3D::FNPF" <<endl<<endl;

    if(p->A10==4)
    cout<<endl<<"REEF3D::NSEWAVE"<<endl<<endl;

    if(p->A10==44)
    cout<<endl<<"REEF3D::NHFLOW"<<endl<<endl;

    if(p->A10==5)
    cout<<endl<<"REEF3D::CFD" <<endl<<endl;
    }

    // 3D Framework
    if(p->A10==3 && p->A300==1)
    {
        p->flagini();
        p->gridini_outflow();
        pgc->flagfield(p);
        pgc->tpflagfield(p);
        makegrid_fnpf(p,pgc);

        pgc->ndflag_update(p);

        pfsg_driver();
    }

    if((p->A10==3 && p->A300==2) || p->A10==4 || p->A10==44 || p->A10==5)
    {
        p->flagini();
        p->gridini_outflow();
        pgc->flagfield(p);
        pgc->tpflagfield(p);
        makegrid(p,pgc);
        makegrid2D(p,pgc);

        pgc->ndflag_update(p);


        if(p->A10==3)
        pffg_driver();

        if(p->A10==4)
        nsewave_driver();

        if(p->A10==44)
        {
        makegrid_nhflow(p,pgc);
        nhflow_driver();
        }

        if(p->A10==5)
        cfd_driver();
    }

    // 2D Framework
    if(p->A10==2)
    {
        p->flagini2D();
        p->gridini2D();
        makegrid2D(p,pgc);
        sf_driver();
    }
}

void driver::cfd_driver()
{
    if(p->mpirank==0)
	cout<<"initialize fdm"<<endl;

	a=new fdm(p);

	aa=a;
    pgc->fdm_update(a);

    logic();
}

void driver::nsewave_driver()
{
    if(p->mpirank==0)
	cout<<"initialize fdm"<<endl;

	a=new fdm(p);

	aa=a;
    pgc->fdm_update(a);

    logic();
}

void driver::nhflow_driver()
{
    if(p->mpirank==0)
	cout<<"initialize fdm"<<endl;

	a=new fdm(p);

	aa=a;
    pgc->fdm_update(a);

    logic();
}

void driver::pfsg_driver()
{
    if(p->mpirank==0)
	cout<<"initialize fdm"<<endl;

    p->grid2Dsize();

    c=new fdm_fnpf(p);

    makegrid_fnpf_cds(p,pgc);

    logic_fnpf_sg();
}

void driver::pffg_driver()
{
    if(p->mpirank==0)
	cout<<"initialize fdm"<<endl;

    a=new fdm(p);

    aa=a;
    pgc->fdm_update(a);

    logic_fnpf_fg();
}

void driver::sf_driver()
{
    if(p->mpirank==0)
	cout<<"initialize fdm"<<endl;

    b=new fdm2D(p);
    bb=b;

    psflow = new sflow_f(p,b,pgc);

    makegrid2D_cds(p,pgc,b);

	psflow->start(p,b,pgc);
}

driver::~driver()
{
}
