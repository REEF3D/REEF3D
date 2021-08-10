/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2021 Hans Bihs

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
#include"patchBC.h"

driver::driver(int& argc, char **argv)
{
	p = new lexer;
	pgc = new ghostcell(argc,argv,p);



    time_t now = time(0);
    char* timenow = ctime(&now);

    if(p->mpirank==0)
    {
    cout<<endl<<"REEF3D (c) 2008-2021 Hans Bihs"<<endl;
    cout<<endl<<":: Open-Source Hydrodynamics" <<endl;
    cout<< timenow << endl;
    cout<<endl<<"v_210802"<<"; "<<BRANCH<<"; "<<VERSION<<endl<<endl;
    }

	p->lexer_read(pgc);
	pgc->gcini(p);
    p->gridini(pgc);
    patchBC_logic();

    if(p->mpirank==0 && p->B90==1)
    cout<<endl<<"!!! IMPORTANT NOTICE: changed input for B91, B93 and B96. please check the manual !!!"<<endl<<endl;

    if(p->mpirank==0)
    {
    if(p->A10==2)
    cout<<endl<<"REEF3D::SFLOW" <<endl<<endl;

    if(p->A10==3)
    cout<<endl<<"REEF3D::FNPF" <<endl<<endl;

    if(p->A10==4)
    cout<<endl<<"REEF3D::PTF" <<endl<<endl;

    if(p->A10==5)
    cout<<endl<<"REEF3D::NSEWAVE"<<endl<<endl;

    if(p->A10==55)
    cout<<endl<<"REEF3D::NHFLOW"<<endl<<endl;

    if(p->A10==6)
    cout<<endl<<"REEF3D::CFD" <<endl<<endl;
    }

// 3D Framework
    // sigma grid
    if(p->A10==3)
    {
        p->flagini();
        p->gridini_patchBC();
        pgc->flagfield(p);
        pgc->tpflagfield(p);
        makegrid_fnpf(p,pgc);

        pgc->ndflag_update(p);

        fnpf_driver();
    }

    // fixed grid
    if(p->A10==4 || p->A10==5 || p->A10==55 || p->A10==6)
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

        if(p->A10==5)
        nsewave_driver();

        if(p->A10==55)
        {
        makegrid_nhflow(p,pgc);
        nhflow_driver();
        }

        if(p->A10==6)
        cfd_driver();
    }

    // 2D Framework
    if(p->A10==2)
    {
        p->flagini2D();
        p->gridini2D();
        makegrid2D(p,pgc);
        pBC->patchBC_ini(p,pgc);
        sflow_driver();
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

void driver::fnpf_driver()
{
    if(p->mpirank==0)
	cout<<"initialize fdm"<<endl;

    p->grid2Dsize();

    c=new fdm_fnpf(p);

    pgc->fdm_fnpf_update(c);

    makegrid_fnpf_cds(p,pgc);

    logic_fnpf();
}

void driver::ptf_driver()
{
    if(p->mpirank==0)
	cout<<"initialize fdm"<<endl;

    a=new fdm(p);

    aa=a;
    pgc->fdm_update(a);

    logic_ptf();
}

void driver::sflow_driver()
{
    if(p->mpirank==0)
	cout<<"initialize fdm"<<endl;

    b=new fdm2D(p);
    bb=b;

    psflow = new sflow_f(p,b,pgc,pBC);

    makegrid2D_cds(p,pgc,b);

	psflow->start(p,b,pgc);
}

driver::~driver()
{
}
