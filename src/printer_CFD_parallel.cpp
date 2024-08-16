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

#include"printer_CFD.h"
#include<string>
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"turbulence.h"
#include"heat.h"
#include"vorticity.h"
#include"data.h"
#include"concentration.h"
#include"multiphase.h"
#include"sediment.h"
#include"print_averaging.h"

void printer_CFD::parallelData(fdm* a, lexer* p, ghostcell* pgc, turbulence *pturb, heat *pheat, data *pdata, concentration *pconc, multiphase *pmp, sediment *psed)
{
    int num=0;

    if(p->P15==1)
    num = p->printcount;

    if(p->P15==2)
    num = p->count;

	outputFormat->parallelFileName(name, "CFD", num);

	ofstream result;
	result.open(name);

	if(result.is_open())
	{
		outputFormat->beginningParallel(p,result);

		result<<"<PPointData>"<<endl;
		result<<"\t<PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>"<<endl;
		
		pmean->name_pvtk(p,a,pgc,result);

		result<<"	<PDataArray type=\"Float32\" Name=\"pressure\"/>"<<endl;

		pturb->name_pvtk(p,a,pgc,result);

		result<<"\t<PDataArray type=\"Float32\" Name=\"eddyv\"/>"<<endl;
		result<<"	<PDataArray type=\"Float32\" Name=\"phi\"/>"<<endl;

		pheat->name_pvtk(p,a,pgc,result);
		
		pmp->name_pvtk(p,a,pgc,result);

		pvort->name_pvtk(p,a,pgc,result);

		pdata->name_pvtk(p,a,pgc,result);

		pconc->name_pvtk(p,a,pgc,result);

		if(p->P24==1 && p->F300==0)
		result<<"	<PDataArray type=\"Float32\" Name=\"rho\"/>"<<endl;

		if(p->P71==1)
		result<<"	<PDataArray type=\"Float32\" Name=\"viscosity\"/>"<<endl;
		
		if(p->P72==1)
		result<<"	<PDataArray type=\"Float32\" Name=\"VOF\"/>"<<endl;
		
		if(p->A10==4)
		result<<"	<PDataArray type=\"Float32\" Name=\"Fi\"/>"<<endl;

		if(p->P26==1)
		{
		result<<"	<PDataArray type=\"Float32\" Name=\"ST_conc\"/>"<<endl;
		}

		if(p->P27==1)
		result<<"	<PDataArray type=\"Float32\" Name=\"topo\"/>"<<endl;
		
		if(p->P76==1)
		psed->name_pvtk_bedload(p,pgc,result);
		
		if(p->P77==1)
		psed->name_pvtk_parameter1(p,pgc,result);

		if(p->P78==1)
		psed->name_pvtk_parameter2(p,pgc,result);

		if(p->P79>=1)
		psed->name_pvtk_bedshear(p,pgc,result);

		if(p->P23==1)
		result<<"	<PDataArray type=\"Float32\" Name=\"test\"/>"<<endl;

		result<<"	<PDataArray type=\"Float32\" Name=\"elevation\"/>"<<endl;

		if(p->P25==1)
		result<<"	<PDataArray type=\"Float32\" Name=\"solid\"/>"<<endl;

		if(p->P28==1)
		result<<"	<PDataArray type=\"Float32\" Name=\"floating\"/>"<<endl;

		if(p->P29==1)
		result<<"	<PDataArray type=\"Float32\" Name=\"walldist\"/>"<<endl;

		result<<"</PPointData>"<<endl;

		outputFormat->endingParallel(result,"CFD",p->M10,num);

		result.close();
	}
	else
	cout<<"Failed to open output file."<<endl;
}