/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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
#include"expdata.h"
#include"concentration.h"
#include"multiphase.h"
#include"sediment.h"
#include"print_averaging.h"

void printer_CFD::parallel(lexer* p, fdm* a, ghostcell* pgc, turbulence *pturb, heat *pheat, expdata *pdata, concentration *pconc, multiphase *pmp, sediment *psed, int num)
{
    outputFormat->parallelFileName(name,sizeof(name),"CFD",num);

    ofstream result;
    result.open(name);

    outputFormat->beginningParallel(p,result);

    result<<"<PPointData>\n";

    result<<"<PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>\n";

    pmean->name_ParaView_parallel(p,result);

    result<<"<PDataArray type=\"Float32\" Name=\"pressure\"/>\n";

    pturb->name_ParaView_parallel(p,result);

    result<<"<PDataArray type=\"Float32\" Name=\"eddyv\"/>\n";

    result<<"<PDataArray type=\"Float32\" Name=\"phi\"/>\n";

    pheat->name_ParaView_parallel(p,result);

    pmp->name_ParaView_parallel(p,result);

    pvort->name_ParaView_parallel(p,result);

    pdata->name_ParaView_parallel(p,result);

    pconc->name_ParaView_parallel(p,result);

    if(p->P24==1 && p->F300==0)
        result<<"<PDataArray type=\"Float32\" Name=\"rho\"/>\n";

    if(p->P71==1)
        result<<"<PDataArray type=\"Float32\" Name=\"viscosity\"/>\n";

    if(p->P72==1)
        result<<"<PDataArray type=\"Float32\" Name=\"VOF\"/>\n";

    if(p->A10==4)
        result<<"<PDataArray type=\"Float32\" Name=\"Fi\"/>\n";

    if(p->P26==1)
        result<<"<PDataArray type=\"Float32\" Name=\"ST_conc\"/>\n";

    if(p->P27==1)
        result<<"<PDataArray type=\"Float32\" Name=\"topo\"/>\n";

    if(p->P76==1)
        psed->name_ParaView_parallel_bedload(p,result);

    if(p->P77==1)
        psed->name_ParaView_parallel_parameter1(p,result);

    if(p->P78==1)
        psed->name_ParaView_parallel_parameter2(p,result);

    if(p->P79>=1)
        psed->name_ParaView_parallel_bedshear(p,result);

    if(p->P23==1)
        result<<"<PDataArray type=\"Float32\" Name=\"test\"/>\n";

    result<<"<PDataArray type=\"Float32\" Name=\"elevation\"/>\n";

    if(p->P25==1)
        result<<"<PDataArray type=\"Float32\" Name=\"solid\"/>\n";

    if(p->P28==1)
        result<<"<PDataArray type=\"Float32\" Name=\"floating\"/>\n";

    if(p->P29==1)
        result<<"<PDataArray type=\"Float32\" Name=\"walldist\"/>\n";

    result<<"</PPointData>\n";

    result<<"<PCellData>\n";

    if(p->P72==1)
        result<<"<PDataArray type=\"Float32\" Name=\"VOF_C\"/>\n";

    result<<"</PCellData>\n";

    outputFormat->endingParallel(result,"CFD",p->M10,num);

    result.close();
}
