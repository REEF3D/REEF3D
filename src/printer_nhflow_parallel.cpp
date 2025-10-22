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

#include"printer_nhflow.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment.h"
#include"nhflow_turbulence.h"

void printer_nhflow::parallel(lexer *p, fdm_nhf *d, ghostcell* pgc, nhflow_turbulence *pnhfturb, sediment *psed, int num)
{
    outputFormat->parallelFileName(name,sizeof(name),"NHFLOW",num);

    ofstream result;
    result.open(name);

    outputFormat->beginningParallel(p,result);

    result<<"<PPointData>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"pressure\"/>\n";
    pnhfturb->name_ParaView_parallel(p,d,pgc,result);
    if(p->P74==1)
        result<<"<PDataArray type=\"Float32\" Name=\"omega_sig\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"elevation\"/>\n";
    if(p->P23==1)
        result<<"<PDataArray type=\"Float32\" Name=\"test\"/>\n";
    if(p->P110==1)
        result<<"<PDataArray type=\"Float32\" Name=\"Hs\"/>\n";
    if(p->P25==1)
        result<<"<PDataArray type=\"Float32\" Name=\"solid\"/>\n";
    if(p->P25==1 || p->P28==1)
        result<<"<PDataArray type=\"Float32\" Name=\"Heaviside\"/>\n";
    if(p->P26==1)
        result<<"<PDataArray type=\"Float32\" Name=\"ST_conc\"/>\n";
    if(p->P28==1)
        result<<"<PDataArray type=\"Float32\" Name=\"floating\"/>\n";
    result<<"</PPointData>\n";

    outputFormat->endingParallel(result,"NHFLOW",p->M10,num);

    result.close();
}
