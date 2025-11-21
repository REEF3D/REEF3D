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

#include"printer_fnpf.h"
#include"lexer.h"

void printer_fnpf::parallel(lexer *p, int num)
{
    outputFormat->parallelFileName(name,sizeof(name),"FNPF",num);

    ofstream result;
    result.open(name);

    outputFormat->beginningParallel(p,result);

    result<<"<PPointData>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"Fi\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"elevation\"/>\n";
    if(p->P23==1)
        result<<"<PDataArray type=\"Float32\" Name=\"test\"/>\n";
    if(p->P110==1)
        result<<"<PDataArray type=\"Float32\" Name=\"Hs\"/>\n";
    if(p->P25==1)
        result<<"<PDataArray type=\"Float32\" Name=\"solid\"/>\n";
    result<<"</PPointData>\n";

    outputFormat->endingParallel(result,"FNPF",p->M10,num);

    result.close();
}
