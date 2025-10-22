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

#include"topo_vtp.h"
#include"lexer.h"
#include"sediment.h"

void topo_vtp::pvtp(lexer* p, sediment *psed, int num)
{
    sprintf(name,"./REEF3D_CFD_Topo/REEF3D-CFD-Topo-%08i.pvtp",num);

    ofstream result;
    result.open(name);

    result<<"<?xml version=\"1.0\"?>\n";
    result<<"<VTKFile type=\"PPolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    result<<"<PPolyData GhostLevel=\"0\">\n";

    result<<"<PPoints>\n";
    result<<"<PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n";
    result<<"</PPoints>\n";

    result<<"<PPointData>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\"/>\n";
    result<<"<PDataArray type=\"Float32\" Name=\"elevation\"/>\n";
    if(p->P76==1)
        psed->name_ParaView_parallel_bedload(p,result);
    if(p->P77==1)
        psed->name_ParaView_parallel_parameter1(p,result);
    if(p->P78==1)
        psed->name_ParaView_parallel_parameter2(p,result);
    if(p->P79>=1)
        psed->name_ParaView_parallel_bedshear(p,result);
    result<<"</PPointData>\n";

    result<<"<Polys>\n";
    result<<"<DataArray type=\"Int32\" Name=\"connectivity\"/>\n";
    result<<"<DataArray type=\"Int32\" Name=\"offsets\"/>\n";
    result<<"</Polys>\n";

    char pname[100];
    for(n=0; n<p->M10; ++n)
    {
        sprintf(pname,"REEF3D-CFD-Topo-%08i-%06i.vtp",num,n+1);
        result<<"<Piece Source=\""<<pname<<"\"/>\n";
    }

    result<<"</PPolyData>\n";
    result<<"</VTKFile>\n";

    result.close();
}
