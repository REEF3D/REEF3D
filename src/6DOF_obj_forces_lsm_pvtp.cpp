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

#include"6DOF_obj.h"
#include"lexer.h"

void sixdof_obj::pvtp(lexer* p, int num)
{
    sprintf(name,"./REEF3D_CFD_6DOF/REEF3D-FB-%08i.pvtp",num);

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
    result<<"<PDataArray type=\"Float32\" Name=\"pressure\"/>\n";
    result<<"</PPointData>\n";

    char pname[100];
    for(n=0; n<p->M10; ++n)
    {
        sprintf(pname,"REEF3D-FB-%08i-%06i.vtp",num,n+1);
        result<<"<Piece Source=\""<<pname<<"\"/>\n";
    }

    result<<"</PPolyData>\n";
    result<<"</VTKFile>\n";

    result.close();
}
