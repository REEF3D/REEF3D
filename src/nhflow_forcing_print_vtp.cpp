/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"nhflow_forcing.h"
#include"lexer.h"
#include<sys/stat.h>
#include<fstream>

void nhflow_forcing::print_vtp(lexer *p)
{
    int offset[100];
    int n=0;
    offset[n]=0;
    ++n;
    offset[n]=offset[n-1]+sizeof(float)*tricount*3*3 + sizeof(int);
    ++n;
    offset[n]=offset[n-1]+sizeof(int)*tricount*3 + sizeof(int);
    ++n;
    //---------------------------------------------

    mkdir("./REEF3D_NHFLOW_FORCING_VTP", 0777);

    char path[300];
    sprintf(path,"./REEF3D_NHFLOW_FORCING_VTP/REEF3D-NHFLOW-FORCING.vtp");

    std::ofstream result;
    result.open(path, std::ios::binary);

    //---------------------------------------------

    vtp3D::beginning(p, result, tricount*3, 0, 0, 0, tricount);

    n=0;
    result<<"<Points>\n";
    result<<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"</Points>\n";

    vtp3D::polys(result, offset, n);

    vtp3D::ending(result);

    //----------------------------------------------------------------------------
    int q,iin,m;
    float ffn;

    //  XYZ
    iin=4*tricount*3*3;
    result.write((char*)&iin, sizeof(int));
    for(m=0;m<tricount;++m)
    for(q=0;q<3;++q)
    {
        ffn=tri_x[m][q];
        result.write((char*)&ffn, sizeof(float));

        ffn=tri_y[m][q];
        result.write((char*)&ffn, sizeof(float));

        ffn=tri_z[m][q];
        result.write((char*)&ffn, sizeof(float));
    }

    //  Connectivity POLYGON
    int count=0;
    iin=4*tricount*3;
    result.write((char*)&iin, sizeof(int));
    for(m=0;m<tricount;++m)
    for(q=0;q<3;++q)
    {
        iin=count;
        result.write((char*)&iin, sizeof(int));
        ++count;
    }

    //  Offset of Connectivity
    iin=4*tricount;
    result.write((char*)&iin, sizeof(int));
    iin=0;
    for(m=0;m<tricount;++m)
    {
        iin+= 3;
        result.write((char*)&iin, sizeof(int));
    }

    vtp3D::footer(result);

    result.close();
}
