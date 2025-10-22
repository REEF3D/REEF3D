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

#include"print_porous.h"
#include"lexer.h"
#include<sys/stat.h>
#include<sys/types.h>

void print_porous::print_vtp(lexer *p)
{
    // Create Folder
    if(p->mpirank==0)
        mkdir("./REEF3D_CFD_Porous",0777);

    char name[100];
    sprintf(name,"./REEF3D_CFD_Porous/REEF3D_Porous-Object.vtp");

    ofstream result;
    result.open(name, ios::binary);

    int offset[200];
    int n = 0;
    offset[n]=0;
    ++n;
    offset[n]=offset[n-1]+sizeof(float)*vertice_num*3 + sizeof(int);
    ++n;
    offset[n]=offset[n-1]+sizeof(int)*polygon_sum + sizeof(int);
    ++n;
    offset[n]=offset[n-1]+sizeof(int)*polygon_num + sizeof(int);
    ++n;
    //---------------------------------------------

    vtp3D::beginning(p,result,vertice_num,0,0,0,polygon_num);

    n=0;
    vtp3D::points(result,offset,n);

    result<<"<Polys>\n";
    result<<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"appended\" offset=\""<<offset[n]<<"\"/>\n";
    ++n;
    result<<"</Polys>\n";

    vtp3D::ending(result);

    //----------------------------------------------------------------------------

    float ffn;
    int iin;
    //  XYZ
    iin=sizeof(float)*vertice_num*3;
    result.write((char*)&iin, sizeof(int));
    for(n=0;n<vertice_num;++n)
    {
        ffn=vertice[n][0];
        result.write((char*)&ffn, sizeof(float));

        ffn=vertice[n][1];
        result.write((char*)&ffn, sizeof(float));

        ffn=vertice[n][2];
        result.write((char*)&ffn, sizeof(float));
    }

    //  Connectivity POLYGON
    iin=sizeof(int)*polygon_sum;
    result.write((char*)&iin, sizeof(int));
    for(n=0;n<polygon_num;++n)
    for(int q=0;q<numvert[n];++q)
    {
        iin=polygon[n][q];
        result.write((char*)&iin, sizeof(int));
    }

    //  Offset of Connectivity
    iin=sizeof(int)*polygon_num;
    result.write((char*)&iin, sizeof(int));
    iin=0;
    for(n=0;n<polygon_num;++n)
    {
        iin+=+ numvert[n];
        result.write((char*)&iin, sizeof(int));
    }

    vtp3D::footer(result);

    result.close();
}
