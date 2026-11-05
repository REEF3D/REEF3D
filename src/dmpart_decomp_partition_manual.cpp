/*--------------------------------------------------------------------
DIVEMesh
Copyright 2008-2026 Hans Bihs

This file is part of REEF3D (port of the DIVEMesh decomposition, dmpart).

DIVEMesh is free software; you can redistribute it and/or modify it
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

#include "dmpart_decomp.h"

namespace dmpart
{

void decomp::partition_manual(lexer* p,dive* a)
{
    a->mx=p->M30_x;
    a->my=p->M30_y;
    a->mz=p->M30_z;

    a->xorig[0]=p->xmin;
    a->yorig[0]=p->ymin;
    a->zorig[0]=p->zmin;
    a->xnode[0]=0;
    a->ynode[0]=0;
    a->znode[0]=0;

    a->xorig[a->mx]=p->xmax;
    a->yorig[a->my]=p->ymax;
    a->zorig[a->mz]=p->zmax;
    a->xnode[a->mx]=p->knox;
    a->ynode[a->my]=p->knoy;
    a->znode[a->mz]=p->knoz;

    // partition
    for(n=1;n<a->mx;n++)
    {
        a->xnode[n]=p->M31[n-1];
        a->xorig[n]=p->XN[a->xnode[n]+marge];
    }

    for(n=1;n<a->my;n++)
    {
        a->ynode[n]=p->M32[n-1];
        a->yorig[n]=p->YN[a->ynode[n]+marge];
    }

    for(n=1;n<a->mz;n++)
    {
        a->znode[n]=p->M33[n-1];
        a->zorig[n]=p->ZN[a->znode[n]+marge];
    }

    // the subdomain of every cell, as partition() sets it: the surfaces, the neighbours and the counts of the
    // decomposition are taken from it (without it the cells kept the subdomain of no partition)
    MALOOP
    {
        a->subgrid(i,j,k)=-1;
        a->subslice(i,j)=-1;
    }

    NLOOP
    SUBLOOP
    {
        a->subgrid(i,j,k)=PARANUM;
        a->subslice(i,j)=PARANUM;
    }
}

void decomp::partition_input_check(lexer* p, dive* a)
{
    if(p->M20<1 || p->M20>5)
    {
        cout<<endl<<"!!! M 20 "<<p->M20<<" is not a decomposition method (1-5) !!!"<<endl<<endl;
        exit(1);
    }

    if(p->M10<1)
    {
        cout<<endl<<"!!! M 10 must be at least 1 !!!"<<endl<<endl;
        exit(1);
    }

    if(p->M20==3 || p->M20==4 || (p->M20==5 && (p->M30_x>0 || p->M30_y>0 || p->M30_z>0)))
    {
        if(p->M30_x<1 || p->M30_y<1 || p->M30_z<1 || p->M30_x*p->M30_y*p->M30_z!=p->M10)
        {
            cout<<endl<<"!!! M 30 "<<p->M30_x<<" "<<p->M30_y<<" "<<p->M30_z
                <<" must be at least 1 in each direction and multiply to M 10 "<<p->M10<<" !!!"<<endl<<endl;
            exit(1);
        }

        if(p->M30_x>a->knox || p->M30_y>a->knoy || p->M30_z>a->knoz)
        {
            cout<<endl<<"!!! M 30 "<<p->M30_x<<" "<<p->M30_y<<" "<<p->M30_z
                <<" asks for more subdomains than cells ("<<a->knox<<" x "<<a->knoy<<" x "<<a->knoz<<") !!!"<<endl<<endl;
            exit(1);
        }
    }
}

void decomp::partition_verify(lexer* p, dive* a)
{
    const char dir[3] = {'x','y','z'};
    int *node[3] = {a->xnode,a->ynode,a->znode};
    const int m[3] = {a->mx,a->my,a->mz};
    const int kn[3] = {a->knox,a->knoy,a->knoz};
    int error=0;

    for(int d=0;d<3;++d)
    {
        if(node[d][0]!=0 || node[d][m[d]]!=kn[d])
        {
            cout<<"!!! "<<dir[d]<<"-partition does not span the grid: nodes "<<node[d][0]<<" ... "<<node[d][m[d]]
                <<", expected 0 ... "<<kn[d]<<" !!!"<<endl;
            error=1;
        }

        for(int n=1;n<=m[d];++n)
        {
            const int width = node[d][n]-node[d][n-1];

            if(width<1)
            {
                cout<<"!!! "<<dir[d]<<"-subdomain "<<n<<" has "<<width<<" cells (nodes "<<node[d][n-1]<<" to "<<node[d][n]<<") !!!"<<endl;
                error=1;
            }
            else if(width<3 && m[d]>1)
            {
                cout<<"warning: "<<dir[d]<<"-subdomain "<<n<<" is only "<<width<<" cells wide"<<endl;
                ddout<<"warning: "<<dir[d]<<"-subdomain "<<n<<" is only "<<width<<" cells wide"<<endl;
            }
        }
    }

    if(error==1)
    {
        cout<<endl<<"!!! invalid domain decomposition, see above !!!"<<endl<<endl;
        exit(1);
    }
}
}
