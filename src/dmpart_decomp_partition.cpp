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

void decomp::partition(lexer* p, dive* a)
{
    int xid,yid,zid;
    int xdiff,ydiff,zdiff;

    xid=yid=zid=0;
    xdiff=ydiff=zdiff=0;

    posx=int(a->knox/a->mx);

    if(a->knox%a->mx!=0 && a->mx>1)
    {
        xid=1;
        xdiff = a->knox - posx*a->mx;
    }

    posy=int(a->knoy/a->my);

    if(a->knoy%a->my!=0 && a->my>1)
    {
        yid=1;
        ydiff = a->knoy - posy*a->my;
    }

    posz=int(a->knoz/a->mz);

    if(a->knoz%a->mz!=0 && a->mz>1)
    {
        zid=1;
        zdiff = a->knoz - posz*a->mz;
    }

    ddout<<endl<<"posx: "<<posx<<"  posy: "<<posy<<"  posz: "<<posz<<endl;
    ddout<<endl<<"xid: "<<xid<<"  yid: "<<yid<<"  zid: "<<zid<<endl;
    ddout<<endl<<"xdiff: "<<xdiff<<"  ydiff: "<<ydiff<<"  zdiff: "<<zdiff<<endl;

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

    for(n=1;n<a->mx;n++)
    {
        if(n>=xdiff)
        xid=0;

        a->xnode[n]=a->xnode[n-1] + (posx+xid);
        a->xorig[n]=p->XN[a->xnode[n]+marge];
    }

    for(n=1;n<a->my;n++)
    {
        if(n>=ydiff)
        yid=0;

        a->ynode[n]=a->ynode[n-1] + (posy+yid);
        a->yorig[n]=p->YN[a->ynode[n]+marge];
    }

    for(n=1;n<a->mz;n++)
    {
        if(n>=zdiff)
        zid=0;

        a->znode[n]=a->znode[n-1] + (posz+zid);
        a->zorig[n]=p->ZN[a->znode[n]+marge];
    }

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

    for(aa=0;aa<=a->mx;++aa)
    ddout<<"xcount"<<aa<<" : "<<xcount[aa]<<"  xnode: "<<a->xnode[aa]<<"  xorig: "<<a->xorig[aa]<<endl;
}

void decomp::partition_voidcheck(lexer* p, dive* a)
{
    int check=1;

    for(q=0;q<p->M10;++q)
    subcell[q]=0;

    q=0;
    NLOOP
    {
        SUBLOOP
        if(a->flag(i,j,k)>0)
        ++subcell[q];

        ++q;
    }

    for(q=0;q<p->M10;++q)
    if(subcell[q]==0)
    check=0;

    if(check==0)
    {
        ddout<<"Void partition detected !!!"<<endl;

        if(maindir==1)
        {
            p->M11=1;
            p->M12=0;
            p->M13=0;
        }
        else if(maindir==2)
        {
            p->M11=0;
            p->M12=1;
            p->M13=0;
        }
        else if(maindir==3)
        {
            p->M11=0;
            p->M12=0;
            p->M13=1;
        }

        costfunc(p,a);
        partition(p,a);
        partition_correction(p,a);
    }
}

void decomp::rank(lexer *p, dive *a, int &dc_i, int &dc_j, int &dc_k, int &rank1, int &rank2, int &rank3)
{
    if(dc_i>=dc_j)
    {
        if(dc_i>=dc_k)
        {
            rank1=1;

            if(dc_j>=dc_k)
            {
                rank2=2;
                rank3=3;
            }
            else if(dc_j<dc_k)
            {
                rank2=3;
                rank3=2;
            }
        }
        else if(dc_i<dc_k)
        {
            rank1=3;
            rank2=1;
            rank3=2;
        }
   }
   else if(dc_i<dc_j)
   {
        if(dc_i>=dc_k)
        {
            rank1=2;
            rank2=1;
            rank3=3;
        }
        else if(dc_i<dc_k)
        {
            rank3=1;

            if(dc_j>=dc_k)
            {
                rank1=2;
                rank2=3;
            }
            else if(dc_j<dc_k)
            {
                rank1=3;
                rank2=2;
            }
        }
   }
}
}
