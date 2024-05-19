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
--------------------------------------------------------------------*/#include"ghostcell.h"
#include"lexer.h"

void ghostcell::gcxslupdate(lexer* p)
{
    //cout<<p->mpirank<<" . "<<p->knox<<" "<<p->origin_i<<" "<<p->gknox<<endl;

    for(n=0;n<p->gcslpara1_count;++n)
    {
    i=p->gcslpara1[n][0];
    j=p->gcslpara1[n][1];

        // 1
        p->gcslpara1[n][3]=1;
        
        //if(i + p->origin_i >= p->gknox-1)
        //p->gcslpara1[n][3]=0;
        
        if(p->flagslice1[Im1J]<0 || p->flagslice1[IJ]<0)
        p->gcslpara1[n][3]=p->Y71;
        
        // 2
        p->gcslpara1[n][4]=1;
        
        if(j + p->origin_j >= p->gknoy-1)
        p->gcslpara1[n][4]=0;
        
        if(p->flagslice2[Im1J]<0 || p->flagslice2[IJ]<0)
        p->gcslpara1[n][4]=p->Y72;
        
        // 4
        p->gcslpara1[n][6]=1;
        
        if(p->flagslice4[Im1J]<0 || p->flagslice4[IJ]<0)
        p->gcslpara1[n][6]=p->Y74;
    }

    for(n=0;n<p->gcslpara2_count;++n)
    {
    i=p->gcslpara2[n][0];
    j=p->gcslpara2[n][1];
        
        // 1
        p->gcslpara2[n][3]=1;
        
        if(i + p->origin_i >= p->gknox-1)
        p->gcslpara2[n][3]=0;
        
        if(p->flagslice1[IJp1]<0 || p->flagslice1[IJ]<0)
        p->gcslpara2[n][3]=p->Y71;
        
        // 2
        p->gcslpara2[n][4]=1;
        
        //if(j + p->origin_j >= p->gknoy-1)
        //p->gcslpara2[n][4]=0;
        
        if(p->flagslice2[IJp1]<0 || p->flagslice2[IJ]<0)
        p->gcslpara2[n][4]=p->Y71;
        
        // 4
        p->gcslpara2[n][6]=1;
        
        if(p->flagslice4[IJp1]<0 || p->flagslice4[IJ]<0)
        p->gcslpara2[n][6]=p->Y74;
    }

    for(n=0;n<p->gcslpara3_count;++n)
    {
    i=p->gcslpara3[n][0];
    j=p->gcslpara3[n][1];

        // 1
        p->gcslpara3[n][3]=1;
        
        if(i + p->origin_i >= p->gknox-1)
        p->gcslpara3[n][3]=0;
        
        if(p->flagslice1[Im1J]<0 || p->flagslice1[IJ]<0)
        p->gcslpara3[n][3]=p->Y71;
        
        // 2
        p->gcslpara3[n][4]=1;
        
        //if(j + p->origin_j >= p->gknoy-1)
        //p->gcslpara3[n][4]=0;
        
        if(p->flagslice2[IJm1]<0 || p->flagslice2[IJ]<0)
        p->gcslpara3[n][4]=p->Y71;
        
        // 4
        p->gcslpara3[n][6]=1;
        
        if(p->flagslice4[IJm1]<0 || p->flagslice4[IJ]<0)
        p->gcslpara3[n][6]=p->Y74;
    }    
    
    for(n=0;n<p->gcslpara4_count;++n)
    {
    i=p->gcslpara4[n][0];
    j=p->gcslpara4[n][1];

        // 1
        p->gcslpara4[n][3]=1;
        
        //if(i + p->origin_i >= p->gknox-1)
        //p->gcslpara4[n][3]=0;
        
        if(p->flagslice1[Ip1J]<0 || p->flagslice1[IJ]<0)
        p->gcslpara4[n][3]=p->Y71;
        
        // 2
        p->gcslpara4[n][4]=1;
        
        if(j + p->origin_j >= p->gknoy-1)
        p->gcslpara4[n][4]=0;
        
        if(p->flagslice2[Ip1J]<0 || p->flagslice2[IJ]<0)
        p->gcslpara4[n][4]=p->Y71;
        
        // 4
        p->gcslpara4[n][6]=1;
        
        if(p->flagslice4[Ip1J]<0 || p->flagslice4[IJ]<0)
        p->gcslpara4[n][6]=p->Y74;
        
    }

   
}

