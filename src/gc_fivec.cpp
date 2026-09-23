/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"ghostcell.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"sliceint.h"
#include<vector>

// (i,j,k) triplets of FLOOP cells that have a solid (flag7<0) neighbour in i
// (list2D) or in i or j (list3D). flag7 is static after the sigma grid is built.
static std::vector<int> fivec_list3D, fivec_list2D;
static bool fivec_listbuilt=false;

// flag7 is fixed after driver_makegrid_sigma, so the cells that can own a
// wall ghost cell are found once. The loops below visit only those cells, in
// the original FLOOP order, instead of sweeping the whole 3D field every call.
void ghostcell::fivec_buildlist(lexer *p)
{
    fivec_list3D.clear();
    fivec_list2D.clear();
    
    FLOOP
    {
        const bool xwall = (p->flag7[FIm1JK]<0 || p->flag7[FIp1JK]<0);
        const bool ywall = (p->flag7[FIJm1K]<0 || p->flag7[FIJp1K]<0);
        
        if(xwall || ywall)
        {
        fivec_list3D.push_back(i);
        fivec_list3D.push_back(j);
        fivec_list3D.push_back(k);
        }
        
        if(xwall)
        {
        fivec_list2D.push_back(i);
        fivec_list2D.push_back(j);
        fivec_list2D.push_back(k);
        }
    }
    
    fivec_listbuilt=true;
}

#define FIVEC_LOOP3D if(!fivec_listbuilt) fivec_buildlist(p); \
    for(size_t qq=0; qq<fivec_list3D.size(); qq+=3) \
    if((i=fivec_list3D[qq], j=fivec_list3D[qq+1], k=fivec_list3D[qq+2], true))

#define FIVEC_LOOP2D if(!fivec_listbuilt) fivec_buildlist(p); \
    for(size_t qq=0; qq<fivec_list2D.size(); qq+=3) \
    if((i=fivec_list2D[qq], j=fivec_list2D[qq+1], k=fivec_list2D[qq+2], true))

void ghostcell::fivec(lexer *p, double *f, sliceint &bc)
{	
    FIVEC_LOOP3D
    {  
        if(p->B98<3||bc(i-1,j)==0)
        if(p->flag7[FIm1JK]<0)
        {
        f[FIm1JK] = f[FIJK];
        f[FIm2JK] = f[FIJK];
        f[FIm3JK] = f[FIJK];
        }
        
        if(p->B98>=3&&bc(i-1,j)==1)
        if(p->flag7[FIm1JK]<0)
        {
        f[FIm1JK] = f[FIJK]-c->Uin[FIm1JK]*1.0*p->DXP[IM1];
        f[FIm2JK] = f[FIJK]-c->Uin[FIm1JK]*2.0*p->DXP[IM1];
        f[FIm3JK] = f[FIJK]-c->Uin[FIm1JK]*3.0*p->DXP[IM1];
        }
        
        if(p->B99<3||bc(i+1,j)==0)
        if(p->flag7[FIp1JK]<0)
        {
        f[FIp1JK] = f[FIJK];
        f[FIp2JK] = f[FIJK];
        f[FIp3JK] = f[FIJK];
        }
        
        if(p->B99>=3||bc(i+1,j)==2)
        if(p->flag7[FIp1JK]<0)
        {
        f[FIp1JK] = f[FIJK]+c->Uin[FIp1JK]*1.0*p->DXP[IM1];
        f[FIp2JK] = f[FIJK]+c->Uin[FIp1JK]*2.0*p->DXP[IM1];
        f[FIp3JK] = f[FIJK]+c->Uin[FIp1JK]*3.0*p->DXP[IM1];
        }
        
        
        if(p->flag7[FIJm1K]<0)
        {
        f[FIJm1K] = f[FIJK];
        f[FIJm2K] = f[FIJK];
        f[FIJm3K] = f[FIJK];
        }
        
        if(p->flag7[FIJp1K]<0)
        {
        f[FIJp1K] = f[FIJK];
        f[FIJp2K] = f[FIJK];
        f[FIJp3K] = f[FIJK];
        }
    }
}

void ghostcell::fivec2D(lexer *p, double *f, sliceint &bc)
{	
    FIVEC_LOOP2D
    {
        if(p->B98<3||bc(i-1,j)==0)
        if(p->flag7[FIm1JK]<0)
        {
        f[FIm1JK] = f[FIJK];
        f[FIm2JK] = f[FIJK];
        f[FIm3JK] = f[FIJK];
        }
        
        if(p->B98>=3&&bc(i-1,j)==1)
        if(p->flag7[FIm1JK]<0)
        {
        f[FIm1JK] = f[FIJK]-c->Uin[FIm1JK]*1.0*p->DXP[IM1];
        f[FIm2JK] = f[FIJK]-c->Uin[FIm1JK]*2.0*p->DXP[IM1];
        f[FIm3JK] = f[FIJK]-c->Uin[FIm1JK]*3.0*p->DXP[IM1];
        }
        
        if(p->B99<3||bc(i+1,j)==0)
        if(p->flag7[FIp1JK]<0)
        {
        f[FIp1JK] = f[FIJK];
        f[FIp2JK] = f[FIJK];
        f[FIp3JK] = f[FIJK];
        }
        
        if(p->B99>=3||bc(i+1,j)==2)
        if(p->flag7[FIp1JK]<0)
        {
        f[FIp1JK] = f[FIJK]+c->Uin[FIp1JK]*1.0*p->DXP[IM1];
        f[FIp2JK] = f[FIJK]+c->Uin[FIp1JK]*2.0*p->DXP[IM1];
        f[FIp3JK] = f[FIJK]+c->Uin[FIp1JK]*3.0*p->DXP[IM1];
        }
    }
}


void ghostcell::fivec_vel(lexer *p, double *f, sliceint &bc)
{	
    FIVEC_LOOP3D
    {  
        if(p->flag7[FIm1JK]<0)
        {
        f[FIm1JK] = f[FIJK];
        f[FIm2JK] = f[FIJK];
        f[FIm3JK] = f[FIJK];
        }
          
        if(p->flag7[FIp1JK]<0)
        {
        f[FIp1JK] = f[FIJK];
        f[FIp2JK] = f[FIJK];
        f[FIp3JK] = f[FIJK];
        }
        
        if(p->flag7[FIJm1K]<0)
        {
        f[FIJm1K] = f[FIJK];
        f[FIJm2K] = f[FIJK];
        f[FIJm3K] = f[FIJK];
        }
        
        if(p->flag7[FIJp1K]<0)
        {
        f[FIJp1K] = f[FIJK];
        f[FIJp2K] = f[FIJK];
        f[FIJp3K] = f[FIJK];
        }
    }
}

void ghostcell::fivec2D_vel(lexer *p, double *f, sliceint &bc)
{	
    FIVEC_LOOP2D
    {
        if(p->flag7[FIm1JK]<0)
        {
        f[FIm1JK] = f[FIJK];
        f[FIm2JK] = f[FIJK];
        f[FIm3JK] = f[FIJK];
        }
        
        if(p->flag7[FIp1JK]<0)
        {
        f[FIp1JK] = f[FIJK];
        f[FIp2JK] = f[FIJK];
        f[FIp3JK] = f[FIJK];
        }
    }
}


