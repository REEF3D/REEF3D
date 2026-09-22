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

#include"position.h"
#include"lexer.h"

int position::posf_sig(int ii, int jj, double zs)
{
    
    i = ii;
    j = jj;
    
    k = 0;
    int IJK_start = IJK;
    
    k = p->knoz-1;
    int IJK_end = IJK;
    
    
    stop=0;

    ks = 0;
    ke = p->knoz;
    
    count=0;
    do{
    kloc = ihalf(ks,ke);
    
    if(count%3==0)
    kloc+=1;
    
    k=kloc;
    
        // out of bounds
        if(zs<p->ZSP[IJK_start])
        {
            kk = -1;
   
         stop=1;
         break;   
        }
        
        // out of bounds
        if(zs>p->ZSP[IJK_end])
        {
            kk = p->knoz;

         stop=1;
         break;   
        }
    
        // matching criterion
        if(zs<p->ZSP[IJK] && zs>=p->ZSP[IJKm1])
        {
            kk = kloc-1;
            
         stop=1;
         break;   
        }
        
        if(zs>=p->ZSP[IJK] && zs<p->ZSP[IJKp1])
        {
            kk = kloc;

         stop=1;
         break;   
        }
        
        // further divksion
        if(zs<p->ZSP[IJK] && zs<p->ZSP[IJKm1])
        ke=kloc;
        
        if(zs>p->ZSP[IJK] && zs>p->ZSP[IJKp1])
        ks=kloc;
        
        
        ++count;
    }while(stop==0 && count<1000);
    
    kk=MAX(kk,-1);
    kk=MIN(kk,p->knoz);
        
    
    return kk;
}

