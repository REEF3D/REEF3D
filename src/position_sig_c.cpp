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

int position::posc_sig(int ii, int jj, double zs)
{
    i = ii;
    j = jj;
    
    i = MAX(i,0);
	i = MIN(i,p->knox-1);
    
	j = MAX(j,0);
	j = MIN(j,p->knoy-1);
    
    k = 0;
    int FIJK_start = FIJK;
    
    k = p->knoz;
    int FIJK_end = FIJK;
    
    stop=0;
    
    ks = 0;
    ke = p->knoz+1;

    count=0;
    do{
    kloc = ihalf(ks,ke);
    
    if(count%3==0)
    kloc+=1;
    
    k=kloc;
    
        
        // out of bounds
        if(zs<p->ZSN[FIJK_start])
        {
            kk = -1;
            
            //cout<<"EXIT 0m"<<endl;
   
         stop=1;
         break;   
        }
        
        // out of bounds
        if(zs>p->ZSN[FIJK_end])
        {
            kk = p->knoz+1;
            
            //cout<<"EXIT 0p "<<" zs: "<<zs<<" p->ZSP[IJK_end]: "<<p->ZSP[IJK_end]<<" i: "<<i<<" j: "<<j<<endl;
   
         stop=1;
         break;   
        }
        
        // matching criterion
        if(zs<p->ZSN[FIJK] && zs>=p->ZSN[FIJKm1] && stop==0)
        {
            kk = kloc-1;
            
            //cout<<"EXIT 1"<<endl;
            
         stop=1;
         break;   
        }
        
        if(zs>=p->ZSN[FIJK] && zs<p->ZSN[FIJKp1] && stop==0)
        {
            kk = kloc;
            
            //cout<<"EXIT 2"<<endl;
   
         stop=1;
         break;   
        }
        
        // further divksion
        if(zs<p->ZSN[FIJK] && zs<p->ZSN[FIJKm1])
        ke=kloc;
        
        if(zs>p->ZSN[FIJK] && zs>p->ZSN[FIJKp1])
        ks=kloc;
        
        
        ++count;
    }while(stop==0 && count<1000);
    
    kk=MAX(kk,-1);
    kk=MIN(kk,p->knoz+1);
    
    
    return kk;
}



