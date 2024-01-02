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
--------------------------------------------------------------------*/
#include"position.h"
#include"lexer.h"

int position::posf_i(double xs)
{
    stop=0;

    is = 0;
    ie = p->knox;
    
    
    count=0;
    do{
    iloc = ihalf(is,ie);
    
    if(count%3==0)
    iloc+=1;
    
        // matching criterion
        if(xs<p->XP[iloc+marge] && xs>=p->XP[iloc-1+marge])
        {
            ii = iloc;
            
         stop=1;
         break;   
        }
        
        if(xs>=p->XP[iloc+marge] && xs<p->XP[iloc+1+marge])
        {
            ii = iloc+1;
   
         stop=1;
         break;   
        }
        
        // out of bounds
        if(xs<p->XP[0+marge])
        {
            ii = -1;

         stop=1;
         break;   
        }
        
        // out of bounds
        if(xs>p->XP[p->knox-1+marge])
        {
            ii = p->knox+1;
            
         stop=1;
         break;   
        }
        
        // further division
        if(xs<p->XP[iloc+marge] && xs<p->XP[iloc-1+marge])
        ie=iloc;
        
        if(xs>p->XP[iloc+marge] && xs>p->XP[iloc+1+marge])
        is=iloc;
        
        ++count;
    }while(stop==0 && count<1000);
    
    ii=MAX(ii,0);
    ii=MIN(ii,p->knox);
    
    return ii;
}

int position::posf_j(double ys)
{
    stop=0;
    
    js = 0;
    je = p->knoy;
    
    count=0;
    do{
    jloc = ihalf(js,je);
    
    if(count%3==0)
    jloc+=1;
    
        // out of bounds
        if(ys<p->YP[0+marge])
        {
            jj = -1;
            
            //cout<<"EXIT 0m"<<endl;
   
         stop=1;
         break;   
        }
        
        // out of bounds
        if(ys>p->YP[p->knoy-1+marge])
        {
            jj = p->knoy+1;
            
            //cout<<"EXIT 0p  "<<p->YP[p->knoy-1]<<endl;
   
         stop=1;
         break;   
        }
    
        // matching criterion
        if(ys<p->YP[jloc+marge] && ys>=p->YP[jloc-1+marge])
        {
            jj = jloc;
            
            //cout<<"EXIT 1"<<endl;
            
         stop=1;
         break;   
        }
        
        if(ys>=p->YP[jloc+marge] && ys<p->YP[jloc+1+marge])
        {
            jj = jloc+1;
            
            //cout<<"EXIT 2"<<endl;
   
         stop=1;
         break;   
        }
        
        // further divjsion
        if(ys<p->YP[jloc+marge] && ys<p->YP[jloc-1+marge])
        je=jloc;
        
        if(ys>p->YP[jloc+marge] && ys>p->YP[jloc+1+marge])
        js=jloc;
        
        
        ++count;
    }while(stop==0 && count<1000);
    
    jj=MAX(jj,0);
    jj=MIN(jj,p->knoy);
    
    return jj;
}



int position::posf_k(double zs)
{
    stop=0;

    ks = 0;
    ke = p->knoz;
    
    count=0;
    do{
    kloc = ihalf(ks,ke);
    
    if(count%3==0)
    kloc+=1;
    
        // out of bounds
        if(zs<p->ZP[0+marge])
        {
            kk = -1;
   
         stop=1;
         break;   
        }
        
        // out of bounds
        if(zs>p->ZP[p->knoz-1+marge])
        {
            kk = p->knoz+1;

         stop=1;
         break;   
        }
    
        // matching criterion
        if(zs<p->ZP[kloc+marge] && zs>=p->ZP[kloc-1+marge])
        {
            kk = kloc;
            
         stop=1;
         break;   
        }
        
        if(zs>=p->ZP[kloc+marge] && zs<p->ZP[kloc+1+marge])
        {
            kk = kloc+1;

         stop=1;
         break;   
        }
        
        // further divksion
        if(zs<p->ZP[kloc+marge] && zs<p->ZP[kloc-1+marge])
        ke=kloc;
        
        if(zs>p->ZP[kloc+marge] && zs>p->ZP[kloc+1+marge])
        ks=kloc;
        
        
        ++count;
    }while(stop==0 && count<1000);
    
    kk=MAX(kk,0);
    kk=MIN(kk,p->knoz);
    
    
    return kk;
}
