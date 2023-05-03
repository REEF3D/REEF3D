/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

int position::posc_i(double xs)
{
    stop=0;
    count=0;

    is = 0;
    ie = p->knox+1;

    count=0;
    do{
    iloc = ihalf(is,ie);
    
    if(count%3==0)
    iloc+=1;
    
        // out of bounds
        if(xs<p->XN[0+marge])
        {
            ii = -1;

         stop=1;
         break;   
        }
        
        // out of bounds
        if(xs>p->XN[p->knox+marge])
        {
            ii = p->knox+1;
            
         stop=1;
         break;   
        }
        
        // matching criterion
        if(xs<p->XN[iloc+marge] && xs>=p->XN[iloc-1+marge] && stop==0)
        {
            ii = iloc-1;
            
         stop=1;
         break;   
        }
        
        if(xs>=p->XN[iloc+marge] && xs<p->XN[iloc+1+marge] && stop==0)
        {
            ii = iloc;
            
         stop=1;
         break;   
        }
        
        // further division
        if(xs<p->XN[iloc+marge] && xs<p->XN[iloc-1+marge])
        ie=iloc;
        
        if(xs>p->XN[iloc+marge] && xs>p->XN[iloc+1+marge])
        is=iloc;
        
        ++count;
    }while(stop==0 && count<1000);
    

    ii=MAX(ii,-1);
    ii=MIN(ii,p->knox+1);

    return ii;
}

int position::posc_j(double ys)
{
    stop=0;
    
    js = 0;
    je = p->knoy+1;

    count=0;
    do{
    jloc = ihalf(js,je);
    
    if(count%3==0)
    jloc+=1;
        
        // out of bounds
        if(ys<p->YN[0+marge])
        {
            jj = -1;
            
            //cout<<"EXIT 0m"<<endl;
   
         stop=1;
         break;   
        }
        
        // out of bounds
        if(ys>p->YN[p->knoy+marge])
        {
            jj = p->knoy+1;
            
            //cout<<"EXIT 0p"<<endl;
   
         stop=1;
         break;   
        }
        
        // matching criterion
        if(ys<p->YN[jloc+marge] && ys>=p->YN[jloc-1+marge] && stop==0)
        {
            jj = jloc-1;
            
    
         stop=1;
         break;   
        }
        
        if(ys>=p->YN[jloc+marge] && ys<p->YN[jloc+1+marge] && stop==0)
        {
            jj = jloc;
            
            //cout<<"EXIT 2"<<endl;
   
         stop=1;
         break;   
        }
        
        
        // further divjsion
        if(ys<p->YN[jloc+marge] && ys<p->YN[jloc-1+marge])
        je=jloc;
        
        if(ys>p->YN[jloc+marge] && ys>p->YN[jloc+1+marge])
        js=jloc;
        
        
        ++count;
    }while(stop==0 && count<1000);
    
    jj=MAX(jj,-1);
    jj=MIN(jj,p->knoy+1);
    
    return jj;
}

int position::posc_k(double zs)
{
    stop=0;
    
    ks = 0;
    ke = p->knoz+1;

    count=0;
    do{
    kloc = ihalf(ks,ke);
    
    if(count%3==0)
    kloc+=1;
        
        // out of bounds
        if(zs<p->ZN[0+marge])
        {
            kk = -1;
            
            //cout<<"EXIT 0m"<<endl;
   
         stop=1;
         break;   
        }
        
        // out of bounds
        if(zs>p->ZN[p->knoz+marge])
        {
            kk = p->knoz+1;
            
            //cout<<"EXIT 0p"<<endl;
   
         stop=1;
         break;   
        }
        
        // matching criterion
        if(zs<p->ZN[kloc+marge] && zs>=p->ZN[kloc-1+marge] && stop==0)
        {
            kk = kloc-1;
            
         stop=1;
         break;   
        }
        
        if(zs>=p->ZN[kloc+marge] && zs<p->ZN[kloc+1+marge] && stop==0)
        {
            kk = kloc;
            
            //cout<<"EXIT 2"<<endl;
   
         stop=1;
         break;   
        }
        
        // further divksion
        if(zs<p->ZN[kloc+marge] && zs<p->ZN[kloc-1+marge])
        ke=kloc;
        
        if(zs>p->ZN[kloc+marge] && zs>p->ZN[kloc+1+marge])
        ks=kloc;
        
        
        ++count;
    }while(stop==0 && count<1000);
    
    kk=MAX(kk,-1);
    kk=MIN(kk,p->knoz+1);
    
    return kk;
}




