/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"wave_lib_wcp.h"
#include"lexer.h"

int wave_lib_wcp::pos_i(double xs)
{
    int is,ie,iloc;
    int stop=0;
    int count=0;
    int ii;
    double XM1,XP,XP1,Xstart,Xend;
    
    is = 0;
    ie = Nx;
    
    
    count=0;
    do{
    iloc = ihalf(is,ie);
    
    if(count%3==0)
    iloc+=1;
    
    Xstart =  X[0];// - 0.5*(X[0] + X[1]);
    Xend   =  X[Nx-1];// + 0.5*(X[Nx-1] + X[Nx]);
    
    XM1 = 0.5*(X[iloc] + X[iloc-1]);
    XP  = 0.5*(X[iloc] + X[iloc+1]);
    XP1 = 0.5*(X[iloc+1] + X[iloc+2]);
    
        // matching criterion
        if(xs<XP && xs>=XM1)
        {
            ii = iloc;
            
         stop=1;
         break;   
        }
        
        if(xs>=XP && xs<XP1)
        {
            ii = iloc+1;
   
         stop=1;
         break;   
        }
        
        // out of bounds
        if(xs<Xstart)
        {
            ii = 0;

         stop=1;
         break;   
        }
        
        // out of bounds
        if(xs>Xend)
        {
            ii = Nx-1;
            
         stop=1;
         break;   
        }
        
        // further division
        if(xs<XP && xs<XM1)
        ie=iloc;
        
        if(xs>XP && xs>XP1)
        is=iloc;
        
        ++count;
    }while(stop==0 && count<1000);
  
  
    cout<<" ii: "<<ii<<endl;
    
    ii=MAX(ii,0);
    ii=MIN(ii,Nx-1);
    
    
    
    return ii;
    
}

int wave_lib_wcp::pos_j(double ys)
{
    int js,je,jloc;
    int stop=0;
    int count=0;
    int jj;
    double YM1,YP,YP1,Ystart,Yend;
    
    js = 0;
    je = Ny;
    
    
    count=0;
    do{
    jloc = ihalf(js,je);
    
    if(count%3==0)
    jloc+=1;
    
    Ystart = 0.5*(Y[0] + Y[1]);
    Yend   = 0.5*(Y[Ny-1] + Y[Ny]);
    
    YM1 = 0.5*(Y[jloc] + Y[jloc-1]);
    YP  = 0.5*(Y[jloc] + Y[jloc+1]);
    YP1 = 0.5*(Y[jloc+1] + Y[jloc+2]);
    
        // matching criterion
        if(ys<YP && ys>=YM1)
        {
            jj = jloc;
            
         stop=1;
         break;   
        }
        
        if(ys>=YP && ys<YP1)
        {
            jj = jloc+1;
   
         stop=1;
         break;   
        }
        
        // out of bounds
        if(ys<Ystart)
        {
            jj = 0;

         stop=1;
         break;   
        }
        
        // out of bounds
        if(ys>Yend)
        {
            jj = Ny-1;
            
         stop=1;
         break;   
        }
        
        // further divjsion
        if(ys<YP && ys<YM1)
        je=jloc;
        
        if(ys>YP && ys>YP1)
        js=jloc;
        
        ++count;
    }while(stop==0 && count<1000);
    
    cout<<" jj: "<<jj<<endl;
    
    jj=MAX(jj,0);
    jj=MIN(jj,Ny);
    
    return jj;
    
}

int wave_lib_wcp::pos_k(double zs, int i, int k)
{
    int ks,ke,kloc;
    int stop=0;
    int count=0;
    int kk;
    
    ks = 0;
    ke = Nz;
    
    
    count=0;
    do{
    kloc = ihalf(ks,ke);
    
    if(count%3==0)
    kloc+=1;
        
        // out of bounds
        if(zs<Z[i][k][0])
        {
            kk = 0;
            
            //cout<<"EXIT 0m"<<endl;
   
         stop=1;
         break;   
        }
        
        // out of bounds
        if(zs>Z[i][k][Nz-1])
        {
            kk = Nz-1;
            
            //cout<<"EXIT 0p"<<endl;
   
         stop=1;
         break;   
        }
        
        // matching criterion
        if(zs<Z[i][k][kloc] && zs>=Z[i][k][kloc-1] && stop==0)
        {
            kk = kloc-1;
            
         stop=1;
         break;   
        }
        
        if(zs>=Z[i][k][kloc] && zs<Z[i][k][kloc+1] && stop==0)
        {
            kk = kloc;
            
            //cout<<"EXIT 2"<<endl;
   
         stop=1;
         break;   
        }
        
        // further divksion
        if(zs<Z[i][k][kloc] && zs<Z[i][k][kloc-1])
        ke=kloc;
        
        if(zs>Z[i][k][kloc] && zs>Z[i][k][kloc+1])
        ks=kloc;
        
        
        ++count;
    }while(stop==0 && count<1000);
    
    cout<<" kk: "<<kk<<endl;
    
    kk=MAX(kk,0);
    kk=MIN(kk,Nz);
    
    return kk;
}


int wave_lib_wcp::ihalf(int a, int b)
{
    int c;
    double d,diff;

    c = b-a;
    
    d = double(c)*0.5;

    c = int(d) + a;
    

    return c;
}
