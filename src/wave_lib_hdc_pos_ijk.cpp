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

#include"wave_lib_hdc.h"
#include"lexer.h"

int wave_lib_hdc::pos_i(lexer *p, double xs)
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
    
    Xstart =  X[0];
    Xend   =  X[Nx-1];
    
    
    XM1 = 0.5*(X[iloc] + X[iloc-1]);
    XP  = 0.5*(X[iloc] + X[iloc+1]);
    XP1 = 0.5*(X[iloc+1] + X[iloc+2]);
    
    
    //if(p->mpirank==0 && xs>9.99)
    //cout<<"POS_x: "<<xs<<" Nx: "<<Nx<<" Xmax: "<<X[Nx-1]<<" iloc: "<<iloc<<" XM1: "<<XM1<<" XP: "<<XP<<" XP1: "<<XP1<<endl;
    
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
        
        if(iloc<=0)
        {
            ii = 0;
            
         stop=1;
         break;   
        }
        
        if(iloc>=Nx-1)
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
  
    ii=MAX(ii,0);
    ii=MIN(ii,Nx-1);
    
    return ii;
}

int wave_lib_hdc::pos_j(lexer *p, double ys)
{
    int js,je,jloc;
    int stop=0;
    int count=0;
    int jj;
    double YM1,YP,YP1,Ystart,Yend;
    
    js = 0;
    je = Ny;
    
    
    count=0;
    if(jdir==1)
    do{
    jloc = ihalf(js,je);
    
    if(count%3==0)
    jloc+=1;
    

    Ystart = Y[0];
    Yend   = Y[Ny-1];
    
    YM1 = 0.5*(Y[jloc] + Y[jloc-1]);
    YP  = 0.5*(Y[jloc] + Y[jloc+1]);
    YP1 = 0.5*(Y[jloc+1] + Y[jloc+2]);
    
    //cout<<"ys: "<<ys<<" YM1: "<<YM1<<" YP: "<<YP<<" YP1: "<<YP1<<endl;
    
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
        
        if(jloc<=0)
        {
            jj = 0;
            
         stop=1;
         break;   
        }
        
        if(jloc>=Ny-1)
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
    
    jj=MAX(jj,0);
    jj=MIN(jj,Ny-1);
    
    
    if(jdir==0)
    jj=0;
    
    return jj;    
}

int wave_lib_hdc::pos_k(lexer *p, double zs, int i, int k)
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

         stop=1;
         break;   
        }
        
        // out of bounds
        if(zs>Z[i][k][Nz-1])
        {
            kk = Nz-1;
            
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

         stop=1;
         break;   
        }
        
        if(kloc<=0)
        {
            kk = 0;
            
         stop=1;
         break;   
        }
        
        if(kloc>=Nz-1)
        {
            kk = Nz-1;
            
         stop=1;
         break;   
        }
        
        // further divsion
        if(zs<Z[i][k][kloc] && zs<Z[i][k][kloc-1])
        ke=kloc;
        
        if(zs>Z[i][k][kloc] && zs>Z[i][k][kloc+1])
        ks=kloc;
        
        
        ++count;
    }while(stop==0 && count<1000);
    
    //cout<<" kk: "<<kk<<endl;
    
    kk=MAX(kk,0);
    kk=MIN(kk,Nz);
    
    
    return kk;
}


int wave_lib_hdc::ihalf(int a, int b)
{
    int c;
    double d,diff;

    c = b-a;
    
    d = double(c)*0.5;

    c = int(d) + a;
    

    return c;
}
