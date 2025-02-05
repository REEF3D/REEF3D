/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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


double wave_lib_hdc::plane_interpol(lexer *p, double **F, double x, double y)
{
    val=0.0;
    
    ii=i;
    jj=j;
    
    xp = x + p->I231;
    yp = y + p->I232;

    if(xp>=Xstart  && xp<Xend && ((yp>=Ystart && yp<Yend)|| jdir==0))
    {
        i = pos_i(p,xp);
        j = pos_j(p,yp);

        val=ccpol2D(p,F,x,y);
    }
    

    i=ii;
    j=jj;
    
    return val;
}

double wave_lib_hdc::ccpol2D(lexer *p, double **F, double x, double y)
{
    
    // wa
    if(xp>X[0] && xp<X[Nx-1])
    {
        wa = (X[i+1]-xp)/(X[i+1]-X[i]);
        
        if(i<Nx-1)
        if((X[i+1]-xp)/(X[i+1]-X[i])<0.0)
        {
        wa = (X[i+2]-xp)/(X[i+2]-X[i+1]);
        ++i;
        }
        
        if(i>0)
        if((X[i+1]-xp)/(X[i+1]-X[i])>1.0)
        {
        wa = (X[i]-xp)/(X[i]-X[i-1]);
        --i;
        }
        
        //cout<<i<<" X[i-2]: "<<X[i-2]<<" X[i-1]: "<<X[i-1]<<" X[i]: "<<X[i]<<" X[i+1]: "<<X[i+1]<<" X[i+2]: "<<X[i+2]<<endl;
    }
    
    if(xp<=X[0] || i<0)
    wa=0.0;
    
    if(xp>=X[Nx-1] || i>=Nx-1)
    wa=1.0;
    

    // wb
    if(Ny==1 || jdir==0)
    wb=1.0;    
        
    if(Ny>1 && jdir==1)
    {
    if(yp>Y[0] && yp<Y[Ny-1])
    {
        wb = (Y[j+1]-yp)/(Y[j+1]-Y[j]);
        
        if((Y[j+1]-yp)/(Y[j+1]-Y[j])<0.0)
        {
        wb = (Y[j+1]-yp)/(Y[j+2]-Y[j+1]);
        ++j;
        }
        
        if((Y[j+1]-yp)/(Y[j+1]-Y[j])>1.0)
        {
        wb = (Y[j]-yp)/(Y[j]-Y[j-1]);
        --j;
        }
    }
    
    if(yp<=Y[0])
    wb=0.0;
    
    if(yp>=Y[Ny-1])
    wb=1.0;
    }
    
    //cout<<"wb2D: "<<wb<<endl;


    v1=v2=v3=v4=0.0;

    
    ip1 = (i+1)<(Nx-1)?(i+1):i;
    jp1 = (j+1)<(Ny-1)?(j+1):j;
    
    iii=i;
    jjj=j;

    i = i<0?0:i;
    j = j<0?0:j;
    
    i = i>(Nx-1)?(Nx-1):i;
    j = j>(Ny-1)?(Ny-1):j;

    
    v1 = F[i][j];
    v2 = F[i][jp1];
    v3 = F[ip1][j];
    v4 = F[ip1][jp1];
    
    //cout<<"i: "<<i<<" j: "<<j<<" ip1: "<<ip1<<" jp1: "<<jp1<<endl;


    x1 = wa*v1 + (1.0-wa)*v3;
    x2 = wa*v2 + (1.0-wa)*v4;

    val = wb*x1 +(1.0-wb)*x2;
    
    //cout<<p->mpirank<<" HDC  i: "<<i<<" j: "<<j<<" xp: "<<xp<<" yp: "<<yp<<" val: "<<val<<" Fi: "<<F[i][j]<<endl;
    
    i=iii;
    j=jjj;
    
    
    return val;
}


double wave_lib_hdc::ccpol2DM(lexer *p, double ***F, double x, double y)
{
    
    // wa
    if(xp>X[0] && xp<X[Nx-1])
    {
        wa = (X[i+1]-xp)/(X[i+1]-X[i]);
        
        if(i<Nx-1)
        if((X[i+1]-xp)/(X[i+1]-X[i])<0.0)
        {
        wa = (X[i+2]-xp)/(X[i+2]-X[i+1]);
        ++i;
        }
        
        if(i>0)
        if((X[i+1]-xp)/(X[i+1]-X[i])>1.0)
        {
        wa = (X[i]-xp)/(X[i]-X[i-1]);
        --i;
        }
        
        //cout<<i<<" X[i-2]: "<<X[i-2]<<" X[i-1]: "<<X[i-1]<<" X[i]: "<<X[i]<<" X[i+1]: "<<X[i+1]<<" X[i+2]: "<<X[i+2]<<endl;
    }
    
    if(xp<=X[0] || i<0)
    wa=0.0;
    
    if(xp>=X[Nx-1] || i>=Nx-1)
    wa=1.0;
    

    // wb
    if(Ny==1 || jdir==0)
    wb=1.0;    
        
    if(Ny>1 && jdir==1)
    {
    if(yp>Y[0] && yp<Y[Ny-1])
    {
        wb = (Y[j+1]-yp)/(Y[j+1]-Y[j]);
        
        if((Y[j+1]-yp)/(Y[j+1]-Y[j])<0.0)
        {
        wb = (Y[j+1]-yp)/(Y[j+2]-Y[j+1]);
        ++j;
        }
        
        if((Y[j+1]-yp)/(Y[j+1]-Y[j])>1.0)
        {
        wb = (Y[j]-yp)/(Y[j]-Y[j-1]);
        --j;
        }
    }
    
    if(yp<=Y[0])
    wb=0.0;
    
    if(yp>=Y[Ny-1])
    wb=1.0;
    }

    v1=v2=v3=v4=0.0;

    
    ip1 = (i+1)<(Nx-1)?(i+1):i;
    jp1 = (j+1)<(Ny-1)?(j+1):j;
    
    iii=i;
    jjj=j;

    i = i<0?0:i;
    j = j<0?0:j;
    
    i = i>(Nx-1)?(Nx-1):i;
    j = j>(Ny-1)?(Ny-1):j;

    
    v1 = F[i][j][0];
    v2 = F[i][jp1][0];
    v3 = F[ip1][j][0];
    v4 = F[ip1][jp1][0];
    
    //cout<<"i: "<<i<<" j: "<<j<<" ip1: "<<ip1<<" jp1: "<<jp1<<endl;


    x1 = wa*v1 + (1.0-wa)*v3;
    x2 = wa*v2 + (1.0-wa)*v4;

    val = wb*x1 +(1.0-wb)*x2;
    
    //cout<<p->mpirank<<" HDC  i: "<<i<<" j: "<<j<<" xp: "<<xp<<" yp: "<<yp<<" val: "<<val<<" Fi: "<<F[i][j]<<endl;
    
    i=iii;
    j=jjj;
    
    
    return val;
}


