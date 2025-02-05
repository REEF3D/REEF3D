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

double wave_lib_hdc::space_interpol(lexer *p, double ***F, double x, double y, double z)
{
    val=0.0;
    
    ii=i;
    jj=j;
    kk=k;
    
    xp = x + p->I231;
    yp = y + p->I232;
    zp = z + p->I233 + p->wd;
    
    
    if(xp>=Xstart  && xp<Xend && ((yp>=Ystart && yp<Yend)|| jdir==0))
    {
        i = pos_i(p,xp);
        j = pos_j(p,yp);
        
        
        // sflow
        if(file_version==2)
        {
        val=ccpol2DM(p,F,xp,yp);
        }
        
        // fnpf
        if(file_version==4)
        {
        k = pos_k(p,zp,i,j);

        val=ccpol3D(p,F,xp,yp,zp);
        }
        
        // nhflow
        if(file_version==5)
        {
        k = pos_k(p,zp,i,j);

        val=ccpol3D(p,F,xp,yp,zp);
        }
    }
    
    i=ii;
    j=jj;
    k=kk;
    
    return val;
}

double wave_lib_hdc::ccpol3D(lexer *p, double ***F, double xp, double yp, double zp)
{
    // wa
    if(xp>X[0] && xp<X[Nx-1])
    {
        wa = (X[i+1]-xp)/(X[i+1]-X[i]);
        
        if(i<Nx-1)
        if((X[i+1]-xp)/(X[i+1]-X[i])<0.0)
        {
        wa = (X[i+1]-xp)/(X[i+2]-X[i+1]);
        ++i;
        }

        if(i>0)
        if((X[i+1]-xp)/(X[i+1]-X[i])>1.0)
        {
        wa = (X[i]-xp)/(X[i]-X[i-1]);
        --i;
        }
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
    
    //wc
    if(zp>Z[i][j][0] && zp<Z[i][j][Nz-1])
    {
        wc = (Z[i][j][k+1]-zp)/(Z[i][j][k+1]-Z[i][j][k]);
        
        if(k<Nz-1)
        if((Z[i][j][k+1]-zp)/p->DZP[KP]<0.0)
        {
        wc = (Z[i][j][k+2]-zp)/(Z[i][j][k+2]-Z[i][j][k+1]);
        ++k;
        }
        
        if(k>0)
        if((Z[i][j][k+1]-zp)/p->DZP[KP]>1.0)
        {
        wc = (Z[i][j][k]-zp)/(Z[i][j][k]-Z[i][j][k-1]);
        --k;
        }
    }
    
    if(zp<=Z[i][j][0])
    wc=0.0;
    
    if(zp>=Z[i][j][Nz-1])
    wc=1.0;


    v1=v2=v3=v4=v5=v6=v7=v8=0.0;

    ip1 = (i+1)<(Nx-1)?(i+1):i;
    jp1 = (j+1)<(Ny-1)?(j+1):j;
    kp1 = (k+1)<(Nz-1)?(k+1):k;
    
    iii=i;
    jjj=j;
    kkk=k;
    
    i = i<0?0:i;
    j = j<0?0:j;
    k = k<0?0:k;
    
    i = i>Nx-1?Nx-1:i;
    j = j>Ny-1?Ny-1:j;
    k = k>Nz-1?Nz-1:k;
    
    
    v1 = F[i][j][k];
    v2 = F[i][jp1][k];
    v3 = F[ip1][j][k];
    v4 = F[ip1][jp1][k];
    
    v5 = F[i][j][kp1];
    v6 = F[i][jp1][kp1];
    v7 = F[ip1][j][kp1];
    v8 = F[ip1][jp1][kp1];
    

    x1 = wa*v1 + (1.0-wa)*v3;
    x2 = wa*v2 + (1.0-wa)*v4;

    x3 = wa*v5 + (1.0-wa)*v7;
    x4 = wa*v6 + (1.0-wa)*v8;

    y1 = wb*x1 +(1.0-wb)*x2;
    y2 = wb*x3 +(1.0-wb)*x4;

    val = wc*y1 +(1.0-wc)*y2;
    
    i=iii;
    j=jjj;
    k=kkk;

    return val;
}


