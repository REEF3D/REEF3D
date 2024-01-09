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
        
        //cout<<"xp: "<<xp<<" zp: "<<zp<<" i: "<<i<<" j: "<<j<<" Xstart: "<<Xstart<<" Xend: "<<Xend<<" Zstart: "<<Z[i][j][0]<<" Zend: "<<Z[i][j][Nz-1]<<endl;
        
        //cout<<i<<" "<<j<<" ZSE: "<<Z[i][j][Nz-3]<<" "<<Z[i][j][Nz-2]<<endl;

        if(file_version!=2)
        {
        k = pos_k(p,zp,i,j);

        val=ccpol3D(p,F,x,y,z);
        
        //cout<<"$#% SPACEPOLVAL "<<val<<endl;
        }
        
        if(file_version==2)
        {
        val=ccpol2DM(p,F,x,y);
        
        //cout<<"$#% SPACEPOLVAL "<<val<<" zp: "<<zp<<" E1: "<<E1[i][j]<<" E2: "<<E2[i][j]<<" E: "<<E[i][j]<<endl;
        }
    }
    
    i=ii;
    j=jj;
    k=kk;
    
    return val;
}

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

double wave_lib_hdc::ccpol3D(lexer *p, double ***F, double x, double y, double z)
{
    // wa
    if(xp>X[0] && xp<X[Nx-1])
    {
        wa = (X[i+1]-xp)/(X[i+1]-X[i]);
        
        if(i<Nx-1)
        if((X[i+1]-xp)/(X[i+1]-X[i])<0.0)
        {
        wa = (X[i+1]-xp)/(X[i+2]-X[i+1]);
        
        //if(p->mpirank==0 && x>9.99)
        //cout<<" wa1a: "<<wa<<" xp: "<<xp<<" X[i+1]: "<<X[i+1]<<" X[i+2]: "<<X[i+2]<<endl;
        ++i;
        }

        if(i>0)
        if((X[i+1]-xp)/(X[i+1]-X[i])>1.0)
        {
        wa = (X[i]-xp)/(X[i]-X[i-1]);
        
        //if(p->mpirank==0 && x>9.99)
        //cout<<" wa1b: "<<wa<<endl;
        --i;
        }
        
        //cout<<i<<" X[i-2]: "<<X[i-2]<<" X[i-1]: "<<X[i-1]<<" X[i]: "<<X[i]<<" X[i+1]: "<<X[i+1]<<" X[i+2]: "<<X[i+2]<<endl;
        
        
    }
    
    if(xp<=X[0] || i<0)
    {
    wa=0.0;
    
    //if(p->mpirank==0 && x>9.99)
    //cout<<" wa2: "<<wa<<endl;
    }
    
    if(xp>=X[Nx-1] || i>=Nx-1)
    {
    wa=1.0;
    
    //if(p->mpirank==0 && x>9.99)
    //cout<<" wa3: "<<wa<<endl;
    }
    
    
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
    
    /*
    cout<<p->mpirank<<" HDC  i: "<<i<<" j: "<<j<<" xp: "<<xp<<" yp: "<<yp<<" val: "<<val<<" Fi: "<<F[i][j][k]<<endl;
    
    cout<<" HDC 3D: "<<v1<<" "<<v2<<" "<<v3<<" "<<v4<<" "<<v5<<" "<<v6<<" "<<v7<<" "<<v7<<" "<<endl;
    cout<<" HDC i: "<<i<<" j: "<<j<<" k: "<<k<<" Z[ijk]: "<<Z[i][j][k]<<" z: "<<z<<endl;*/
    
    /*
    if(p->mpirank==0 && x>9.99)
    cout<<" Nx: "<<Nx<<" i: "<<i<<" j: "<<j<<" k: "<<k<<" X[i]: "<<X[i]<<" | U[i]: "<<U[i][0][2]
        <<" x: "<<x<<" z: "<<z<<" val: "<<val<<" | wa: "<<wa<<" wb: "<<wb<<" wc: "<<wc<<endl;*/
    

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


