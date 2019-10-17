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

double wave_lib_wcp::space_interpol(lexer *p, double ***F, double x, double y, double z)
{
    
    ii=i;
    jj=j;
    kk=k;
    
    xp = x + p->I231;
    yp = y + p->I232;
    zp = z + p->I233;
    
    i = pos_i(xp);
    j = pos_j(yp);
    k = pos_k(zp,i,j);
    
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
    
    if(xp<=X[0])
    wa=0.0;
    
    if(xp>=X[Nx-1])
    wa=1.0;
    
    
    // wb
    if(Ny==1)
    wb=1.0;    
        
    if(Ny>1)
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
        wc = (p->ZN[KP1]-zp)/p->DZP[KP];
        
        if(k<Nz-1)
        if((p->ZN[KP1]-zp)/p->DZP[KP]<0.0)
        {
        wc = (p->ZN[KP2]-zp)/p->DZP[KP1];
        ++k;
        }
        
        if(k>0)
        if((p->ZN[KP1]-zp)/p->DZP[KP]>1.0)
        {
        wc = (p->ZN[KP]-zp)/p->DZP[KM1];
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

    i=ii;
    j=jj;
    k=kk;

    return val;
}

double wave_lib_wcp::plane_interpol(lexer *p, double **F, double x, double y)
{
    ii=i;
    jj=j;
    
    xp = x + p->I231;
    yp = y + p->I232;
    
    i = pos_i(xp);
    j = pos_j(yp);
    
    /*cout<<"xs_xe: "<<X[0]<<" "<<X[Nx-1]<<endl;
    cout<<"pos_xy: "<<xp<<" "<<yp<<endl;
    cout<<"pos_ij: "<<i<<" "<<j<<endl;
    cout<<"Nxy: "<<Nx<<" "<<Ny<<endl;*/
		
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
    
    if(xp<=X[0])
    wa=0.0;
    
    if(xp>=X[Nx-1])
    wa=1.0;
    
    
    
    // wb
    if(Ny==1)
    wb=1.0;    
        
    if(Ny>1)
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
    
    cout<<"wa: "<<wa<<" wb: "<<wb<<endl;


    v1=v2=v3=v4=0.0;

    
    ip1 = (i+1)<(Nx-1)?(i+1):i;
    jp1 = (j+1)<(Ny-1)?(j+1):j;
    
    v1 = F[i][j];
    v2 = F[i][jp1];
    v3 = F[ip1][j];
    v4 = F[ip1][jp1];



    x1 = wa*v1 + (1.0-wa)*v3;
    x2 = wa*v2 + (1.0-wa)*v4;

    val = wb*x1 +(1.0-wb)*x2;

    i=ii;
    j=jj;
    
    cout<<p->mpirank<<" WCP  i: "<<i<<" j: "<<j<<" xp: "<<xp<<" yp: "<<yp<<" val: "<<val<<" Fi: "<<F[i][j]<<endl;
    
    return val;
}


