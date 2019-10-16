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
    
    xp += p->I231;
    yp += p->I232;
    zp += p->I233;
    
    i = p->posc_i(xp);
    j = p->posc_j(yp);
    k = p->posf_k(zp);
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXN[IP];
    
    if((p->XP[IP1]-xp)/p->DXN[IP]<0.0)
    {
    wa = (p->XP[IP2]-xp)/p->DXN[IP1];
    ++i;
    }
    
    if((p->XP[IP1]-xp)/p->DXN[IP]>1.0)
    {
    wa = (p->XP[IP]-xp)/p->DXN[IM1];
    --i;
    }
    
    
    // wb
    wb = (p->YP[JP1]-yp)/p->DYN[JP];
    
    if((p->YP[JP1]-yp)/p->DYN[JP]<0.0)
    {
    wb = (p->YP[JP2]-yp)/p->DYN[JP1];
    ++j;
    }
    
    if((p->YP[JP1]-yp)/p->DYN[JP]>1.0)
    {
    wb = (p->YP[JP]-yp)/p->DYN[JM1];
    --j;
    }
    
    
    //wc
    wc = (p->ZN[KP1]-zp)/p->DZP[KP];
    
    if((p->ZN[KP1]-zp)/p->DZP[KP]<0.0)
    {
    wc = (p->ZN[KP2]-zp)/p->DZP[KP1];
    ++k;
    }
    
    if((p->ZN[KP1]-zp)/p->DZP[KP]>1.0)
    {
    wc = (p->ZN[KP]-zp)/p->DZP[KM1];
    --k;
    }


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
    
    xp += p->I231;
    yp += p->I232;
    
    i = p->posc_i(xp);
    j = p->posc_j(yp);
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXN[IP];
    
    if((p->XP[IP1]-xp)/p->DXN[IP]<0.0)
    {
    wa = (p->XP[IP2]-xp)/p->DXN[IP1];
    ++i;
    }
    
    if((p->XP[IP1]-xp)/p->DXN[IP]>1.0)
    {
    wa = (p->XP[IP]-xp)/p->DXN[IM1];
    --i;
    }
    
    
    // wb
    wb = (p->YP[JP1]-yp)/p->DYN[JP];
    
    if((p->YP[JP1]-yp)/p->DYN[JP]<0.0)
    {
    wb = (p->YP[JP2]-yp)/p->DYN[JP1];
    ++j;
    }
    
    if((p->YP[JP1]-yp)/p->DYN[JP]>1.0)
    {
    wb = (p->YP[JP]-yp)/p->DYN[JM1];
    --j;
    }


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
    
    return val;
}


