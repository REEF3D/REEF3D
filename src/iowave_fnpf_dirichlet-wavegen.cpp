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

#include"iowave.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"

void iowave::dirichlet_wavegen_fnpf(lexer *p, fdm_fnpf *c, ghostcell* pgc, double *Fi, double *Uin, slice &Fifsf, slice &etaf)
{
    double etax;
    
    // 
    count=0;
    for(n=0;n<p->gcslin_count;n++)
    {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        if(h_switch==1)
        {
        //etaf(i,j)   = eta(i,j);
        etaf(i-1,j) = eta(i,j);
        etaf(i-2,j) = eta(i,j);
        etaf(i-3,j) = eta(i,j);
        }
        
        
        if(h_switch==0)
        {
        if(p->A329==1 || p->count<=2)
        etax = -(1.0/9.81) * (Fifsfval[count]-Fifsfval0[count])/p->dt;
        
        if(p->A329==2 && p->count>2)
        etax = -(1.0/9.81) * (-1.5*Fifsfval[count] + 2.0*Fifsfval0[count] - 0.5*Fifsfval1[count])/(-1.5*time_n + 2.0*time_0 - 0.5*time_1);

        etaf(i-1,j) = etaf(i,j) + etax*1.0*p->DXP[IM1];
        etaf(i-2,j) = etaf(i,j) + etax*2.0*p->DXP[IM1];
        etaf(i-3,j) = etaf(i,j) + etax*3.0*p->DXP[IM1];
        }
        
        if(p->A329==1)
        {
        Fifsf(i-1,j) = Fifsf(i,j) - Fifsfval[count]*1.0*p->DXP[IM1];
        Fifsf(i-2,j) = Fifsf(i,j) - Fifsfval[count]*2.0*p->DXP[IM1];
        Fifsf(i-3,j) = Fifsf(i,j) - Fifsfval[count]*3.0*p->DXP[IM1];
        }
        
        if(p->A329>=2)
        {
        Fifsf(i-1,j) = (4.0/3.0)*Fifsf(i,j) - (1.0/3.0)*Fifsf(i+1,j) - (2.0/3.0)*Fifsfval[count]*(-1.5*p->XP[IM1] + 2.0*p->XP[IP] - 0.5*p->XP[IP1]);
        Fifsf(i-2,j) = (4.0/3.0)*Fifsf(i,j) - (1.0/3.0)*Fifsf(i+1,j) - (2.0/3.0)*Fifsfval[count]*(-1.5*p->XP[IM2] + 2.0*p->XP[IP] - 0.5*p->XP[IP1]);
        Fifsf(i-3,j) = (4.0/3.0)*Fifsf(i,j) - (1.0/3.0)*Fifsf(i+1,j) - (2.0/3.0)*Fifsfval[count]*(-1.5*p->XP[IM3] + 2.0*p->XP[IP] - 0.5*p->XP[IP1]);
        }
    
        ++count;
    }
    
    /*
    count=0;
    for(n=0;n<p->gcslin_count;n++)
    {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        FKLOOP
        FPCHECK
        {
        Fi[FIm1JK] = Fi[FIJK] - Uinval[count]*1.0*p->DXP[IM1];
        Fi[FIm2JK] = Fi[FIJK] - Uinval[count]*2.0*p->DXP[IM1];
        Fi[FIm3JK] = Fi[FIJK] - Uinval[count]*3.0*p->DXP[IM1];
        
        ++count;
        }
    }*/
    
    // Uin
    count=0;
    for(n=0;n<p->gcslin_count;n++)
    {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        FKLOOP
        FPCHECK
        {
        Uin[FIm1JK] = Uinval[count]; 
        
        // moving paddle: 2nd-order Taylor terms of the paddle BC about x=0
        // phi_x = X_t + X_z*phi_z - X*phi_xx,  phi_xx = -(phi_yy + phi_zz)
        if(p->B119==1 && (p->B92==21 || p->B92==22) && p->knoz>=2)
        Uin[FIm1JK] += paddle_taylor_fnpf(p,pgc,Fi);
        
        ++count;
        }
    }
}

double iowave::paddle_taylor_fnpf(lexer *p, ghostcell *pgc, double *Fi)
{
    double zz = p->ZSN[FIJK]-p->phimean;
    
    double X  = wave_paddle_X(p,pgc,zz);
    double Xz = wave_paddle_Xz(p,pgc,zz);
    
    if(fabs(X)<1.0e-20 && fabs(Xz)<1.0e-20)
    return 0.0;
    
    // three-point stencil in the vertical column at the first interior node
    int k0 = k-1;
    
    if(k==0)
    k0 = 0;
    
    if(k==p->knoz)
    k0 = k-2;
    
    const int id = FIJK;
    const int d0 = k0-k;
    
    double z0 = p->ZSN[id+d0], z1 = p->ZSN[id+d0+1], z2 = p->ZSN[id+d0+2];
    double f0 = Fi[id+d0],     f1 = Fi[id+d0+1],     f2 = Fi[id+d0+2];
    double zc = p->ZSN[id];
    
    double c0 = 1.0/((z0-z1)*(z0-z2));
    double c1 = 1.0/((z1-z0)*(z1-z2));
    double c2 = 1.0/((z2-z0)*(z2-z1));
    
    double fzz = 2.0*(f0*c0 + f1*c1 + f2*c2);
    double fz  = f0*c0*(2.0*zc-z1-z2) + f1*c1*(2.0*zc-z0-z2) + f2*c2*(2.0*zc-z0-z1);
    
    if(k==0)
    fz = 0.0;   // bed
    
    double fyy = 0.0;
    
    if(p->y_dir>0.5)
    fyy = 2.0*((Fi[FIJp1K]-Fi[FIJK])/p->DYP[JP] - (Fi[FIJK]-Fi[FIJm1K])/p->DYP[JM1])/(p->DYP[JP]+p->DYP[JM1]);
    
    return Xz*fz + X*(fzz + fyy);
}
