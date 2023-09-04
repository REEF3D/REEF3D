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

#include"iowave.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"


void iowave::dirichlet_wavegen_fnpf(lexer *p, fdm_fnpf *c, ghostcell* pgc, double *Fi, double *Uin, slice &Fifsf, slice &eta)
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
        eta(i,j)   = etaval[count];
        eta(i-1,j) = etaval[count];
        eta(i-2,j) = etaval[count];
        eta(i-3,j) = etaval[count];
        }
        
        
        if(h_switch==0)
        {
        if(p->A329==1 || p->count<=2)
        etax = -(1.0/9.81) * (Fifsfval[count]-Fifsfval0[count])/p->dt;
        
        if(p->A329==2 && p->count>2)
        etax = -(1.0/9.81) * (-1.5*Fifsfval[count] + 2.0*Fifsfval0[count] - 0.5*Fifsfval1[count])/(-1.5*time_n + 2.0*time_0 - 0.5*time_1);

        eta(i-1,j) = eta(i,j) + etax*1.0*p->DXP[IM1];
        eta(i-2,j) = eta(i,j) + etax*2.0*p->DXP[IM1];
        eta(i-3,j) = eta(i,j) + etax*3.0*p->DXP[IM1];
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
        
        
        
    }
    
    // Uin
    count=0;
    for(n=0;n<p->gcslin_count;n++)
    {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        FKLOOP
        FPCHECK
        {// add eta guard
        Uin[FIm1JK] = Uinval[count]; 
        
        ++count;
        }
    }
}
