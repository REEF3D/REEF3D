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
    
    // wave theories (h_switch 1) are ramped as the relaxation zone (B 101); wavemakers carry their own signal
    const double rmp = h_switch==1 ? ramp(p) : 1.0;
    
    // B 98 4: active absorption at the wavemaker. The deviation of the computed surface from the target
    // (bound waves of the 2nd-order theories included) leaves the domain as a shallow-water wave,
    // u = -eta_R sqrt(g/h) on top of the target velocity, the exact condition for long waves
    // (free infragravity waves reflected in the domain). Before, B 98 4 ran as B 98 3 (no absorption).
    // B 98 6 (gen_lp) runs as B 98 4 with the edge targets of fnpf_precalc_edge (the relaxation zone
    // in front keeps its own)
    const bool absorb = p->B98==4 && h_switch==1;
    
    if((int)ua_edge.size()<p->gcslin_count)
    ua_edge.assign(p->gcslin_count,0.0);
    
    // B 98 7: short part of the deviation per column (two-band absorption)
    if(ab2_in && (int)ab2_sin.size()<p->gcslin_count)
    ab2_sin.assign(p->gcslin_count,0.0);
    
    count=0;
    for(n=0;n<p->gcslin_count;n++)
    {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        ua_edge[n] = 0.0;
        
        if(absorb && !ab2_in)
        ua_edge[n] = -(etaf(i,j) - rmp*(gen_lp ? edge_eta[n] : eta(i,j)))*sqrt(9.81/c->depth(i,j));
        
        // B 98 7: long part with the shallow-water velocity (here), short part with the peak-wave
        // profile (per level below)
        if(absorb && ab2_in)
        {
        const double eR = etaf(i,j) - rmp*(gen_lp ? edge_eta[n] : eta(i,j));
        const double eL = ab2_lowpass(p,ab2_lin,ab2_cin,n,eR,ab2_tau_in);
        ua_edge[n] = -eL*sqrt(9.81/c->depth(i,j));
        ab2_sin[n] = eR - eL;
        }

        if(h_switch==1)
        {
        // absorbing: zero gradient of the computed surface, else the (ramped) target
        const double eg = absorb ? etaf(i,j) : rmp*eta(i,j);
        etaf(i-1,j) = eg;
        etaf(i-2,j) = eg;
        etaf(i-3,j) = eg;
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
        
        // surface velocity at the edge
        double us = rmp*(gen_lp ? edge_us[count] : Fifsfval[count]) + ua_edge[n];
        
        if(absorb && ab2_in)
        us -= ab2_prof(c->depth(i,j),etaf(i,j))*ab2_sin[n];
        
        if(p->A329==1)
        {
        Fifsf(i-1,j) = Fifsf(i,j) - us*1.0*p->DXP[IM1];
        Fifsf(i-2,j) = Fifsf(i,j) - us*2.0*p->DXP[IM1];
        Fifsf(i-3,j) = Fifsf(i,j) - us*3.0*p->DXP[IM1];
        }
        
        if(p->A329>=2)
        {
        Fifsf(i-1,j) = (4.0/3.0)*Fifsf(i,j) - (1.0/3.0)*Fifsf(i+1,j) - (2.0/3.0)*us*(-1.5*p->XP[IM1] + 2.0*p->XP[IP] - 0.5*p->XP[IP1]);
        Fifsf(i-2,j) = (4.0/3.0)*Fifsf(i,j) - (1.0/3.0)*Fifsf(i+1,j) - (2.0/3.0)*us*(-1.5*p->XP[IM2] + 2.0*p->XP[IP] - 0.5*p->XP[IP1]);
        Fifsf(i-3,j) = (4.0/3.0)*Fifsf(i,j) - (1.0/3.0)*Fifsf(i+1,j) - (2.0/3.0)*us*(-1.5*p->XP[IM3] + 2.0*p->XP[IP] - 0.5*p->XP[IP1]);
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
        Uin[FIm1JK] = rmp*(gen_lp ? edge_uin[count] : Uinval[count]) + ua_edge[n]; 
        
        if(ab2_in && h_switch==1)
        Uin[FIm1JK] -= ab2_prof(c->depth(i,j),p->ZSN[FIJK]-p->phimean)*ab2_sin[n];
        
        // moving paddle: 2nd-order Taylor terms of the paddle BC about x=0,
        // phi_x = X_t + X_z*phi_z - X*phi_xx, with phi from the linear paddle
        // solution (prescribed, no feedback from the computed field)
        if(p->B119==1 && (p->B92==21 || p->B92==22 || p->B92==23))
        Uin[FIm1JK] += wave_paddle_Q(p,pgc,p->ZSN[FIJK]-p->phimean);
        
        ++count;
        }
    }
}

