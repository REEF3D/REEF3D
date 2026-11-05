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
    const bool absorb_th = p->B98==4 && h_switch==1;
    
    // B 98 4/7 with the piston wavemaker (B 92 21): active wave absorption at the paddle (ARC). The deviation
    // of the computed surface at the first cell from the linear surface elevation of the paddle signal at
    // that distance (progressive and evanescent modes, wave_lib_wavemaker2nd::make_etatable) is the reflected
    // wave; it is absorbed with the same paddle velocity correction as for the wave theories (shallow-water
    // u = -eta_R sqrt(g/h), Schaffer & Klopman 2000, or the two-band version with B 98 7)
    const bool absorb_wm = p->B98==4 && h_switch==0 && wave_arc();
    
    // B 98 8 (SIRW, arc_sirw): the reflected wave of the band from two gauges in front of the paddle, plus the
    // shallow-water absorption of the long band of the deviation from the target
    const bool absorb_sw = sirw_in.on;
    const bool absorb = absorb_th || absorb_wm || absorb_sw;
    
    if(absorb_sw)
    sirw_in.update(p,pgc,etaf,p->simtime);
    
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
        
        // target: wave theory (ramped) or linear elevation of the paddle signal at the first cell
        const double etaT = absorb_wm ? wave_eta_arc(p,pgc,p->XP[IP]-p->global_xmin)
                                      : rmp*(gen_lp ? edge_eta[n] : eta(i,j));
        
        if(absorb && !ab2_in && !absorb_sw)
        ua_edge[n] = -(etaf(i,j) - etaT)*sqrt(9.81/c->depth(i,j));
        
        // B 98 8: the band of the gauges per level below; what it does not explain (the deviation from the target
        // minus the predicted reflected wave: long waves below the band, errors of the prediction) is absorbed
        // with the shallow-water velocity after a first-order low pass (removes the near field and bound waves
        // of the generated waves; its phase lag stays below 90 degrees, so the absorption damps at all
        // frequencies; with a third-order low pass the 2nd and 4th seiche modes of the tank grew)
        if(absorb_sw && (absorb_th || absorb_wm))
        ua_edge[n] = -sw_lowpass(p,n,etaf(i,j)-etaT-sirw_in.eta(p->simtime))*sqrt(9.81/c->depth(i,j));
        
        // B 98 7: long part with the shallow-water velocity (here), short part with the peak-wave
        // profile (per level below)
        if(absorb && ab2_in)
        {
        const double eR = etaf(i,j) - etaT;
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
        
        // B 98 8: band of the gauges, velocity at the surface
        if(absorb_sw)
        us += sirw_in.velocity(p->simtime,etaf(i,j));
        
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
        
        // B 98 8: band of the gauges, velocity of this level
        if(absorb_sw)
        Uin[FIm1JK] += sirw_in.velocity(p->simtime,p->ZSN[FIJK]-p->phimean);
        
        if(ab2_in && absorb)
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

