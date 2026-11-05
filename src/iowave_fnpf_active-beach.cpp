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

void iowave::active_beach_fnpf(lexer *p, fdm_fnpf *c, ghostcell* pgc, double *Fi, double *Uin, slice &Fifsf, slice &eta)
{
    // active beach at the x+ outflow: the outgoing wave leaves with u = eta sqrt(g/h)
    // B 99 3: shallow-water velocity, uniform over the depth (long waves)
    // B 99 4: linear-theory profile of the B 92 wave, eta omega cosh(k(h+z))/sinh(kh)
    // B 99 5: shallow-water flux, carried linearly between z = B 123 and the surface
    // B 99 7: two-band, shallow-water velocity for the low-passed deviation, B 99 4 profile for the rest
    double eta_R,Uc,eta_T,eta_M,hd,z,zabs,wsf,fac,multiplier;
    
    // B 99 8: gauges and filter once per step on all ranks (global sums)
    if(sirw_out.on)
    sirw_out.update(p,pgc,eta,p->simtime);
    
		for(n=0;n<p->gcslout_count;++n)
		{
		i=p->gcslout[n][0];
		j=p->gcslout[n][1];
		
		eta_T = 0.0;
		eta_M = eta(i,j); 
		eta_R = eta_M-eta_T;
        
        hd = c->depth(i,j);
        
        Uc=eta_R*sqrt(9.81/p->wd);
        
        // B 99 8 (SIRW): shallow-water velocity of what the predicted outgoing wave does not explain (long waves below
        // the band, errors of the prediction), the band of the gauges per level below
        const double eLs = sirw_out.on ? eta_R - sirw_out.eta(p->simtime) : 0.0;
        
        // B 99 7: long and short part of the deviation
        double eL=0.0, eS=0.0;
        if(ab2_out)
        {
        eL = ab2_lowpass(p,ab2_lout,ab2_cout,n,eta_R,ab2_tau_out);
        eS = eta_R - eL;
        }
        
			FKLOOP 
			{
				z=p->ZSN[FIJK]-p->phimean;
				
				if(p->B99==3 && !sirw_out.on)
				Uc=eta_R*sqrt(9.81/hd);
				
				// was sinh(depth*depth) and the still-water depth in cosh
				if(p->B99==4)
				Uc=eta_R*p->ww*(cosh(p->wk*(hd+z))/sinh(p->wk*hd));
                
                // was the CFD position pos_z() and an uninitialised surface level
                if(p->B99==5)
                {
                zabs = p->ZSN[FIJK];
                wsf  = p->phimean + eta_R;
                
                    if(zabs>p->B123 && wsf>p->B123)
                    {
                    fac = (zabs-p->B123)/(wsf-p->B123);
                    multiplier = 2.0*((wsf)/(wsf-p->B123));
                    Uc =   multiplier*fac*eta_R*sqrt(9.81/hd);
                    }
                    
                    else
                    Uc=0.0;
                }
				
                // B 99 8: band of the gauges with its profile, shallow water for the rest
                if(sirw_out.on)
                Uc = eLs*sqrt(9.81/hd) + sirw_out.velocity(p->simtime,z);
                
                // B 99 7: shallow-water velocity of the long part, peak-wave profile for the short part
                if(ab2_out)
                Uc = eL*sqrt(9.81/hd) + ab2_prof(hd,z)*eS;
				
				Uin[FIp1JK]=Uc;
			}
        
        // free surface: the surface velocity (top level)
        if(p->A329==1)
        { 
        Fifsf(i+1,j) = Fifsf(i,j) + Uc*1.0*p->DXP[IP1];
        Fifsf(i+2,j) = Fifsf(i,j) + Uc*2.0*p->DXP[IP1];
        Fifsf(i+3,j) = Fifsf(i,j) + Uc*3.0*p->DXP[IP1];
        }
        
        if(p->A329>=2)
        {
        Fifsf(i+1,j) = (4.0/3.0)*Fifsf(i,j) - (1.0/3.0)*Fifsf(i-1,j) - (2.0/3.0)*Uc*(-0.5*p->XP[IM1] + 2.0*p->XP[IP] - 1.5*p->XP[IP1]);
        Fifsf(i+2,j) = (4.0/3.0)*Fifsf(i,j) - (1.0/3.0)*Fifsf(i-1,j) - (2.0/3.0)*Uc*(-0.5*p->XP[IM1] + 2.0*p->XP[IP] - 1.5*p->XP[IP2]);
        Fifsf(i+3,j) = (4.0/3.0)*Fifsf(i,j) - (1.0/3.0)*Fifsf(i-1,j) - (2.0/3.0)*Uc*(-0.5*p->XP[IM1] + 2.0*p->XP[IP] - 1.5*p->XP[IP3]);
        }
        
        // ghost potential with the velocity of its level (was the surface velocity at all levels)
        FKLOOP
        FPCHECK
        {
        Fi[FIp1JK] = Fi[FIJK] + Uin[FIp1JK]*1.0*p->DXP[IP1];
        Fi[FIp2JK] = Fi[FIJK] + Uin[FIp1JK]*2.0*p->DXP[IP1];
        Fi[FIp3JK] = Fi[FIJK] + Uin[FIp1JK]*3.0*p->DXP[IP1];
        }
		}
}

// two-band active absorption (B 98 7, B 99 7): frequency and wave number of the reference period
// (B 528, default Tp of the B 92 wave), low-pass time scales (B 527 at the wavemaker, B 526 at the
// beach; default a quarter of the reference period)
void iowave::ab2_ini(lexer *p)
{
    const double g = 9.81;
    const double Tr = p->B528>0.0 ? p->B528 : (p->wTp>0.0 ? p->wTp : p->wT);
    const double hd = p->wd;
    
    ab2_w = Tr>0.0 ? 2.0*PI/Tr : 0.0;
    
    // linear dispersion, Newton from the shallow-water wave number
    double k = ab2_w/sqrt(g*hd);
    for(int q=0; q<50; ++q)
    {
    const double th = tanh(k*hd);
    k -= (g*k*th - ab2_w*ab2_w)/(g*th + g*k*hd*(1.0-th*th));
    }
    ab2_k = k;
    
    ab2_tau_in  = p->B527>0.0 ? p->B527 : 0.25*Tr;
    ab2_tau_out = p->B526>0.0 ? p->B526 : 0.25*Tr;
    
    if(p->mpirank==0)
    {
    if(ab2_in)
    cout<<"iowave: two-band active absorption at the wavemaker (B 98 7): low pass over "<<ab2_tau_in<<" s, short part with the profile of omega_p "<<ab2_w<<" k_p "<<ab2_k<<endl;
    if(ab2_out)
    cout<<"iowave: two-band active beach (B 99 7): low pass over "<<ab2_tau_out<<" s, short part with the profile of omega_p "<<ab2_w<<" k_p "<<ab2_k<<endl;
    }
}

// first-order low pass of the deviation of column n, x_L <- x_L + dt/tau (x - x_L), once per step
double iowave::ab2_lowpass(lexer *p, std::vector<double> &f, std::vector<int> &c, int n, double x, double tau)
{
    if((int)f.size()<=n)
    {
    f.resize(n+1,0.0);
    c.resize(n+1,-1);
    }
    
    if(c[n]!=p->count)
    {
    f[n] += MIN(p->dt/tau,1.0)*(x - f[n]);
    c[n] = p->count;
    }
    
    return f[n];
}

// velocity per unit elevation of the linear peak wave at depth hd, level z (still water 0)
double iowave::ab2_prof(double hd, double z)
{
    const double kh = ab2_k*hd;
    
    if(kh<1.0e-6)
    return sqrt(9.81/MAX(hd,1.0e-6));
    
    return ab2_w*cosh(ab2_k*MAX(hd+z,0.0))/sinh(kh);
}
