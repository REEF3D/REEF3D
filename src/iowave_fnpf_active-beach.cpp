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
    double eta_R,Uc,eta_T,eta_M,hd,z,zabs,wsf,fac,multiplier;
    
		for(n=0;n<p->gcslout_count;++n)
		{
		i=p->gcslout[n][0];
		j=p->gcslout[n][1];
		
		eta_T = 0.0;
		eta_M = eta(i,j); 
		eta_R = eta_M-eta_T;
        
        hd = c->depth(i,j);
        
        Uc=eta_R*sqrt(9.81/p->wd);
        
			FKLOOP 
			{
				z=p->ZSN[FIJK]-p->phimean;
				
				if(p->B99==3)
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
