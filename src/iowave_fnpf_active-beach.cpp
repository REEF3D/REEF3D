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

void iowave::active_beach_fnpf(lexer *p, fdm_fnpf *c, ghostcell* pgc, double *Fi, double *Uin, slice &Fifsf, slice &eta)
{
        double eta_R,Uc,Un,Vc,Wc,eta_T,eta_M,wsf;
		double posx,posy,posz,uvel,vvel,uabs,fx,fy,pval,fp;
        double fxdir,fydir;
		double x=0.0;
		double z=0.0;
		double fac1,fac,multiplier;
		int aa,bb;

		// U / FI
		for(n=0;n<p->gcslout_count;++n)
		{
		i=p->gcslout[n][0];
		j=p->gcslout[n][1];

		eta_T = 0.0;
		eta_M = eta(i,j); 
		eta_R = eta_M-eta_T;

        aa=bb=0;
		if(p->gcslout[n][3]==1)
		aa=-1;
		
		if(p->gcslout[n][3]==4)
		aa=1;
		
		if(p->gcslout[n][3]==3)
		bb=-1;
		
		if(p->gcslout[n][3]==2)
		bb=1;

        fx=1.0;
        Uc=eta_R*sqrt(9.81/p->wd);
        
        
            
			FKLOOP 
			{
				z=p->ZSN[FIJK]-p->phimean;
				
				if(p->B99==3)
				Uc=eta_R*sqrt(9.81/c->depth(i,j));
				
				if(p->B99==4)
				Uc=eta_R*p->ww*(cosh(p->wk*(p->wd+z))/sinh(c->depth(i,j)*c->depth(i,j)));
                
               if(p->B99==5)
               {
                   if(p->pos_z()>p->B123)
                   {
                   fac = (p->pos_z()-p->B123)/(wsf-p->B123);
                   multiplier = 2.0*((wsf)/(wsf-p->B123));
                   Uc =   multiplier*fac*eta_R*sqrt(9.81/c->depth(i,j));
                   }
                   
                   if(p->pos_z()<=p->B123)
                   Uc=0.0;
               }
                   

				Uin[FIp1JK]=Uc*fx;
			}
          
        
        if(p->A329==1)
        { 
        Fifsf(i+1,j) = Fifsf(i,j) + Uc*fx*1.0*p->DXP[IP1];
        Fifsf(i+2,j) = Fifsf(i,j) + Uc*fx*2.0*p->DXP[IP1];
        Fifsf(i+3,j) = Fifsf(i,j) + Uc*fx*3.0*p->DXP[IP1];
        }
        
        if(p->A329>=2)
        {
        Fifsf(i+1,j) = (4.0/3.0)*Fifsf(i,j) - (1.0/3.0)*Fifsf(i-1,j) - (2.0/3.0)*Uc*fx*(-0.5*p->XP[IM1] + 2.0*p->XP[IP] - 1.5*p->XP[IP1]);
        Fifsf(i+2,j) = (4.0/3.0)*Fifsf(i,j) - (1.0/3.0)*Fifsf(i-1,j) - (2.0/3.0)*Uc*fx*(-0.5*p->XP[IM1] + 2.0*p->XP[IP] - 1.5*p->XP[IP2]);
        Fifsf(i+3,j) = (4.0/3.0)*Fifsf(i,j) - (1.0/3.0)*Fifsf(i-1,j) - (2.0/3.0)*Uc*fx*(-0.5*p->XP[IM1] + 2.0*p->XP[IP] - 1.5*p->XP[IP3]);
        }
        
        FKLOOP
        FPCHECK
        {
        Fi[FIp1JK] = Fi[FIJK] + Uc*fx*1.0*p->DXP[IP1];
        Fi[FIp2JK] = Fi[FIJK] + Uc*fx*2.0*p->DXP[IP1];
        Fi[FIp3JK] = Fi[FIJK] + Uc*fx*3.0*p->DXP[IP1];
        }

		}
              
}
