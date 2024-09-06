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
Author: Fabian Knoblauch
--------------------------------------------------------------------*/

#include"VOF_PLIC.h"
#include"gradient.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"convection.h"
#include"solver.h"
#include"ghostcell.h"
#include"ioflow.h"
#include"fluid_update_vof.h"
#include"heat.h"
#include"hires.h"
#include"weno_hj.h"
#include"hric.h"

void VOF_PLIC::transportPhi_Bonn
(
    fdm* a,
    lexer* p,
    int nSweep,
    int sweep
)
{
    double Gp, Gm, uip, uim, dtdxi;
    if(sweep==0)
    {
        uip=a->u(i,j,k);
        uim=a->u(i-1,j,k);
        dtdxi=p->dt/p->DXN[IP];
        if(uip>=0.0)
            Gp=uip*(phistep(i,j,k)+0.5*p->DXP[IP]*(1.0-uip*(p->dt/p->DXP[IP]))*(phistep(i+1,j,k)-phistep(i,j,k))/(p->DXP[IP]));
        else
            Gp=uip*(phistep(i,j,k)+0.5*p->DXP[IP]*(1.0+uip*(p->dt/p->DXP[IP]))*(phistep(i+1,j,k)-phistep(i,j,k))/(p->DXP[IP]));
            
        if(uim<=0.0)
            Gm=uim*(phistep(i,j,k)+0.5*p->DXP[IM1]*(-1.0-uim*(p->dt/p->DXP[IM1]))*(phistep(i,j,k)-phistep(i-1,j,k))/(p->DXP[IM1]));
        else
            Gm=uim*(phistep(i,j,k)+0.5*p->DXP[IM1]*(-1.0+uim*(p->dt/p->DXP[IM1]))*(phistep(i,j,k)-phistep(i-1,j,k))/(p->DXP[IM1]));
    }
    else if(sweep==1)
    {
        uip=a->v(i,j,k);
        uim=a->v(i-1,j,k);
        dtdxi=p->dt/p->DYN[JP];
        if(uip>=0.0)
            Gp=uip*(phistep(i,j,k)+0.5*p->DYP[JP]*(1.0-uip*(p->dt/p->DYP[JP]))*(phistep(i,j+1,k)-phistep(i,j,k))/(p->DYP[JP]));
        else
            Gp=uip*(phistep(i,j,k)+0.5*p->DYP[JP]*(1.0+uip*(p->dt/p->DYP[JP]))*(phistep(i,j+1,k)-phistep(i,j,k))/(p->DYP[JP]));
            
        if(uim<=0.0)
            Gm=uim*(phistep(i,j,k)+0.5*p->DYP[JM1]*(-1.0-uim*(p->dt/p->DYP[JM1]))*(phistep(i,j,k)-phistep(i,j-1,k))/(p->DYP[JM1]));
        else
            Gm=uim*(phistep(i,j,k)+0.5*p->DYP[JM1]*(-1.0+uim*(p->dt/p->DYP[JM1]))*(phistep(i,j,k)-phistep(i,j-1,k))/(p->DYP[JM1]));
            
    }
    else
    {
        uip=a->w(i,j,k);
        uim=a->w(i,j,k-1);
        dtdxi=p->dt/p->DZN[KP];
        if(uip>=0.0)
            Gp=uip*(phistep(i,j,k)+0.5*p->DZP[KP]*(1.0-uip*(p->dt/p->DZP[IP]))*(phistep(i,j,k+1)-phistep(i,j,k))/(p->DZP[IP]));
        else
            Gp=uip*(phistep(i,j,k)+0.5*p->DZP[KP]*(1.0+uip*(p->dt/p->DZP[IP]))*(phistep(i,j,k+1)-phistep(i,j,k))/(p->DZP[IP]));
            
        if(uim<=0.0)
            Gm=uim*(phistep(i,j,k)+0.5*p->DZP[IM1]*(-1.0-uim*(p->dt/p->DZP[IM1]))*(phistep(i,j,k)-phistep(i,j,k-1))/(p->DZP[IM1]));
        else
            Gm=uim*(phistep(i,j,k)+0.5*p->DZP[IM1]*(-1.0+uim*(p->dt/p->DZP[IM1]))*(phistep(i,j,k)-phistep(i,j,k-1))/(p->DZP[IM1]));
    }
    
    if(nSweep==0)
        phiS0(i,j,k)=(a->phi(i,j,k)-dtdxi*(Gp-Gm))/(1.0-dtdxi*(uip-uim));
    else if(nSweep==1)
        phiS1(i,j,k)=phiS0(i,j,k)*(1.0+dtdxi*(uip-uim))-dtdxi*(Gp-Gm);
    else
        phiS2(i,j,k)=phiS1(i,j,k)+phiS0(i,j,k)*dtdxi*(uip-uim)-dtdxi*(Gp-Gm);
}

void VOF_PLIC::transportVOF_Bonn
(
    fdm* a,
    lexer* p,
    int nSweep,
    int sweep
)
{
    double Gp, Gm, uip, uim, dtdxi;
    LOOP
    {
        if(sweep==0)
        {
        
            uip=a->u(i,j,k);
            uim=a->u(i-1,j,k);
            dtdxi=p->dt/p->DXN[IP];
            if(uip>=0.0)
                Gp=V_w_p(i,j,k)/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
            else
                Gp=V_w_m(i+1,j,k)/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
                
            if(uim<=0.0)
                Gm=V_w_m(i,j,k)/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
            else
                Gm=V_w_p(i-1,j,k)/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
        }
        else if(sweep==1)
        {
            uip=a->v(i,j,k);
            uim=a->v(i,j-1,k);
            dtdxi=p->dt/p->DYN[JP];
            
            if(uip>=0.0)
                Gp=V_w_p(i,j,k)/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
            else
                Gp=V_w_m(i,j+1,k)/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
                
            if(uim<=0.0)
                Gm=V_w_m(i,j,k)/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
            else
                Gm=V_w_p(i,j-1,k)/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
           
        }
        else
        {
            uip=a->w(i,j,k);
            uim=a->w(i,j,k-1);
            dtdxi=p->dt/p->DZN[KP];
            if(uip>=0.0)
                Gp=V_w_p(i,j,k)/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
            else
                Gp=V_w_m(i,j,k+1)/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
            
            if(uim<=0.0)
                Gm=V_w_m(i,j,k)/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
            else
                Gm=V_w_p(i,j,k-1)/(p->DXN[IP]*p->DYN[JP]*p->DZN[KP]);
                
        }
        
        if(nSweep==0)
        {
            vofS0(i,j,k)=(vofstep(i,j,k)-(Gp-Gm))/(1.0-dtdxi*(uip-uim));
        }
        else if(nSweep==1)
        {
            vofS1(i,j,k)=vofS0(i,j,k)*(1.0+dtdxi*(uip-uim))-(Gp-Gm);
        }
        else 
        {
            vofS2(i,j,k)=vofS1(i,j,k)+vofS0(i,j,k)*dtdxi*(uip-uim)-(Gp-Gm);
        }
    }
}