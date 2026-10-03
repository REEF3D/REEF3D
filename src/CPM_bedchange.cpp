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
Authors: Hans Bihs, Alexander Hanke
--------------------------------------------------------------------*/

#include"CPM.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"sediment_fdm.h"

// bed interface from the particles: iso-surface of the solid volume fraction theta_bed = Q 26 (1-S 24)
// first order level set estimate, reinitialised afterwards by reinitopo
void CPM::topo_update(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s)
{
    double gx,gy,gz,grad,h,val;
    double Tm,Tp;
    
    BASELOOP
    {
        h = p->j_dir==1 ? (1.0/3.0)*(p->DXN[IP]+p->DYN[JP]+p->DZN[KP]) : 0.5*(p->DXN[IP]+p->DZN[KP]);
        
        Tm = wallcell(p,i-1,j,k) ? Ts(i,j,k) : Ts(i-1,j,k);
        Tp = wallcell(p,i+1,j,k) ? Ts(i,j,k) : Ts(i+1,j,k);
        gx = (Tp-Tm)/(p->DXP[IM1]+p->DXP[IP]);
        
        gy=0.0;
        if(p->j_dir==1)
        {
        Tm = wallcell(p,i,j-1,k) ? Ts(i,j,k) : Ts(i,j-1,k);
        Tp = wallcell(p,i,j+1,k) ? Ts(i,j,k) : Ts(i,j+1,k);
        gy = (Tp-Tm)/(p->DYP[JM1]+p->DYP[JP]);
        }
        
        Tm = wallcell(p,i,j,k-1) ? (1.0-p->S24) : Ts(i,j,k-1);
        Tp = wallcell(p,i,j,k+1) ? Ts(i,j,k) : Ts(i,j,k+1);
        gz = (Tp-Tm)/(p->DZP[KM1]+p->DZP[KP]);
        
        grad = sqrt(gx*gx + gy*gy + gz*gz);
        grad = MAX(grad, theta_bed/(2.0*h));
        
        val = (theta_bed - Ts(i,j,k))/grad;
        
        val = MAX(val,-3.0*h);
        val = MIN(val, 3.0*h);
        
        a->topo(i,j,k) = val;
    }
    
    pgc->start4a(p,a->topo,150);
}

// bed elevation from the topo level set
void CPM::bedzh_update(lexer *p, fdm *a, ghostcell *pgc, sediment_fdm *s)
{
    double h;
    
    ILOOP
    JLOOP
    {
        h = s->bedzh(i,j);
        
        if(a->topo(i,j,0)>=0.0 && p->nb5<0)
        h = p->ZN[0+marge];
        
        KLOOP
        PBASECHECK
        if(k>0 && a->topo(i,j,k-1)<0.0 && a->topo(i,j,k)>=0.0)
            h = -(a->topo(i,j,k-1)*p->DZP[KM1])/(a->topo(i,j,k)-a->topo(i,j,k-1)) + p->pos_z()-p->DZP[KM1];
            
        s->bedzh(i,j)=h;
        a->bed(i,j)=h;
    }
    
    pgc->gcsl_start4(p,s->bedzh,1);
    pgc->gcsl_start4(p,a->bed,50);
}
