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
#include"fdm_ptf.h"
#include"ghostcell.h"


void iowave::wavegen_precalc_decomp_space_dirichlet_ptf(lexer* p ,ghostcell* pgc)
{
    double fsfloc;
    int qn;

    count=0;
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];
    
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);
		
		// Wave Generation
        if(p->B98==3 && h_switch==1)
        {
            // Zone 1
                for(qn=0;qn<wave_comp;++qn)
                {
                etaval_S_sin[count][qn] = wave_eta_space_sin(p,pgc,xg,yg,qn);
                etaval_S_cos[count][qn] = wave_eta_space_cos(p,pgc,xg,yg,qn);
                }
            ++count;
		}
    }
    pgc->gcsl_start4(p,eta,50);
    
    // Fifsf
    count=0;
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];
    
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);
        
        z = eta(i,j);
		
		// Wave Generation
        if(p->B98==3 && h_switch==1)
        {
            // Zone 1
            if(dg<dist1)
            {
                for(qn=0;qn<wave_comp;++qn)
                {
                Fifsfval_S_sin[count][qn] = wave_fi_space_sin(p,pgc,xg,yg,z,qn);
                Fifsfval_S_cos[count][qn] = wave_fi_space_cos(p,pgc,xg,yg,z,qn);
                }
            ++count;
            }
		}
    }
    
    
    // Uin
        count=0;
		for(n=0;n<p->gcslin_count;n++)
        {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        xg=xgen(p);
        yg=ygen(p);
        x1=xgen1(p);
        y2=ygen2(p);
        
            KLOOP
            PCHECK
            {
            z=p->ZSN[FIJK]-p->phimean;

                for(qn=0;qn<wave_comp;++qn)
                {
                uval_S_sin[count][qn] = wave_u_space_sin(p,pgc,xg,yg,z,qn);
                uval_S_cos[count][qn] = wave_u_space_cos(p,pgc,xg,yg,z,qn);
                }
            ++count;
            }
        }



}


void iowave::ptf_precalc_dirichlet_ini(lexer *p, ghostcell *pgc)
{    
    int dbcount=0;
    
    upt_count=vpt_count=wpt_count=ppt_count=ept_count=dbcount=0;
    
    
    // FI ------------------------------------------------
    LOOP
    {
        db = distbeach(p); 
        
        if(p->B99==1||p->B99==2)
		{
            if(db<1.0e20)
            ++dbcount;
        }
    }	
    
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];
        KLOOP
        ++ept_count;
    }
    
    upt_count=vpt_count=wpt_count=ppt_count=ept_count;
    
    if(p->B89==1)
    {
        if(p->B92==5)
        wave_comp = 5;
        
        if(p->B92==31 || p->B92==41 || p->B92==51)
        wave_comp = p->wN;
    }
      
    // precalc array alloc
    p->Darray(etaval,ept_count);
    p->Darray(Fival,ppt_count);
    p->Darray(Uinval,ppt_count);
    p->Darray(Fifsfval,ept_count);
     p->Darray(Fifsfval0,ept_count);
     p->Darray(Fifsfval1,ept_count);
    p->Darray(uval,upt_count);
    
    p->Darray(rb3val,dbcount);
    
    if(p->B89==1) 
    {
    p->Darray(etaval_S_sin,ept_count,wave_comp);
    p->Darray(Fival_S_sin,ppt_count,wave_comp);
    p->Darray(Fifsfval_S_sin,ept_count,wave_comp);
    p->Darray(uval_S_sin,ppt_count,wave_comp);
    
    p->Darray(etaval_S_cos,ept_count,wave_comp);
    p->Darray(Fival_S_cos,ppt_count,wave_comp);
    p->Darray(Fifsfval_S_cos,ept_count,wave_comp);
    p->Darray(uval_S_cos,ppt_count,wave_comp);
    
    p->Darray(etaval_T_sin,wave_comp);
    p->Darray(Fival_T_sin,wave_comp);
    p->Darray(Fifsfval_T_sin,wave_comp);
    p->Darray(uval_T_sin,wave_comp);
    
    p->Darray(etaval_T_cos,wave_comp);
    p->Darray(Fival_T_cos,wave_comp);
    p->Darray(Fifsfval_T_cos,wave_comp);
    p->Darray(uval_T_cos,wave_comp);
    }
}


void iowave::ptf_precalc_dirichlet(lexer *p, ghostcell *pgc)
{
        count=0;
		for(n=0;n<p->gcslin_count;n++)
        {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        xg=xgen(p);
        yg=ygen(p);
        x1=xgen1(p);
        y2=ygen2(p);
        

        eta(i,j) = wave_eta(p,pgc,xg,yg);
        etaval[count] = eta(i,j);
        
        z = eta(i,j);
        
        time_1=time_0;
        time_0=p->simtime;
        time_n=p->simtime+p->dt;
        Fifsfval1[count] = Fifsfval0[count];
        Fifsfval0[count] = Fifsfval[count];
        
        Fifsfval[count] = wave_u(p,pgc,xg,yg,z);
        ++count;
        }
        
        
        // Uin
        count=0;
		for(n=0;n<p->gcslin_count;n++)
        {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        xg=xgen(p);
        yg=ygen(p);
        x1=xgen1(p);
        y2=ygen2(p);
        
            KLOOP
            PCHECK
            {
            z=p->ZSN[FIJK]-p->phimean;

            Uinval[count] = wave_u(p,pgc,xg,yg,z);
            ++count;
            }
        }
        
        
    count=0;
    ILOOP 
    JLOOP 
    {

		db = distbeach(p);
        
        KLOOP 
        PCHECK
        {
                    
            if(p->B99==1||p->B99==2)
            {
                // Zone 2
                if(db<dist2)
                {
                rb3val[count] = rb3(p,db);
                ++count;
                }
            }
        }
    }
        
}

void iowave::dirichlet_wavegen_ptf(lexer *p, fdm_ptf *a, ghostcell* pgc, double *Fi_, double *Uin_, slice &Fifsf, slice &eta)
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
        
        KLOOP
        PCHECK
        {
        Fi_[Im1JK] = Fi_[IJK] - Uinval[count]*1.0*p->DXP[IM1];
        Fi_[Im2JK] = Fi_[IJK] - Uinval[count]*2.0*p->DXP[IM1];
        Fi_[Im3JK] = Fi_[IJK] - Uinval[count]*3.0*p->DXP[IM1];
        
        ++count;
        }
        
        
        
    }
    
    // Uin
    
    count=0;
    for(n=0;n<p->gcslin_count;n++)
    {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        KLOOP
        PCHECK
        {// add eta guard
        Uin_[Im1JK] = Uinval[count]; 
        
        ++count;
        }
    }
}
