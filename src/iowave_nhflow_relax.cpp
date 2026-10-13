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
#include"fdm_nhf.h"
#include"ghostcell.h"

void iowave::WL_relax(lexer *p, ghostcell *pgc, slice &WL, slice &depth)
{
    starttime=pgc->timer();
    
    // beach zones given as B 520 method 2 relax the water level too, like B 99 1 / 2
    // (and like FNPF); before, only their velocities were relaxed
    const bool beach_wl = p->B99==1 || p->B99==2 || zones.user_beach() || beach_lp;
    
	count=0;
    if(!zc_built) zonecol_build(p);
    ZSLICELOOP4
    {
		dg = distgen(p);
		db = distbeach(p);
        

		// Wave Generation
        if(p->B98==2 && h_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            { 
            WETDRYDEEP
            {
            // with a background (B 523): background + ramped waves
            const int b = bg_on ? gen_bg(p) : -1;
            
            if(b<0 && !gen_lp)
            WL(i,j) = (1.0-relax4_wg(i,j))*ramp(p)*(eta(i,j) + depth(i,j)) + relax4_wg(i,j) * WL(i,j);
            
            // combined generation (B 98 6): target plus the low-passed deviation
            if(b<0 && gen_lp)
            {
            const double tg = ramp(p)*(eta(i,j) + depth(i,j));
            WL(i,j) = (1.0-relax4_wg(i,j))*(tg + gen_dev(p,lpg_wl,p->imax*p->jmax,IJ,WL(i,j)-tg)) + relax4_wg(i,j) * WL(i,j);
            }
            
            if(b>=0)
            WL(i,j) = (1.0-relax4_wg(i,j))*(depth(i,j) + bgs.eta(b,p->pos_x(),p->pos_y()) + ramp(p)*eta(i,j)) + relax4_wg(i,j) * WL(i,j);
            }
            ++count;
            }
		}
        
		
		// Numerical Beach
		if(beach_wl)
		{
            // Zone 2
            if(db<1.0e20)
            {
            if(p->wet[IJ]==1)
            {
            // with a background (B 523): relax to the background level instead of still water
            const int b = bg_on ? beach_bg(p) : -1;
            
            if(b<0 && beach_lp)
            WL(i,j) = (1.0-relax4_nb(i,j))*beach_target(p,lp_wl,p->imax*p->jmax,IJ,depth(i,j),WL(i,j)) + relax4_nb(i,j)*WL(i,j);
            else
            if(b<0)
            WL(i,j) = (1.0-relax4_nb(i,j))*depth(i,j) + relax4_nb(i,j)*WL(i,j);
            
            if(b>=0)
            WL(i,j) = (1.0-relax4_nb(i,j))*(depth(i,j) + bgs.eta(b,p->pos_x(),p->pos_y())) + relax4_nb(i,j)*WL(i,j);
            }
            }
        }
    }
    p->wavecalctime+=pgc->timer()-starttime;
}

void iowave::U_relax(lexer *p, ghostcell *pgc, double *U, double *UH)
{
    starttime=pgc->timer();
    
    count=0;
    if(!zc_built) zonecol_build(p);
    ZLOOP
    {
         dg = distgen(p);
		db = distbeach(p);
        
		// Wave Generation
		if(p->B98==2 && u_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            WETDRYDEEP
            {
            const int b = bg_on ? gen_bg(p) : -1;
            
            if(b<0 && !gen_lp)
            {
            U[IJK]  = (1.0-relax4_wg(i,j))*ramp(p)*uval[count] + relax4_wg(i,j)*U[IJK];
            UH[IJK] = (1.0-relax4_wg(i,j))*ramp(p)*UHval[count] + relax4_wg(i,j)*UH[IJK];
            }
            
            // combined generation (B 98 6): target plus the low-passed deviation
            if(b<0 && gen_lp)
            {
            const int nn = p->imax*p->jmax*(p->kmax+2);
            const double tu  = ramp(p)*uval[count];
            const double tuh = ramp(p)*UHval[count];
            U[IJK]  = (1.0-relax4_wg(i,j))*(tu + gen_dev(p,lpg_u,nn,IJK,U[IJK]-tu)) + relax4_wg(i,j)*U[IJK];
            UH[IJK] = (1.0-relax4_wg(i,j))*(tuh + gen_dev(p,lpg_uh,nn,IJK,UH[IJK]-tuh)) + relax4_wg(i,j)*UH[IJK];
            }
            
            if(b>=0)
            {
            double ub,vb;
            bgs.vel(b,col_h0[IJ],p->pos_x(),p->pos_y(),ub,vb);
            if(bgs.profiled(b))
            {
            const double fp = bg_prof(p,b,col_h0[IJ]+bgs.eta(b,p->pos_x(),p->pos_y()));
            ub *= fp;
            vb *= fp;
            }
            const double ut = ub + ramp(p)*(uval[count]-p->Ui);
            const double ht = col_h0[IJ] + bgs.eta(b,p->pos_x(),p->pos_y()) + ramp(p)*eta(i,j);
            U[IJK]  = (1.0-relax4_wg(i,j))*ut + relax4_wg(i,j)*U[IJK];
            UH[IJK] = (1.0-relax4_wg(i,j))*ht*ut + relax4_wg(i,j)*UH[IJK];
            }
            }
            ++count;
            }
		}
		
		// Numerical Beach
        if(p->B99==1||p->B99==2||beach_relax==1)
		{
            // Zone 2
            if(db<1.0e20)
            {
            // B 97 1: waves on a current, relax to the inflow current instead of to rest
            if(p->B97==1)
            {
            U[IJK]  = relax4_nb(i,j)*U[IJK]  + (1.0-relax4_nb(i,j))*ramp(p)*p->Ui;
            UH[IJK] = relax4_nb(i,j)*UH[IJK] + (1.0-relax4_nb(i,j))*ramp(p)*p->Ui*p->WL[IJ];
            }
            
            const int b = bg_on ? beach_bg(p) : -1;
            
            if(p->B97==0 && b<0)
            {
            const int nn = p->imax*p->jmax*(p->kmax+2);
            if(beach_lp)
            U[IJK] = relax4_nb(i,j)*U[IJK] + (1.0-relax4_nb(i,j))*beach_target(p,lp_u,nn,IJK,0.0,U[IJK]);
            else
            U[IJK] = relax4_nb(i,j)*U[IJK];
            if(beach_lp)
            UH[IJK] = relax4_nb(i,j)*UH[IJK] + (1.0-relax4_nb(i,j))*beach_target(p,lp_uh,nn,IJK,0.0,UH[IJK]);
            else
            UH[IJK] = relax4_nb(i,j)*UH[IJK];
            }
            
            // background (B 523): relax to its current
            if(p->B97==0 && b>=0)
            {
            double ub,vb;
            bgs.vel(b,col_h0[IJ],p->pos_x(),p->pos_y(),ub,vb);
            if(bgs.profiled(b))
            {
            const double fp = bg_prof(p,b,col_h0[IJ]+bgs.eta(b,p->pos_x(),p->pos_y()));
            ub *= fp;
            vb *= fp;
            }
            U[IJK]  = relax4_nb(i,j)*U[IJK]  + (1.0-relax4_nb(i,j))*ub;
            UH[IJK] = relax4_nb(i,j)*UH[IJK] + (1.0-relax4_nb(i,j))*ub*(col_h0[IJ]+bgs.eta(b,p->pos_x(),p->pos_y()));
            }
            }
        }
    }
    p->wavecalctime+=pgc->timer()-starttime;
}

void iowave::V_relax(lexer *p, ghostcell *pgc, double *V, double *VH)
{ 
    starttime=pgc->timer();
    
    count=0;
    if(!zc_built) zonecol_build(p);
    if(p->j_dir==1)
    ZLOOP
    {
         dg = distgen(p);
		db = distbeach(p);
        
		// Wave Generation
		if(p->B98==2 && v_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            WETDRYDEEP
            {
            const int b = bg_on ? gen_bg(p) : -1;
            
            if(b<0 && !gen_lp)
            {
            V[IJK]  = (1.0-relax4_wg(i,j))*ramp(p)*vval[count] + relax4_wg(i,j)*V[IJK];
            VH[IJK] = (1.0-relax4_wg(i,j))*ramp(p)*VHval[count] + relax4_wg(i,j)*VH[IJK];
            }
            
            // combined generation (B 98 6): target plus the low-passed deviation
            if(b<0 && gen_lp)
            {
            const int nn = p->imax*p->jmax*(p->kmax+2);
            const double tv  = ramp(p)*vval[count];
            const double tvh = ramp(p)*VHval[count];
            V[IJK]  = (1.0-relax4_wg(i,j))*(tv + gen_dev(p,lpg_v,nn,IJK,V[IJK]-tv)) + relax4_wg(i,j)*V[IJK];
            VH[IJK] = (1.0-relax4_wg(i,j))*(tvh + gen_dev(p,lpg_vh,nn,IJK,VH[IJK]-tvh)) + relax4_wg(i,j)*VH[IJK];
            }
            
            if(b>=0)
            {
            double ub,vb;
            bgs.vel(b,col_h0[IJ],p->pos_x(),p->pos_y(),ub,vb);
            if(bgs.profiled(b))
            {
            const double fp = bg_prof(p,b,col_h0[IJ]+bgs.eta(b,p->pos_x(),p->pos_y()));
            ub *= fp;
            vb *= fp;
            }
            const double vt = vb + ramp(p)*vval[count];
            const double ht = col_h0[IJ] + bgs.eta(b,p->pos_x(),p->pos_y()) + ramp(p)*eta(i,j);
            V[IJK]  = (1.0-relax4_wg(i,j))*vt + relax4_wg(i,j)*V[IJK];
            VH[IJK] = (1.0-relax4_wg(i,j))*ht*vt + relax4_wg(i,j)*VH[IJK];
            }
            }
            ++count;
            }
		}
		
		// Numerical Beach
		if(p->B99==1||p->B99==2||beach_relax==1)
		{	
            // Zone 2
            if(db<1.0e20)
            {
            const int b = bg_on ? beach_bg(p) : -1;
            
            if(b<0)
            {
            const int nn = p->imax*p->jmax*(p->kmax+2);
            if(beach_lp)
            V[IJK] = relax4_nb(i,j)*V[IJK] + (1.0-relax4_nb(i,j))*beach_target(p,lp_v,nn,IJK,0.0,V[IJK]);
            else
            V[IJK] = relax4_nb(i,j)*V[IJK];
            if(beach_lp)
            VH[IJK] = relax4_nb(i,j)*VH[IJK] + (1.0-relax4_nb(i,j))*beach_target(p,lp_vh,nn,IJK,0.0,VH[IJK]);
            else
            VH[IJK] = relax4_nb(i,j)*VH[IJK];
            }
            
            if(b>=0)
            {
            double ub,vb;
            bgs.vel(b,col_h0[IJ],p->pos_x(),p->pos_y(),ub,vb);
            if(bgs.profiled(b))
            {
            const double fp = bg_prof(p,b,col_h0[IJ]+bgs.eta(b,p->pos_x(),p->pos_y()));
            ub *= fp;
            vb *= fp;
            }
            V[IJK]  = relax4_nb(i,j)*V[IJK]  + (1.0-relax4_nb(i,j))*vb;
            VH[IJK] = relax4_nb(i,j)*VH[IJK] + (1.0-relax4_nb(i,j))*vb*(col_h0[IJ]+bgs.eta(b,p->pos_x(),p->pos_y()));
            }
            }
            
        }
    }
    p->wavecalctime+=pgc->timer()-starttime;
}

void iowave::W_relax(lexer *p, ghostcell *pgc, double *W, double *WH)
{   
    starttime=pgc->timer();
    
    count=0;
    if(!zc_built) zonecol_build(p);
    ZLOOP
    {
        dg = distgen(p);
        db = distbeach(p);
        
		// Wave Generation
		if(p->B98==2 && w_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            WETDRYDEEP
            {
            const int b = bg_on ? gen_bg(p) : -1;
            
            if(b<0 && !gen_lp)
            {
            W[IJK]  = (1.0-relax4_wg(i,j))*ramp(p)*wval[count] + relax4_wg(i,j)*W[IJK];
            WH[IJK] = (1.0-relax4_wg(i,j))*ramp(p)*WHval[count] + relax4_wg(i,j)*WH[IJK];
            }
            
            // combined generation (B 98 6): target plus the low-passed deviation
            if(b<0 && gen_lp)
            {
            const int nn = p->imax*p->jmax*(p->kmax+2);
            const double tv  = ramp(p)*wval[count];
            const double tvh = ramp(p)*WHval[count];
            W[IJK]  = (1.0-relax4_wg(i,j))*(tv + gen_dev(p,lpg_w,nn,IJK,W[IJK]-tv)) + relax4_wg(i,j)*W[IJK];
            WH[IJK] = (1.0-relax4_wg(i,j))*(tvh + gen_dev(p,lpg_wh,nn,IJK,WH[IJK]-tvh)) + relax4_wg(i,j)*WH[IJK];
            }
            
            if(b>=0)
            {
            const double ht = col_h0[IJ] + bgs.eta(b,p->pos_x(),p->pos_y()) + ramp(p)*eta(i,j);
            W[IJK]  = (1.0-relax4_wg(i,j))*ramp(p)*wval[count] + relax4_wg(i,j)*W[IJK];
            WH[IJK] = (1.0-relax4_wg(i,j))*ht*ramp(p)*wval[count] + relax4_wg(i,j)*WH[IJK];
            }
            }
            ++count;
            }
		}
		
		// Numerical Beach
        if(p->B99==1||p->B99==2||beach_relax==1)
		{
            // Zone 2
            if(db<1.0e20)
            {
            const int nn = p->imax*p->jmax*(p->kmax+2);
            if(beach_lp)
            W[IJK] = relax4_nb(i,j)*W[IJK] + (1.0-relax4_nb(i,j))*beach_target(p,lp_w,nn,IJK,0.0,W[IJK]);
            else
            W[IJK] = relax4_nb(i,j)*W[IJK];
            if(beach_lp)
            WH[IJK] = relax4_nb(i,j)*WH[IJK] + (1.0-relax4_nb(i,j))*beach_target(p,lp_wh,nn,IJK,0.0,WH[IJK]);
            else
            WH[IJK] = relax4_nb(i,j)*WH[IJK];
            }
        }
    }
    p->wavecalctime+=pgc->timer()-starttime;		
}

void iowave::P_relax(lexer *p, ghostcell *pgc, double *P)
{
    starttime=pgc->timer();
    if(!zc_built) zonecol_build(p);
    ZFLOOP
    {
        dg = distgen(p);
        db = distbeach(p);
        
        // Numerical Beach
        if(p->B99==1||p->B99==2||beach_relax==1)
        {            
            // Zone 2
            if(db<1.0e20)

            P[FIJK] = relax4_nb(i,j)*P[FIJK];
        }
    }	
    p->wavecalctime+=pgc->timer()-starttime;
}

void iowave::turb_relax_nhflow(lexer *p, fdm_nhf *d, ghostcell *pgc, double *F)
{
    starttime=pgc->timer();
    
    if(!zc_built) zonecol_build(p);
    ZLOOP
    {
        dg = distgen(p);    
        db = distbeach(p);

		// Wave Generation
		if(p->B98==2 && u_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            F[IJK] = relax4_wg(i,j)*F[IJK];

		}
        
        // Numerical Beach
        if(p->B99==1||p->B99==2||beach_relax==1)
		{
            // Zone 2
            if(db<1.0e20)
            {
            F[IJK] = relax4_nb(i,j)*F[IJK];
            }
        }
    }
    
    p->wavecalctime+=pgc->timer()-starttime;
}

// combined beach (B 99 6): low-passed state of a cell, updated once per step,
// f <- f + dt/tau (x - f), starting from the still-water value f0
double iowave::beach_target(lexer *p, lowpass &lp, int n, int idx, double f0, double x)
{
    if(lp.f.empty())
    {
        lp.f.assign(n,0.0);
        lp.c.assign(n,-1);
    }
    
    if(lp.c[idx]<0)
    lp.f[idx] = f0;
    
    if(lp.c[idx]!=p->count)
    {
        lp.f[idx] += MIN(p->dt/beach_tau,1.0)*(x - lp.f[idx]);
        lp.c[idx] = p->count;
    }
    
    return lp.f[idx];
}

// combined generation (B 98 6): deviation of the state x from the target T, low-passed over gen_tau and
// updated once per step and cell, d <- d + dt/tau ((x - T) - d), starting from 0. The zone relaxes to
// T + d: components of the deviation much longer than tau (free long waves) pass the zone to the
// absorbing edge, shorter ones are relaxed to the target.
double iowave::gen_dev(lexer *p, lowpass &lp, int n, int idx, double dev)
{
    if(lp.f.empty())
    {
        lp.f.assign(n,0.0);
        lp.c.assign(n,-1);
    }
    
    if(lp.c[idx]!=p->count)
    {
        lp.f[idx] += MIN(p->dt/gen_tau,1.0)*(dev - lp.f[idx]);
        lp.c[idx] = p->count;
    }
    
    return lp.f[idx];
}
