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
#include"ghostcell.h"
#include"slice.h"

void iowave::eta_relax(lexer *p, ghostcell *pgc, slice &f)
{
    starttime=pgc->timer();
    
	count=0;
    // only relaxation-zone cells, SLICELOOP4 order (cached geometry)
    if(!rz4_built) relaxzone4_build(p);
    for(size_t rzq=0; rzq<rz4_i.size(); ++rzq)
    {
        i = rz4_i[rzq];
        j = rz4_j[rzq];
        dg = rz4_dg[rzq];
        db = rz4_db[rzq];
        

		// Wave Generation
        if(p->B98==2 && h_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            { 
            // with a background (B 523): background + ramped waves
            const int bi = bg_on ? gen_bg(p) : -1;
            const double eb = bi>=0 ? bgs.eta(bi,p->XP[IP],p->YP[JP]) : 0.0;
            
            if(bi<0)
            f(i,j) = (1.0-relax4_wg(i,j))*ramp(p)*eta(i,j) + relax4_wg(i,j) * f(i,j);
            else
            f(i,j) = (1.0-relax4_wg(i,j))*(eb + ramp(p)*eta(i,j)) + relax4_wg(i,j) * f(i,j);
            ++count;
            }
		}
        
		
		// Numerical Beach (B 99 1 / 2; B 520 beach zones in SFLOW, FNPF as before)
		if(p->B99==1 || p->B99==2 || (beach_relax==1 && p->A10==2))
		{
            // Zone 2
            if(p->A10!=3 || p->A348==1 || p->A348==2)
            if(db<1.0e20)
            {
            if(p->wet[IJ]==1)
            {
            // with a background (B 523): relax to the background level instead of still water
            const int bi = bg_on ? beach_bg(p) : -1;
            const double eb = bi>=0 ? bgs.eta(bi,p->XP[IP],p->YP[JP]) : 0.0;
            
            // combined beach (B 99 6): towards the low-passed level
            const double tg = beach_lp ? beach_target(p,lp_wl,p->imax*p->jmax,IJ,eb,f(i,j)) : eb;
            
            f(i,j) = relax4_nb(i,j)*f(i,j) + (1.0-relax4_nb(i,j))*tg;
            }
            }
        }
    }
    
    p->wavecalctime+=pgc->timer()-starttime;
}

// SFLOW (cell-centred HLL): um/vm/wm_relax(U,UH,WL) relax the depth-averaged velocity
// and the conserved momentum towards the precalculated wave velocities at the cell
// centres (wavegen_2D_precalc); WL is the relaxed water depth of the stage.

void iowave::um_relax(lexer *p, ghostcell *pgc, slice &U, slice &UH, slice &WL)
{
    starttime=pgc->timer();
    
    count=0;
    
    // only relaxation-zone cells, SLICELOOP4 order (cached geometry)
    if(!rz4_built) relaxzone4_build(p);

    for(size_t rzq=0; rzq<rz4_i.size(); ++rzq)
    {
        i = rz4_i[rzq];
        j = rz4_j[rzq];
        dg = rz4_dg[rzq];
        db = rz4_db[rzq];
        
        // Wave Generation
		if(p->B98==2 && u_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            if(p->wet[IJ]==1)
            {
            // with a background (B 523): background current + ramped waves
            double ub=0.0, vb=0.0;
            const int bi = bg_on ? gen_bg(p) : -1;
            if(bi>=0)
            bgs.vel(bi,col_h0[IJ],p->XP[IP],p->YP[JP],ub,vb);
            
            if(bi<0)
            {
            U(i,j)  = (1.0-relax4_wg(i,j))*ramp(p)*uval[count] + relax4_wg(i,j)*U(i,j);
            UH(i,j) = (1.0-relax4_wg(i,j))*ramp(p)*uval[count]*WL(i,j) + relax4_wg(i,j)*UH(i,j);
            }
            else
            {
            const double tg = ub + ramp(p)*uval[count];
            U(i,j)  = (1.0-relax4_wg(i,j))*tg + relax4_wg(i,j)*U(i,j);
            UH(i,j) = (1.0-relax4_wg(i,j))*tg*WL(i,j) + relax4_wg(i,j)*UH(i,j);
            }
            }
            ++count;
            }
		}
		
		// Numerical Beach
        if(p->B99==1 || p->B99==2 || beach_relax==1)
		{
            // Zone 2: with a background (B 523) towards the background current
            if(db<1.0e20)
            {
            double ub=0.0, vb=0.0;
            const int bi = bg_on ? beach_bg(p) : -1;
            if(bi>=0)
            bgs.vel(bi,col_h0[IJ],p->XP[IP],p->YP[JP],ub,vb);
            
            // combined beach (B 99 6): towards the low-passed state
            const int nn = p->imax*p->jmax;
            const double tg  = beach_lp ? beach_target(p,lp_u,nn,IJ,ub,U(i,j)) : ub;
            const double tgh = beach_lp ? beach_target(p,lp_uh,nn,IJ,ub*WL(i,j),UH(i,j)) : ub*WL(i,j);
            
            U(i,j)  = relax4_nb(i,j)*U(i,j) + (1.0-relax4_nb(i,j))*tg;
            UH(i,j) = relax4_nb(i,j)*UH(i,j) + (1.0-relax4_nb(i,j))*tgh;
            }
        }
    }
    
    p->wavecalctime+=pgc->timer()-starttime;
}

void iowave::vm_relax(lexer *p, ghostcell *pgc, slice &V, slice &VH, slice &WL)
{
    starttime=pgc->timer();
    
    count=0;
    
    if(p->j_dir==1)
    {
    // only relaxation-zone cells, SLICELOOP4 order (cached geometry)
    if(!rz4_built) relaxzone4_build(p);

    for(size_t rzq=0; rzq<rz4_i.size(); ++rzq)
    {
        i = rz4_i[rzq];
        j = rz4_j[rzq];
        dg = rz4_dg[rzq];
        db = rz4_db[rzq];
        
        // Wave Generation
		if(p->B98==2 && v_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            if(p->wet[IJ]==1)
            {
            // with a background (B 523): background current + ramped waves
            double ub=0.0, vb=0.0;
            const int bi = bg_on ? gen_bg(p) : -1;
            if(bi>=0)
            bgs.vel(bi,col_h0[IJ],p->XP[IP],p->YP[JP],ub,vb);
            
            if(bi<0)
            {
            V(i,j)  = (1.0-relax4_wg(i,j))*ramp(p)*vval[count] + relax4_wg(i,j)*V(i,j);
            VH(i,j) = (1.0-relax4_wg(i,j))*ramp(p)*vval[count]*WL(i,j) + relax4_wg(i,j)*VH(i,j);
            }
            else
            {
            const double tg = vb + ramp(p)*vval[count];
            V(i,j)  = (1.0-relax4_wg(i,j))*tg + relax4_wg(i,j)*V(i,j);
            VH(i,j) = (1.0-relax4_wg(i,j))*tg*WL(i,j) + relax4_wg(i,j)*VH(i,j);
            }
            }
            ++count;
            }
		}
		
		// Numerical Beach
        if(p->B99==1 || p->B99==2 || beach_relax==1)
		{
            // Zone 2: with a background (B 523) towards the background current
            if(db<1.0e20)
            {
            double ub=0.0, vb=0.0;
            const int bi = bg_on ? beach_bg(p) : -1;
            if(bi>=0)
            bgs.vel(bi,col_h0[IJ],p->XP[IP],p->YP[JP],ub,vb);
            
            // combined beach (B 99 6): towards the low-passed state
            const int nn = p->imax*p->jmax;
            const double tg  = beach_lp ? beach_target(p,lp_v,nn,IJ,vb,V(i,j)) : vb;
            const double tgh = beach_lp ? beach_target(p,lp_vh,nn,IJ,vb*WL(i,j),VH(i,j)) : vb*WL(i,j);
            
            V(i,j)  = relax4_nb(i,j)*V(i,j) + (1.0-relax4_nb(i,j))*tg;
            VH(i,j) = relax4_nb(i,j)*VH(i,j) + (1.0-relax4_nb(i,j))*tgh;
            }
        }
    }
    }
    
    p->wavecalctime+=pgc->timer()-starttime;
}

void iowave::wm_relax(lexer *p, ghostcell *pgc, slice &W, slice &WH, slice &WL)
{
    starttime=pgc->timer();
    
    count=0;
    
    // only relaxation-zone cells, SLICELOOP4 order (cached geometry)
    if(!rz4_built) relaxzone4_build(p);

    for(size_t rzq=0; rzq<rz4_i.size(); ++rzq)
    {
        i = rz4_i[rzq];
        j = rz4_j[rzq];
        dg = rz4_dg[rzq];
        db = rz4_db[rzq];
        
        // Wave Generation
		if(p->B98==2 && w_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            if(p->wet[IJ]==1)
            {
            W(i,j)  = (1.0-relax4_wg(i,j))*ramp(p)*wval[count] + relax4_wg(i,j)*W(i,j);
            WH(i,j) = (1.0-relax4_wg(i,j))*ramp(p)*wval[count]*WL(i,j) + relax4_wg(i,j)*WH(i,j);
            }
            ++count;
            }
		}
		
		// Numerical Beach
        if(p->B99==1 || p->B99==2 || beach_relax==1)
		{
            // Zone 2
            if(db<1.0e20)
            {
            if(beach_lp)
            {
            const int nn = p->imax*p->jmax;
            W(i,j)  = relax4_nb(i,j)*W(i,j)  + (1.0-relax4_nb(i,j))*beach_target(p,lp_w,nn,IJ,0.0,W(i,j));
            WH(i,j) = relax4_nb(i,j)*WH(i,j) + (1.0-relax4_nb(i,j))*beach_target(p,lp_wh,nn,IJ,0.0,WH(i,j));
            }
            else
            {
            W(i,j)  = relax4_nb(i,j)*W(i,j);
            WH(i,j) = relax4_nb(i,j)*WH(i,j);
            }
            }
        }
    }
    
    p->wavecalctime+=pgc->timer()-starttime;
}

void iowave::pm_relax(lexer *p, ghostcell *pgc, slice &f)
{
	starttime=pgc->timer();
    
    // only relaxation-zone cells, SLICELOOP4 order (cached geometry)
    if(!rz4_built) relaxzone4_build(p);
    for(size_t rzq=0; rzq<rz4_i.size(); ++rzq)
    {
        i = rz4_i[rzq];
        j = rz4_j[rzq];
        xg = rz4_xg[rzq];
        yg = rz4_yg[rzq];
        dg = rz4_dg[rzq];
        db = rz4_db[rzq];
		
		// Wave Generation
        if(p->B98==2)
        {
            // Zone 1
            //if(dg<1.0e20)
            //f(i,j) = relax4_wg(i,j) * f(i,j);
		}
		
		// Numerical Beach
		if(p->B99==1 || p->B99==2 || beach_relax==1)
		{
            // Zone 2
            if(db<1.0e20)
            f(i,j) = relax4_nb(i,j)*f(i,j);
        }
    }
    
    p->wavecalctime+=pgc->timer()-starttime;
}
