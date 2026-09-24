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
            f(i,j) = (1.0-relax4_wg(i,j))*ramp(p)*eta(i,j) + relax4_wg(i,j) * f(i,j);
            ++count;
            }
		}
        
		
		// Numerical Beach
		if(p->B99==1 || p->B99==2)
		{
            // Zone 2
            if(p->A10!=3 || p->A348==1 || p->A348==2)
            if(db<1.0e20)
            {
            if(p->wet[IJ]==1)
            f(i,j) = relax4_nb(i,j)*f(i,j);
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
            U(i,j)  = (1.0-relax4_wg(i,j))*ramp(p)*uval[count] + relax4_wg(i,j)*U(i,j);
            UH(i,j) = (1.0-relax4_wg(i,j))*ramp(p)*uval[count]*WL(i,j) + relax4_wg(i,j)*UH(i,j);
            }
            ++count;
            }
		}
		
		// Numerical Beach
        if(p->B99==1 || p->B99==2)
		{
            // Zone 2
            if(db<1.0e20)
            {
            U(i,j)  = relax4_nb(i,j)*U(i,j);
            UH(i,j) = relax4_nb(i,j)*UH(i,j);
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
            V(i,j)  = (1.0-relax4_wg(i,j))*ramp(p)*vval[count] + relax4_wg(i,j)*V(i,j);
            VH(i,j) = (1.0-relax4_wg(i,j))*ramp(p)*vval[count]*WL(i,j) + relax4_wg(i,j)*VH(i,j);
            }
            ++count;
            }
		}
		
		// Numerical Beach
        if(p->B99==1 || p->B99==2)
		{
            // Zone 2
            if(db<1.0e20)
            {
            V(i,j)  = relax4_nb(i,j)*V(i,j);
            VH(i,j) = relax4_nb(i,j)*VH(i,j);
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
        if(p->B99==1 || p->B99==2)
		{
            // Zone 2
            if(db<1.0e20)
            {
            W(i,j)  = relax4_nb(i,j)*W(i,j);
            WH(i,j) = relax4_nb(i,j)*WH(i,j);
            }
        }
    }
    
    p->wavecalctime+=pgc->timer()-starttime;
}

void iowave::ws_relax(lexer *p, ghostcell *pgc, slice &W, slice &bed, slice &eta)
{
    starttime=pgc->timer();
    
	double wval=0.0;
    
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

        z=eta(i,j);

        wval = wave_w(p,pgc,xg,yg,z);
        
        
        // Wave Generation
		if(p->B98==2 && w_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            W(i,j) = (1.0-relax4_wg(i,j))*ramp(p)*wval + relax4_wg(i,j)*W(i,j);
		}
		
		// Numerical Beach
        if(p->B99==1 || p->B99==2)
		{
            // Zone 2
            if(db<1.0e20)
            W(i,j) = relax4_nb(i,j)*W(i,j);
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
		if(p->B99==1 || p->B99==2)
		{
            // Zone 2
            if(db<1.0e20)
            f(i,j) = relax4_nb(i,j)*f(i,j);
        }
    }
    
    p->wavecalctime+=pgc->timer()-starttime;
}
