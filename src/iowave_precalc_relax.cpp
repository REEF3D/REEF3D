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
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"iowave.h"
#include"lexer.h"
#include"ghostcell.h"

void iowave::wavegen_precalc_relax(lexer *p, ghostcell *pgc)
{
    double fsfloc;
    
    // pre-calc every iteration
    count=0;
    SLICELOOP4
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);
		
		// Wave Generation
        if(p->B98==2 && h_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            eta(i,j) = wave_eta(p,pgc,xg,yg);
            etaval[count] = eta(i,j);
            ++count;
            }
		}
    }
    pgc->gcsl_start4(p,eta,50);
    
    count=0;
    ULOOP
    {
		xg = xgen1(p);
        yg = ygen1(p);
        dg = distgen(p);
		db = distbeach(p);
        
        zloc1 = p->pos1_z();
        fsfloc = 0.5*(eta(i,j)+eta(i+1,j)) + p->phimean;
    

        if(zloc1<=fsfloc)
        {
        if(zloc1<=p->phimean)
        z=-(fabs(p->phimean-zloc1));
		
		if(zloc1>p->phimean)
        z=(fabs(p->phimean-zloc1));
        }
        
        if(zloc1>fsfloc)
        z = 0.5*(eta(i,j)+eta(i+1,j));
		
		// Wave Generation
		if(p->B98==2 && u_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            if(zloc1<=fsfloc+epsi)
            uval[count] = wave_u(p,pgc,xg,yg,z) + p->Ui;
            
            if(zloc1>fsfloc+epsi)
            uval[count] = 0.0;
            
            ++count;
            }
		}
    }
		
    count=0;
    VLOOP
    {
        xg = xgen2(p);
        yg = ygen2(p);
        dg = distgen(p);
		db = distbeach(p);
        
        zloc2 = p->pos2_z();
        fsfloc = 0.5*(eta(i,j)+eta(i,j+1)) + p->phimean;
    

        if(zloc2<=fsfloc)
        {
        if(zloc2<=p->phimean)
        z=-(fabs(p->phimean-zloc2));
		
		if(zloc2>p->phimean)
        z=(fabs(p->phimean-zloc2));
        }
        
        if(zloc2>fsfloc)
        z = 0.5*(eta(i,j)+eta(i,j+1));

		// Wave Generation
		if(p->B98==2 && v_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            if(zloc2<=fsfloc+epsi)
            vval[count] = wave_v(p,pgc,xg,yg,z);
            
            if(zloc2>fsfloc+epsi)
            vval[count] = 0.0;
            
            ++count;
            }
		}
    }

    count=0;
    WLOOP
    {
        xg = xgen(p);
        yg = ygen(p);
        dg = distgen(p);
		db = distbeach(p);
        
        zloc3 = p->pos3_z();
        fsfloc = eta(i,j) + p->phimean;
    

        if(zloc3<=fsfloc)
        {
        if(zloc3<=p->phimean)
        z=-(fabs(p->phimean-zloc3));
		
		if(zloc3>p->phimean)
        z=(fabs(p->phimean-zloc3));
        }
        
        if(zloc3>fsfloc)
        z = eta(i,j);


		// Wave Generation		
		if(p->B98==2 && w_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            if(zloc3<=fsfloc+epsi)
            wval[count] = wave_w(p,pgc,xg,yg,z);
            
            if(zloc3>fsfloc+epsi)
            wval[count] = 0.0;
            
            ++count;
            }
		}
    }	

    count=0;
    FLUIDLOOP
    {
        xg = xgen(p);
        yg = ygen(p);
		dg = distgen(p);
		db = distbeach(p);
        
		if(p->pos_z()<=p->phimean)
        z=-(fabs(p->phimean-p->pos_z()));
		
		if(p->pos_z()>p->phimean)
        z=(fabs(p->phimean-p->pos_z()));
		
		// Wave Generation
        if(p->B98==2 && h_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            lsval[count] = eta(i,j)+p->phimean-p->pos_z();
            
            ++count;
            }
		}
    }
    
    count=0;
    if(p->A10==3)
    FLOOP
    {
        xg = xgen(p);
        yg = ygen(p);
        dg = distgen(p);
		db = distbeach(p);
        
        zloc4 = p->pos_z();
        fsfloc = eta(i,j) + p->phimean;
        
        z=p->ZSN[FIJK]-p->phimean;
        
		// Wave Generation
		if(p->B98==2 && f_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            { 
            Fival[count] = wave_fi(p,pgc,xg,yg,z);
            ++count;
            }
		}
    }
    
    count=0;
    if(p->A10==3)
    FLUIDLOOP
    {
		
        xg = xgen(p);
        yg = ygen(p);
        dg = distgen(p);
		db = distbeach(p);
        
        zloc4 = p->pos_z();
        fsfloc = eta(i,j) + p->phimean;
    
        if(zloc4<=fsfloc)
        {
        if(zloc4<=p->phimean)
        z=-(fabs(p->phimean-zloc4));
		
		if(zloc4>p->phimean)
        z=(fabs(p->phimean-zloc4));
  
        if(zloc4>fsfloc)
        z = eta(i,j);
        }
		
		// Wave Generation		
		if(p->B98==2 && f_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            { 
            if(zloc4<=fsfloc+epsi)
            Fival[count] = wave_fi(p,pgc,xg,yg,z);
            
            if(zloc4>fsfloc+epsi)
            Fival[count] = 0.0;
            
            ++count;
            }
		}
    }
    

    count=0;
    SLICELOOP4
    {
		
        xg = xgen(p);
        yg = ygen(p);
        dg = distgen(p);
		db = distbeach(p);
        
        z = eta(i,j);
		
		// Wave Generation
		if(p->B98==2 && f_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            { 
            if(zloc4<=fsfloc+epsi || p->A10==3)
            Fifsfval[count] = wave_fi(p,pgc,xg,yg,z);
            
            if(zloc4>fsfloc+epsi && p->A10==4)
            Fifsfval[count] = 0.0;
            
            ++count;
            }
		}
    }

}
    
