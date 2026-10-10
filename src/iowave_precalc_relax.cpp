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

void iowave::wavegen_precalc_relax(lexer *p, ghostcell *pgc)
{
    double fsfloc;
    
    p->wavetime = p->simtime;
    
    // generation-zone columns, registered for the cached evaluation: per column the cell
    // centre (q = 3g), the u point (3g+1) and the v point (3g+2); before, the theory was
    // evaluated directly in every zone cell
    if(!gen_built) cfd_genzone_build(p,pgc);
    const bool zsel = zones.has_sources();
    
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
            if(zsel)
            select_sources(gen_src[gen_idx[IJ]]);   // zone sources (B 524)
            eta(i,j) = wave_eta_c(p,pgc,3*gen_idx[IJ]);
            }
		}
    }
    pgc->gcsl_start4(p,eta,50);
    
    count=0;
    if(!zc_built) zonecol_build(p);
    ZULOOP
    {
        xg = xgen1(p);
        yg = ygen1(p);
        dg = distgen(p);
        db = distbeach(p);

        zloc1 = p->pos1_z();
        fsfloc = 0.5*(eta(i,j)+eta(i+1,j)) + p->phimean;
    
        if(zloc1<=fsfloc)
        z = zloc1-p->phimean;
        
        if(zloc1>fsfloc)
        z = 0.5*(eta(i,j)+eta(i+1,j));
		
		// Wave Generation
		if(p->B98==2 && u_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            if(zsel)
            select_sources(gen_src[gen_idx[IJ]]);   // zone sources (B 524)
            if(zloc1<=fsfloc+epsi)
            {
            double uw,vw,ww;
            wave_uvw_c(p,pgc,3*gen_idx[IJ]+1,z,uw,vw,ww);
            uval[count] = uw + p->Ui;
            }
            
            if(zloc1>fsfloc+epsi)
            uval[count] = 0.0;
            
            ++count;
            }
		}
    }
		
    count=0;
    if(!zc_built) zonecol_build(p);
    ZVLOOP
    {
        xg = xgen2(p);
        yg = ygen2(p);
        dg = distgen(p);
        db = distbeach(p);

        zloc2 = p->pos2_z();
        fsfloc = 0.5*(eta(i,j)+eta(i,j+1)) + p->phimean;
    

        if(zloc2<=fsfloc)
        z = zloc2-p->phimean;
        
        if(zloc2>fsfloc)
        z = 0.5*(eta(i,j)+eta(i,j+1));

		// Wave Generation
		if(p->B98==2 && v_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            if(zsel)
            select_sources(gen_src[gen_idx[IJ]]);   // zone sources (B 524)
            if(zloc2<=fsfloc+epsi)
            {
            double uw,vw,ww;
            wave_uvw_c(p,pgc,3*gen_idx[IJ]+2,z,uw,vw,ww);
            vval[count] = vw;
            }
            
            if(zloc2>fsfloc+epsi)
            vval[count] = 0.0;
            
            ++count;
            }
		}
    }

    count=0;
    if(!zc_built) zonecol_build(p);
    ZWLOOP
    {
        xg = xgen(p);
        yg = ygen(p);
        dg = distgen(p);
        db = distbeach(p);

        zloc3 = p->pos3_z();
        fsfloc = eta(i,j) + p->phimean;
    
        if(zloc3<=fsfloc)
        z = zloc3-p->phimean;
        
        if(zloc3>fsfloc)
        z = eta(i,j);


		// Wave Generation		
		if(p->B98==2 && w_switch==1)
        {
            // Zone 1
            if(dg<1.0e20)
            {
            if(zsel)
            select_sources(gen_src[gen_idx[IJ]]);   // zone sources (B 524)
            if(zloc3<=fsfloc+epsi)
            {
            double uw,vw,ww;
            wave_uvw_c(p,pgc,3*gen_idx[IJ],z,uw,vw,ww);
            wval[count] = ww;
            }
            
            if(zloc3>fsfloc+epsi)
            wval[count] = 0.0;
            
            ++count;
            }
		}
    }	

    count=0;
    if(!zc_built) zonecol_build(p);
    ZLOOP
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
            select_zone_at(p);   // zone sources (B 524)
            lsval[count] = eta(i,j)+p->phimean-p->pos_z();
            
            ++count;
            }
		}
    }
    
    count=0;
    
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
            select_zone_at(p);   // zone sources (B 524)
            Fival[count] = wave_fi(p,pgc,xg,yg,z);
            ++count;
            }
		}
    }
    
    count=0;
    if(p->A10==3)
    LOOP
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
            select_zone_at(p);   // zone sources (B 524)
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
            select_zone_at(p);   // zone sources (B 524)
            if(zloc4<=fsfloc+epsi || p->A10==3)
            Fifsfval[count] = wave_fi(p,pgc,xg,yg,z);
            
            
            ++count;
            }
		}
    }
    
    if(p->F80==4)
    {
        field4 &vg = *vofgen;

    LOOP
        {
        if((eta(i,j)+p->phimean)>=(p->pos_z()+0.5*p->DZN[KP]))
            vg(i,j,k)=1.0;
        else if((eta(i,j)+p->phimean)<=(p->pos_z()-0.5*p->DZN[KP]))
            vg(i,j,k)=0.0;
        else
            vg(i,j,k)=(eta(i,j)+p->phimean-(p->pos_z()-0.5*p->DZN[KP]))/p->DZN[KP];
                    
        if(vg(i,j,k)>1.0)
            vg(i,j,k)=1.0;
        else if(vg(i,j,k)<0.0)
            vg(i,j,k)=0.0;
        }
        
    }
    
    if(zones.has_sources())
    select_sources(nullptr);
}
    
