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
#include"ghostcell.h"

void iowave::wavegen_precalc_space_dirichlet(lexer *p, ghostcell *pgc)
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
            if(dg<1.0e20)
            {
                for(qn=0;qn<wave_comp;++qn)
                {
                etaval_S_sin[count][qn] = wave_eta_space_sin(p,pgc,xg,yg,qn);
                etaval_S_cos[count][qn] = wave_eta_space_cos(p,pgc,xg,yg,qn);
                }
            ++count;
            }
		}
    }
    pgc->gcsl_start4(p,eta,50);
    


    count=0;
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];
    
    xg = xgen1(p);
    yg = ygen1(p);
        
        KLOOP
        UCHECK
        {
        zloc1 = p->pos1_z();

        if(zloc1<=p->phimean)
        z=-(fabs(p->phimean-zloc1));
		
		if(zloc1>p->phimean)
        z=(fabs(p->phimean-zloc1));
        
        if(zloc1>p->phimean+p->wA)
        z = 0.5*(eta(i,j)+eta(i+1,j));

		// Wave Generation
		if(p->B98>=3 && u_switch==1)
        {
                for(qn=0;qn<wave_comp;++qn)
                {
                uval_S_sin[count][qn] = wave_u_space_sin(p,pgc,xg,yg,z,qn);
                uval_S_cos[count][qn] = wave_u_space_cos(p,pgc,xg,yg,z,qn);
                }
            ++count;
		}
        }
    }


    count=0;
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];
    
    xg = xgen2(p);
    yg = ygen2(p);
        
        KLOOP
        VCHECK
        {
        zloc2 = p->pos2_z();

        if(zloc2<=p->phimean)
        z=-(fabs(p->phimean-zloc2));
		
		if(zloc2>p->phimean)
        z=(fabs(p->phimean-zloc2));
        
        if(zloc2>p->phimean+p->wA)
        z = 0.5*(eta(i,j)+eta(i,j+1));
        
		// Wave Generation		
		if(p->B98>=3 && v_switch==1)
        {
                for(qn=0;qn<wave_comp;++qn)
                {
                vval_S_sin[count][qn] = wave_v_space_sin(p,pgc,xg,yg,z,qn);
                vval_S_cos[count][qn] = wave_v_space_cos(p,pgc,xg,yg,z,qn);
                }
            ++count;
		}
        }
    }


    
    count=0;
    count=0;
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];
    
    xg = xgen(p);
    yg = ygen(p);
        
        KWLOOP
        WCHECK
        {
        zloc3 = p->pos3_z();

        if(zloc3<=p->phimean)
        z=-(fabs(p->phimean-zloc3));
		
		if(zloc3>p->phimean)
        z=(fabs(p->phimean-zloc3));
        
        if(zloc3>p->phimean+p->wA)
        z = eta(i,j);
        
		// Wave Generation
		if(p->B98>=3 && w_switch==1)
        {
                for(qn=0;qn<wave_comp;++qn)
                {
                wval_S_sin[count][qn] = wave_w_space_sin(p,pgc,xg,yg,z,qn);
                wval_S_cos[count][qn] = wave_w_space_cos(p,pgc,xg,yg,z,qn);
                }
            ++count;
		}
        }
    }	


    
    count=0;
    if(f_switch==1)
    for(n=0;n<p->gcslin_count;n++)
    {
    i=p->gcslin[n][0];
    j=p->gcslin[n][1];
    
    xg = xgen1(p);
    yg = ygen1(p);
    dg = distgen(p);
    db = distbeach(p);
        
        KLOOP
        PCHECK
        {
        zloc4 = p->pos_z();

        if(zloc4<=p->phimean)
        z=-(fabs(p->phimean-zloc4));
		
		if(zloc4>p->phimean)
        z=(fabs(p->phimean-zloc4));
		
		// Wave Generation
		if(p->B98==3)
        {
            // Zone 1
            if(dg<dist1)
            {
                for(qn=0;qn<wave_comp;++qn)
                {
                Fival_S_sin[count][qn] = wave_fi_space_sin(p,pgc,xg,yg,z,qn);
                Fival_S_cos[count][qn] = wave_fi_space_cos(p,pgc,xg,yg,z,qn);
                }
            ++count;
            }
		}
        }
    }

}
