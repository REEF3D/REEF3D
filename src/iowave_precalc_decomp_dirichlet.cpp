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

void iowave::wavegen_precalc_decomp_dirichlet(lexer *p, ghostcell *pgc)
{
    
double fsfloc;
int qn;

        count=0;
		for(n=0;n<p->gcin_count;n++)
		{
		i=p->gcin[n][0];
		j=p->gcin[n][1];
        
        x=xgen(p);
        y=ygen(p);
        x1=xgen1(p);
        y2=ygen2(p);
        

    // ETA
        eta(i,j) = 0.0;
        etaval[count] = 0.0;
            
            for(qn=0;qn<wave_comp;++qn)
            {
            eta(i,j) += etaval_S_cos[count][qn]*etaval_T_cos[qn] - etaval_S_sin[count][qn]*etaval_T_sin[qn];
            etaval[count] = eta(i,j);
            }
        }
    
        count=0;
		for(n=0;n<p->gcin_count;n++)
		{
		i=p->gcin[n][0];
		j=p->gcin[n][1];
		k=p->gcin[n][2];
        
        x=xgen(p);
        y=ygen(p);
        x1=xgen1(p);
        y2=ygen2(p);
        

            
        zloc3 = p->pos3_z();
        zloc4 = p->pos_z();
        
        fsfloc = eta(i,j) + p->phimean;
    
        // Z3
        if(zloc3<=fsfloc)
        {
        if(zloc3<=p->phimean)
        z3 = -(fabs(p->phimean-zloc3));
		
		if(zloc3>p->phimean)
        z3 = (fabs(p->phimean-zloc3));
        }
        
        if(zloc3>fsfloc)
        z3 = eta(i,j);
        
        // Z
        if(zloc4<=fsfloc)
        {
        if(zloc4<=p->phimean)
        z=-(fabs(p->phimean-zloc4));
		
		if(zloc4>p->phimean)
        z=(fabs(p->phimean-zloc4));
        }
        
        if(zloc4>fsfloc)
        z = eta(i,j);
        
        uval[count] = 0.0;
        vval[count] = 0.0;
        wval[count] = 0.0;
        
    // U
        if(zloc4<=fsfloc+epsi)
        for(qn=0;qn<wave_comp;++qn)
        uval[count] += uval_S_cos[count][qn]*uval_T_cos[qn] - uval_S_sin[count][qn]*uval_T_sin[qn];
            
        if(zloc4>fsfloc+epsi)
        uval[count] = 0.0;

    // V
        if(zloc4<=fsfloc+epsi)
        for(qn=0;qn<wave_comp;++qn)
            vval[count] += vval_S_cos[count][qn]*vval_T_cos[qn] - vval_S_sin[count][qn]*vval_T_sin[qn];
            
        if(zloc4>fsfloc+epsi)
        vval[count] = 0.0;
        
    // W
        if(zloc4<=fsfloc+epsi)
        for(qn=0;qn<wave_comp;++qn)
        wval[count] += wval_S_cos[count][qn]*wval_T_sin[qn] + wval_S_sin[count][qn]*wval_T_cos[qn];
            
        if(zloc4>fsfloc+epsi)
        wval[count] = 0.0;
        
    // LS
        lsval[count] = eta(i,j) + p->phimean - p->pos_z();

        ++count;
        }

}
