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

#include"fnpf_breaking.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"

void fnpf_breaking::filter(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &f, int outer_iter, int inner_iter)
{
    double he,hw,hn,hs,hp;
    double dhe, dhw, dhn, dhs,dhp;
    
    // A343 2/3: never filter a dry cell, and treat dry neighbours as mirror
    // cells (zero gradient) instead of averaging with their held values
    if(p->A343>=2)
    {
        if(p->wet[IJ]==0)
        return;
        
        const bool ws = (p->wet[Im1J]==1);
        const bool wn = (p->wet[Ip1J]==1);
        const bool we = (p->wet[IJm1]==1 || p->j_dir==0);
        const bool ww = (p->wet[IJp1]==1 || p->j_dir==0);
        
        if(!(ws && wn && we && ww))
        {
            const double cw = (p->j_dir==1) ? 0.125 : 0.25;
            
            for(int qn=0;qn<outer_iter;++qn)
            {
                hp = f(i,j);
                hs = ws ? f(i-1,j) : hp;
                hn = wn ? f(i+1,j) : hp;
                he = (p->j_dir==1 && we) ? f(i,j-1) : hp;
                hw = (p->j_dir==1 && ww) ? f(i,j+1) : hp;
                
                // predictor
                if(p->j_dir==1)
                f(i,j) = 0.5*hp + cw*(hs + hn + he + hw);
                else
                f(i,j) = 0.5*hp + cw*(hs + hn);
                
                // corrector
                for(int qqn=0;qqn<inner_iter;++qqn)
                {
                    const double d0 = hp - f(i,j);
                    
                    dhs = ws ? hs - f(i-1,j) : d0;
                    dhn = wn ? hn - f(i+1,j) : d0;
                    dhe = (p->j_dir==1 && we) ? he - f(i,j-1) : d0;
                    dhw = (p->j_dir==1 && ww) ? hw - f(i,j+1) : d0;
                    
                    if(p->j_dir==1)
                    dhp = 0.5*d0 + cw*(dhs + dhn + dhe + dhw);
                    else
                    dhp = 0.5*d0 + cw*(dhs + dhn);
                    
                    f(i,j) += dhp;
                }
            }
            return;
        }
    }
    
    
    if(p->j_dir==0)
	for(int qn=0;qn<outer_iter;++qn)
	{
		hp = f(i,j);
        hs = f(i-1,j);
        hn = f(i+1,j);

        // predictor

		f(i,j) = 0.5*hp + 0.25*(hs + hn);
		
        // corrector
		for(int qqn=0;qqn<inner_iter;++qqn)
		{
            dhp = hp - f(i,j);
            dhs = hs - f(i-1,j);
            dhn = hn - f(i+1,j);
            
            dhp = 0.5*dhp+ 0.25*(dhs + dhn);
            f(i,j) += dhp;
		}
    }
    
    
    if(p->j_dir==1)
	for(int qn=0;qn<outer_iter;++qn)
	{
		hp = f(i,j);
        hs = f(i-1,j);
        hn = f(i+1,j);
        he = f(i,j-1);
        hw = f(i,j+1);
		
        // predictor

		f(i,j) = 0.5*hp + 0.125*(hs + hn + he + hw);
		
        // corrector
		for(int qqn=0;qqn<inner_iter;++qqn)
		{
            dhp = hp - f(i,j);
            dhs = hs - f(i-1,j);
            dhn = hn - f(i+1,j);
            dhe = he - f(i,j-1);
            dhw = hw - f(i,j+1);
            
            dhp = 0.5*dhp+ 0.125*(dhs + dhn + dhe + dhw);
            f(i,j) += dhp;
		}
    }
}
