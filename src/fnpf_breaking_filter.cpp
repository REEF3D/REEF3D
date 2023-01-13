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

#include"fnpf_breaking.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"

void fnpf_breaking::filter(lexer *p, fdm_fnpf *c,ghostcell *pgc, slice &f)
{
    double he,hw,hn,hs,hp;
    double dhe, dhw, dhn, dhs,dhp;
    
    int outer_iter = p->A361;
    int inner_iter = p->A362;
    
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
