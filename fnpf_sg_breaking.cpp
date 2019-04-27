/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"fnpf_sg_fsfbc.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"

void fnpf_sg_fsfbc::breaking(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &eta_n, slice &Fifsf, double alpha)
{
    int ii,jj;
    
    if(p->A346>=1)
    SLICELOOP4
    c->breaking(i,j)=0;
    
    if((p->A346==1 || p->A346==3) && p->count>1)
    SLICELOOP4
    {
            
            if( (eta(i,j)-eta_n(i,j))/(alpha*p->dt) > p->A347*sqrt(9.81*c->WL(i,j)))
            {
            c->breaking(i,j)=1;
            }
    }
    
    if((p->A346==2 || p->A346==3) && p->count>1)
    SLICELOOP4
    {
            
            if( (eta(i+1,j)-eta(i-1,j))/(p->DXP[IM1] + p->DXP[IP])   < -p->A348)
            {
                c->breaking(i,j)=1;
                c->breaking(i-1,j)=1;
                c->breaking(i-2,j)=1;
                c->breaking(i-3,j)=1;
            }
            
            if( (eta(i+1,j)-eta(i-1,j))/(p->DXP[IM1] + p->DXP[IP])   > p->A348)
            {
                c->breaking(i,j)=1;
                c->breaking(i+1,j)=1;
                c->breaking(i+2,j)=1;
                c->breaking(i+3,j)=1;
            }
            
            if( (eta(i,j+1)-eta(i,j-1))/(p->DYP[JM1] + p->DYP[JP])   < -p->A348)
            {
                c->breaking(i,j)=1;
                c->breaking(i,j-1)=1;
                c->breaking(i,j-2)=1;
                c->breaking(i,j-3)=1;
            }
            
            if( (eta(i,j+1)-eta(i,j-1))/(p->DYP[JM1] + p->DYP[JP])    > p->A348)
            {
                c->breaking(i,j)=1;
                c->breaking(i,j+1)=1;
                c->breaking(i,j+2)=1;
                c->breaking(i,j+3)=1;
            }
    }

    
    if(p->A346>=1)
    SLICELOOP4
    {
        if(c->breaking(i,j)==1 || c->breaking(i-1,j)==1 || c->breaking(i+1,j)==1 || c->breaking(i,j-1)==1 || c->breaking(i,j+1)==1)
        {
         filter(p,c,pgc,eta);
         filter(p,c,pgc,Fifsf);
        }   
    }
}

void fnpf_sg_fsfbc::filter(lexer *p, fdm_fnpf *c,ghostcell *pgc, slice &f)
{
    double he,hw,hn,hs,hp;
    double dhe, dhw, dhn, dhs,dhp;
    
    int outer_iter = p->A350;
    int inner_iter = p->A351;
    
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