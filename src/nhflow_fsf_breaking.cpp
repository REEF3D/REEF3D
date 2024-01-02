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

#include"nhflow_fsf_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void nhflow_fsf_f::breaking(lexer* p, fdm_nhf* d, ghostcell* pgc, slice& eta, slice &eta_n, double alpha)
{   
    
    SLICELOOP4
    d->breaking(i,j)=0;
        
    if(p->A550>=1)
    SLICELOOP4
    if(p->wet[IJ]==1)
    {       
            if(p->A551==1 || p->A551==3)
            if( (eta(i,j)-eta_n(i,j))/(alpha*p->dt) > p->A247*sqrt(9.81*d->WL(i,j)))
            d->breaking(i,j)=1;
            
            if(p->A551==1 || p->A551==3)
            {
            if((eta(i+1,j)-eta(i-1,j))/(2.0*p->DXM)   < -p->A355 || (eta(i+1,j)-eta(i-1,j))/(2.0*p->DXM)   > p->A355)
            d->breaking(i,j)=1;
            
            if((eta(i,j+1)-eta(i,j+1))/(2.0*p->DXM)   < -p->A355 || (eta(i,j+1)-eta(i,j-1))/(2.0*p->DXM)   > p->A355)
            d->breaking(i,j)=1;
            }
    }
    
    
    SLICELOOP4
    if(p->wet[Ip1J]==0 || p->wet[Im1J]==0 || p->wet[IJp1]==0 || p->wet[IJm1]==0)
    d->breaking(i,j)=1;
    
    if(p->A553==1)
    SLICELOOP4
    if(p->wet[IJ]==1)
    {
        if(p->wet[Ip1J]==0 || p->wet[Ip2J]==0)
        d->breaking(i,j)=1;
            
        if(p->wet[Im1J]==0 || p->wet[Im2J]==0)
        d->breaking(i,j)=1;
        
        
        if(p->wet[IJp1]==0 || p->wet[IJp2]==0)
        d->breaking(i,j)=1;
            
        if(p->wet[IJm1]==0 || p->wet[IJm2]==0)
        d->breaking(i,j)=1;
            
        //if(p->wet[IJp1]==0)
       // d->breaking(i,j)=1;

        //if(p->wet[IJm1]==0)
        //d->breaking(i,j)=1;
        
        if(eta(i,j) + p->wd - d->bed(i,j) < 0.01)
        {
            d->breaking(i-1,j)=1;
            d->breaking(i,j)=1;
            d->breaking(i+1,j)=1;
            
            d->breaking(i,j-1)=1;
            d->breaking(i,j+1)=1;
        }
    }
    
    pgc->gcsl_start4int(p,d->breaking,50);
    
    // filter
        if(p->A552==1)
        SLICELOOP4
        if(d->breaking(i,j)==2)
        {
         filter(p,d,pgc,d->WL);
        }   
        
        if(p->A552==2)
        SLICELOOP4
        if(d->breaking(i,j)==1)
        {
         filter(p,d,pgc,d->WL);
        }   
        
        if(p->A552==3)
        SLICELOOP4
        if(d->breaking(i,j)>=1)
        {
         filter(p,d,pgc,d->WL);
        } 
}

void nhflow_fsf_f::filter(lexer *p, fdm_nhf *d, ghostcell *pgc, slice &f)
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