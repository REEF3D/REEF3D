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

#include"nhflow_fsf_rk.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"

void nhflow_fsf_rk::breaking(lexer* p, fdm_nhf* d, ghostcell* pgc, slice& eta, slice &eta_n, double alpha)
{   
    
    SLICELOOP4
    d->breaking(i,j)=0;
        
    if(p->A550>=1)
    SLICELOOP4
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
    
    
    
}