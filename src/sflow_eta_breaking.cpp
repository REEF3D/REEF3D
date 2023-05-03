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
#include"sflow_eta.h"
#include"lexer.h"
#include"fdm2D.h"
#include"ghostcell.h"

void sflow_eta::breaking(lexer* p, fdm2D* b, ghostcell* pgc, slice &eta, slice &eta_n, double alpha)
{    
    if(p->A246>=1)
    SLICELOOP4
    {
            if( (eta(i,j)-eta_n(i,j))/(alpha*p->dt) > p->A247*sqrt(9.81*b->hp(i,j)))
            b->breaking(i,j)=1;
            
            if(p->A246==2)
            {
            if((eta(i+1,j)-eta(i-1,j))/(2.0*p->DXM)   < -p->A355 || (eta(i+1,j)-eta(i-1,j))/(2.0*p->DXM)   > p->A355)
            b->breaking(i,j)=1;
            
            if((eta(i,j+1)-eta(i,j+1))/(2.0*p->DXM)   < -p->A355 || (eta(i,j+1)-eta(i,j-1))/(2.0*p->DXM)   > p->A355)
            b->breaking(i,j)=1;
            }
    }
    
    if(p->A242==1)
    SLICELOOP4
    if(p->wet[IJ]==1)
    {
        if(p->wet[Ip1J]==0 || p->wet[Ip2J]==0)
        b->breaking(i,j)=1;
            
        if(p->wet[Im1J]==0 || p->wet[Im2J]==0)
        b->breaking(i,j)=1;
        
        
        if(p->wet[IJp1]==0 || p->wet[IJp2]==0)
        b->breaking(i,j)=1;
            
        if(p->wet[IJm1]==0 || p->wet[IJm2]==0)
        b->breaking(i,j)=1;
            
        //if(p->wet[IJp1]==0)
       // b->breaking(i,j)=1;

        //if(p->wet[IJm1]==0)
        //b->breaking(i,j)=1;
        
        if(b->hp(i,j) < 0.01)
            {
            b->breaking(i-1,j)=1;
            b->breaking(i,j)=1;
            b->breaking(i+1,j)=1;
            
            b->breaking(i,j-1)=1;
            b->breaking(i,j+1)=1;
            }
    }
    
    
    if(p->A242>=2)
    SLICELOOP4
    if(p->wet[IJ]==1)
    {
        if(p->wet[Ip1J]==0 || p->wet[Ip2J]==0)
        b->breaking(i,j)=1;
            
        if(p->wet[Im1J]==0 || p->wet[Im2J]==0)
        b->breaking(i,j)=1;
        
        
        
        if(fabs((b->depth(i+1,j) - b->depth(i-1,j))/p->DXM) > 2.0)// || fabs((b->hp(i+1,j) - b->hp(i-1,j))/p->DXM) > 2.0 )
        {
        b->breaking(i-1,j)=1;
        b->breaking(i,j)=1;
        b->breaking(i+1,j)=1;
        
        b->breaking(i-2,j)=1;
        b->breaking(i+2,j)=1;
        
        b->breaking(i-3,j)=1;
        b->breaking(i+3,j)=1;
        }
            
        if(fabs((b->depth(i,j+1) - b->depth(i,j-1))/p->DXM) > 2.0)// || fabs((b->hp(i,j+1) - b->hp(i,j-1))/p->DXM) > 2.0 )
        {
        b->breaking(i,j-1)=1;
        b->breaking(i,j)=1;
        b->breaking(i,j+1)=1;
        
        b->breaking(i,j-2)=1;
        b->breaking(i,j+2)=1;
        
        b->breaking(i,j-3)=1;
        b->breaking(i,j+3)=1;
        }
        
            if(p->A242==3)
            if(b->hp(i,j) < 0.02)
            {
            b->breaking(i-1,j)=1;
            b->breaking(i,j)=1;
            b->breaking(i+1,j)=1;
            
            b->breaking(i,j-1)=1;
            b->breaking(i,j+1)=1;
            }

            
        //if(p->wet(i,j+1)==0)
       // b->breaking(i,j)=1;

        //if(p->wet(i,j-1)==0)
        //b->breaking(i,j)=1;
    }
    
    if(p->B77==2)
    for(n=0;n<p->gcslout_count;++n)
    {
		i=p->gcslout[n][0];
		j=p->gcslout[n][1];
        
        b->breaking(i-1,j)=1;
        b->breaking(i,j)=1;
    }
    
    pgc->gcsl_start4int(p,b->breaking,1);
}


void sflow_eta::breaking_persist(lexer* p, fdm2D* b, ghostcell* pgc, slice &eta, slice &eta_n, double alpha)
{    
    SLICELOOP4
    b->breaking_print(i,j)=double(b->breaking(i,j));
    
    SLICELOOP4
    {
        if(p->A248==1)
        {
            if(b->breaking(i,j)>=1)
            {
            if( (eta(i,j)-eta_n(i,j))/(alpha*p->dt) > p->A249*sqrt(9.81*b->hp(i,j)))
            b->breaking(i,j)=1;
            
            if( (eta(i,j)-eta_n(i,j))/(alpha*p->dt) <= p->A249*sqrt(9.81*b->hp(i,j)))
            b->breaking(i,j)=0;
            }
        }
        
        
    }    
}
