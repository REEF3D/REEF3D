/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

void sflow_eta::wetdry_nb(lexer* p, fdm2D* b, ghostcell* pgc, slice &eta, slice &P, slice &Q, slice &ws)
{
    double eps=1.0e-6;
    
    SLICELOOP4
	b->hp(i,j) = eta(i,j) + b->depth(i,j);
    
    pgc->gcsl_start4(p,b->hp,gcval_eta);
    
        if(p->count==0)
        {
            SLICELOOP4
            if(eta(i,j) <= wd_criterion - b->depth(i,j) + eps)
            {
            temp[IJ]=0;
            eta(i,j) = wd_criterion - b->depth(i,j);
            b->hp(i,j) = wd_criterion;
            }
            
            SLICEBASELOOP
            if(p->flagslice4[IJ]<0)
            {
            p->wet[IJ]=0;
            temp[IJ]=0;
            }
            
        pgc->gcsl_start4(p,eta,gcval_eta);
        pgc->gcsl_start4(p,b->hp,gcval_eta);
        }
        
    SLICELOOP4
    {
    p->wet_n[IJ] = p->wet[IJ];
    temp[IJ] = p->wet[IJ];
    }
    
    pgc->gcsl_start4Vint(p,p->wet,50);
    
    // ----------------
    // ----------------

    SLICELOOP4
    {
        if(p->wet[IJ]==0)
        {
            if(p->wet[Ip1J]==1 && eta(i,j)<eta(i+1,j) && b->hp(i+1,j)>wd_criterion+eps)
            temp[IJ]=1;
            
            if(p->wet[Im1J]==1 && eta(i,j)<eta(i-1,j) && b->hp(i-1,j)>wd_criterion+eps)
            temp[IJ]=1;
            
            if(p->wet[IJp1]==1 && eta(i,j)<eta(i,j+1) && b->hp(i,j+1)>wd_criterion+eps && p->j_dir==1)
            temp[IJ]=1;
            
            if(p->wet[IJm1]==1 && eta(i,j)<eta(i,j-1) && b->hp(i,j-1)>wd_criterion+eps && p->j_dir==1)
            temp[IJ]=1;
        }
        
        else              
        if(eta(i,j) < wd_criterion - b->depth(i,j) + eps)
        {
        temp[IJ]=0;
        eta(i,j) = wd_criterion - b->depth(i,j);
        b->hp(i,j) = wd_criterion;
        }
    }
    
    SLICELOOP4
    p->wet[IJ] = temp[IJ];
    
    pgc->gcsl_start4Vint(p,p->wet,50);

    
    SLICELOOP1
      {
          if(p->wet[IJ]==1 && p->wet[Ip1J]==1)
           b->wet1(i,j)=1;
           
          if(p->wet[IJ]==0 || p->wet[Ip1J]==0)// || b->hx(i,j)<wd_criterion)
          {
           b->P(i,j)=0.0; 
           P(i,j)=0.0; 
           b->wet1(i,j)=0;
          }
      }
      
      SLICELOOP2
      {
          if(p->wet[IJ]==1 && p->wet[IJp1]==1)
           b->wet2(i,j)=1;
           
          if(p->wet[IJ]==0 || p->wet[IJp1]==0)// || b->hy(i,j)<wd_criterion)
          {
           b->Q(i,j)=0.0; 
           Q(i,j)=0.0; 
           b->wet2(i,j)=0;
          }
      }
    pgc->gcsl_start1int(p,b->wet1,50);
    pgc->gcsl_start2int(p,b->wet2,50);
    /*
    SLICELOOP1
    b->hx(i,j) = MAX(b->hx(i,j), wd_criterion);
    
    SLICELOOP2
    b->hy(i,j) = MAX(b->hy(i,j), wd_criterion);
    */

}
