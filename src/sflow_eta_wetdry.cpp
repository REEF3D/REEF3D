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

void sflow_eta::wetdry(lexer* p, fdm2D* b, ghostcell* pgc, slice &P, slice &Q, slice &ws)
{
    if(p->A243>=1)
    {
      SLICELOOP4
      {
          
          if(b->hp(i,j)>=wd_criterion)
          p->wet[IJ]=1;
              
          if(b->hp(i,j)<wd_criterion)
          {
           ws(i,j)=0.0; 
           b->ws(i,j)=0.0;
           p->wet[IJ]=0;
          }
      }
      
      pgc->gcsl_start4Vint(p,p->wet,50);
      
      
      SLICELOOP1
      {
          if(b->hx(i,j)>=wd_criterion)
           b->wet1(i,j)=1;
           
          if(b->hx(i,j)<wd_criterion || (p->wet[IJ]==0 && p->wet[Ip1J]==0))
          {
           b->P(i,j)=0.0; 
           P(i,j)=0.0; 
           b->wet1(i,j)=0;
          }
      }
      
      SLICELOOP2
      {
          if(b->hy(i,j)>=wd_criterion)
           b->wet2(i,j)=1;
           
          if(b->hy(i,j)<wd_criterion || (p->wet[IJ]==0 && p->wet[IJp1]==0))
          {
           b->Q(i,j)=0.0; 
           Q(i,j)=0.0; 
           b->wet2(i,j)=0;
          }
      }
    pgc->gcsl_start1int(p,b->wet1,50);
    pgc->gcsl_start2int(p,b->wet2,50);
      
    // gcslin update
    if(p->count<=1)
    {
        for(n=0;n<p->gcslin_count;n++)
        {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];
        
        if(p->wet[IJ]==0)
        p->gcslin[n][5]=0;
        }
    }
    
    }
    
}
