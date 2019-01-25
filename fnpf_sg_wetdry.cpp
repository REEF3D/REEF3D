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

void fnpf_sg_fsfbc::wetdry(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf) 
{  
    SLICELOOP4
    c->wet(i,j)=1;
    
    
    
    /*
    if(p->A243>=1)
    {
      SLICELOOP4
      {
          if(b->hp(i,j)>=wd_criterion)
          b->wet4(i,j)=1;
              
          if(b->hp(i,j)<wd_criterion)
          {
           ws(i,j)=0.0; 
           b->ws(i,j)=0.0;
           //b->wb(i,j)=0.0;
           b->wet4(i,j)=0;
          }
      }
      
      pgc->gcsl_start4int(p,b->wet4,50);
      
      
      SLICELOOP1
      //if(b->wet4(i,j)==0)
      if(b->hx(i,j)<wd_criterion)
      {
       b->P(i,j)=0.0; 
       P(i,j)=0.0; 
      }
      
      SLICELOOP2
      //if(b->wet4(i,j)==0)
      if(b->hy(i,j)<wd_criterion)
      {
       b->Q(i,j)=0.0; 
       Q(i,j)=0.0; 
      }
    }*/
    
}