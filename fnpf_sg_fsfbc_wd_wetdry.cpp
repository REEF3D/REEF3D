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

#include"fnpf_sg_fsfbc_wd.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"

void fnpf_sg_fsfbc_wd::wetdry(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf) 
{   
      SLICELOOP4
      {     
          c->wet(i,j)=1;
          
          if(p->A343>=1)
          {  
              if(eta(i,j) + p->wd - c->bed(i,j) < c->wd_criterion)
              {
              c->wet(i,j)=0;
              //eta(i,j)=0.5*c->wd_criterion+c->bed(i,j)-p->wd;
              }
          
              if(p->A343==2 || p->A343==4)
              if(c->wet(i,j)==0)
              if((eta(i-1,j)>eta(i,j)  && c->wet(i-1,j)==1)
              || (eta(i+1,j)>eta(i,j)  && c->wet(i+1,j)==1)
              || (eta(i,j-1)>eta(i,j)  && c->wet(i,j-1)==1)
              || (eta(i,j+1)>eta(i,j)  && c->wet(i,j+1)==1))
              {
              c->wet(i,j)=1;  
              eta(i,j)=c->wd_criterion;
              }      
          }    
      }
      
      pgc->gcsl_start4int(p,c->wet,50);
      
      if(p->A343==3 || p->A343==4)
      SLICELOOP4
      {     
          if(c->wet(i,j)==1)
          if((eta(i,j)<eta(i-1,j)  && c->wet(i-1,j)==0) 
          || (eta(i,j)<eta(i+1,j)  && c->wet(i+1,j)==0) 
          || (eta(i,j)<eta(i,j-1)  && c->wet(i,j-1)==0) 
          || (eta(i,j)<eta(i,j+1)  && c->wet(i,j+1)==0) )
          {
          c->wet(i,j)=0;  
          eta(i,j)=c->wd_criterion*0.5;
          }   
      }

      
      pgc->gcsl_start4int(p,c->wet,50);

}
