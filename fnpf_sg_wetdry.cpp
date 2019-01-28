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
    if(p->A343>=1)
    {
      SLICELOOP4
      {
          if(p->count<1)
          if(eta(i,j) + p->wd - c->bed(i,j) >= wd_criterion)
          {
            if(c->wet(i,j)==0)
            {
            cout<<"WETTING !!!!!!!!!!!!!!!!!!!!!"<<endl;
            
            /*Fifsf(i,j) = Fifsf(i-1,j);
            FKLOOP
            c->Fi[IJK] = c->Fi[Im1JK];*/
            }
          c->wet(i,j)=1;
          
          }
              
          if(eta(i,j) + p->wd - c->bed(i,j) < wd_criterion)
          c->wet(i,j)=0;
      }
      
      pgc->gcsl_start4int(p,c->wet,50);
    }
}