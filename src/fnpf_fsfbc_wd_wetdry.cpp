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

#include"fnpf_fsfbc_wd.h"
#include"lexer.h"
#include"fdm_fnpf.h"
#include"ghostcell.h"
#include"fnpf_coastline.h"

void fnpf_fsfbc_wd::wetdry(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &eta, slice &Fifsf) 
{   
      SLICELOOP4
      c->wet_n(i,j)=p->wet[IJ];
      
      SLICELOOP4
      {     
          p->wet[IJ]=1;
          
          if(p->A343>=1)
          if(eta(i,j) + p->wd - c->bed(i,j) < c->wd_criterion)
          p->wet[IJ]=0;

      } 
      
      pgc->gcsl_start4Vint(p,p->wet,50);
      
      pcoast->start(p,pgc,c->coastline,p->wet,c->wet_n);
}

void fnpf_fsfbc_wd::coastline_eta(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &f) 
{
    SLICELOOP4
    {
        if(c->coastline(i,j)>=0.0)
        {
            db = c->coastline(i,j);
            
            if(db<dist3)
            {
            f(i,j) = rb3(p,db)*f(i,j);
        
            }
        }
        
        if(c->coastline(i,j)<0.0 && p->A343==1)
        f(i,j)=0.0;
    }
}

void fnpf_fsfbc_wd::coastline_fi(lexer *p, fdm_fnpf *c, ghostcell *pgc, slice &f) 
{
    SLICELOOP4
    {
        if(c->coastline(i,j)>=0.0)
        {
            db = c->coastline(i,j);
            
            if(db<dist4)
            {
            f(i,j) = rb4(p,db)*f(i,j);
        
            }
        }
        
        if(c->coastline(i,j)<0.0 && p->A343==1)
        f(i,j)=0.0;
    }
}

double fnpf_fsfbc_wd::rb3(lexer *p, double x)
{
    double r=0.0;

    x=(dist3-fabs(x))/(dist3);
    x=MAX(x,0.0);
    
    r = 1.0 - (exp(pow(x,3.5))-1.0)/(EE-1.0);

	return r;
}

double fnpf_fsfbc_wd::rb4(lexer *p, double x)
{
    double r=0.0;

    x=(dist4-fabs(x))/(dist4);
    x=MAX(x,0.0);
    
    r = 1.0 - (exp(pow(x,3.5))-1.0)/(EE-1.0);

	return r;
}
