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

#include"ptf_fsfbc.h"
#include"lexer.h"
#include"fdm_ptf.h"
#include"ghostcell.h"
#include"ptf_coastline.h"

void ptf_fsfbc::wetdry(lexer *p, fdm_ptf *a, ghostcell *pgc, slice &eta, slice &Fifsf) 
{   
    /*
    SLICELOOP4
    a->WL(i,j) = eta(i,j) + p->wd - a->bed(i,j);

    
    pgc->gcsl_start4(p,a->WL,50);
      
    
    SLICELOOP4
    {
    p->wet_n[IJ] = p->wet[IJ];
    temp[IJ] = p->wet[IJ];
    }
     
    SLICELOOP4
    {
        if(p->wet[IJ]==0)
        {
            if(p->wet[Ip1J]==1 && eta(i,j)<eta(i+1,j) && a->WL(i+1,j)>a->wd_criterion+eps)
            temp[IJ]=1;
            
            if(p->wet[Im1J]==1 && eta(i,j)<eta(i-1,j) && a->WL(i-1,j)>a->wd_criterion+eps)
            temp[IJ]=1;
            
            if(p->wet[IJp1]==1 && eta(i,j)<eta(i,j+1) && a->WL(i,j+1)>a->wd_criterion+eps && p->j_dir==1)
            temp[IJ]=1;
            
            if(p->wet[IJm1]==1 && eta(i,j)<eta(i,j-1) && a->WL(i,j-1)>a->wd_criterion+eps && p->j_dir==1)
            temp[IJ]=1;
        }
        
        else              
        if(a->WL(i,j)<=a->wd_criterion)
        {
        temp[IJ]=0;
        eta(i,j) = a->wd_criterion - a->depth(i,j);
        a->WL(i,j) = eta(i,j) + a->depth(i,j);
        }
    }
    
    SLICELOOP4
    p->wet[IJ] = temp[IJ];
    

    pgc->gcsl_start4Vint(p,p->wet,50);
    pgc->gcsl_start4(p,eta,gcval_eta);
    pgc->gcsl_start4(p,a->WL,gcval_eta);*/
    
    // wetdry old
    
      SLICELOOP4
      a->wet_n(i,j)=p->wet[IJ];
      
      if(p->count<2)
      SLICELOOP4
      {     
          p->wet[IJ]=1;
          
          if(p->A343>=1)
          if(p->wd - a->bed(i,j) < a->wd_criterion)
          p->wet[IJ]=0;
          
        //if(p->wet[IJ]==0)
        //Fifsf(i,j) = 0.0;

      } 
      
      pgc->gcsl_start4Vint(p,p->wet,50);
      
      pcoast->start(p,a,pgc,a->coastline,p->wet,a->wet_n);
      
      
    
    // check
    /*SLICELOOP4
    {
    eta(i,j) = MAX(eta(i,j), -p->wd + a->bed(i,j) + a->wd_criterion);

    a->WL(i,j) = MAX(a->wd_criterion, a->eta(i,j) + p->wd - a->bed(i,j));
    
    //p->wet[IJ]=1;
    
    //if(p->wet[IJ]==0)
    //Fifsf(i,j) = 0.0;
    }*/
    
    
    SLICELOOP4
    a->test2D(i,j) = double (p->wet[IJ]);
    
}

