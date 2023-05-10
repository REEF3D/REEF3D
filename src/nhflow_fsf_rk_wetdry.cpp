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

void nhflow_fsf_rk::wetdry(lexer* p, fdm_nhf* d, ghostcell* pgc, double *U, double *V, double *W, slice &eta)
{
    double wl;
    
    SLICELOOP4
    {
          
    wl =  eta(i,j) + p->wd - d->bed(i,j);
    
          if(wl>=wd_criterion)
          p->wet[IJ]=1;
              
          if(wl<wd_criterion)
          {
              p->wet[IJ]=0;
              
              KLOOP
              PCHECK
              {
               U[IJK] = 0.0;
               V[IJK] = 0.0;
               W[IJK] = 0.0;
              }
          }
    }
    
    SLICELOOP4
    if(eta(i,j)< -p->wd  + d->bed(i,j) - wd_criterion+1.0e-20)
    eta(i,j) = -p->wd  + d->bed(i,j) - wd_criterion - 1.0e-15;
    
    pgc->gcsl_start4Vint(p,p->wet,50);
    
}