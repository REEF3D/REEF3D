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

#include"ioflow_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"patchBC_interface.h"

void ioflow_f::fsfinflow_nhflow(lexer *p, fdm_nhf* d, ghostcell* pgc, slice &WL)
{
     // fsf
    p->phiout=0.0;
    count=0;
    for(n=0;n<p->gcslout_count;n++)
    {
        i=p->gcslout[n][0];
        j=p->gcslout[n][1];
        
        p->phiout+=d->eta(i,j);
        
        ++count;
    }
    
    p->phiout=pgc->globalsum(p->phiout);
    count=pgc->globalisum(count);
    
    p->phiout=p->phiout/(int(count)>0?int(count):1.0e20);
    
    p->phiout += p->phimean;
    
    //cout<<"phiout: "<<p->phiout<<" phimean: "<<p->phimean<<endl;
    
    
    // set fsf 
    double wsfout=p->phimean;
    
    if(p->F62>1.0e-20)
    wsfout=p->F62;
    
    for(n=0;n<p->gcslout_count;n++)
    {
        i=p->gcslout[n][0];
        j=p->gcslout[n][1];
        
        if(p->F50==2 || p->F50==4)
        {
        WL(i,j) = wsfout-d->bed(i,j);
        WL(i+1,j) = wsfout-d->bed(i,j);
        WL(i+2,j) = wsfout-d->bed(i,j);
        WL(i+3,j) = wsfout-d->bed(i,j);
        
        d->eta(i,j) = WL(i,j) - d->depth(i,j);
        d->eta(i+1,j) = WL(i+1,j) - d->depth(i,j);
        d->eta(i+2,j) = WL(i+2,j) - d->depth(i,j);
        d->eta(i+3,j) = WL(i+3,j) - d->depth(i,j);
        }
        
    }
    
    
    
}