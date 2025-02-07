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

#include"ioflow_f.h"
#include"lexer.h"
#include"fdm_nhf.h"
#include"ghostcell.h"
#include"patchBC_interface.h"

void ioflow_f::fsfinflow_nhflow(lexer *p, fdm_nhf* d, ghostcell* pgc, slice &WL)
{
    // -------------------------------------
     // fsf out
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
    
    
    // -------------------------------------
    // Find Hi
    double eta_in;
    count=0;
    zval=0.0;
    for(n=0;n<p->gcslin_count;n++)
    {
        i=p->gcslin[n][0];
        j=p->gcslin[n][1];

        zval+=d->eta(i,j);
        ++count;

    }

    count=pgc->globalisum(count);
    zval=pgc->globalsum(zval);
    
    eta_in=zval/double(count);
    
    p->Hi = p->wd+eta_in;
    
    
    // -------------------------------------
    // set fsf 
    double wsfout=p->phimean;
    double f=1.0;
    
    if(p->F62>1.0e-20)
    {
        if(p->F64==0)
        wsfout=p->F62;
        
        if(p->F64>0)
        {
        if(p->count<p->F64)
        f = 0.5*cos(PI + PI*double(p->count)/double(p->F64)) + 0.5;
        
        if(p->count>=p->F64)
        f = 1.0;
        
        wsfout = f*p->F62 + (1.0-f)*p->F60;
        //cout<<"wsfout: "<<wsfout<<" f: "<<f<<endl;
        }
    }
    
    //cout<<"wsfout: "<<wsfout<<" phiout: "<< p->phiout<<endl;
    
    for(n=0;n<p->gcslout_count;n++)
    {
        i=p->gcslout[n][0];
        j=p->gcslout[n][1];
        
        if(p->F50==2 || p->F50==4)
        if(p->wet[IJ]==1)
        {
        
    //cout<<"wsfout: "<<wsfout<<" phiout: "<< wsfout-d->bed(i,j)<<endl;
        WL(i,j)   = wsfout - d->bed(i,j);
        WL(i+1,j) = wsfout - d->bed(i,j);
        WL(i+2,j) = wsfout - d->bed(i,j);
        WL(i+3,j) = wsfout - d->bed(i,j);
        
        d->eta(i,j)   = WL(i,j)   - d->depth(i,j);
        d->eta(i+1,j) = WL(i+1,j) - d->depth(i,j);
        d->eta(i+2,j) = WL(i+2,j) - d->depth(i,j);
        d->eta(i+3,j) = WL(i+3,j) - d->depth(i,j);
        }
        
    }
}