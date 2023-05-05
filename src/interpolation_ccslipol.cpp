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

#include"interpolation.h"
#include"slice.h"
#include"lexer.h"

double interpolation::ccslipol1(slice& f, double xp, double yp)
{
    ii=i;
    jj=j;
    
    i = p->posf_i(xp);
    j = p->posc_j(yp);
		
    // wa
    wa = (p->XN[IP1]-xp)/p->DXP[IP];
    
    if((p->XN[IP1]-xp)/p->DXP[IP]<0.0)
    {
    wa = (p->XN[IP2]-xp)/p->DXP[IP1];
    ++i;
    }
    
    if((p->XN[IP1]-xp)/p->DXP[IP]>1.0)
    {
    wa = (p->XN[IP]-xp)/p->DXP[IM1];
    --i;
    }
    
    
    // wb
    wb = (p->YP[JP1]-yp)/p->DYN[JP];
    
    if((p->YP[JP1]-yp)/p->DYN[JP]<0.0)
    {
    wb = (p->YP[JP2]-yp)/p->DYN[JP1];
    ++j;
    }
    
    if((p->YP[JP1]-yp)/p->DYN[JP]>1.0)
    {
    wb = (p->YP[JP]-yp)/p->DYN[JM1];
    --j;
    }
    
    value = lintsl1(f,i,j,wa,wb);

    i=ii;
    j=jj;

    return value;
}

double interpolation::ccslipol2(slice& f, double xp, double yp)
{
    ii=i;
    jj=j;
    
    i = p->posc_i(xp);
    j = p->posf_j(yp);
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXN[IP];
    
    if((p->XP[IP1]-xp)/p->DXN[IP]<0.0)
    {
    wa = (p->XP[IP2]-xp)/p->DXN[IP1];
    ++i;
    }
    
    if((p->XP[IP1]-xp)/p->DXN[IP]>1.0)
    {
    wa = (p->XP[IP]-xp)/p->DXN[IM1];
    --i;
    }
    
    
    // wb
    wb = (p->YN[JP1]-yp)/p->DYP[JP];
    
    if((p->YN[JP1]-yp)/p->DYP[JP]<0.0)
    {
    wb = (p->YN[JP2]-yp)/p->DYP[JP1];
    ++j;
    }
    
    if((p->YN[JP1]-yp)/p->DYP[JP]>1.0)
    {
    wb = (p->YN[JP]-yp)/p->DYP[JM1];
    --j;
    }

    value = lintsl2(f,i,j,wa,wb);

    i=ii;
    j=jj;

    return value;
}

double interpolation::ccslipol4(slice& f, double xp, double yp)
{
    ii=i;
    jj=j;
    
    i = p->posc_i(xp);
    j = p->posc_j(yp);
    
    
		
    // wa
    wa = (p->XP[IP1]-xp)/p->DXN[IP];
    
    if((p->XP[IP1]-xp)/p->DXN[IP]<0.0)
    {
    wa = (p->XP[IP2]-xp)/p->DXN[IP1];
    ++i;
    }
    
    if((p->XP[IP1]-xp)/p->DXN[IP]>1.0)
    {
    wa = (p->XP[IP]-xp)/p->DXN[IM1];
    --i;
    }
    

    // wb
    wb = (p->YP[JP1]-yp)/p->DYN[JP];
    
    if((p->YP[JP1]-yp)/p->DYN[JP]<0.0)
    {
    wb = (p->YP[JP2]-yp)/p->DYN[JP1];
    ++j;
    }
    
    if((p->YP[JP1]-yp)/p->DYN[JP]>1.0)
    {
    wb = (p->YP[JP]-yp)/p->DYN[JM1];
    --j;
    }
    
    
    value =  lintsl4(f,i,j,wa,wb);

    i=ii;
    j=jj;
    
    //cout<<i<<" "<<j<<" "<<wa<<" "<<wb<<" | "<<value<<endl;
    

    return value;
}


