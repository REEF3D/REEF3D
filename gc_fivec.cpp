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

#include"ghostcell.h"
#include"lexer.h"
#include"sliceint.h"

void ghostcell::fivec(lexer *p, double *f, sliceint &bc)
{	
    FLOOP
    {  
        if(p->B98!=3||bc(i-1,j)==0)
        if(p->flag7[FIm1JK]<0)
        {
        f[FIm1JK] = f[FIJK];
        f[FIm2JK] = f[FIJK];
        f[FIm3JK] = f[FIJK];
        }
          
        if(p->B99!=3||bc(i-1,j)==0)
        if(p->flag7[FIp1JK]<0)
        {
        f[FIp1JK] = f[FIJK];
        f[FIp2JK] = f[FIJK];
        f[FIp3JK] = f[FIJK];
        }
        
        if(p->flag7[FIJm1K]<0)
        {
        f[FIJm1K] = f[FIJK];
        f[FIJm2K] = f[FIJK];
        f[FIJm3K] = f[FIJK];
        }
        
        if(p->flag7[FIJp1K]<0)
        {
        f[FIJp1K] = f[FIJK];
        f[FIJp2K] = f[FIJK];
        f[FIJp3K] = f[FIJK];
        }
    }
}

void ghostcell::fivec2D(lexer *p, double *f)
{	
    FLOOP
    {
        if(p->flag7[FIm1JK]<0)
        {
        f[FIm1JK] = f[FIJK];
        f[FIm2JK] = f[FIJK];
        f[FIm3JK] = f[FIJK];
        }
        
        if(p->flag7[FIp1JK]<0)
        {
        f[FIp1JK] = f[FIJK];
        f[FIp2JK] = f[FIJK];
        f[FIp3JK] = f[FIJK];
        }
    }
}


