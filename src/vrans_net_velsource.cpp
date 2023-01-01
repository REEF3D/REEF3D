/*--------------------------------------------------------------------
REEF3D
Copyright 2018-2023 Tobias Martin

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTIBILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Hans Bihs
--------------------------------------------------------------------*/

#include"vrans_net.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void vrans_net::u_source(lexer *p, fdm *a)
{
    count=0;
    ULOOP
	{
        a->rhsvec.V[count++] -= Fx_net(i,j,k);
	}
	
}

void vrans_net::v_source(lexer *p, fdm *a)
{  
    count=0;
    VLOOP
	{
       a->rhsvec.V[count++] -= Fy_net(i,j,k);        
	}   
}

void vrans_net::w_source(lexer *p, fdm *a)
{   
    count=0;
    WLOOP
	{
        a->rhsvec.V[count++] -= Fz_net(i,j,k);        
	}		
}
