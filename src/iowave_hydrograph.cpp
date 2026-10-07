/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2026 Hans Bihs

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

#include"iowave.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include<fstream>


double iowave::hydrograph_ipol(lexer *p, ghostcell* pgc, double ** hydro, int hydrocount)
{
	double val;
    
    for(n=0;n<hydrocount-1;++n)
    if(p->simtime>=hydro[n][0] && p->simtime<hydro[n+1][0])
	{
    val = ((hydro[n+1][1]-hydro[n][1])/(hydro[n+1][0]-hydro[n][0]))*(p->simtime-hydro[n][0]) + hydro[n][1];
	}
    
    if(p->count==0 )
    val = hydro[0][1];
	
	if(p->simtime>=hydro[hydrocount-1][0])
	val=hydro[hydrocount-1][1];
	
	return val;
	
}
