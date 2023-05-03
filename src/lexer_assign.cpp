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

#include"lexer.h"

void lexer::control_calc()
{
	dt=dt_old=0.0;
	simtime=I50;
	sedtime=0.0;
	dtsed=0.0;
	presstime=veltime=lsmtime=reinitime=reinitime=turbtime=0.0;
    fsitime=fbtime=0.0;
    fbdt=fbmax=0.0;
	printouttime=0.0;
	xtime=0.0;
	gctime=0.0;
	totaltime=0.0;
	meantime=0.0;
	Xmeantime=Xtotaltime=0.0;
	gcmeantime=gctotaltime=0.0;
}

void lexer::assign_margin()
{	
    margin=3; 
    
    if(A311==7)
	margin=4;
    
	imax=knox+2*margin;
	jmax=knoy+2*margin;
	kmax=knoz+2*margin;
    kmaxF=knoz+1+2*margin;
	
	imin=-margin;
	jmin=-margin;
	kmin=-margin;
}

int lexer::maxparacount()
{
        maxpara=0;

        maxpara=MAX(gcpara1_count,gcpara2_count);

        maxpara=MAX(maxpara,gcpara3_count);
        maxpara=MAX(maxpara,gcpara4_count);
        maxpara=MAX(maxpara,gcpara5_count);
        maxpara=MAX(maxpara,gcpara6_count);

        return maxpara;
}



