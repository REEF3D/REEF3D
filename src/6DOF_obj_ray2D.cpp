/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"6DOF_obj.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"fieldint.h"

void sixdof_obj::ray_cast_2D(lexer *p, ghostcell *pgc)
{
	SLICELOOP4
	{
        fsio(i,j)=1;
	    fs(i,j)=1.0e9;
        Rxmin(i,j)=1.0e9;
        Rxmax(i,j)=-1.0e9;
        Rymin(i,j)=1.0e9;
        Rymax(i,j)=-1.0e9;
        Ls(i,j) = 1.0e-6;
        Bs(i,j) = 1.0e-6;
        draft(i,j) = 0.0;
	}
	
    for(int rayiter=0; rayiter<2; ++rayiter)
    {
        if(rayiter==0)
        {
            ray_cast_2D_io_x(p,pgc,0,tricount);
            ray_cast_2D_io_ycorr(p,pgc,0,tricount);
        }
    
        if(rayiter==1)
        {
            pgc->gcslparax_int(p,fsio,1);
            
            ray_cast_2D_x(p,pgc,0,tricount);
            ray_cast_2D_y(p,pgc,0,tricount);
            ray_cast_2D_z(p,pgc,0,tricount);
        }
    }
    
	SLICELOOP4
    {
        if(fsio(i,j)==-1)
        fs(i,j) = -fabs(fs(i,j));
        
        if(fsio(i,j)==1)
        fs(i,j) = fabs(fs(i,j));
    }
    
	SLICELOOP4
	{
		if(fs(i,j) > 10.0*p->DXM)
		fs(i,j) = 10.0*p->DXM;
		
		if(fs(i,j) < -10.0*p->DXM)
		fs(i,j) = -10.0*p->DXM;
	}
    
    SLICELOOP4
	{
    if(Rxmax(i,j)>-1.0e9 && Rxmin(i,j)<1.0e9)
    Ls(i,j) = Rxmax(i,j)-Rxmin(i,j);
    
    if(Rymax(i,j)>-1.0e8 && Rymin(i,j)<1.0e8)
    Bs(i,j) = Rymax(i,j)-Rymin(i,j);
    }
    
	pgc->gcsl_start4(p,fs,50);
    pgc->gcsl_start4(p,draft,50);
}





