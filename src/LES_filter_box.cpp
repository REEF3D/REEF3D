/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2024 Hans Bihs

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

#include"LES_filter_box.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"strain.h"

LES_filter_box::LES_filter_box(lexer* p, fdm* a) : strain(p,a)
{

}

LES_filter_box::~LES_filter_box()
{
}

void LES_filter_box::start(lexer *p, fdm *a, ghostcell *pgc, field &uprime, field &vprime, field &wprime, int gcval)
{
	//    vel_label=veleval(p,gcv);


	if(gcval==10)
	{
        
        ULOOP
        uprime(i,j,k) = a->u(i,j,k);
        
        pgc->start1(p,uprime,gcval);

        
	}

	if(gcval==11)
	{

        VLOOP
        vprime(i,j,k) = a->v(i,j,k);
        
        pgc->start2(p,vprime,gcval);

        
	}

	if(gcval==12)
	{

        WLOOP
        wprime(i,j,k) = a->w(i,j,k);
        
        pgc->start3(p,wprime,gcval);


	}
 
}




