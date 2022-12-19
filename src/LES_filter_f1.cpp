/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2022 Hans Bihs

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

#include"LES_filter_f1.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"
#include"strain.h"

LES_filter_f1::LES_filter_f1(lexer* p, fdm* a) : strain(p,a), ubar(p), vbar(p), wbar(p), etta_u1(p), etta_v1(p), etta_w1(p), etta_u2(p), etta_v2(p), etta_w2(p)
{

}

LES_filter_f1::~LES_filter_f1()
{
}

void LES_filter_f1::start(lexer *p, fdm *a, ghostcell *pgc, field &uprime, field &vprime, field &wprime, int gcval)
{
//    vel_label=veleval(p,gcv);


	if(gcval==10)
	{

        ULOOP
        etta_u1(i,j,k) = -0.25*(a->u(i+1,j,k) -2.0*a->u(i,j,k) + a->u(i-1,j,k));
                
        ULOOP
        etta_u1(i,j,k) = a->u(i,j,k) - etta_u1(i,j,k);
        
        pgc->start1(p,etta_u1,gcval);
        
        ULOOP
        etta_u2(i,j,k) = -0.25*(etta_u1(i,j+1,k) -2.0*etta_u1(i,j,k) + etta_u1(i,j-1,k));
                
        ULOOP
        etta_u2(i,j,k) = etta_u1(i,j,k) - etta_u2(i,j,k);
        
        pgc->start1(p,etta_u2,gcval);
        
        ULOOP
        ubar(i,j,k) = -0.25*(etta_u2(i,j,k+1) -2.0*etta_u2(i,j,k) + etta_u2(i,j,k-1));
                
        ULOOP
        ubar(i,j,k) = etta_u2(i,j,k) - ubar(i,j,k);
        
        pgc->start1(p,ubar,gcval);
        
        ULOOP
        uprime(i,j,k) = a->u(i,j,k) - ubar(i,j,k);
        
        pgc->start1(p,uprime,gcval);

        
	}

	if(gcval==11)
	{
        
        VLOOP
        etta_v1(i,j,k) = -0.25*(a->v(i+1,j,k) -2.0*a->v(i,j,k) + a->v(i-1,j,k));
                
        VLOOP
        etta_v1(i,j,k) = a->v(i,j,k) - etta_v1(i,j,k);
        
        pgc->start2(p,etta_v1,gcval);
        
        VLOOP
        etta_v2(i,j,k) = -0.25*(etta_v1(i,j+1,k) -2.0*etta_v1(i,j,k) + etta_v1(i,j-1,k));
                
        VLOOP
        etta_v2(i,j,k) = etta_v1(i,j,k) - etta_v2(i,j,k);
        
        pgc->start2(p,etta_v2,gcval);
        
        VLOOP
        vbar(i,j,k) = -0.25*(etta_v2(i,j,k+1) -2.0*etta_v2(i,j,k) + etta_v2(i,j,k-1));
                
        VLOOP
        vbar(i,j,k) = etta_v2(i,j,k) - vbar(i,j,k);
        
        pgc->start2(p,vbar,gcval);
        
        VLOOP
        vprime(i,j,k) = a->v(i,j,k) - vbar(i,j,k);
        
        pgc->start2(p,vprime,gcval);

        
	}

	if(gcval==12)
	{
        
        WLOOP
        etta_w1(i,j,k) = -0.25*(a->w(i+1,j,k) -2.0*a->w(i,j,k) + a->w(i-1,j,k));
                
        WLOOP
        etta_w1(i,j,k) = a->w(i,j,k) - etta_w1(i,j,k);
        
        pgc->start3(p,etta_w1,gcval);
        
        WLOOP
        etta_w2(i,j,k) = -0.25*(etta_w1(i,j+1,k) -2.0*etta_w1(i,j,k) + etta_w1(i,j-1,k));
                
        WLOOP
        etta_w2(i,j,k) = etta_w1(i,j,k) - etta_w2(i,j,k);
        
        pgc->start3(p,etta_w2,gcval);
        
        WLOOP
        wbar(i,j,k) = -0.25*(etta_w2(i,j,k+1) -2.0*etta_w2(i,j,k) + etta_w2(i,j,k-1));
                
        WLOOP
        wbar(i,j,k) = etta_w2(i,j,k) - wbar(i,j,k);
        
        pgc->start3(p,wbar,gcval);
        
        WLOOP
        wprime(i,j,k) = a->w(i,j,k) - wbar(i,j,k);
        
        pgc->start3(p,wprime,gcval);


	}
 
}

/*int LES_filter_f1::veleval(lexer *p, int gcv)
{
//	Velocities


	if(gcv==10)
	return 1;
	
	if(gcv==11)
	return 2;
	
	if(gcv==12)
	return 3;
    

	else
	return 0;
}

*/


