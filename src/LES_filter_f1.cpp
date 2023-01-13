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
        etta_u1(i,j,k) = -0.5*(p->DXN[IP]*(a->u(i+1,j,k) - a->u(i,j,k)) - p->DXN[IP1]*(a->u(i,j,k) - a->u(i-1,j,k)))/(p->DXN[IP]+p->DXN[IP1]);
                
        ULOOP
        etta_u1(i,j,k) = a->u(i,j,k) - etta_u1(i,j,k);
        
        pgc->start1(p,etta_u1,gcval);
        
        ULOOP
        etta_u2(i,j,k) = -0.5*(p->DYP[JM1]*(etta_u1(i,j+1,k) - etta_u1(i,j,k)) - p->DYP[JP]*(etta_u1(i,j,k) - etta_u1(i,j-1,k)))/(p->DYP[JP]+p->DYP[JM1]);
                
        ULOOP
        etta_u2(i,j,k) = etta_u1(i,j,k) - etta_u2(i,j,k);
        
        pgc->start1(p,etta_u2,gcval);
        
        ULOOP
        ubar(i,j,k) = -0.5*(p->DZP[KM1]*(etta_u2(i,j,k+1) - etta_u2(i,j,k)) - p->DZP[KP]*(etta_u2(i,j,k) - etta_u2(i,j,k-1)))/(p->DZP[KP]+p->DZP[KM1]);
                
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
        etta_v1(i,j,k) = -0.5*(p->DXP[IM1]*(a->v(i+1,j,k) - a->v(i,j,k)) - p->DXP[IP]*(a->v(i,j,k) - a->v(i-1,j,k)))/(p->DXP[IP]+p->DXP[IM1]);
                
        VLOOP
        etta_v1(i,j,k) = a->v(i,j,k) - etta_v1(i,j,k);
        
        pgc->start2(p,etta_v1,gcval);
        
        VLOOP
        etta_v2(i,j,k) = -0.5*(p->DYN[JP]*(etta_v1(i,j+1,k) - etta_v1(i,j,k)) - p->DYN[JP1]*(etta_v1(i,j,k) - etta_v1(i,j-1,k)))/(p->DYN[JP]+p->DYN[JP1]);
                
        VLOOP
        etta_v2(i,j,k) = etta_v1(i,j,k) - etta_v2(i,j,k);
        
        pgc->start2(p,etta_v2,gcval);
        
        VLOOP
        vbar(i,j,k) = -0.5*(p->DZP[KM1]*(etta_v2(i,j,k+1) - etta_v2(i,j,k)) - p->DZP[KP]*(etta_v2(i,j,k) - etta_v2(i,j,k-1)))/(p->DZP[KP]+p->DZP[KM1]);
                
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
        etta_w1(i,j,k) = -0.5*(p->DXP[IM1]*(a->w(i+1,j,k) - a->w(i,j,k)) - p->DXP[IP]*(a->w(i,j,k) - a->w(i-1,j,k)))/(p->DXP[IP]+p->DXP[IM1]);
                
        WLOOP
        etta_w1(i,j,k) = a->w(i,j,k) - etta_w1(i,j,k);
        
        pgc->start3(p,etta_w1,gcval);
        
        WLOOP
        etta_w2(i,j,k) = -0.5*(p->DYP[JM1]*(etta_w1(i,j+1,k) - etta_w1(i,j,k)) - p->DYP[JP]*(etta_w1(i,j,k) - etta_w1(i,j-1,k)))/(p->DYP[JP]+p->DYP[JM1]);
                
        WLOOP
        etta_w2(i,j,k) = etta_w1(i,j,k) - etta_w2(i,j,k);
        
        pgc->start3(p,etta_w2,gcval);
        
        WLOOP
        wbar(i,j,k) = -0.5*(p->DZN[KP]*(etta_w2(i,j,k+1) - etta_w2(i,j,k)) - p->DZN[KP1]*(etta_w2(i,j,k) - etta_w2(i,j,k-1)))/(p->DZN[KP]+p->DZN[KP1]);
                
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


