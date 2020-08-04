/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2020 Hans Bihs

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

#include"icds2.h"
#include"lexer.h"
#include"fdm.h"
#include"flux_face_CDS2.h"
#include"flux_face_CDS4.h"
#include"flux_face_CDS2_vrans.h"
#include"flux_face_FOU.h"
#include"flux_face_FOU_vrans.h"
#include"flux_face_QOU.h"

icds2::icds2 (lexer *p)
{
    if(p->B269==0)
    {
        if(p->D11==1)
        pflux = new flux_face_FOU(p);
        
        if(p->D11==2)
        pflux = new flux_face_CDS2(p);
        
        if(p->D11==3)
        pflux = new flux_face_QOU(p);
        
        if(p->D11==4)
        pflux = new flux_face_CDS4(p);
    }
    
    if(p->B269>=1 || p->S10==2)
    {
        if(p->D11==1)
        pflux = new flux_face_FOU_vrans(p);
        
        if(p->D11==2)
        pflux = new flux_face_CDS2_vrans(p);
        
        if(p->D11==3)
        pflux = new flux_face_FOU_vrans(p);
        
        if(p->D11==4)
        pflux = new flux_face_CDS2(p);
    }
}

icds2::~icds2()
{
}

void icds2::start(lexer* p, fdm* a, field& b, int ipol, field& uvel, field& vvel, field& wvel)
{
    count=0;

    if(ipol==1)
    ULOOP
    aij(p,a,b,1,uvel,vvel,wvel);

    if(ipol==2)
    VLOOP
    aij(p,a,b,2,uvel,vvel,wvel);

    if(ipol==3)
    WLOOP
    aij(p,a,b,3,uvel,vvel,wvel);

    if(ipol==4)
    LOOP
    aij(p,a,b,4,uvel,vvel,wvel);

    if(ipol==5)
    LOOP
    aij(p,a,b,5,uvel,vvel,wvel);
}

void icds2::aij(lexer* p,fdm* a,field& b,int ipol, field& uvel, field& vvel, field& wvel)
{		
	pflux->u_flux(a,ipol,uvel,ivel1,ivel2);
    pflux->v_flux(a,ipol,vvel,jvel1,jvel2);
    pflux->w_flux(a,ipol,wvel,kvel1,kvel2);

	 
	 a->M.p[count] = ( ivel2 - ivel1
					 + jvel2 - jvel1
					 + kvel2 - kvel1)/(p->DXM*2.0);
	 
	 a->M.s[count] = -ivel1/(p->DXM*2.0);
	 a->M.n[count] =  ivel2/(p->DXM*2.0);
	 
	 a->M.e[count] = -jvel1/(p->DXM*2.0);
	 a->M.w[count] =  jvel2/(p->DXM*2.0);
	 
	 a->M.b[count] = -kvel1/(p->DXM*2.0);
	 a->M.t[count] =  kvel2/(p->DXM*2.0);
	 
	 ++count;
}





