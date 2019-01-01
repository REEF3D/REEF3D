/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2019 Hans Bihs

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
--------------------------------------------------------------------*/

#include"6DOF_f.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_f::motion_ext(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->X210==1 || p->X211==1)
	motion_fixed(p,a,pgc);
    
    if(p->X221==1)
	motion_vec(p,a,pgc);
    
    if(p->X210==1 || p->X211==1 || p->X221==1)
    {
	if(p->X11_u==2)
	dxg = p->dt*Uext;
	
	if(p->X11_v==2)
	dyg = p->dt*Vext;

	if(p->X11_w==2)
	dzg = p->dt*Wext;
	
	if(p->X11_p==2)
	dphi = p->dt*Pext;
	
	if(p->X11_q==2)
	dtheta = p->dt*Qext;
	
	if(p->X11_r==2)
	dpsi = p->dt*Rext;
    }
	
}


void sixdof_f::motion_ext_quaternion(lexer *p, fdm *a, ghostcell *pgc)
{
    if(p->X210==1 || p->X211==1)
	motion_fixed(p,a,pgc);
    
    if(p->X221==1)
	motion_vec(p,a,pgc);
    
    if(p->X210==1 || p->X211==1 || p->X221==1)
    {
		if(p->X11_u==2)
		dxg = p->dt*Uext;
		
		if(p->X11_v==2)
		dyg = p->dt*Vext;

		if(p->X11_w==2)
		dzg = p->dt*Wext;
		
		if(p->X11_p==2)
		cout<<"not implemented yet"<<endl;
		
		if(p->X11_q==2)
		cout<<"not implemented yet"<<endl;
		
		if(p->X11_r==2)
		cout<<"not implemented yet"<<endl;
    }
	
}


void sixdof_f::preventMotion(lexer *p)
{
	if(p->X11_u != 1) Xe = 0.0;
	if(p->X11_v != 1) Ye = 0.0;
	if(p->X11_w != 1) Ze = 0.0;
	if(p->X11_p != 1) Ke = 0.0;
	if(p->X11_q != 1) Me = 0.0;
	if(p->X11_r != 1) Ne = 0.0;
}