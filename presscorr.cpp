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

#include"presscorr.h"
#include"lexer.h"
#include"fdm.h"

presscorr::presscorr (lexer * p) : density(p)
{
}

presscorr::~presscorr()
{
}

void presscorr::istart(lexer* p, fdm* a, field &apu, field &apv, field &apw, field &pcorr)
{

	n=0;
	LOOP
	{	 
	
         a->M.s[n]    = -1.0/(roface(p,a,-1,0,0)*(apu(i-1,j,k)>1.0e-8?apu(i-1,j,k):(1.0/p->dt))*p->DXP[IM1]*p->DXN[IP])*p->x_dir;

         a->M.e[n]    = -1.0/(roface(p,a,0,-1,0)*(apv(i,j-1,k)>1.0e-8?apv(i,j-1,k):(1.0/p->dt))*p->DYP[JM1]*p->DYN[JP])*p->y_dir;

         a->M.b[n]    = -1.0/(roface(p,a,0,0,-1)*(apw(i,j,k-1)>1.0e-8?apw(i,j,k-1):(1.0/p->dt))*p->DZP[KM1]*p->DZN[KP])*p->z_dir;


         a->M.p[n]    = 1.0/(roface(p,a,-1,0,0)*(apu(i-1,j,k)>1.0e-8?apu(i-1,j,k):(1.0/p->dt))*p->DXP[IM1]*p->DXN[IP])*p->x_dir
					   + 1.0/(roface(p,a,1,0,0)*(apu(i,j,k)>1.0e-8?apu(i,j,k):(1.0/p->dt))*p->DXP[IP]*p->DXN[IP])*p->x_dir

                       + 1.0/(roface(p,a,0,-1,0)*(apv(i,j-1,k)>1.0e-8?apv(i,j-1,k):(1.0/p->dt))*p->DYP[JM1]*p->DYN[JP])*p->y_dir
                       + 1.0/(roface(p,a,0,1,0)*(apv(i,j,k)>1.0e-8?apv(i,j,k):(1.0/p->dt))*p->DYP[JP]*p->DYN[JP])*p->y_dir

				       + 1.0/(roface(p,a,0,0,-1)*(apw(i,j,k-1)>1.0e-8?apw(i,j,k-1):(1.0/p->dt))*p->DZP[KM1]*p->DZN[KP])*p->z_dir
				       + 1.0/(roface(p,a,0,0,1)*(apw(i,j,k)>1.0e-8?apw(i,j,k):(1.0/p->dt))*p->DZP[KP]*p->DZN[KP])*p->z_dir;


         a->M.t[n]    = -1.0/(roface(p,a,0,0,1)*(apw(i,j,k)>1.0e-8?apw(i,j,k):(1.0/p->dt))*p->DZP[KP]*p->DZN[KP])*p->z_dir;

         a->M.w[n]    = -1.0/(roface(p,a,0,1,0)*(apv(i,j,k)>1.0e-8?apv(i,j,k):(1.0/p->dt))*p->DYP[JP]*p->DYN[JP])*p->y_dir;

         a->M.n[n]    = -1.0/(roface(p,a,1,0,0)*(apu(i,j,k)>1.0e-8?apu(i,j,k):(1.0/p->dt))*p->DXP[IP]*p->DXN[IP])*p->x_dir;
     ++n;
	}
}

void presscorr::estart(lexer* p,fdm* a, field &press)
{
}
