/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2023 Hans Bihs

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

#include<iostream>
#include<fstream>
#include"6DOF_gc.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

void sixdof_gc::print_E_position(lexer *p, fdm *a, ghostcell *pgc)
{
	if(p->count%p->X19==0)
	eposout<<p->simtime<<" \t "<<xg<<" \t "<<yg<<" \t "<<zg<<" \t "<<phi*(180.0/PI)<<" \t "<<theta*(180.0/PI)<<" \t "<<psi*(180.0/PI)<<endl;
}

void sixdof_gc::print_E_velocity(lexer *p, fdm *a, ghostcell *pgc)
{
	if(p->count%p->X19==0)
	evelout<<p->simtime<<" \t "<<Ue<<" \t "<<Ve<<" \t "<<We<<" \t "<<Pe*(180.0/PI)<<" \t "<<Qe*(180.0/PI)<<" \t "<<Re*(180.0/PI)<<endl;
}

void sixdof_gc::print_E_force(lexer *p, fdm *a, ghostcell *pgc)
{
	
	if(p->count%p->X19==0)
	eforceout<<p->simtime<<" \t "<<Xe<<" \t "<<Ye<<" \t "<<Ze<<" \t "<<Ke<<" \t "<<Me<<" \t "<<Ne<<endl;
}

void sixdof_gc::print_S_force(lexer *p, fdm *a, ghostcell *pgc)
{
	
	if(p->count%p->X19==0)
	sforceout<<p->simtime<<" \t "<<Xs<<" \t "<<Ys<<" \t "<<Zs<<" \t "<<Ks<<" \t "<<Ms<<" \t "<<Ne<<endl;
}
