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

#include"roughness.h"
#include"lexer.h"
#include"fdm.h"
#include"ghostcell.h"

roughness::roughness(lexer *p) 
{
}

roughness::~roughness()
{
}

double roughness::ks_val(lexer *p, fdm* a,int ii,int jj, int kk, int cs, int bc)
{
	ks=p->B50;
	
	if(cs==1 && p->B51>0.0)
	ks=p->B51;	
	
	if(cs==2 && p->B52>0.0)
	ks=p->B52;	
	
	if(cs==3 && p->B53>0.0)
	ks=p->B53;	
	
	if(cs==4 && p->B54>0.0)
	ks=p->B54;	
	
	if(cs==5 && p->B55>0.0)
	ks=p->B55;	
	
	if(cs==6 && p->B56>0.0)
	ks=p->B56;	
	
	if(bc==5)
	ks=p->S21*p->S20;
    
    if(p->G3==1 && p->S10>0 && (a->topo(i-1,j,k)<0.0 || a->topo(i+1,j,k-1)<0.0 || a->topo(i,j-1,k)<0.0 || a->topo(i,j+1,k)<0.0 || a->topo(i,j,k-1)<0.0))
    ks=p->S21*p->S20;
			
		
	return ks;
}
