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

#include"ghostcell.h"
#include"lexer.h"
#include"slice.h"

int ghostcell::gcsleval1(lexer *p, int gcv, int bc, int cs)
{
    // general Neuman
    if(gcv==40 || gcv==50 || gcv==1)
	return 4;
    
//Wall
	// Parallel
	if((bc==21||bc==22||bc==7||bc==5)&&(cs==2||cs==3||cs==5||cs==6)&&(gcv==10||gcv==1||gcv==20))
	return gclabel_u;
	
    // Orthogonal
	else
	if((bc==21||bc==22||bc==5||(bc==7&&awa_label==0))&&(cs==1||cs==4)&&(gcv==10||gcv==20||gcv==1))
	return 5;

//Inflow: none

	
//Outflow
	else
	if((bc==2)&&(cs==1||cs==4) && (gcv==10||gcv==20||gcv==1))
	return 4;

//Symmetry
	else
	if(bc==3 && (cs==2||cs==3||cs==5||cs==6) && (gcv==10||gcv==20||gcv==1))
	return 4;

	else
	if(bc==3 && (cs==1||cs==4)&&(gcv==10||gcv==20||gcv==1))
	return 5;
    /*
    else
	if(bc==8 && (gcv==10||gcv==20||gcv==1) && p->B99==4)
	return 8;*/
	
//Hx
    else
    if((bc==1||bc==6)&&(gcv==51||gcv==52||gcv==54))
	return 4;
    
    else
    if((bc==2||bc==7)&&(gcv==51||gcv==54))
	return 4;
    
    else
    if(bc==2)
	return 4;
    
    else
    if(bc==8 && p->B99==3)
	return 4;
    
    /*
    else
    if(bc==8 && p->B99==4)
	return 4;*/
    
    else
    if((bc==21||bc==3)&&(gcv==51||gcv==52||gcv==53||gcv==54))
	return 4;
    
    
    else
    return -1;
}


void ghostcell::gcsldistro1(lexer *p, slice &f, int ii, int jj, int nn, double dist,  int gcv, int bc, int cs)
{
    i=ii;
	j=jj;
	n=nn;

	bc_label=gcsleval1(p,gcv,bc,cs);

	if(bc_label==4)
	gcsl_neumann(f,gcv,bc,cs);
	
	if(bc_label==5)
	gcsl_noslip(f,gcv,bc,cs);
    
    if(bc_label==7)
	gcsl_outflow(p,f,gcv,bc,cs);
    
    if(bc_label==8)
	gcsl_sommerfeld(p,f,gcv,bc,cs);
    
    if(bc_label==9)
	gcsl_outflow_fsf(p,f,gcv,bc,cs);
    
}
