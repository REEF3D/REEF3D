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

#include"ghostcell.h"
#include"lexer.h"
#include"slice.h"

int ghostcell::gcsleval2(lexer *p, int gcv, int bc, int cs)
{
    // general Neuman
    if(gcv==40 || gcv==50 || gcv==1)
	return 4;
    
//Wall
    // Parallel	
	if((bc==21||bc==22||bc==7||bc==6||bc==5||bc==1)&&(cs==1||cs==4||cs==5||cs==6)&&(gcv==11||gcv==21||gcv==2))
	return gclabel_v;
	
    // Orthogonal
	else
	if((bc==21||bc==22||bc==5||bc==7)&&(cs==2||cs==3)&&(gcv==11||gcv==21||gcv==2))
	return 5;
	
//Inflow 
    else
	if((bc==1) && (gcv==11||gcv==21||gcv==2) && (cs==1||cs==4||cs==5||cs==6))
	return 4;
    
//Patch    
    else
	if((bc==111 || bc==112 || bc==121 || bc==122) && (gcv==11||gcv==2||gcv==21||gcv==8))
	return 4;
    
//Outflow
	else
	if((bc==2)&&(cs==2||cs==3) && (gcv==11||gcv==21||gcv==2))
	return 4;
	
// Symmetry
	else
	if(bc==3 && (cs==1||cs==4||cs==5||cs==6) && (gcv==11||gcv==21|| gcv==2))
	return 4;

	else
	if(bc==3 && (cs==2||cs==3)&&(gcv==11||gcv==21|| gcv==2))
	return 5;
	
//Hy
    else
    if((bc==1||bc==6)&&(gcv==52||gcv==54))
	return 4;
    
    else
    if((bc==2||bc==7)&&(gcv==51||gcv==54))
	return 4;
    
    else
    if(bc==8 && p->B99==3)
	return 4;
    
    else
    if(bc==8 && p->B99==4)
	return 4;
    
    else
    if((bc==21||bc==3)&&(gcv==51||gcv==52||gcv==53||gcv==54))
	return 4;
    
    //Patch Hy
    else
	if((bc==221 || bc==211 || bc==121 || bc==111) && (gcv==55||gcv==51||gcv==52||gcv==53||gcv==54))
	return 42;
    
    else
	if((bc==222 || bc==212 || bc==122 || bc==112) && (gcv==55||gcv==51||gcv==52||gcv==53||gcv==54))
	return 4;
    
    
    else
    return -1;
}


void ghostcell::gcsldistro2(lexer *p, slice &f, int ii, int jj, int nn, double dist,  int gcv, int bc, int cs)
{
    i=ii;
	j=jj;
	n=nn;

	bc_label=gcsleval2(p,gcv,bc,cs);

	if(bc_label==4)
	gcsl_neumann(f,gcv,bc,cs);
    
    if(bc_label==42)
	gcsl_neumann_hy(f,gcv,bc,cs);
	
	if(bc_label==5)
	gcsl_noslip(f,gcv,bc,cs);
    
    if(bc_label==7)
	gcsl_outflow(p,f,gcv,bc,cs);
    
    if(bc_label==8)
	gcsl_sommerfeld(p,f,gcv,bc,cs);
 
}

void ghostcell::gcsldistro2int(lexer *p, sliceint &f, int ii, int jj, int nn, double dist,  int gcv, int bc, int cs)
{
    i=ii;
	j=jj;
	n=nn;

	gcsl_neumann_int(f,gcv,bc,cs);    
}
