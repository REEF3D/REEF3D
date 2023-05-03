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
--------------------------------------------------------------------*/

#include"ghostcell.h"
#include<math.h>

int ghostcell::gceval4a(lexer *p, int gcv, int bc, int cs)
{

	//topo
	if((bc==21||bc==22||bc==5||bc==3||bc==6||bc==7||bc==8)&&(cs==5||cs==6)&&(gcv==151 || gcv==152 || gcv==153))
	return 75;

	else
	if((bc==21||bc==22||bc==5||bc==3||bc==6||bc==7||bc==8)&&(cs!=5&&cs!=6)&&(gcv==151 || gcv==152 || gcv==153))
	return 74;

	else
	if((bc==2&&gcv==151) || (bc==1&&gcv==152))
	return 74;

	else
	if(gcv==150 || gcv==154)
	return 74;
    
    else
	if(gcv==159)
	return 79;

	//topo for bedload
	else
	if((bc==21||bc==22||bc==5||bc==3||bc==6||bc==7||bc==8)&&(cs==5||cs==6)&&(gcv==161 || gcv==162 || gcv==163))
	return 75;

	else
	if((bc==21||bc==22||bc==5||bc==3||bc==6||bc==7||bc==8)&&(cs!=5&&cs!=6)&&(gcv==161 || gcv==162 || gcv==163))
	return 74;

	else
	if((bc==2&&gcv==161) || (bc==1&&gcv==162))
	return 74;
	
	// fb
	else
	if(gcv==50)
	return 75;
    
    
//Level Set	
    
    else
	if((bc==21||bc==22||bc==5||bc==41||bc==6||bc==7||bc==8||bc==9) && (gcv==51 || gcv==52 || gcv==53 || gcv==54))
	return 74;
    
	else
	if((bc==3) && (gcv==51 || gcv==52 || gcv==53 || gcv==54))
	return 74;

	else
	if(bc==1&&(gcv==52 || gcv==54))
	return 74;

	else
	if((bc==2)&&(gcv==51 || gcv==54))
	return 74;

	else
	if(gcv==50)
	return 74;
	
	// porosity
	else
	if(gcv==1)
	return 75;

	else
	return 0;
}

void ghostcell::gcdistro4a(lexer *p,field& f, int ii, int jj, int kk, int nn, double dist,  int gcv, int bc, int cs)
{
    i=ii;
	j=jj;
	k=kk;
	n=nn;
	
	cs = fabs(cs);
    

	bc_label=gceval4a(p,gcv,bc,cs);

	if(bc_label==74 || bc_label==75)
	neumann_all(f,gcv,bc,cs);
    
    if(bc_label==79)
    extend(p,f,dist,gcv,bc,cs);
}

void ghostcell::gcdistro4avec(lexer *p, fdm* a, vec &vec, int ii, int jj, int kk, double dist,  int gcv, int bc, int cs, int id)
{
    i=ii;
	j=jj;
	k=kk;
	
	cs = fabs(cs);

	bc_label=gceval4a(p,gcv,bc,cs);
    
	if(bc_label==74 || bc_label==75)
	gcV_neumann_all(vec,gcv,bc,cs,id);
}

void ghostcell::gcdistro6vec(lexer *p, fdm* a, vec &vec, int ii, int jj, int kk, double dist,  int gcv, int bc, int cs, int id)
{
    i=ii;
	j=jj;
	k=kk;
	
	cs = fabs(cs);

	bc_label=gceval4a(p,gcv,bc,cs);
    
	if(bc_label==74 || bc_label==75)
	gcV_neumann_6V(vec,gcv,bc,cs,id);
}

