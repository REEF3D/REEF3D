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

#include"ghostcell.h"
#include"lexer.h"
#include"slice.h"

int ghostcell::gcsleval4(lexer *p, int gcv, int bc, int cs)
{   
    // general Neuman
    if(gcv==40 || gcv==50 || gcv==1 )
	return 4;
    
    // vertical w
    else
	if(gcv==12 && bc!=1)
	return 4;
    
    else
	if(gcv==12 && bc==1)
	return 5;
    
    
    // pressure 40
    else
	if((bc==3||bc==21||bc==6) && (gcv==41 || gcv==42 || gcv==43 || gcv==44))
	return 4;
    
    else
	if((bc==1||bc==2) && (gcv==41 || gcv==42 || gcv==43 || gcv==44))
	return 5;
    
    else
	if((bc==8) && (gcv==41 || gcv==42 || gcv==43 || gcv==44))
	return 4;
    
    /*
    else
	if((bc==3||bc==21||bc==2) && (gcv==41 || gcv==42 || gcv==43 || gcv==44))
	return 4;
    
    else
	if((bc==1||bc==6) && (gcv==42  || gcv==44))
	return 4;
    
    else
	if(bc==2 && (gcv==41  || gcv==44))
	return 4;*/
    
    // Fifsf 60
    else
    if(gcv==60 && (p->B98==2 || (p->B98==3 && bc!=6)))
    return 4;
    
    // eta 55
    else
    if(gcv==55)// && (p->B98==2 || (p->B98==3 && bc!=6)))
    return 4;
    
    
    else
    if((bc==1||bc==6) && (gcv==52||gcv==54))
	return 4;
    
    else
    if((bc==2||bc==7)&&(gcv==51||gcv==54))
	return 4;
    
    else
    if(bc==8 && (gcv==51||gcv==52||gcv==53||gcv==54) &&p->B99==3)
	return 4;
    
    else
    if(bc==8 && (gcv==51||gcv==52||gcv==53||gcv==54) &&p->B99==4)
	return 8;
    
    else
    if((bc==21||bc==3) && (gcv==51||gcv==52||gcv==53||gcv==54))
	return 4;
    
    // Fifsf
    else
    if(gcv==160 && (p->B98==2 || (p->B98==3 && bc!=6)))
    return 14;
    
    
    // eta 150
    else
    if(gcv==155)// && (p->B98==2 || (p->B98==3 && bc!=6)))
    return 14;
    
    
    else
    if((bc==1||bc==6) && (gcv==152||gcv==154))
	return 14;
    
    
    else
    if((bc==2||bc==7)&&(gcv==151||gcv==154))
	return 14;
    
    else
    if(bc==8 && (gcv==151||gcv==152||gcv==153||gcv==154) &&p->B99==3)
	return 14;
    
    else
    if(bc==8 && (gcv==151||gcv==152||gcv==153||gcv==154) &&p->B99==4)
	return 8;
    
    else
    if((bc==21||bc==3) && (gcv==151||gcv==152||gcv==153||gcv==154))
	return 14;
    
    // Turbulence
	else
	if(gcv==20)
	return 4;
    
    else
	if(gcv==24)
	return 4;
    
    else
	if(gcv==30)
	return 4;
    
    else
    return -1;
}

void ghostcell::gcsldistro4(lexer *p, slice &f, int ii, int jj, int nn, double dist,  int gcv, int bc, int cs)
{
    i=ii;
	j=jj;
	n=nn;
	
	bc_label=gcsleval4(p,gcv,bc,cs);

	if(bc_label==4)
	gcsl_neumann(f,gcv,bc,cs);
    
    if(bc_label==14)
	gcsl_neumann_x(f,gcv,bc,cs);
    
    if(bc_label==5)
	gcsl_noslip(f,gcv,bc,cs);
    
    if(bc_label==8)
	gcsl_sommerfeld(p,f,gcv,bc,cs);
}

void ghostcell::gcsldistro4V(lexer *p, vec2D &f, cpt2D &C, int ii, int jj, int nn, double dist,  int gcv, int bc, int cs, int id)
{
    i=ii;
	j=jj;
	n=nn;
	
	bc_label=gcsleval4(p,gcv,bc,cs);

	if(bc_label==4)
	gcsl_neumannV(f,C,gcv,bc,cs,id);
    
    if(bc_label==14)
	gcsl_neumannV_x(f,C,gcv,bc,cs,id);
    
    if(bc_label==5)
	gcsl_noslipV(f,C,gcv,bc,cs,id);
}

void ghostcell::gcsldistro4int(lexer *p, sliceint &f, int ii, int jj, int nn, double dist,  int gcv, int bc, int cs)
{
    i=ii;
	j=jj;
	n=nn;

	gcsl_neumann_int(f,gcv,bc,cs);    
}
