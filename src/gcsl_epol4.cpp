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
	return 4;
    
    else
	if((bc==8) && (gcv==41 || gcv==42 || gcv==43 || gcv==44))
	return 4;
    
    else
	if((bc==211 || bc==212 || bc==112 || bc==111) && (gcv==41 || gcv==42 || gcv==43 || gcv==44))
	return 4;
    
    // Fifsf 60
    else
    if((cs==2 || cs==3) && gcv==60)
    return 4;
    
    else
    if(cs==1 && p->B98<=2 && gcv==60)
    return 4;
    
    else
    if(cs==4 && p->B99<=2 && gcv==60)
    return 4;
    
    // eta 55
    else
    if(gcv==55)
    return 4;
    
    //Patch eta / Hx / Hy
    else
	if((bc==221 || bc==211 || bc==121 || bc==111) && (gcv==50||gcv==51||gcv==52||gcv==53||gcv==54))
	return 4;
    
    //Hp / eta
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
    if((bc==21||bc==3)&&(gcv==51||gcv==52||gcv==53||gcv==54))
	return 4;
    
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
    
    
    
    // eta
    else
    if((bc==1||bc==6) && (gcv==52||gcv==54))
	return 4;
    
    else
    if((bc==2||bc==7) && (gcv==51||gcv==54))
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
    if(((bc!=1 && bc!=6) || p->B98<=2) && cs==1 && gcv==160)
    return 14;
    
    else
    if(((bc!=7 && bc!=8) || p->B99<=2) && cs==4 && gcv==160)
    return 14;
    
    
    // eta 150
    else
    if(((bc!=1 && bc!=6 && p->B98<=2))  &&gcv==155)
    return 14;
    
    else
    if((bc==2||bc==7) && gcv==155)
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
	if(gcv==24 && bc!=1)
	return 4;
    
    else
	if(gcv==30)
	return 4;
    
    // Potential Ini
	else
	if((bc==21||bc==22||bc==5||bc==41||bc==42||bc==43||bc==9)&&(gcv==49))
	return 4;

	else
	if((bc==2||bc==1||bc==6||bc==7||bc==8)&&(gcv==49))
	return 7;
	
	else
	if(bc==3 && gcv==49)
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
    
    if(bc_label==7)
	gcsl_potentialbc(p,f,bc,cs);
    
    if(bc_label==8)
	gcsl_sommerfeld(p,f,gcv,bc,cs);
}


void ghostcell::gcsldistro4int(lexer *p, sliceint &f, int ii, int jj, int nn, double dist,  int gcv, int bc, int cs)
{
    i=ii;
	j=jj;
	n=nn;

	gcsl_neumann_int(f,gcv,bc,cs);    
}

void ghostcell::gcsldistro4Vint(lexer *p, int *f, int ii, int jj, int nn, double dist,  int gcv, int bc, int cs)
{
    i=ii;
	j=jj;
	n=nn;

	gcsl_neumann_V_int(p,f,gcv,bc,cs);    
}
