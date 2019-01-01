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
-------------------------------------------------------------------
--------------------------------------------------------------------*/

#include"ghostcell.h"

int ghostcell::gceval1(lexer *p, int gcv, int bc, int cs)
{
//	Velocities

    // Parallel
	//Wall
	if((bc==21||bc==22||bc==7)&&(cs==2||cs==3||cs==5||cs==6)&&(gcv==10||gcv==1||gcv==20))
	return gclabel_u;
	
	if((bc==21||bc==22||bc==7)&&(cs==2||cs==3||cs==5||cs==6)&&(gcv==110))
	return 5;
	
	if((bc==21||bc==22||bc==7)&&(cs==2||cs==3||cs==5||cs==6)&&(gcv==114))
	return gclabel_u;
    
    // Topo
    if((bc==5)&&(cs==2||cs==3||cs==5||cs==6)&&(gcv==10||gcv==1||gcv==20))
	return gclabel_utopo;
	
	if((bc==5)&&(cs==2||cs==3||cs==5||cs==6)&&(gcv==110))
	return 5;
	
	if((bc==5)&&(cs==2||cs==3||cs==5||cs==6)&&(gcv==114))
	return gclabel_utopo;
	
	else
	if((bc==21||bc==22||bc==5)&& gcv==14)
	return 4;
	
    // Orthogonal
	else
	if((bc==21||bc==22||bc==5||(bc==7&&awa_label==0))&&(cs==1||cs==4)&&(gcv==10||gcv==20||gcv==1))
	return gclabel_u_orth;

	else
	if((bc==21||bc==22||bc==5||bc==7)&&(cs==1||cs==4)&&gcv==7)
	return gclabel_vel;

//Inflow
	else
	if((bc==1 && (gcv==20)))
	return 4;
    
    else
	if((bc==6 && (cs==1||cs==4) && (gcv==10||gcv==20||gcv==1||gcv==7)))
	return gclabel_u_in;
    
    
	
//Outflow
	else
	if((bc==2 && gclabel_outflow==1) && (gcv==10||gcv==20||gcv==1) && (cs==2||cs==3||cs==5||cs==6))
	return 4;
	
	else
	if((bc==2 && gclabel_outflow==1) && (gcv==10||gcv==20||gcv==1) && (cs==1||cs==4))
	return gclabel_u_out;

//Free Surface
	else
	if(bc==3 && (cs==2||cs==3||cs==5||cs==6) && (gcv==10||gcv==20||gcv==17||gcv==1))
	return 4;

	else
	if(bc==3 && (cs==1||cs==4)&&(gcv==10||gcv==20||gcv==17||gcv==1))
	return gclabel_u_orth;
	
	else
	if((bc==9) && cs==6 && (gcv==10||gcv==20||gcv==17||gcv==1))
	return 4;
	
// 6DOF
	else
	if(bc==41||bc==42||bc==43)
	return 11;

//PISO Velcorr
//Wall
	else
	if((bc==21||bc==22||bc==5||bc==7)&&(cs==2||cs==3||cs==5||cs==6)&&(gcv==17))
	return gclabel_u;

	else
	if((bc==21||bc==22||bc==5||bc==7)&&(cs==1||cs==4)&&(gcv==17))
	return 5;

	else
	if((bc==2||bc==1||bc==6) && (gcv==17))
	return 4;

	else
	return 0;
}

void ghostcell::gcdistro1(lexer *p,field& f, int ii, int jj, int kk, int nn, double dist,  int gcv, int bc, int cs)
{
    i=ii;
	j=jj;
	k=kk;
	n=nn;

	bc_label=gceval1(p,gcv,bc,cs);

	if(bc_label==1)
	dirichlet_ortho(p,f,dist,gcv,bc,cs);

	if(bc_label==2)
	dirichlet_para(f,dist,gcv,bc,cs);

	if(bc_label==3)
	extend(f,dist,gcv,bc,cs);

	if(bc_label==4)
	neumann(f,gcv,bc,cs);

	if(bc_label==5)
	noslip(f,dist,gcv,bc,cs);
	
	if(bc_label==6)
	outflow(p,f,gcv,bc,cs);
    
    if(bc_label==7)
	sommerfeld(p,f,gcv,bc,cs);
	
	if(bc_label==11)
	fbvel1(p,f,dist,gcv,bc,cs);
}
