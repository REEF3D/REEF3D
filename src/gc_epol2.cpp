/*--------------------------------------------------------------------
REEF3D
Copyright 2008-202 Hans Bihs

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

int ghostcell::gceval2(lexer *p, int gcv, int bc, int cs)
{
//	Velocities    if(gcv==50)	return 4;

    // Parallel
	//Wall
	if((bc==21||bc==22||bc==7||bc==6)&&(cs==1||cs==4||cs==5||cs==6)&&(gcv==11||gcv==2))
	return gclabel_v;
	
	if((bc==21||bc==22||bc==7||bc==6)&&(cs==1||cs==4||cs==5||cs==6)&&(gcv==111))
	return 5;
	
	if((bc==21||bc==22||bc==7||bc==6)&&(cs==1||cs==4||cs==5||cs==6)&&(gcv==115))
	return gclabel_v;
    
    if((bc==21||bc==22||bc==7||bc==6)&&(cs==1||cs==4||cs==5||cs==6)&&(gcv==118))
	return 4;
    
    // Topo
    if((bc==5)&&(cs==1||cs==4||cs==5||cs==6)&&(gcv==11||gcv==2))
	return gclabel_vtopo;
	
	if((bc==5)&&(cs==1||cs==4||cs==5||cs==6)&&(gcv==111))
	return 5;
	
	if((bc==5)&&(cs==1||cs==4||cs==5||cs==6)&&(gcv==115))
	return gclabel_vtopo;
    
    if((bc==5)&&(cs==1||cs==4||cs==5||cs==6)&&(gcv==118))
	return 4;
	
	else
	if((bc==21||bc==22||bc==5) && gcv==15)
	return 4;
	
    // Orthogonal
	else
	if((bc==21||bc==22||bc==5||bc==7)&&(cs==2||cs==3)&&(gcv==11||gcv==2))
	return gclabel_v_orth;

	else
	if((bc==21||bc==22||bc==5||bc==7)&&(cs==2||cs==3)&&gcv==8)
	return gclabel_vel;
	

//Inflow
    else
	if((bc==6  && (gcv==11||gcv==2||gcv==8)))
	return gclabel_v_in;
	
//Outflow
	else
	if((bc==2 && gclabel_outflow==1) && (gcv==11||gcv==2) && (cs==1||cs==4||cs==5||cs==6))
	return 4;
	
	else
	if((bc==2 && gclabel_outflow==1) && (gcv==11||gcv==2) && (cs==2||cs==3))
	return gclabel_v_out;
    
//Patch    
    else
	if((bc==111 || bc==112 || bc==121 || bc==122) && (gcv==11||gcv==2||gcv==8))
	return 4;


//Free Surface
	else
	if(bc==3 && (cs==1||cs==4||cs==5||cs==6) && (gcv==11||gcv==18 || gcv==2))
	return 4;

	else
	if(bc==3 && (cs==2||cs==3)&&(gcv==11||gcv==18 || gcv==2))
	return 1;
	
	else
	if((bc==9) && cs==6 && (gcv==11||gcv==18 || gcv==2))
	return 4;
	
// 6DOF
	else
	if(bc==41||bc==42||bc==43)
	return 9;


     else
	if(gcv==999)
	return 99;

	else
	return 0;
}


void ghostcell::gcdistro2(lexer *p,field& f, int ii, int jj, int kk, int nn, double dist,  int gcv, int bc, int cs)
{
    i=ii;
	j=jj;
	k=kk;
	n=nn;

	bc_label=gceval2(p,gcv,bc,cs);


	if(bc_label==1)
	dirichlet_ortho(p,f,dist,gcv,bc,cs);

	if(bc_label==2)
	dirichlet_para(p,f,dist,gcv,bc,cs);

	if(bc_label==3)
	extend(p,f,dist,gcv,bc,cs);

	if(bc_label==4)
	neumann(f,gcv,bc,cs);

	if(bc_label==5)
	noslip(f,dist,gcv,bc,cs);
	
	if(bc_label==6)
	outflow(p,f,gcv,bc,cs);
    
    if(bc_label==7)
	sommerfeld(p,f,gcv,bc,cs);
	
	if(bc_label==9)
	fbvel2(p,f,dist,gcv,bc,cs);
    
    if(bc_label==11)
	dirichlet_ortho_reflect(p,f,dist,gcv,bc,cs);

	if(bc_label==12)
	dirichlet_para_reflect(p,f,dist,gcv,bc,cs);
    
    if(bc_label==99)
	gcb_debug(f,gcv,bc,cs);
}

void ghostcell::gcdistro2V(lexer *p, double *f, int ii, int jj, int kk, int nn, double dist,  int gcv, int bc, int cs){}