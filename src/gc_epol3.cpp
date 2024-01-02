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

int ghostcell::gceval3(lexer *p, int gcv, int bc, int cs)
{
//	Velocities
    if(gcv==50)
	return 4;
    
    // Parallel
	// Wall
	if((bc==21||bc==22||(bc==7&&awa_lable==0))&&(cs==2||cs==3||cs==1||cs==4)&&(gcv==12||gcv==3))
	return gclabel_w;
	
	if((bc==21||bc==22||(bc==7&&awa_lable==0))&&(cs==2||cs==3||cs==1||cs==4)&&(gcv==112))
	return 5;
	
	if((bc==21||bc==22||(bc==7&&awa_lable==0))&&(cs==2||cs==3||cs==1||cs==4)&&(gcv==116))
	return gclabel_w;
    
    if((bc==21||bc==22||(bc==7&&awa_lable==0))&&(cs==2||cs==3||cs==1||cs==4)&&(gcv==119))
	return 4;
    
    // Topo
    if((bc==5)&&(cs==2||cs==3||cs==1||cs==4)&&(gcv==12||gcv==3))
	return gclabel_wtopo;
	
	if((bc==5)&&(cs==2||cs==3||cs==1||cs==4)&&(gcv==112))
	return 5;
	
	if((bc==5)&&(cs==2||cs==3||cs==1||cs==4)&&(gcv==116))
	return gclabel_wtopo;
    
    if((bc==5)&&(cs==2||cs==3||cs==1||cs==4)&&(gcv==119))
	return 4;
	
	else
	if((bc==21||bc==22||bc==5) && gcv==16)
	return 4;
	
    // Othogonal
	else
	if((bc==21||bc==22||bc==5||(bc==7&&awa_lable==0))&&(cs==6)&&(gcv==12||gcv==3))
	return gclabel_w_orth;
    
    else
	if((bc==21||bc==22||bc==5||(bc==7&&awa_lable==0))&&(cs==5)&&(gcv==12||gcv==3)&&p->A10==6)
	return gclabel_w_orth;
    
    else
	if((bc==21||bc==22||bc==5||(bc==7&&awa_lable==0))&&(cs==5)&&(gcv==12||gcv==3)&&p->A10==5)
	return gclabel_w_orth;

	else
	if((bc==21||bc==22||bc==5||(bc==7&&awa_lable==0))&&(cs==5||cs==6)&&gcv==9)
	return gclabel_vel;

//Inflow	
    else
	if((bc==6 && (gcv==12||gcv==3||gcv==9)))
	return gclabel_w_in;
	
//Outflow
	else
	if((bc==2 && gclabel_outflow==1) && (gcv==12||gcv==3) && (cs==2||cs==3||cs==1||cs==4))
	return 4;
	
	else
	if((bc==2 && gclabel_outflow==1) && (gcv==12||gcv==3) && (cs==5||cs==6))
	return gclabel_w_out;
    
//Patch    
    else
	if((bc==111 || bc==112 || bc==121 || bc==122) && (gcv==12||gcv==3||gcv==9))
	return 4;

//Free Surface

	else
	if((bc==3) && (cs==2||cs==3||cs==1||cs==4) && (gcv==12||gcv==19 || gcv==3))
	return 4;

	else
	if(bc==3 && (cs==5||cs==6)&&(gcv==12||gcv==19 || gcv==3) && p->A10!=3 && p->A10!=5)
	return 5;
    
    else
	if(bc==3 && (cs==5||cs==6)&&(gcv==12||gcv==19 || gcv==3) && p->A10==3)
	return 4;
    
    else
	if(bc==3 && (cs==5||cs==6)&&(gcv==12||gcv==19||gcv==3||gcv==112) && p->A10==5)
	return 9;
	
	else
	if(bc==9 && cs==6 && (gcv==12||gcv==19 || gcv==3))
	return 4;
    
//Omega_sig
    //else
	//if(bc==3 && cs==6 && gcv==17)
	//return 5;
    
    else
	if(bc==21 && cs==5 && gcv==17)
	return 5;

	else
	if((bc!=3 || cs!=6) && gcv==17)
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

void ghostcell::gcdistro3(lexer *p,field& f, int ii, int jj, int kk, int nn, double dist,  int gcv, int bc, int cs)
{
    i=ii;
	j=jj;
	k=kk;
	n=nn;

	bc_label=gceval3(p,gcv,bc,cs);

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
    
    if(bc_label==8)
	kinematic_bed(p,f,dist,gcv,bc,cs);
    
	if(bc_label==9)
	fbvel3(p,f,dist,gcv,bc,cs);
    
    if(bc_label==11)
	dirichlet_ortho_reflect(p,f,dist,gcv,bc,cs);

	if(bc_label==12)
	dirichlet_para_reflect(p,f,dist,gcv,bc,cs);
    
    if(bc_label==99)
	gcb_debug(f,gcv,bc,cs);
}

void ghostcell::gcdistro3V(lexer *p, double *f, int ii, int jj, int kk, int nn, double dist,  int gcv, int bc, int cs)
{
    
}
