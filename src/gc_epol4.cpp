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

int ghostcell::gceval4(lexer *p, int gcv, int bc, int cs)
{

//Level Set

	if((bc==21||bc==22||bc==5||bc==41||bc==42||bc==43||bc==7||bc==8||bc==9||bc==41||bc==221||bc==211||bc==121||bc==111) 
        && (gcv==51 || gcv==52 || gcv==53 || gcv==54))
	return gclabel_lsm;
	
	else
	if((bc==3||bc==221||bc==211||bc==121||bc==111) && (gcv==51 || gcv==52 || gcv==53 || gcv==54))
	return 4;

	else
	if((bc==2||bc==221||bc==211||bc==121||bc==111) && (gcv==51 || gcv==54))
	return 4;
    
    else
	if((bc==1||bc==6||bc==221||bc==211||bc==121||bc==111) && (gcv==52 || gcv==54))
	return 4;

	else
	if(gcv==50)
	return 4;
    
    // inflow
    else
	if((bc==1||bc==221||bc==211||bc==121||bc==111) && (gcv==52 || gcv==54))
	return gclabel_lsm_in;
    
    if((bc==6 ) && (gcv==51 || gcv==52 || gcv==53 || gcv==54) )
	return gclabel_lsm_in;

// Pressure    
	else
	if((bc==21||bc==22||bc==5||bc==3||bc==211||bc==212||bc==112||bc==111) && gcv==40)
	return gclabel_press;
    
    // wavegen
    else
	if(((bc==6&&pressin_lable==0)||bc==211||bc==212||bc==112||bc==111) && gcv==40)
	return gclabel_press;
    
    // awa beach
    else
	if(((bc==7&&awa_lable==0)||bc==211||bc==212||bc==112||bc==111) && gcv==40)
	return gclabel_press;
    
    // inflow
    else
	if(((bc==1&&pressin_lable==0)||bc==211||bc==212||bc==112||bc==111) && gcv==40)
	return gclabel_press_in;
    
    // outflow
    else
	if(( (bc==2&&pressout_lable==0) ||bc==211||bc==212||bc==112||bc==111) && gcv==40)
	return gclabel_press;
    
    // amtosphere
    else
	if(bc==9 && gcv==40)
	return 21;
	
// ro
    else
	if(gcv==1)
	return 4;
	
	// ro
	if(gcv==2 && (cs!=5 && bc!=5 && bc!=21))
	return 4;
	
// Turbulence kin
	else
	if((bc==21||bc==22||bc==5||bc==41||bc==42||bc==43||bc==9) && gcv==20)
	return gclabel_k;

	else
	if((bc==3||bc==2) && (cs!=6||bc!=3)  && gcv==20)
	return 4;
	
	else
	if((cs==6 && bc==3) && gcv==20)
	return 5;
	
	else
	if((bc==6 || bc==7 || bc==8) && gcv==20)
	return 5;

// Turbulence eps
	else
	if((bc==21||bc==22||bc==5||bc==41||bc==42||bc==43||bc==6||bc==7||bc==8||bc==9) && gcv==30)
	return gclabel_e;

	else
	if((bc==3||bc==2) && gcv==30)
	return 4;

	else
	if(bc==1 && gcv==30)
	return 4;


// Turbulence eddyv
    else
	if((bc==21||bc==22||bc==5||bc==41||bc==42||bc==43||bc==9)&&(gcv==24))
	return 4;

	else
	if((bc==3||bc==2||bc==1)&&(gcv==24))
	return 4;
    
    else
	if(bc==1 && gcv==24)
	return 5;
	
	else
	if((cs==6 && bc==3) && (gcv==24))
	return 5;
	
	else
	if((cs!=6 || bc!=3) && (gcv==24))
	return 4;
	
	else
	if((bc==6 || bc==7 || bc==8) && gcv==24)
	return 5;
    
// omega (sigma coordinate)
    // Parallel
	// Wall
	if((bc==21||bc==22||bc==5||(bc==7&&awa_lable==0))&&(cs==2||cs==3||cs==1||cs==4)&&(gcv==12))
	return 4;

    // Othogonal
	else
	if((bc==21||bc==22||bc==5||(bc==7&&awa_lable==0))&&(cs==6)&&(gcv==12))
	return 5;

    //Inflow	
    else
	if((bc==6 && gcv==12))
	return 4;
	
    //Outflow
	else
	if((bc==2 && gclabel_outflow==1) && (gcv==12||gcv==3) && (cs==2||cs==3||cs==1||cs==4))
	return 4;
	
	else
	if((bc==2 && gclabel_outflow==1) && (gcv==12) && (cs==5||cs==6))
	return 5;
    
    //Patch    
    else
	if((bc==111 || bc==112 || bc==121 || bc==122) && (gcv==12))
	return 4;

    //Free Surface
	else
	if((bc==3) && (cs==2||cs==3||cs==1||cs==4) && (gcv==12))
	return 4;

	else
	if(bc==3 && (cs==5||cs==6)&&(gcv==12) && p->A10==5)
	return 4;
	
// VOF
	else
	if((bc==21||bc==22||bc==5||bc==41||bc==42||bc==43||bc==3||bc==6||bc==7||bc==8||bc==9) && (gcv==71 || gcv==72 || gcv==73 || gcv==74))
	return 4;

	else
	if(bc==1&&(gcv==72 || gcv==74))
	return 4;

	else
	if((bc==2)&&(gcv==71 || gcv==74))
	return 4;

	else
	if((bc==21||bc==22||bc==5||bc==41||bc==42||bc==43||bc==1||bc==2||bc==3||bc==6||bc==7||bc==8||bc==9) && gcv==70)
	return 4;

	else
	if((bc==21||bc==22||bc==5||bc==41||bc==42||bc==43||bc==3||bc==6||bc==7||bc==8||bc==9) && gcv==75)
	return 3;
	
// Pk Velocity
	else
	if((bc==21||bc==22||bc==5||bc==41)&&(gcv==101||gcv==102||gcv==103))
	return 5;
	
	//Outflow, Inflow
	else
	if((bc==2||bc==1||bc==7||bc==8||bc==6) && (gcv==101||gcv==102||gcv==103))
	return 4;
	
	// Free Surface Uvel
	else
	if(bc==3 && (cs==2||cs==3||cs==5||cs==6) && gcv==101)
	return 4;

	else
	if(bc==3 && (cs==1||cs==4) && gcv==101)
	return 5;
	
	// Free Surface Vvel
	else
	if(bc==3 && (cs==1||cs==4||cs==5||cs==6) && gcv==102)
	return 4;

	else
	if(bc==3 && (cs==2||cs==3) && gcv==102)
	return 5;
	
	// Free Surface Wvel
	else
	if(bc==3 && (cs==1||cs==4||cs==2||cs==3) && gcv==103)
	return 4;

	else
	if(bc==3 && (cs==5||cs==6) && gcv==103)
	return 5;
	
// Suspended Sediment
	else
	if(gcv==60)
	return 4;

// Heat
	else
	if(gcv==80 && ((p->H61==1 && cs==1) || (p->H62==1 && cs==2) || (p->H63==1 && cs==3) 
                || (p->H64==1 && cs==4) || (p->H65==1 && cs==5) || (p->H66==1 && cs==6)))
	return 61;
    
    else
	if(gcv==80)
	return 4;

	else
	if(gcv==81)
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
    
// Potential Waves
	else
	if((bc==21||bc==22||bc==5||bc==41||bc==42||bc==43||bc==7||bc==8||bc==9)&&(cs!=5)&&(gcv==250))
	return 4;

	else
	if((bc==2||bc==1||bc==7||bc==6)&&(gcv==250))
	return 4;
	
	else
	if(bc==3 && (cs!=6) && gcv==250)
	return 4;
    
    else
	if(gcv==999)
	return 99;
    
// NHFLOW
    else
	if((bc==21||bc==22||bc==5||bc==3||(bc==2&&pressout_lable==0)||bc==6||(bc==7&&awa_lable==0)||bc==211||bc==212||bc==112||bc==111) && cs!=6 && gcv==540)
	return 4;
    
    else
	if(bc==3 && cs==6 && gcv==540)
	return 11;
    
    
	else
	return 0;
}


void ghostcell::gcdistro4(lexer *p, field &f, int ii, int jj, int kk, int nn, double dist,  int gcv, int bc, int cs)
{
    i=ii;
	j=jj;
	k=kk;
	n=nn;

	bc_label=gceval4(p,gcv,bc,cs);

    if(bc_label==22)
	lsm(p,f,dist,gcv,bc,cs);
    
	if(bc_label==3)
	extend(p,f,dist,gcv,bc,cs);

	if(bc_label==4)
	neumann(f,gcv,bc,cs);
    
    if(bc_label==5)
	noslip(f,dist,gcv,bc,cs);

	if(bc_label==6)
	extend(p,f,dist,gcv,bc,cs);

	if(bc_label==7)
	potentialbc(p,f,bc,cs);

	if(bc_label==8)
	neumann_press(p,f,dist,gcv,bc,cs);
	
	if(bc_label==9)
	fbpress(p,f,dist,gcv,bc,cs);
	
	if(bc_label==10)
	gravity_press(p,f,dist,gcv,bc,cs);
    
    if(bc_label==11)
    nhpress(p,f,dist,gcv,bc,cs);

	if(bc_label==21)
	atmosphere(p,f,gcv,bc,cs);
    
    if(bc_label==61)
	heatbc(p,f,gcv,bc,cs);
    
    if(bc_label==99)
	gcb_debug(f,gcv,bc,cs);
}

void ghostcell::gcdistro4V(lexer *p, double *f, int ii, int jj, int kk, int nn, double dist,  int gcv, int bc, int cs)
{
    
}

void ghostcell::gcdistro4vec(lexer *p, fdm* a, vec &vec, int ii, int jj, int kk, double dist,  int gcv, int bc, int cs, int id)
{
    i=ii;
	j=jj;
	k=kk;

	bc_label=gceval4(p,gcv,bc,cs);
	
	if(bc_label==22)
	gcV_lsm(p,vec,dist,gcv,bc,cs,id);
    
    if(bc_label==3)
	extendV(p,a,vec,dist,gcv,bc,cs);
    
	if(bc_label==4)
	gcV_neumann(vec,gcv,bc,cs,id);

}



