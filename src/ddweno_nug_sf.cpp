/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

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

#include"ddweno_nug_sf.h"
#include"fdm.h"
#include"lexer.h"
#include"ghostcell.h"
#include"vec.h"
#include"cpt.h"

ddweno_nug_sf::ddweno_nug_sf(lexer* pp):weno_nug_func(pp)
{
    p=pp;
    
    modus=0;
}

ddweno_nug_sf::~ddweno_nug_sf()
{
}

double ddweno_nug_sf::ddwenox(fdm* a, field &f, double uw)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    uf=0;
    
    int check=0;
    
	grad=0.0;

	if(uw>0.0)
	{
    if(modus==0)
	iqmin0(a,f);
    
    if(modus==1)
	iqmax1(a,f);
    
    if(modus==2)
	iqmax2(a,f);
    
    if(modus==3)
	iqmax3(a,f);
    
	is_min_x();
	weight_min_x();
    
	grad = w1x*(q4 + qfx[IP][uf][0][0]*(q3-q4) - qfx[IP][uf][0][1]*(q5-q4))
    
         + w2x*(q3 + qfx[IP][uf][1][0]*(q4-q3) - qfx[IP][uf][1][1]*(q2-q3))
          
         + w3x*(q2 + qfx[IP][uf][2][0]*(q1-q2) + qfx[IP][uf][2][1]*(q3-q2));
	}

	if(uw<0.0)
	{
    if(modus==0)
	iqmax0(a,f);
    
    if(modus==1)
	iqmax1(a,f);
    
    if(modus==2)
	iqmax2(a,f);
    
    if(modus==3)
	iqmax3(a,f);
    
	is_max_x();
    //weight_max_sfcheck_x(a);
	weight_max_x();
    
	grad = w1x*(q4 + qfx[IP][uf][3][0]*(q3-q4) + qfx[IP][uf][3][1]*(q5-q4))
    
         + w2x*(q3 + qfx[IP][uf][4][0]*(q2-q3) - qfx[IP][uf][4][1]*(q4-q3))
          
         + w3x*(q2 + qfx[IP][uf][5][0]*(q3-q2) - qfx[IP][uf][5][1]*(q1-q2));
	}
    
	return grad;
}

double ddweno_nug_sf::ddwenoy(fdm* a, field &f, double uw)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    vf=0;
    
	grad=0.0;

	if(uw>0.0)
	{
    if(modus==0)
	jqmin0(a,f);
    
    if(modus==1)
	jqmin1(a,f);
    
    if(modus==2)
	jqmin2(a,f);
    
    if(modus==3)
	jqmin3(a,f);
    
	is_min_y();
	weight_min_y();
	
	grad = w1y*(q4 + qfy[JP][vf][0][0]*(q3-q4) - qfy[JP][vf][0][1]*(q5-q4))
    
         + w2y*(q3 + qfy[JP][vf][1][0]*(q4-q3) - qfy[JP][vf][1][1]*(q2-q3))
          
         + w3y*(q2 + qfy[JP][vf][2][0]*(q1-q2) + qfy[JP][vf][2][1]*(q3-q2));
	}

	if(uw<0.0)
	{
    if(modus==0)
	jqmax0(a,f);
    
    if(modus==1)
	jqmax1(a,f);
    
    if(modus==2)
	jqmax2(a,f);
    
    if(modus==3)
	jqmax3(a,f);
    
	is_max_y();
	weight_max_y();
	
	grad = w1y*(q4 + qfy[JP][vf][3][0]*(q3-q4) + qfy[JP][vf][3][1]*(q5-q4))
    
         + w2y*(q3 + qfy[JP][vf][4][0]*(q2-q3) - qfy[JP][vf][4][1]*(q4-q3))
          
         + w3y*(q2 + qfy[JP][vf][5][0]*(q3-q2) - qfy[JP][vf][5][1]*(q1-q2));
	}

	return grad;
}

double ddweno_nug_sf::ddwenoz(fdm* a, field &f, double uw)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    wf=0;

	grad=0.0;

	if(uw>0.0)
	{
    if(modus==0)
	kqmin0(a,f);
    
    if(modus==1)
	kqmin1(a,f);
    
    if(modus==2)
	kqmin2(a,f);
    
    if(modus==3)
	kqmin3(a,f);
    
	is_min_z();
	weight_min_z();

	grad = w1z*(q4 + qfz[KP][wf][0][0]*(q3-q4) - qfz[KP][wf][0][1]*(q5-q4))
    
         + w2z*(q3 + qfz[KP][wf][1][0]*(q4-q3) - qfz[KP][wf][1][1]*(q2-q3))
          
         + w3z*(q2 + qfz[KP][wf][2][0]*(q1-q2) + qfz[KP][wf][2][1]*(q3-q2));
	}


	if(uw<0.0)
	{
    if(modus==0)
	kqmax0(a,f);
    
    if(modus==1)
	kqmax1(a,f);
    
    if(modus==2)
	kqmax2(a,f);
    
    if(modus==3)
	kqmax3(a,f);
    
	is_max_z();
	weight_max_z();
    
	grad = w1z*(q4 + qfz[KP][wf][3][0]*(q3-q4) + qfz[KP][wf][3][1]*(q5-q4))
    
         + w2z*(q3 + qfz[KP][wf][4][0]*(q2-q3) - qfz[KP][wf][4][1]*(q4-q3))
          
         + w3z*(q2 + qfz[KP][wf][5][0]*(q3-q2) - qfz[KP][wf][5][1]*(q1-q2));
	}
    
	return grad;
}

void ddweno_nug_sf::iqmin0(fdm* a, field &f)
{   
    q1 = (f.V[Im2JK] - f.V[Im3JK])/DX[IM3];
	q2 = (f.V[Im1JK] - f.V[Im2JK])/DX[IM2];
	q3 = (f.V[IJK]   - f.V[Im1JK])/DX[IM1];
	q4 = (f.V[Ip1JK] - f.V[IJK]  )/DX[IP];
	q5 = (f.V[Ip2JK] - f.V[Ip1JK])/DX[IP1];
}

void ddweno_nug_sf::jqmin0(fdm* a, field &f)
{
    q1 = (f.V[IJm2K] - f.V[IJm3K])/DY[JM3];
	q2 = (f.V[IJm1K] - f.V[IJm2K])/DY[JM2];
	q3 = (f.V[IJK]   - f.V[IJm1K])/DY[JM1];
	q4 = (f.V[IJp1K] - f.V[IJK]  )/DY[JP];
	q5 = (f.V[IJp2K] - f.V[IJp1K])/DY[JP1];
}

void ddweno_nug_sf::kqmin0(fdm* a, field &f)
{
    q1 = (f.V[IJKm2] - f.V[IJKm3])/DZ[KM3];
	q2 = (f.V[IJKm1] - f.V[IJKm2])/DZ[KM2];
	q3 = (f.V[IJK]   - f.V[IJKm1])/DZ[KM1];
	q4 = (f.V[IJKp1] - f.V[IJK]  )/DZ[KP];
	q5 = (f.V[IJKp2] - f.V[IJKp1])/DZ[KP1];
}

void ddweno_nug_sf::iqmax0(fdm* a, field &f)
{
    q1 = (f.V[Im1JK] - f.V[Im2JK])/DX[IM2];
    q2 = (f.V[IJK]   - f.V[Im1JK])/DX[IM1];
    q3 = (f.V[Ip1JK] - f.V[IJK]  )/DX[IP];
    q4 = (f.V[Ip2JK] - f.V[Ip1JK])/DX[IP1];
    q5 = (f.V[Ip3JK] - f.V[Ip2JK])/DX[IP2];
}

void ddweno_nug_sf::jqmax0(fdm* a, field &f)
{
    q1 = (f.V[IJm1K] - f.V[IJm2K])/DY[JM2];
    q2 = (f.V[IJK]   - f.V[IJm1K])/DY[JM1];
    q3 = (f.V[IJp1K] - f.V[IJK]  )/DY[JP];
    q4 = (f.V[IJp2K] - f.V[IJp1K])/DY[JP1];
    q5 = (f.V[IJp3K] - f.V[IJp2K])/DY[JP2];
}

void ddweno_nug_sf::kqmax0(fdm* a, field &f)
{
    q1 = (f.V[IJKm1] - f.V[IJKm2])/DZ[KM2];
    q2 = (f.V[IJK]   - f.V[IJKm1])/DZ[KM1];
    q3 = (f.V[IJKp1] - f.V[IJK]  )/DZ[KP];
    q4 = (f.V[IJKp2] - f.V[IJKp1])/DZ[KP1];
    q5 = (f.V[IJKp3] - f.V[IJKp2])/DZ[KP2];
}
    
void ddweno_nug_sf::iqmin1(fdm* a, field &f)
{   
    q1=q2=q3=q4=q5=0.0;
    
    if(a->solid(i-2,j,k)>0.0 && a->topo(i-2,j,k)>0.0 && a->solid(i-3,j,k)>0.0 && a->topo(i-3,j,k)>0.0)
	q1 = (f.V[Im2JK] - f.V[Im3JK])/DX[IM3];
    
    if(a->solid(i-1,j,k)>0.0 && a->topo(i-1,j,k)>0.0 && a->solid(i-2,j,k)>0.0 && a->topo(i-2,j,k)>0.0)
	q2 = (f.V[Im1JK] - f.V[Im2JK])/DX[IM2];
    
    if(a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0 && a->solid(i-1,j,k)>0.0 && a->topo(i-1,j,k)>0.0)
	q3 = (f.V[IJK]   - f.V[Im1JK])/DX[IM1];
    
    if(a->solid(i+1,j,k)>0.0 && a->topo(i+1,j,k)>0.0 && a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0)
	q4 = (f.V[Ip1JK] - f.V[IJK]  )/DX[IP];
    
    if(a->solid(i+2,j,k)>0.0 && a->topo(i+2,j,k)>0.0 && a->solid(i+1,j,k)>0.0 && a->topo(i+1,j,k)>0.0)
	q5 = (f.V[Ip2JK] - f.V[Ip1JK])/DX[IP1];

}

void ddweno_nug_sf::jqmin1(fdm* a, field &f)
{
    q1=q2=q3=q4=q5=0.0;
    
    if(a->solid(i,j-2,k)>0.0 && a->topo(i,j-2,k)>0.0 && a->solid(i,j-3,k)>0.0 && a->topo(i,j-3,k)>0.0)
	q1 = (f.V[IJm2K] - f.V[IJm3K])/DY[JM3];
    
    if(a->solid(i,j-1,k)>0.0 && a->topo(i,j-1,k)>0.0 && a->solid(i,j-2,k)>0.0 && a->topo(i,j-2,k)>0.0)
	q2 = (f.V[IJm1K] - f.V[IJm2K])/DY[JM2];
    
    if(a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0 && a->solid(i-1,j,k)>0.0 && a->topo(i-1,j,k)>0.0)
	q3 = (f.V[IJK]   - f.V[IJm1K])/DY[JM1];
    
    if(a->solid(i,j+1,k)>0.0 && a->topo(i,j+1,k)>0.0 && a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0)
	q4 = (f.V[IJp1K] - f.V[IJK]  )/DY[JP];
    
    if(a->solid(i,j+2,k)>0.0 && a->topo(i,j+2,k)>0.0 && a->solid(i,j+1,k)>0.0 && a->topo(i,j+1,k)>0.0)
	q5 = (f.V[IJp2K] - f.V[IJp1K])/DY[JP1];
}

void ddweno_nug_sf::kqmin1(fdm* a, field &f)
{
    q1=q2=q3=q4=q5=0.0;
    
    if(a->solid(i,j,k-2)>0.0 && a->topo(i,j,k-2)>0.0 && a->solid(i,j,k-3)>0.0 && a->topo(i,j,k-3)>0.0)
	q1 = (f.V[IJKm2] - f.V[IJKm3])/DZ[KM3];
    
    if(a->solid(i,j,k-1)>0.0 && a->topo(i,j,k-1)>0.0 && a->solid(i,j,k-2)>0.0 && a->topo(i,j,k-2)>0.0)
	q2 = (f.V[IJKm1] - f.V[IJKm2])/DZ[KM2];
    
    if(a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0 && a->solid(i,j,k-1)>0.0 && a->topo(i,j,k-1)>0.0)
	q3 = (f.V[IJK]   - f.V[IJKm1])/DZ[KM1];
    
    if(a->solid(i,j,k+1)>0.0 && a->topo(i,j,k+1)>0.0 && a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0)
	q4 = (f.V[IJKp1] - f.V[IJK]  )/DZ[KP];
    
    if(a->solid(i,j,k+2)>0.0 && a->topo(i,j,k+2)>0.0 && a->solid(i,j,k+1)>0.0 && a->topo(i,j,k+1)>0.0)
	q5 = (f.V[IJKp2] - f.V[IJKp1])/DZ[KP1];
}

void ddweno_nug_sf::iqmax1(fdm* a, field &f)
{
    q1=q2=q3=q4=q5=0.0;
    
    if(a->solid(i-1,j,k)>0.0 && a->topo(i-1,j,k)>0.0 && a->solid(i-2,j,k)>0.0 && a->topo(i-2,j,k)>0.0)
	q1 = (f.V[Im1JK] - f.V[Im2JK])/DX[IM2];
    
    if(a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0 && a->solid(i-1,j,k)>0.0 && a->topo(i-1,j,k)>0.0)
    q2 = (f.V[IJK]   - f.V[Im1JK])/DX[IM1];
    
    if(a->solid(i+1,j,k)>0.0 && a->topo(i+1,j,k)>0.0 && a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0)
    q3 = (f.V[Ip1JK] - f.V[IJK]  )/DX[IP];
    
    if(a->solid(i+2,j,k)>0.0 && a->topo(i+2,j,k)>0.0 && a->solid(i+1,j,k)>0.0 && a->topo(i+1,j,k)>0.0)
    q4 = (f.V[Ip2JK] - f.V[Ip1JK])/DX[IP1];
    
    if(a->solid(i+3,j,k)>0.0 && a->topo(i+3,j,k)>0.0 && a->solid(i+2,j,k)>0.0 && a->topo(i+2,j,k)>0.0)
    q5 = (f.V[Ip3JK] - f.V[Ip2JK])/DX[IP2];
}

void ddweno_nug_sf::jqmax1(fdm* a, field &f)
{
    q1=q2=q3=q4=q5=0.0;
    
    if(a->solid(i,j-1,k)>0.0 && a->topo(i,j-1,k)>0.0 && a->solid(i,j-2,k)>0.0 && a->topo(i,j-2,k)>0.0)
	q1 = (f.V[IJm1K] - f.V[IJm2K])/DY[JM2];
    
    if(a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0 && a->solid(i,j-1,k)>0.0 && a->topo(i,j-1,k)>0.0)
    q2 = (f.V[IJK]   - f.V[IJm1K])/DY[JM1];
    
    if(a->solid(i,j+1,k)>0.0 && a->topo(i,j+1,k)>0.0 && a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0)
    q3 = (f.V[IJp1K] - f.V[IJK]  )/DY[JP];
    
    if(a->solid(i,j+2,k)>0.0 && a->topo(i,j+2,k)>0.0 && a->solid(i,j+1,k)>0.0 && a->topo(i,j+1,k)>0.0)
    q4 = (f.V[IJp2K] - f.V[IJp1K])/DY[JP1];
    
    if(a->solid(i,j+3,k)>0.0 && a->topo(i,j+3,k)>0.0 && a->solid(i,j+2,k)>0.0 && a->topo(i,j+2,k)>0.0)
    q5 = (f.V[IJp3K] - f.V[IJp2K])/DY[JP2];
}

void ddweno_nug_sf::kqmax1(fdm* a, field &f)
{
    q1=q2=q3=q4=q5=0.0;
    
    if(a->solid(i,j,k-1)>0.0 && a->topo(i,j,k-1)>0.0 && a->solid(i,j,k-2)>0.0 && a->topo(i,j,k-2)>0.0)
	q1 = (f.V[IJKm1] - f.V[IJKm2])/DZ[KM2];
    
    if(a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0 && a->solid(i,j,k-1)>0.0 && a->topo(i,j,k-1)>0.0)
    q2 = (f.V[IJK]   - f.V[IJKm1])/DZ[KM1];
    
    if(a->solid(i,j,k+1)>0.0 && a->topo(i,j,k+1)>0.0 && a->solid(i,j,k)>0.0 && a->topo(i,j,k)>0.0)
    q3 = (f.V[IJKp1] - f.V[IJK]  )/DZ[KP];
    
    if(a->solid(i,j,k+2)>0.0 && a->topo(i,j,k+2)>0.0 && a->solid(i,j,k+1)>0.0 && a->topo(i,j,k+1)>0.0)
    q4 = (f.V[IJKp2] - f.V[IJKp1])/DZ[KP1];
    
    if(a->solid(i,j,k+3)>0.0 && a->topo(i,j,k+3)>0.0 && a->solid(i,j,k+2)>0.0 && a->topo(i,j,k+2)>0.0)
    q5 = (f.V[IJKp3] - f.V[IJKp2])/DZ[KP2];
}

void ddweno_nug_sf::iqmin2(fdm* a, field &f)
{
    q1=q2=q3=q4=q5=0.0;
   
    if(a->fbh4(i-2,j,k)  < 0.5 && a->fbh4(i-3,j,k)  < 0.5)
    q1 = (f.V[Im2JK] - f.V[Im3JK])/DX[IM3];
    
    if(a->fbh4(i-1,j,k)  < 0.5 && a->fbh4(i-2,j,k)  < 0.5)
	q2 = (f.V[Im1JK] - f.V[Im2JK])/DX[IM2];
    
    if(a->fbh4(i,j,k)  < 0.5 && a->fbh4(i-1,j,k)  < 0.5)
	q3 = (f.V[IJK]   - f.V[Im1JK])/DX[IM1];
    
    if(a->fbh4(i+1,j,k)  < 0.5 && a->fbh4(i,j,k)  < 0.5)
	q4 = (f.V[Ip1JK] - f.V[IJK]  )/DX[IP];
    
    if(a->fbh4(i+2,j,k)  < 0.5 && a->fbh4(i+1,j,k)  < 0.5)
	q5 = (f.V[Ip2JK] - f.V[Ip1JK])/DX[IP1];
}

void ddweno_nug_sf::jqmin2(fdm* a, field &f)
{
    q1=q2=q3=q4=q5=0.0;
    
    if(a->fbh4(i,j-2,k)  < 0.5 && a->fbh4(i,j-3,k)  < 0.5)
    q1 = (f.V[IJm2K] - f.V[IJm3K])/DY[JM3];
    
    if(a->fbh4(i,j-1,k)  < 0.5 && a->fbh4(i,j-2,k)  < 0.5)
	q2 = (f.V[IJm1K] - f.V[IJm2K])/DY[JM2];
    
    if(a->fbh4(i,j,k)  < 0.5 && a->fbh4(i,j-1,k)  < 0.5)
	q3 = (f.V[IJK]   - f.V[IJm1K])/DY[JM1];
    
    if(a->fbh4(i,j+1,k)  < 0.5 && a->fbh4(i,j,k)  < 0.5)
	q4 = (f.V[IJp1K] - f.V[IJK]  )/DY[JP];
    
    if(a->fbh4(i,j+2,k)  < 0.5 && a->fbh4(i,j+1,k)  < 0.5)
	q5 = (f.V[IJp2K] - f.V[IJp1K])/DY[JP1];
}

void ddweno_nug_sf::kqmin2(fdm* a, field &f)
{
    q1=q2=q3=q4=q5=0.0;
    
    if(a->fbh4(i,j,k-2)  < 0.5 && a->fbh4(i,j,k-3)  < 0.5)
    q1 = (f.V[IJKm2] - f.V[IJKm3])/DZ[KM3];
    
    if(a->fbh4(i,j,k-1)  < 0.5 && a->fbh4(i,j,k-2)  < 0.5)
	q2 = (f.V[IJKm1] - f.V[IJKm2])/DZ[KM2];
    
    if(a->fbh4(i,j,k)  < 0.5 && a->fbh4(i,j,k-1)  < 0.5)
	q3 = (f.V[IJK]   - f.V[IJKm1])/DZ[KM1];
    
    if(a->fbh4(i,j,k+1)  < 0.5 && a->fbh4(i,j,k)  < 0.5)
	q4 = (f.V[IJKp1] - f.V[IJK]  )/DZ[KP];
    
    if(a->fbh4(i,j,k+2)  < 0.5 && a->fbh4(i,j,k+1)  < 0.5)
	q5 = (f.V[IJKp2] - f.V[IJKp1])/DZ[KP1];
}

void ddweno_nug_sf::iqmax2(fdm* a, field &f)
{
    q1=q2=q3=q4=q5=0.0;
    
    if(a->fbh4(i-1,j,k)  < 0.5 && a->fbh4(i-2,j,k)  < 0.5)
    q1 = (f.V[Im1JK] - f.V[Im2JK])/DX[IM2];
    
    if(a->fbh4(i,j,k)  < 0.5 && a->fbh4(i-1,j,k)  < 0.5)
    q2 = (f.V[IJK]   - f.V[Im1JK])/DX[IM1];
    
    if(a->fbh4(i+1,j,k)  < 0.5 && a->fbh4(i,j,k)  < 0.5)
    q3 = (f.V[Ip1JK] - f.V[IJK]  )/DX[IP];
    
    if(a->fbh4(i+2,j,k)  < 0.5 && a->fbh4(i+1,j,k)  < 0.5)
    q4 = (f.V[Ip2JK] - f.V[Ip1JK])/DX[IP1];
    
    if(a->fbh4(i+3,j,k)  < 0.5 && a->fbh4(i+2,j,k)  < 0.5)
    q5 = (f.V[Ip3JK] - f.V[Ip2JK])/DX[IP2];
}

void ddweno_nug_sf::jqmax2(fdm* a, field &f)
{
    q1=q2=q3=q4=q5=0.0;
    
    if(fabs(a->fb(i,j-1,k))  > 0.6*p->DXM && fabs(a->fb(i,j-2,k))  > 0.6*p->DXM)
    q1 = (f.V[IJm1K] - f.V[IJm2K])/DY[JM2];
    
    if(fabs(a->fb(i,j,k))  > 0.6*p->DXM && fabs(a->fb(i,j-1,k))  > 0.6*p->DXM)
    q2 = (f.V[IJK]   - f.V[IJm1K])/DY[JM1];
    
    if(fabs(a->fb(i,j+1,k))  > 0.6*p->DXM && fabs(a->fb(i,j,k))  > 0.6*p->DXM)
    q3 = (f.V[IJp1K] - f.V[IJK]  )/DY[JP];
    
    if(fabs(a->fb(i,j+2,k))  > 0.6*p->DXM && fabs(a->fb(i,j+1,k))  > 0.6*p->DXM)
    q4 = (f.V[IJp2K] - f.V[IJp1K])/DY[JP1];
    
    if(fabs(a->fb(i,j+3,k))  > 0.6*p->DXM && fabs(a->fb(i,j+2,k))  > 0.6*p->DXM)
    q5 = (f.V[IJp3K] - f.V[IJp2K])/DY[JP2];
}

void ddweno_nug_sf::kqmax2(fdm* a, field &f)
{
    q1=q2=q3=q4=q5=0.0;
    
    if(a->fbh4(i,j-1,k)  < 0.5 && a->fbh4(i,j-2,k)  < 0.5)
    q1 = (f.V[IJKm1] - f.V[IJKm2])/DZ[KM2];
    
    if(a->fbh4(i,j,k)  < 0.5 && a->fbh4(i,j-1,k)  < 0.5)
    q2 = (f.V[IJK]   - f.V[IJKm1])/DZ[KM1];
    
    if(a->fbh4(i,j+1,k)  < 0.5 && a->fbh4(i,j,k)  < 0.5)
    q3 = (f.V[IJKp1] - f.V[IJK]  )/DZ[KP];
    
    if(a->fbh4(i,j+2,k)  < 0.5 && a->fbh4(i,j+1,k)  < 0.5)
    q4 = (f.V[IJKp2] - f.V[IJKp1])/DZ[KP1];
    
    if(a->fbh4(i,j+3,k)  < 0.5 && a->fbh4(i,j+2,k)  < 0.5)
    q5 = (f.V[IJKp3] - f.V[IJKp2])/DZ[KP2];
}


void ddweno_nug_sf::iqmin3(fdm* a, field &f)
{   
    q1 = a->fbh5(i-3,j,k)*(f.V[Im2JK] - f.V[Im3JK])/DX[IM3];
	q2 = a->fbh5(i-2,j,k)*(f.V[Im1JK] - f.V[Im2JK])/DX[IM2];
	q3 = a->fbh5(i-1,j,k)*(f.V[IJK]   - f.V[Im1JK])/DX[IM1];
	q4 = a->fbh5(i,j,k)*(f.V[Ip1JK] - f.V[IJK]  )/DX[IP];
	q5 = a->fbh5(i+1,j,k)*(f.V[Ip2JK] - f.V[Ip1JK])/DX[IP1];
}

void ddweno_nug_sf::jqmin3(fdm* a, field &f)
{
    q1 = a->fbh5(i,j-3,k)*(f.V[IJm2K] - f.V[IJm3K])/DY[JM3];
	q2 = a->fbh5(i,j-2,k)*(f.V[IJm1K] - f.V[IJm2K])/DY[JM2];
	q3 = a->fbh5(i,j-1,k)*(f.V[IJK]   - f.V[IJm1K])/DY[JM1];
	q4 = a->fbh5(i,j,k)*(f.V[IJp1K] - f.V[IJK]  )/DY[JP];
	q5 = a->fbh5(i,j+1,k)*(f.V[IJp2K] - f.V[IJp1K])/DY[JP1];
}

void ddweno_nug_sf::kqmin3(fdm* a, field &f)
{
    q1 = a->fbh5(i,j,k-3)*(f.V[IJKm2] - f.V[IJKm3])/DZ[KM3];
	q2 = a->fbh5(i,j,k-2)*(f.V[IJKm1] - f.V[IJKm2])/DZ[KM2];
	q3 = a->fbh5(i,j,k-1)*(f.V[IJK]   - f.V[IJKm1])/DZ[KM1];
	q4 = a->fbh5(i,j,k)*(f.V[IJKp1] - f.V[IJK]  )/DZ[KP];
	q5 = a->fbh5(i,j,k+1)*(f.V[IJKp2] - f.V[IJKp1])/DZ[KP1];
}

void ddweno_nug_sf::iqmax3(fdm* a, field &f)
{
    q1 = a->fbh5(i-1,j,k)*(f.V[Im1JK] - f.V[Im2JK])/DX[IM2];
    q2 = a->fbh5(i,j,k)*(f.V[IJK]   - f.V[Im1JK])/DX[IM1];
    q3 = a->fbh5(i+1,j,k)*(f.V[Ip1JK] - f.V[IJK]  )/DX[IP];
    q4 = a->fbh5(i+2,j,k)*(f.V[Ip2JK] - f.V[Ip1JK])/DX[IP1];
    q5 = a->fbh5(i+3,j,k)*(f.V[Ip3JK] - f.V[Ip2JK])/DX[IP2];
}

void ddweno_nug_sf::jqmax3(fdm* a, field &f)
{
    q1 = a->fbh5(i,j-1,k)*(f.V[IJm1K] - f.V[IJm2K])/DY[JM2];
    q2 = a->fbh5(i,j,k)*(f.V[IJK]   - f.V[IJm1K])/DY[JM1];
    q3 = a->fbh5(i,j+1,k)*(f.V[IJp1K] - f.V[IJK]  )/DY[JP];
    q4 = a->fbh5(i,j+2,k)*(f.V[IJp2K] - f.V[IJp1K])/DY[JP1];
    q5 = a->fbh5(i,j+3,k)*(f.V[IJp3K] - f.V[IJp2K])/DY[JP2];
}

void ddweno_nug_sf::kqmax3(fdm* a, field &f)
{
    q1 = a->fbh5(i,j,k-1)*(f.V[IJKm1] - f.V[IJKm2])/DZ[KM2];
    q2 = a->fbh5(i,j,k)*(f.V[IJK]   - f.V[IJKm1])/DZ[KM1];
    q3 = a->fbh5(i,j,k+1)*(f.V[IJKp1] - f.V[IJK]  )/DZ[KP];
    q4 = a->fbh5(i,j,k+2)*(f.V[IJKp2] - f.V[IJKp1])/DZ[KP1];
    q5 = a->fbh5(i,j,k+3)*(f.V[IJKp3] - f.V[IJKp2])/DZ[KP2];
}
