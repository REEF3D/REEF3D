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

#include"ddweno3_f_nug.h"
#include"lexer.h"
#include"field.h"
#include"slice.h"
#include"ghostcell.h"
#include"vec.h"
#include"cpt.h"

ddweno3_f_nug::ddweno3_f_nug(lexer* pp):weno3_nug_func(pp)
{
    p=pp;
}

ddweno3_f_nug::~ddweno3_f_nug()
{
}

double ddweno3_f_nug::ddwenox(field& f, double uw)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    uf=0;
    
	grad=0.0;

	if(uw>=0.0)
	{
	iqmin(p,f);
	is_min_x();
	weight_min_x();

	grad = w1x*(qfx[IP][uf][0][0]*q2 + qfx[IP][uf][0][1]*q3)
    
         + w2x*(qfx[IP][uf][1][0]*q2 - qfx[IP][uf][1][1]*q1);
	}

	if(uw<0.0)
	{
	iqmax(p,f);
	is_max_x();
	weight_max_x();
    
    
	grad = w1x*(qfx[IP][uf][2][0]*q2 - qfx[IP][uf][2][1]*q3)
    
         + w2x*(qfx[IP][uf][3][0]*q1 + qfx[IP][uf][3][1]*q2);
	}

	return grad;
}

double ddweno3_f_nug::ddwenoy(field& f, double uw)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    vf=0;
    
	grad=0.0;

	if(uw>=0.0)
	{
	jqmin(p,f);
	is_min_y();
	weight_min_y();
	
	grad = w1y*(qfy[JP][vf][0][0]*q2 + qfy[JP][vf][0][1]*q3)
    
         + w2y*(qfy[JP][vf][1][0]*q2 - qfy[JP][vf][1][1]*q1);
	}

	if(uw<0.0)
	{
	jqmax(p,f);
	is_max_y();
	weight_max_y();
	
	grad = w1y*(qfy[JP][vf][2][0]*q2 - qfy[JP][vf][2][1]*q3)
    
         + w2y*(qfy[JP][vf][3][0]*q1 + qfy[JP][vf][3][1]*q2);
	}

	return grad;
}

double ddweno3_f_nug::ddwenoz(field& f, double uw)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    wf=0;
    
    
	grad=0.0;

	if(uw>=0.0)
	{
	kqmin(p,f);
	is_min_z();
	weight_min_z();

	grad = w1z*(qfz[KP][wf][0][0]*q2 + qfz[KP][wf][0][1]*q3)
    
         + w2z*(qfz[KP][wf][1][0]*q2 - qfz[KP][wf][1][1]*q1);
	}


	if(uw<0.0)
	{
	kqmax(p,f);
	is_max_z();
	weight_max_z();
    
	grad = w1z*(qfz[KP][wf][2][0]*q2 - qfz[KP][wf][2][1]*q3)
    
         + w2z*(qfz[KP][wf][3][0]*q1 + qfz[KP][wf][3][1]*q2);
	}
    
	return grad;
}



double ddweno3_f_nug::dswenox(slice& f, double uw)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    uf=0;
    
	grad=0.0;

	if(uw>=0.0)
	{
	isqmin(p,f);
	is_min_x();
	weight_min_x();

	grad = w1x*(qfx[IP][uf][0][0]*q2 + qfx[IP][uf][0][1]*q3)
    
         + w2x*(qfx[IP][uf][1][0]*q2 - qfx[IP][uf][1][1]*q1);
	}

	if(uw<0.0)
	{
	isqmax(p,f);
	is_max_x();
	weight_max_x();

	grad = w1x*(qfx[IP][uf][2][0]*q2 - qfx[IP][uf][2][1]*q3)
    
         + w2x*(qfx[IP][uf][3][0]*q1 + qfx[IP][uf][3][1]*q2);
	}

	return grad;
}

double ddweno3_f_nug::dswenoy(slice& f, double uw)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    vf=0;
    
	grad=0.0;

	if(uw>=0.0)
	{
	jsqmin(p,f);
	is_min_y();
	weight_min_y();
	
	grad = w1y*(qfy[JP][vf][0][0]*q2 + qfy[JP][vf][0][1]*q3)
    
         + w2y*(qfy[JP][vf][1][0]*q2 - qfy[JP][vf][1][1]*q1);
	}

	if(uw<0.0)
	{
	jsqmax(p,f);
	is_max_y();
	weight_max_y();
	
	grad = w1y*(qfy[JP][vf][2][0]*q2 - qfy[JP][vf][2][1]*q3)
    
         + w2y*(qfy[JP][vf][3][0]*q1 + qfy[JP][vf][3][1]*q2);
	}

	return grad;
}

void ddweno3_f_nug::iqmin(lexer *p,field& f)
{	
	q1 = (f(i-1,j,k)-f(i-2,j,k))/DX[IM2];
	q2 = (f(i,j,k)-f(i-1,j,k))/DX[IM1];
	q3 = (f(i+1,j,k)-f(i,j,k))/DX[IP];
}

void ddweno3_f_nug::jqmin(lexer *p,field& f)
{
	q1 = (f(i,j-1,k)-f(i,j-2,k))/DY[JM2];
	q2 = (f(i,j,k)-f(i,j-1,k))/DY[JM1];
	q3 = (f(i,j+1,k)-f(i,j,k))/DY[JP];
}

void ddweno3_f_nug::kqmin(lexer *p,field& f)
{
	q1 = (f(i,j,k-1)-f(i,j,k-2))/DZ[KM2];
	q2 = (f(i,j,k)-f(i,j,k-1))/DZ[KM1];
	q3 = (f(i,j,k+1)-f(i,j,k))/DZ[KP];
}

void ddweno3_f_nug::iqmax(lexer *p,field& f)
{
    q1 = (f(i,j,k)-f(i-1,j,k))/DX[IM1];
	q2 = (f(i+1,j,k)-f(i,j,k))/DX[IP];
	q3 = (f(i+2,j,k)-f(i+1,j,k))/DX[IP1];
}

void ddweno3_f_nug::jqmax(lexer *p,field& f)
{
	q1 = (f(i,j,k)-f(i,j-1,k))/DY[JM1];
	q2 = (f(i,j+1,k)-f(i,j,k))/DY[JP];
	q3 = (f(i,j+2,k)-f(i,j+1,k))/DY[JP1];
}

void ddweno3_f_nug::kqmax(lexer *p,field& f)
{
	q1 = (f(i,j,k)-f(i,j,k-1))/DZ[KM1];
	q2 = (f(i,j,k+1)-f(i,j,k))/DZ[KP];
	q3 = (f(i,j,k+2)-f(i,j,k+1))/DZ[KP1];
}



void ddweno3_f_nug::isqmin(lexer *p,slice& f)
{	
	q1 = (f(i-1,j)-f(i-2,j))/DX[IM2];
	q2 = (f(i,j)-f(i-1,j))/DX[IM1];
	q3 = (f(i+1,j)-f(i,j))/DX[IP];
}

void ddweno3_f_nug::jsqmin(lexer *p,slice& f)
{
	q1 = (f(i,j-1)-f(i,j-2))/DY[JM2];
	q2 = (f(i,j)-f(i,j-1))/DY[JM1];
	q3 = (f(i,j+1)-f(i,j))/DY[JP];
}

void ddweno3_f_nug::isqmax(lexer *p,slice& f)
{
    q1 = (f(i,j)-f(i-1,j))/DX[IM1];
	q2 = (f(i+1,j)-f(i,j))/DX[IP];
	q3 = (f(i+2,j)-f(i+1,j))/DX[IP1];
}

void ddweno3_f_nug::jsqmax(lexer *p,slice& f)
{
	q1 = (f(i,j)-f(i,j-1))/DY[JM1];
	q2 = (f(i,j+1)-f(i,j))/DY[JP];
	q3 = (f(i,j+2)-f(i,j+1))/DY[JP1];
}
