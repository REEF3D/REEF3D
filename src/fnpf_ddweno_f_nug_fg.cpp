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

#include"fnpf_ddweno_f_nug.h"
#include"fdm_fnpf.h"
#include"lexer.h"
#include"field.h"
#include"slice.h"
#include"ghostcell.h"
#include"vec.h"
#include"cpt.h"

double fnpf_ddweno_f_nug::ddwenox(field& f, double uw)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    uf=0;
    
	grad=0.0;

	if(uw>0.0)
	{
	iqmin(f);
	is_min_x();
	weight_min_x();

	grad = w1x*(q4 + qfx[IP][uf][0][0]*(q3-q4) - qfx[IP][uf][0][1]*(q5-q4))
    
         + w2x*(q3 + qfx[IP][uf][1][0]*(q4-q3) - qfx[IP][uf][1][1]*(q2-q3))
          
         + w3x*(q2 + qfx[IP][uf][2][0]*(q1-q2) + qfx[IP][uf][2][1]*(q3-q2));
	}

	if(uw<0.0)
	{
	iqmax(f);
	is_max_x();
	weight_max_x();
    
    
	grad = w1x*(q4 + qfx[IP][uf][3][0]*(q3-q4) + qfx[IP][uf][3][1]*(q5-q4))
    
         + w2x*(q3 + qfx[IP][uf][4][0]*(q2-q3) - qfx[IP][uf][4][1]*(q4-q3))
          
         + w3x*(q2 + qfx[IP][uf][5][0]*(q3-q2) - qfx[IP][uf][5][1]*(q1-q2));
	}

	return grad;
}

double fnpf_ddweno_f_nug::ddwenoy(field& f, double uw)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    vf=0;
    
	grad=0.0;

	if(uw>0.0)
	{
	jqmin(f);
	is_min_y();
	weight_min_y();
	
	grad = w1y*(q4 + qfy[JP][vf][0][0]*(q3-q4) - qfy[JP][vf][0][1]*(q5-q4))
    
         + w2y*(q3 + qfy[JP][vf][1][0]*(q4-q3) - qfy[JP][vf][1][1]*(q2-q3))
          
         + w3y*(q2 + qfy[JP][vf][2][0]*(q1-q2) + qfy[JP][vf][2][1]*(q3-q2));
	}

	if(uw<0.0)
	{
	jqmax(f);
	is_max_y();
	weight_max_y();
	
	grad = w1y*(q4 + qfy[JP][vf][3][0]*(q3-q4) + qfy[JP][vf][3][1]*(q5-q4))
    
         + w2y*(q3 + qfy[JP][vf][4][0]*(q2-q3) - qfy[JP][vf][4][1]*(q4-q3))
          
         + w3y*(q2 + qfy[JP][vf][5][0]*(q3-q2) - qfy[JP][vf][5][1]*(q1-q2));
	}

	return grad;
}

double fnpf_ddweno_f_nug::ddwenoz(field& f, double uw)
{
    DX = p->DXP;
    DY = p->DYP;
    DZ = p->DZP;
    wf=0;
    
    
	grad=0.0;

	if(uw>0.0)
	{
	kqmin(f);
	is_min_z();
	weight_min_z();

	grad = w1z*(q4 + qfz[KP][wf][0][0]*(q3-q4) - qfz[KP][wf][0][1]*(q5-q4))
    
         + w2z*(q3 + qfz[KP][wf][1][0]*(q4-q3) - qfz[KP][wf][1][1]*(q2-q3))
          
         + w3z*(q2 + qfz[KP][wf][2][0]*(q1-q2) + qfz[KP][wf][2][1]*(q3-q2));
	}


	if(uw<0.0)
	{
	kqmax(f);
	is_max_z();
	weight_max_z();
    
	grad = w1z*(q4 + qfz[KP][wf][3][0]*(q3-q4) + qfz[KP][wf][3][1]*(q5-q4))
    
         + w2z*(q3 + qfz[KP][wf][4][0]*(q2-q3) - qfz[KP][wf][4][1]*(q4-q3))
          
         + w3z*(q2 + qfz[KP][wf][5][0]*(q3-q2) - qfz[KP][wf][5][1]*(q1-q2));
	}
    
	return grad;
}

void fnpf_ddweno_f_nug::iqmin(field& f)
{	
	q1 = (f(i-2,j,k)-f(i-3,j,k))/DX[IM3];
	q2 = (f(i-1,j,k)-f(i-2,j,k))/DX[IM2];
	q3 = (f(i,j,k)-f(i-1,j,k))/DX[IM1];
	q4 = (f(i+1,j,k)-f(i,j,k))/DX[IP];
	q5 = (f(i+2,j,k)-f(i+1,j,k))/DX[IP1];
}

void fnpf_ddweno_f_nug::jqmin(field& f)
{
	q1 = (f(i,j-2,k)-f(i,j-3,k))/DY[JM3];
	q2 = (f(i,j-1,k)-f(i,j-2,k))/DY[JM2];
	q3 = (f(i,j,k)-f(i,j-1,k))/DY[JM1];
	q4 = (f(i,j+1,k)-f(i,j,k))/DY[JP];
	q5 = (f(i,j+2,k)-f(i,j+1,k))/DY[JP1];
}

void fnpf_ddweno_f_nug::kqmin(field& f)
{
	q1 = (f(i,j,k-2)-f(i,j,k-3))/DZ[KM3];
	q2 = (f(i,j,k-1)-f(i,j,k-2))/DZ[KM2];
	q3 = (f(i,j,k)-f(i,j,k-1))/DZ[KM1];
	q4 = (f(i,j,k+1)-f(i,j,k))/DZ[KP];
	q5 = (f(i,j,k+2)-f(i,j,k+1))/DZ[KP1];
}

void fnpf_ddweno_f_nug::iqmax(field& f)
{
    q1 = (f(i-1,j,k)-f(i-2,j,k))/DX[IM2];
	q2 = (f(i,j,k)-f(i-1,j,k))/DX[IM1];
	q3 = (f(i+1,j,k)-f(i,j,k))/DX[IP];
	q4 = (f(i+2,j,k)-f(i+1,j,k))/DX[IP1];
	q5 = (f(i+3,j,k)-f(i+2,j,k))/DX[IP2];
}

void fnpf_ddweno_f_nug::jqmax(field& f)
{
	q1 = (f(i,j-1,k)-f(i,j-2,k))/DY[JM2];
	q2 = (f(i,j,k)-f(i,j-1,k))/DY[JM1];
	q3 = (f(i,j+1,k)-f(i,j,k))/DY[JP];
	q4 = (f(i,j+2,k)-f(i,j+1,k))/DY[JP1];
	q5 = (f(i,j+3,k)-f(i,j+2,k))/DY[JP2];
}

void fnpf_ddweno_f_nug::kqmax(field& f)
{
	q1 = (f(i,j,k-1)-f(i,j,k-2))/DZ[KM2];
	q2 = (f(i,j,k)-f(i,j,k-1))/DZ[KM1];
	q3 = (f(i,j,k+1)-f(i,j,k))/DZ[KP];
	q4 = (f(i,j,k+2)-f(i,j,k+1))/DZ[KP1];
	q5 = (f(i,j,k+3)-f(i,j,k+2))/DZ[KP2];
}

